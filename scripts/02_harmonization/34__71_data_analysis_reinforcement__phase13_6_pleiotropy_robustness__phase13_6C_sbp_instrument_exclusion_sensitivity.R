options(stringsAsFactors = FALSE)

message("Running Phase13.6C SBP instrument-level exclusion sensitivity...")

write_tsv <- function(x, file) {
  write.table(
    x,
    file = file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = ""
  )
}

safe_read <- function(path) {
  if (!file.exists(path)) return(NULL)
  if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    read.delim(gzfile(path), check.names = FALSE, stringsAsFactors = FALSE)
  } else {
    read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  }
}

as_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

fmt <- function(x, digits = 4) {
  ifelse(is.na(x), NA, formatC(as.numeric(x), digits = digits, format = "fg", flag = "#"))
}

root <- normalizePath(".", winslash = "/", mustWork = TRUE)

out_dir <- file.path(
  root,
  "71_data_analysis_reinforcement",
  "phase13_6_pleiotropy_robustness"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

master_status_file <- file.path(
  root,
  "71_data_analysis_reinforcement",
  "phase13_master_status.tsv"
)

paired_file <- file.path(
  root,
  "71_data_analysis_reinforcement",
  "phase13_4_empirical_covariance",
  "phase13_4B_paired_internal_sbp_component_input.tsv"
)

radial_file <- file.path(
  root,
  "22_sbp_robustness",
  "phase4_2_sbp_radial_outlier_snps.tsv"
)

measured_iop_file <- file.path(
  root,
  "36_iop_external_validation_inputs",
  "mr_input",
  "SBP__IOP.external_mr_input.tsv.gz"
)

stopifnot(file.exists(paired_file))

paired <- safe_read(paired_file)
radial <- safe_read(radial_file)
measured_iop <- safe_read(measured_iop_file)

required_cols <- c(
  "SNP",
  "chr",
  "pos",
  "beta_exposure",
  "se_exposure",
  "beta_outcome_nonIOP",
  "se_outcome_nonIOP",
  "beta_outcome_IOP",
  "se_outcome_IOP"
)

missing_required <- setdiff(required_cols, names(paired))
if (length(missing_required) > 0) {
  stop("Missing required columns in paired input: ", paste(missing_required, collapse = ", "))
}

paired$beta_exposure <- as_num(paired$beta_exposure)
paired$se_exposure <- as_num(paired$se_exposure)
paired$beta_outcome_nonIOP <- as_num(paired$beta_outcome_nonIOP)
paired$se_outcome_nonIOP <- as_num(paired$se_outcome_nonIOP)
paired$beta_outcome_IOP <- as_num(paired$beta_outcome_IOP)
paired$se_outcome_IOP <- as_num(paired$se_outcome_IOP)

paired <- paired[
  is.finite(paired$beta_exposure) &
    is.finite(paired$se_exposure) &
    is.finite(paired$beta_outcome_nonIOP) &
    is.finite(paired$se_outcome_nonIOP) &
    is.finite(paired$beta_outcome_IOP) &
    is.finite(paired$se_outcome_IOP) &
    paired$se_outcome_nonIOP > 0 &
    paired$se_outcome_IOP > 0,
]

ivw_fit <- function(dat, outcome = c("nonIOP", "IOP")) {
  outcome <- match.arg(outcome)
  
  bx <- as_num(dat$beta_exposure)
  
  if (outcome == "nonIOP") {
    by <- as_num(dat$beta_outcome_nonIOP)
    sy <- as_num(dat$se_outcome_nonIOP)
  } else {
    by <- as_num(dat$beta_outcome_IOP)
    sy <- as_num(dat$se_outcome_IOP)
  }
  
  keep <- is.finite(bx) & is.finite(by) & is.finite(sy) & sy > 0
  bx <- bx[keep]
  by <- by[keep]
  sy <- sy[keep]
  
  w <- 1 / (sy^2)
  n <- length(bx)
  
  beta <- sum(w * bx * by) / sum(w * bx^2)
  residual <- by - beta * bx
  Q <- sum(w * residual^2)
  Q_df <- n - 1
  phi <- ifelse(Q_df > 0, max(1, Q / Q_df), NA_real_)
  se <- sqrt(phi / sum(w * bx^2))
  z <- beta / se
  p <- 2 * pnorm(-abs(z))
  
  data.frame(
    n_instruments = n,
    beta = beta,
    se = se,
    z = z,
    p = p,
    Q = Q,
    Q_df = Q_df,
    Q_p = ifelse(Q_df > 0, pchisq(Q, df = Q_df, lower.tail = FALSE), NA_real_),
    phi = phi,
    stringsAsFactors = FALSE
  )
}

contrast_fit <- function(dat) {
  non <- ivw_fit(dat, "nonIOP")
  iop <- ivw_fit(dat, "IOP")
  
  beta_diff <- iop$beta - non$beta
  se_diff_r0 <- sqrt(iop$se^2 + non$se^2)
  z <- beta_diff / se_diff_r0
  p <- 2 * pnorm(-abs(z))
  
  data.frame(
    n_instruments = min(non$n_instruments, iop$n_instruments),
    beta_nonIOP = non$beta,
    se_nonIOP = non$se,
    p_nonIOP = non$p,
    beta_IOP = iop$beta,
    se_IOP = iop$se,
    p_IOP = iop$p,
    beta_difference_IOP_minus_nonIOP = beta_diff,
    se_difference_working_r0 = se_diff_r0,
    z_difference = z,
    p_difference_working_r0 = p,
    Q_nonIOP = non$Q,
    Q_p_nonIOP = non$Q_p,
    phi_nonIOP = non$phi,
    Q_IOP = iop$Q,
    Q_p_IOP = iop$Q_p,
    phi_IOP = iop$phi,
    direction_preserved = beta_diff > 0,
    stringsAsFactors = FALSE
  )
}

# Full-data diagnostics for Q contribution and leverage.
full_non <- ivw_fit(paired, "nonIOP")
full_iop <- ivw_fit(paired, "IOP")
full_contrast <- contrast_fit(paired)

compute_snp_diagnostics <- function(dat) {
  bx <- as_num(dat$beta_exposure)
  
  w_non <- 1 / (as_num(dat$se_outcome_nonIOP)^2)
  res_non <- as_num(dat$beta_outcome_nonIOP) - full_non$beta * bx
  Q_non <- w_non * res_non^2
  
  w_iop <- 1 / (as_num(dat$se_outcome_IOP)^2)
  res_iop <- as_num(dat$beta_outcome_IOP) - full_iop$beta * bx
  Q_iop <- w_iop * res_iop^2
  
  leverage_non <- w_non * bx^2 / sum(w_non * bx^2)
  leverage_iop <- w_iop * bx^2 / sum(w_iop * bx^2)
  
  # Leave-one-SNP contrast influence.
  full_diff <- full_contrast$beta_difference_IOP_minus_nonIOP
  loo_diff <- rep(NA_real_, nrow(dat))
  
  for (i in seq_len(nrow(dat))) {
    d <- dat[-i, , drop = FALSE]
    loo_diff[i] <- contrast_fit(d)$beta_difference_IOP_minus_nonIOP
  }
  
  data.frame(
    SNP = dat$SNP,
    chr = dat$chr,
    pos = dat$pos,
    Q_contribution_nonIOP = Q_non,
    Q_contribution_IOP = Q_iop,
    Q_contribution_max = pmax(Q_non, Q_iop, na.rm = TRUE),
    leverage_nonIOP = leverage_non,
    leverage_IOP = leverage_iop,
    leverage_max = pmax(leverage_non, leverage_iop, na.rm = TRUE),
    leave_one_contrast_difference = loo_diff,
    leave_one_delta_from_full_contrast = loo_diff - full_diff,
    abs_leave_one_delta_from_full_contrast = abs(loo_diff - full_diff),
    stringsAsFactors = FALSE
  )
}

diag <- compute_snp_diagnostics(paired)

# Radial outliers.
radial_outlier_snps <- character(0)

if (!is.null(radial) && nrow(radial) > 0 && "SNP" %in% names(radial)) {
  if ("radial_outlier_flag" %in% names(radial)) {
    flag <- as.character(radial$radial_outlier_flag)
    radial_outlier_snps <- unique(radial$SNP[toupper(flag) %in% c("TRUE", "T", "1", "YES")])
  } else {
    radial_outlier_snps <- unique(radial$SNP)
  }
}

# Measured IOP association annotation.
if (!is.null(measured_iop) && nrow(measured_iop) > 0 && "SNP" %in% names(measured_iop)) {
  iop_cols <- names(measured_iop)
  
  beta_iop_col <- if ("beta_outcome" %in% iop_cols) "beta_outcome" else NA_character_
  se_iop_col <- if ("se_outcome" %in% iop_cols) "se_outcome" else NA_character_
  p_iop_col <- if ("pval_outcome" %in% iop_cols) "pval_outcome" else NA_character_
  
  m <- measured_iop[, c("SNP", beta_iop_col, se_iop_col, p_iop_col), drop = FALSE]
  names(m) <- c("SNP", "beta_measured_IOP", "se_measured_IOP", "p_measured_IOP")
  m$beta_measured_IOP <- as_num(m$beta_measured_IOP)
  m$se_measured_IOP <- as_num(m$se_measured_IOP)
  m$p_measured_IOP <- as_num(m$p_measured_IOP)
  
  diag <- merge(diag, m, by = "SNP", all.x = TRUE, sort = FALSE)
} else {
  diag$beta_measured_IOP <- NA_real_
  diag$se_measured_IOP <- NA_real_
  diag$p_measured_IOP <- NA_real_
}

# Scenario definitions.
top_n <- function(x, prop) {
  n <- max(1, ceiling(length(x) * prop))
  names(sort(x, decreasing = TRUE))[seq_len(n)]
}

diag_named_Q <- diag$Q_contribution_max
names(diag_named_Q) <- diag$SNP

diag_named_lev <- diag$leverage_max
names(diag_named_lev) <- diag$SNP

diag_named_inf <- diag$abs_leave_one_delta_from_full_contrast
names(diag_named_inf) <- diag$SNP

diag_named_iop_abs <- abs(diag$beta_measured_IOP)
names(diag_named_iop_abs) <- diag$SNP

scenario_list <- list(
  full_391_reference = character(0),
  exclude_radial_outliers_any_component = radial_outlier_snps,
  exclude_top_1pct_Q_contribution_any_component = top_n(diag_named_Q, 0.01),
  exclude_top_5pct_Q_contribution_any_component = top_n(diag_named_Q, 0.05),
  exclude_top_1pct_leverage_any_component = top_n(diag_named_lev, 0.01),
  exclude_top_1pct_contrast_influence = top_n(diag_named_inf, 0.01)
)

if (any(is.finite(diag$p_measured_IOP))) {
  scenario_list$exclude_measured_IOP_p_lt_0p05 <- unique(diag$SNP[is.finite(diag$p_measured_IOP) & diag$p_measured_IOP < 0.05])
  scenario_list$exclude_measured_IOP_p_lt_0p01 <- unique(diag$SNP[is.finite(diag$p_measured_IOP) & diag$p_measured_IOP < 0.01])
  scenario_list$exclude_top_5pct_abs_measured_IOP_association <- top_n(diag_named_iop_abs[is.finite(diag_named_iop_abs)], 0.05)
}

scenario_list$conservative_local_outlier_union <- unique(c(
  scenario_list$exclude_radial_outliers_any_component,
  scenario_list$exclude_top_1pct_Q_contribution_any_component,
  scenario_list$exclude_top_1pct_leverage_any_component,
  scenario_list$exclude_top_1pct_contrast_influence
))

# Run scenarios.
scenario_results <- list()
removed_rows <- list()

for (sc in names(scenario_list)) {
  remove_snps <- unique(scenario_list[[sc]])
  keep <- !(paired$SNP %in% remove_snps)
  dat <- paired[keep, , drop = FALSE]
  
  fit <- contrast_fit(dat)
  
  scenario_results[[length(scenario_results) + 1]] <- data.frame(
    scenario = sc,
    n_removed = length(unique(remove_snps)),
    removed_fraction = length(unique(remove_snps)) / nrow(paired),
    fit,
    delta_difference_from_full = fit$beta_difference_IOP_minus_nonIOP - full_contrast$beta_difference_IOP_minus_nonIOP,
    relative_change_from_full = (
      fit$beta_difference_IOP_minus_nonIOP - full_contrast$beta_difference_IOP_minus_nonIOP
    ) / full_contrast$beta_difference_IOP_minus_nonIOP,
    manuscript_safe_interpretation = ifelse(
      fit$direction_preserved,
      "Positive IOP-minus-nonIOP contrast direction preserved after exclusion.",
      "Contrast direction not preserved after exclusion; treat as robustness caveat."
    ),
    claim_boundary = "Instrument-exclusion sensitivity only; does not upgrade inference beyond hypothesis-generating interpretation.",
    stringsAsFactors = FALSE
  )
  
  if (length(remove_snps) > 0) {
    tmp <- data.frame(
      scenario = sc,
      SNP = remove_snps,
      stringsAsFactors = FALSE
    )
    tmp <- merge(tmp, diag, by = "SNP", all.x = TRUE, sort = FALSE)
    removed_rows[[length(removed_rows) + 1]] <- tmp
  }
}

scenario_results <- do.call(rbind, scenario_results)
removed_instruments <- if (length(removed_rows) > 0) do.call(rbind, removed_rows) else data.frame()

scenario_file <- file.path(out_dir, "phase13_6C_sbp_instrument_exclusion_sensitivity.tsv")
removed_file <- file.path(out_dir, "phase13_6C_removed_instruments_by_scenario.tsv")
diag_file <- file.path(out_dir, "phase13_6C_sbp_instrument_diagnostics.tsv")

write_tsv(scenario_results, scenario_file)
write_tsv(removed_instruments, removed_file)
write_tsv(diag, diag_file)

# Summary.
non_reference <- scenario_results[scenario_results$scenario != "full_391_reference", , drop = FALSE]

summary <- data.frame(
  metric = c(
    "reference_contrast",
    "reference_p_working_r0",
    "n_scenarios_excluding_reference",
    "all_exclusion_scenarios_direction_preserved",
    "minimum_exclusion_contrast",
    "maximum_exclusion_contrast",
    "maximum_absolute_delta_from_full",
    "largest_change_scenario",
    "measured_IOP_association_available",
    "radial_outlier_snps_removed",
    "new_primary_MR_estimates_created",
    "exploratory_sensitivity_created",
    "claim_level",
    "claim_upgrade_allowed"
  ),
  value = c(
    full_contrast$beta_difference_IOP_minus_nonIOP,
    full_contrast$p_difference_working_r0,
    nrow(non_reference),
    all(non_reference$direction_preserved),
    min(non_reference$beta_difference_IOP_minus_nonIOP, na.rm = TRUE),
    max(non_reference$beta_difference_IOP_minus_nonIOP, na.rm = TRUE),
    max(abs(non_reference$delta_difference_from_full), na.rm = TRUE),
    non_reference$scenario[which.max(abs(non_reference$delta_difference_from_full))],
    any(is.finite(diag$p_measured_IOP)),
    length(radial_outlier_snps),
    FALSE,
    TRUE,
    "HYPOTHESIS_GENERATING_NOT_CONFIRMATORY",
    "NO"
  ),
  stringsAsFactors = FALSE
)

summary_file <- file.path(out_dir, "phase13_6C_exclusion_sensitivity_summary.tsv")
write_tsv(summary, summary_file)

# Manuscript inserts.
methods_insert <- c(
  "# SBP instrument-level exclusion sensitivity",
  "",
  "To assess whether the SBP component-divergence pattern was driven by a small number of influential or potentially pleiotropic instruments, we performed instrument-level exclusion sensitivity analyses. Exclusion sets included radial outliers, instruments with high component-level Q contribution, high leverage, high leave-one-SNP contrast influence, and instruments showing association with measured IOP. For each exclusion scenario, IVW random-effects component estimates and the IOP-dependent minus IOP-independent contrast were recomputed.",
  "",
  "These analyses were interpreted as robustness and influence diagnostics only, not as confirmatory evidence."
)

methods_file <- file.path(out_dir, "phase13_6C_methods_insert_instrument_exclusion_sensitivity.md")
writeLines(methods_insert, methods_file, useBytes = TRUE)

all_preserved <- all(non_reference$direction_preserved)

results_insert <- c(
  "# SBP instrument-level exclusion sensitivity",
  "",
  paste0(
    "The full 391-instrument SBP contrast was ",
    fmt(full_contrast$beta_difference_IOP_minus_nonIOP, 5),
    ". Across ",
    nrow(non_reference),
    " exclusion scenarios, the contrast ranged from ",
    fmt(min(non_reference$beta_difference_IOP_minus_nonIOP, na.rm = TRUE), 5),
    " to ",
    fmt(max(non_reference$beta_difference_IOP_minus_nonIOP, na.rm = TRUE), 5),
    "."
  ),
  "",
  paste0(
    "The positive IOP-minus-nonIOP contrast direction was ",
    ifelse(all_preserved, "preserved in all exclusion scenarios.", "not preserved in at least one exclusion scenario and should be treated as a robustness caveat.")
  ),
  "",
  "These analyses support an instrument-level robustness characterization of the SBP component-divergence pattern but do not upgrade the inference beyond a hypothesis-generating interpretation."
)

results_file <- file.path(out_dir, "phase13_6C_results_insert_instrument_exclusion_sensitivity.md")
writeLines(results_insert, results_file, useBytes = TRUE)

# Claim-safety audit.
texts_to_scan <- paste(
  c(
    methods_insert,
    results_insert,
    scenario_results$manuscript_safe_interpretation,
    scenario_results$claim_boundary
  ),
  collapse = "\n"
)

danger_patterns <- data.frame(
  pattern = c(
    "confirmed",
    "validated the mechanism",
    "SBP acts through IOP",
    "mediation proof",
    "pathway proof",
    "SBP protects against NTG",
    "SBP increases IOP-dependent glaucoma risk",
    "HTG validation",
    "hypertension has no effect",
    "definitive proof",
    "establishes causality",
    "causal effect of SBP",
    "proof of mediation",
    "mechanism confirmation"
  ),
  risk = c(
    "confirmatory_overclaim",
    "mechanism_overclaim",
    "pathway_overclaim",
    "mediation_overclaim",
    "pathway_overclaim",
    "protective_overclaim",
    "directional_causal_overclaim",
    "HTG_overclaim",
    "hypertension_null_overclaim",
    "proof_overclaim",
    "causal_overclaim",
    "causal_overclaim",
    "mediation_overclaim",
    "mechanism_overclaim"
  ),
  stringsAsFactors = FALSE
)

count_fixed_hits <- function(pattern, text) {
  hit <- gregexpr(pattern, text, ignore.case = TRUE, fixed = TRUE)[[1]]
  if (length(hit) == 1 && hit[1] == -1) return(0L)
  length(hit)
}

danger_patterns$hits <- vapply(
  danger_patterns$pattern,
  count_fixed_hits,
  integer(1),
  text = texts_to_scan
)

audit_file <- file.path(out_dir, "phase13_6C_claim_safety_audit.tsv")
write_tsv(danger_patterns, audit_file)

high_risk_hits <- sum(danger_patterns$hits)

status <- data.frame(
  field = c(
    "phase",
    "paired_input_found",
    "scenario_results_created",
    "removed_instruments_created",
    "instrument_diagnostics_created",
    "summary_created",
    "methods_insert_created",
    "results_insert_created",
    "claim_safety_audit_created",
    "high_risk_phrase_hits",
    "n_reference_instruments",
    "n_exclusion_scenarios",
    "all_exclusion_scenarios_direction_preserved",
    "minimum_exclusion_contrast",
    "maximum_exclusion_contrast",
    "maximum_absolute_delta_from_full",
    "new_primary_MR_estimates_created",
    "exploratory_sensitivity_created",
    "claim_level",
    "claim_upgrade_allowed",
    "phase13_6C_passed"
  ),
  value = c(
    "Phase13.6C",
    file.exists(paired_file),
    file.exists(scenario_file),
    file.exists(removed_file),
    file.exists(diag_file),
    file.exists(summary_file),
    file.exists(methods_file),
    file.exists(results_file),
    file.exists(audit_file),
    high_risk_hits,
    nrow(paired),
    nrow(non_reference),
    all(non_reference$direction_preserved),
    min(non_reference$beta_difference_IOP_minus_nonIOP, na.rm = TRUE),
    max(non_reference$beta_difference_IOP_minus_nonIOP, na.rm = TRUE),
    max(abs(non_reference$delta_difference_from_full), na.rm = TRUE),
    FALSE,
    TRUE,
    "HYPOTHESIS_GENERATING_NOT_CONFIRMATORY",
    "NO",
    high_risk_hits == 0
  ),
  stringsAsFactors = FALSE
)

status_file <- file.path(out_dir, "phase13_6C_status.tsv")
write_tsv(status, status_file)

if (file.exists(master_status_file)) {
  master_status <- read.delim(master_status_file, check.names = FALSE)
  idx <- which(master_status$phase == "Phase13.6")
  if (length(idx) == 1) {
    master_status$status[idx] <- ifelse(
      high_risk_hits == 0,
      "PASSED_ROBUSTNESS_MATRIX_AND_INSTRUMENT_EXCLUSION_SENSITIVITY",
      "FAILED_PHASE13_6C_REVIEW_REQUIRED"
    )
    master_status$qc_status[idx] <- ifelse(high_risk_hits == 0, "PASSED", "REVIEW_REQUIRED")
    master_status$primary_output[idx] <- "phase13_6B_SBP_robustness_matrix.tsv; phase13_6C_sbp_instrument_exclusion_sensitivity.tsv"
  }
  write_tsv(master_status, master_status_file)
}

message("Phase13.6C completed.")
message("Scenario results: ", scenario_file)
message("Removed instruments: ", removed_file)
message("Diagnostics: ", diag_file)
message("Summary: ", summary_file)
message("Status: ", status_file)
message("High-risk phrase hits: ", high_risk_hits)
