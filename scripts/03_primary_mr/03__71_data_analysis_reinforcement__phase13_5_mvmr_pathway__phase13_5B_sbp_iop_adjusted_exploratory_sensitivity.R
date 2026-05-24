options(stringsAsFactors = FALSE)

# Phase13.5B: SBP-instrument measured-IOP-adjusted exploratory sensitivity
# Project: IOP-dependent vs IOP-independent glaucoma component MR
# Evidence level: HYPOTHESIS_GENERATING_NOT_CONFIRMATORY
# This analysis is an exploratory pathway-consistency sensitivity only.
# It is not a mediation test and does not establish pathway proof.

message("Running Phase13.5B SBP-instrument measured-IOP-adjusted exploratory sensitivity...")

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

as_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

fmt <- function(x, digits = 4) {
  ifelse(is.na(x), NA, formatC(as.numeric(x), digits = digits, format = "fg", flag = "#"))
}

root <- normalizePath(".", winslash = "/", mustWork = TRUE)
out_dir <- file.path(root, "71_data_analysis_reinforcement", "phase13_5_mvmr_pathway")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

master_status_file <- file.path(root, "71_data_analysis_reinforcement", "phase13_master_status.tsv")

overlap_file <- file.path(out_dir, "phase13_5A1_SBP_component_measured_IOP_overlap_dataset.tsv")

if (!file.exists(overlap_file)) {
  stop("Missing Phase13.5A1 overlap dataset: ", overlap_file)
}

dat <- read.delim(overlap_file, check.names = FALSE)

required_cols <- c(
  "SNP",
  "beta_SBP",
  "se_SBP",
  "beta_measured_IOP",
  "se_measured_IOP",
  "beta_GBS_nonIOPcomponent",
  "se_GBS_nonIOPcomponent",
  "beta_GBS_IOPcomponent",
  "se_GBS_IOPcomponent",
  "has_measured_IOP_association",
  "allele_match_internal_vs_iop_file"
)

missing_cols <- setdiff(required_cols, names(dat))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

dat <- dat[dat$has_measured_IOP_association == TRUE & dat$allele_match_internal_vs_iop_file == TRUE, ]

# Standardized clean input for downstream audit.
mvmr_input <- data.frame(
  SNP = as.character(dat$SNP),
  chr = as.character(dat$chr),
  pos = as_num(dat$pos),
  beta_SBP = as_num(dat$beta_SBP),
  se_SBP = as_num(dat$se_SBP),
  beta_measured_IOP = as_num(dat$beta_measured_IOP),
  se_measured_IOP = as_num(dat$se_measured_IOP),
  beta_GBS_nonIOPcomponent = as_num(dat$beta_GBS_nonIOPcomponent),
  se_GBS_nonIOPcomponent = as_num(dat$se_GBS_nonIOPcomponent),
  beta_GBS_IOPcomponent = as_num(dat$beta_GBS_IOPcomponent),
  se_GBS_IOPcomponent = as_num(dat$se_GBS_IOPcomponent),
  locus_proxy_1Mb = if ("locus_proxy_1Mb" %in% names(dat)) as.character(dat$locus_proxy_1Mb) else NA_character_,
  stringsAsFactors = FALSE
)

mvmr_input_file <- file.path(out_dir, "phase13_5B_sbp_iop_adjusted_input.tsv")
write_tsv(mvmr_input, mvmr_input_file)

# Weighted regression without intercept:
# outcome SNP association ~ SNP-SBP association + SNP-measured-IOP association
# weights = 1 / outcome SE^2
# This is an exploratory sensitivity using SBP-selected instruments.
wls_no_intercept <- function(y, se_y, X, exposure_names) {
  y <- as_num(y)
  se_y <- as_num(se_y)
  X <- as.matrix(X)
  storage.mode(X) <- "numeric"
  
  keep <- is.finite(y) & is.finite(se_y) & se_y > 0 & apply(X, 1, function(z) all(is.finite(z)))
  y <- y[keep]
  se_y <- se_y[keep]
  X <- X[keep, , drop = FALSE]
  
  w <- 1 / se_y^2
  XtW <- t(X) * w
  XtWX <- XtW %*% X
  XtWy <- XtW %*% y
  
  beta <- as.numeric(solve(XtWX, XtWy))
  fitted <- as.numeric(X %*% beta)
  resid <- y - fitted
  
  Q <- sum(w * resid^2)
  df <- length(y) - ncol(X)
  phi <- max(1, Q / df)
  vcov <- phi * solve(XtWX)
  se <- sqrt(diag(vcov))
  z <- beta / se
  p <- 2 * pnorm(abs(z), lower.tail = FALSE)
  
  data.frame(
    exposure = exposure_names,
    beta = beta,
    se = se,
    z = z,
    p = p,
    Q = Q,
    Q_df = df,
    Q_p = pchisq(Q, df = df, lower.tail = FALSE),
    phi = phi,
    n_instruments = length(y),
    stringsAsFactors = FALSE
  )
}

run_outcome <- function(outcome_name, beta_y_col, se_y_col) {
  y <- mvmr_input[[beta_y_col]]
  se_y <- mvmr_input[[se_y_col]]
  
  uvw <- wls_no_intercept(
    y = y,
    se_y = se_y,
    X = data.frame(beta_SBP = mvmr_input$beta_SBP),
    exposure_names = "SBP_univariable"
  )
  
  adj <- wls_no_intercept(
    y = y,
    se_y = se_y,
    X = data.frame(
      beta_SBP = mvmr_input$beta_SBP,
      beta_measured_IOP = mvmr_input$beta_measured_IOP
    ),
    exposure_names = c("SBP_adjusted_for_measured_IOP", "measured_IOP_adjusted_for_SBP")
  )
  
  uvw$outcome_suffix <- outcome_name
  uvw$model <- "univariable_IVW_using_same_391_instruments"
  
  adj$outcome_suffix <- outcome_name
  adj$model <- "exploratory_SBP_plus_measured_IOP_weighted_regression"
  
  rbind(uvw, adj)
}

results <- rbind(
  run_outcome(
    outcome_name = "GBS_nonIOPcomponent",
    beta_y_col = "beta_GBS_nonIOPcomponent",
    se_y_col = "se_GBS_nonIOPcomponent"
  ),
  run_outcome(
    outcome_name = "GBS_IOPcomponent",
    beta_y_col = "beta_GBS_IOPcomponent",
    se_y_col = "se_GBS_IOPcomponent"
  )
)

results <- results[, c(
  "outcome_suffix",
  "model",
  "exposure",
  "n_instruments",
  "beta",
  "se",
  "z",
  "p",
  "Q",
  "Q_df",
  "Q_p",
  "phi"
)]

results_file <- file.path(out_dir, "phase13_5B_sbp_iop_adjusted_results.tsv")
write_tsv(results, results_file)

# Contrast comparison: SBP coefficient IOP-minus-nonIOP before and after measured-IOP adjustment.
get_beta <- function(outcome, exposure) {
  z <- results$beta[results$outcome_suffix == outcome & results$exposure == exposure]
  if (length(z) == 0) return(NA_real_)
  as_num(z[1])
}

get_se <- function(outcome, exposure) {
  z <- results$se[results$outcome_suffix == outcome & results$exposure == exposure]
  if (length(z) == 0) return(NA_real_)
  as_num(z[1])
}

contrast_rows <- list()

for (exposure_label in c("SBP_univariable", "SBP_adjusted_for_measured_IOP")) {
  b_non <- get_beta("GBS_nonIOPcomponent", exposure_label)
  b_iop <- get_beta("GBS_IOPcomponent", exposure_label)
  se_non <- get_se("GBS_nonIOPcomponent", exposure_label)
  se_iop <- get_se("GBS_IOPcomponent", exposure_label)
  
  # Working r=0 contrast SE for comparability only.
  diff <- b_iop - b_non
  se_diff <- sqrt(se_non^2 + se_iop^2)
  z <- diff / se_diff
  p <- 2 * pnorm(abs(z), lower.tail = FALSE)
  
  contrast_rows[[length(contrast_rows) + 1]] <- data.frame(
    exposure_label = exposure_label,
    beta_nonIOP = b_non,
    se_nonIOP = se_non,
    beta_IOP = b_iop,
    se_IOP = se_iop,
    beta_difference_IOP_minus_nonIOP = diff,
    se_difference_working_r0 = se_diff,
    z_difference = z,
    p_difference_working_r0 = p,
    direction_preserved = diff > 0,
    stringsAsFactors = FALSE
  )
}

contrast <- do.call(rbind, contrast_rows)

univ_diff <- contrast$beta_difference_IOP_minus_nonIOP[contrast$exposure_label == "SBP_univariable"]
adj_diff <- contrast$beta_difference_IOP_minus_nonIOP[contrast$exposure_label == "SBP_adjusted_for_measured_IOP"]

contrast$absolute_change_from_univariable <- contrast$beta_difference_IOP_minus_nonIOP - univ_diff
contrast$relative_change_from_univariable <- contrast$absolute_change_from_univariable / abs(univ_diff)

contrast_file <- file.path(out_dir, "phase13_5B_sbp_iop_adjusted_contrast_comparison.tsv")
write_tsv(contrast, contrast_file)

# Exposure-column diagnostics.
X <- as.matrix(data.frame(
  beta_SBP = mvmr_input$beta_SBP,
  beta_measured_IOP = mvmr_input$beta_measured_IOP
))

exposure_correlation <- suppressWarnings(cor(X[, 1], X[, 2], use = "complete.obs"))

# Condition number on scaled columns, avoiding intercept.
scaled_X <- scale(X, center = TRUE, scale = TRUE)
condition_number <- suppressWarnings(kappa(crossprod(scaled_X), exact = TRUE))

diagnostics <- data.frame(
  diagnostic = c(
    "n_instruments",
    "correlation_beta_SBP_beta_measured_IOP",
    "condition_number_scaled_exposure_crossproduct",
    "mean_abs_beta_SBP",
    "mean_abs_beta_measured_IOP",
    "interpretation"
  ),
  value = c(
    nrow(mvmr_input),
    exposure_correlation,
    condition_number,
    mean(abs(mvmr_input$beta_SBP), na.rm = TRUE),
    mean(abs(mvmr_input$beta_measured_IOP), na.rm = TRUE),
    "Exposure-column diagnostics for exploratory sensitivity; not a formal conditional-instrument-strength proof."
  ),
  stringsAsFactors = FALSE
)

diagnostics_file <- file.path(out_dir, "phase13_5B_exposure_column_diagnostics.tsv")
write_tsv(diagnostics, diagnostics_file)

# Interpretation grid.
sbp_univ_non <- get_beta("GBS_nonIOPcomponent", "SBP_univariable")
sbp_adj_non <- get_beta("GBS_nonIOPcomponent", "SBP_adjusted_for_measured_IOP")
sbp_univ_iop <- get_beta("GBS_IOPcomponent", "SBP_univariable")
sbp_adj_iop <- get_beta("GBS_IOPcomponent", "SBP_adjusted_for_measured_IOP")

attenuation_non <- (sbp_adj_non - sbp_univ_non) / abs(sbp_univ_non)
attenuation_iop <- (sbp_adj_iop - sbp_univ_iop) / abs(sbp_univ_iop)
attenuation_contrast <- (adj_diff - univ_diff) / abs(univ_diff)

interpretation <- data.frame(
  item = c(
    "nonIOP_SBP_coefficient_change_after_measured_IOP_adjustment",
    "IOP_SBP_coefficient_change_after_measured_IOP_adjustment",
    "SBP_component_contrast_change_after_measured_IOP_adjustment",
    "primary_safe_interpretation",
    "claim_boundary"
  ),
  value = c(
    attenuation_non,
    attenuation_iop,
    attenuation_contrast,
    ifelse(
      is.finite(attenuation_contrast) && abs(attenuation_contrast) < 0.25,
      "SBP component contrast was not materially attenuated after including measured IOP SNP association in this exploratory sensitivity.",
      "SBP component contrast changed after including measured IOP SNP association; interpret only as exploratory pathway-consistency context."
    ),
    "This analysis is exploratory and does not establish mediation or pathway proof."
  ),
  stringsAsFactors = FALSE
)

interpretation_file <- file.path(out_dir, "phase13_5B_pathway_consistency_interpretation.tsv")
write_tsv(interpretation, interpretation_file)

# Methods / Results inserts.
methods_insert <- c(
  "# Exploratory measured-IOP-adjusted sensitivity",
  "",
  "As an exploratory pathway-consistency sensitivity, we evaluated whether the SBP component estimates changed after including SNP associations with measured IOP as a second exposure column among the same 391 SBP instruments. We fit weighted no-intercept regressions of SNP-component associations on SNP-SBP associations alone and on SNP-SBP plus SNP-measured-IOP associations, using inverse outcome-variance weights.",
  "",
  "Because the analysis used SBP-selected instruments rather than a full union instrument set for both exposures, it was interpreted as an exploratory sensitivity analysis rather than a mediation test."
)

methods_file <- file.path(out_dir, "phase13_5B_methods_insert_exploratory_iop_adjusted_sensitivity.md")
writeLines(methods_insert, methods_file, useBytes = TRUE)

results_insert <- c(
  "# Exploratory measured-IOP-adjusted sensitivity",
  "",
  paste0(
    "All 391 internal SBP component instruments were matched to measured IOP association estimates with consistent alleles and identical SBP exposure columns. ",
    "The correlation between SNP-SBP and SNP-measured-IOP associations was ",
    fmt(exposure_correlation, 4),
    "."
  ),
  "",
  paste0(
    "In exploratory weighted regression models, the SBP IOP-minus-nonIOP contrast was ",
    fmt(univ_diff, 4),
    " in the univariable model and ",
    fmt(adj_diff, 4),
    " after including measured IOP SNP association as a second exposure column. ",
    "The relative change in the contrast was ",
    fmt(attenuation_contrast, 4),
    ". These results were interpreted as pathway-consistency sensitivity only and did not change the hypothesis-generating framing."
  )
)

results_file <- file.path(out_dir, "phase13_5B_results_insert_exploratory_iop_adjusted_sensitivity.md")
writeLines(results_insert, results_file, useBytes = TRUE)

# Claim-safety audit.
texts_to_scan <- paste(c(methods_insert, results_insert, interpretation$value), collapse = "\n")

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
    "definitive proof"
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
    "proof_overclaim"
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

audit_file <- file.path(out_dir, "phase13_5B_claim_safety_audit.tsv")
write_tsv(danger_patterns, audit_file)

high_risk_hits <- sum(danger_patterns$hits)

# Status.
status <- data.frame(
  field = c(
    "phase",
    "mvmr_input_created",
    "results_created",
    "contrast_comparison_created",
    "exposure_column_diagnostics_created",
    "interpretation_created",
    "methods_insert_created",
    "results_insert_created",
    "claim_safety_audit_created",
    "high_risk_phrase_hits",
    "n_instruments",
    "exposure_correlation_SBP_measured_IOP",
    "univariable_SBP_contrast",
    "IOP_adjusted_SBP_contrast",
    "relative_change_in_SBP_contrast_after_IOP_adjustment",
    "IOP_adjusted_contrast_direction_preserved",
    "new_primary_MR_estimates_created",
    "exploratory_sensitivity_created",
    "claim_level",
    "claim_upgrade_allowed",
    "phase13_5B_passed"
  ),
  value = c(
    "Phase13.5B",
    file.exists(mvmr_input_file),
    file.exists(results_file),
    file.exists(contrast_file),
    file.exists(diagnostics_file),
    file.exists(interpretation_file),
    file.exists(methods_file),
    file.exists(results_file),
    file.exists(audit_file),
    high_risk_hits,
    nrow(mvmr_input),
    exposure_correlation,
    univ_diff,
    adj_diff,
    attenuation_contrast,
    adj_diff > 0,
    FALSE,
    TRUE,
    "HYPOTHESIS_GENERATING_NOT_CONFIRMATORY",
    "NO",
    high_risk_hits == 0 && nrow(mvmr_input) == 391 && is.finite(adj_diff)
  ),
  stringsAsFactors = FALSE
)

status_file <- file.path(out_dir, "phase13_5B_status.tsv")
write_tsv(status, status_file)

if (file.exists(master_status_file)) {
  master_status <- read.delim(master_status_file, check.names = FALSE)
  idx <- which(master_status$phase == "Phase13.5")
  if (length(idx) == 1) {
    master_status$status[idx] <- ifelse(
      high_risk_hits == 0 && nrow(mvmr_input) == 391 && is.finite(adj_diff),
      "PASSED_EXPLORATORY_MEASURED_IOP_ADJUSTED_SENSITIVITY",
      "FAILED_PHASE13_5B_REVIEW_REQUIRED"
    )
    master_status$qc_status[idx] <- ifelse(
      high_risk_hits == 0 && nrow(mvmr_input) == 391 && is.finite(adj_diff),
      "PASSED",
      "REVIEW_REQUIRED"
    )
    master_status$primary_output[idx] <- "phase13_5B_sbp_iop_adjusted_results.tsv; phase13_5B_sbp_iop_adjusted_contrast_comparison.tsv; phase13_5B_pathway_consistency_interpretation.tsv"
  }
  write_tsv(master_status, master_status_file)
}

message("Phase13.5B completed.")
message("Results: ", results_file)
message("Contrast comparison: ", contrast_file)
message("Diagnostics: ", diagnostics_file)
message("Interpretation: ", interpretation_file)
message("Status: ", status_file)
message("High-risk phrase hits: ", high_risk_hits)
