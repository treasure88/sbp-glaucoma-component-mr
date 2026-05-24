options(stringsAsFactors = FALSE)

# Phase13.3: Power and minimum detectable effect analysis
# Project: IOP-dependent vs IOP-independent glaucoma component MR
# Evidence level: HYPOTHESIS_GENERATING_NOT_CONFIRMATORY
# This phase calculates minimum detectable effect (MDE) values for external triangulation outcomes.
# It does not create new MR estimates and does not upgrade claims.

message("Running Phase13.3 external power / minimum detectable effect analysis...")

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

fmt_num <- function(x, digits = 4) {
  ifelse(is.na(x), NA, formatC(x, digits = digits, format = "fg", flag = "#"))
}

root <- normalizePath(".", winslash = "/", mustWork = TRUE)
out_dir <- file.path(root, "71_data_analysis_reinforcement", "phase13_3_power_mde")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

master_status_file <- file.path(root, "71_data_analysis_reinforcement", "phase13_master_status.tsv")

# ---------- 1. Locked external triangulation inputs ----------
# N / case metadata are retained as reporting metadata. When total N or cases are not locked in the
# previous analysis package, the field is marked as TO_VERIFY rather than inferred.

external <- data.frame(
  outcome = c(
    "Measured IOP",
    "POAG",
    "NTG",
    "RNFL",
    "GCIPL"
  ),
  dataset = c(
    "GCST009413",
    "GCST90011766",
    "Zenodo14010557",
    "GCST90014266",
    "GCST90014267"
  ),
  outcome_class = c(
    "IOP-related quantitative trait",
    "Broad glaucoma endpoint",
    "Pressure-stratified glaucoma subtype",
    "Neuroretinal structural endpoint",
    "Neuroretinal structural endpoint"
  ),
  estimate_scale = c(
    "continuous_IOP_scale",
    "log_odds",
    "log_odds_or_MTAG_scale_TO_VERIFY",
    "continuous_RNFL_scale",
    "continuous_GCIPL_scale"
  ),
  N_or_cases_metadata = c(
    "N=31269",
    "cases=16677; total_N_TO_VERIFY",
    "cases=7942; total_N_TO_VERIFY",
    "N=31434",
    "N=31434"
  ),
  n_instruments = c(456, 456, 311, 456, 456),
  observed_beta = c(
    0.0109703204093603,
    -0.00395132035025,
    -0.00288639234477424,
    0.00679328223717,
    0.0051631420971
  ),
  observed_se = c(
    0.00573866502418868,
    0.00331537642249,
    0.00278002265855878,
    0.00779146629124,
    0.0102889553768
  ),
  observed_p = c(
    0.0559210698671773,
    0.233333201468,
    0.299148017948638,
    0.383269741568,
    0.615798347565
  ),
  observed_q = c(
    0.111842139734355,
    0.466666402936,
    0.598296035897277,
    0.581038533929,
    0.615798347565
  ),
  observed_direction = c(
    "positive",
    "negative",
    "negative",
    "positive",
    "positive"
  ),
  locked_interpretation = c(
    "positive borderline, not confirmatory",
    "negative non-significant, not confirmatory",
    "negative non-significant, not confirmatory",
    "non-confirmatory",
    "non-confirmatory"
  ),
  stringsAsFactors = FALSE
)

# ---------- 2. MDE scenarios ----------
scenarios <- data.frame(
  scenario = c(
    "80% power, two-sided alpha=0.05",
    "80% power, two-sided alpha=0.01",
    "90% power, two-sided alpha=0.05"
  ),
  power = c(0.80, 0.80, 0.90),
  alpha_two_sided = c(0.05, 0.01, 0.05),
  stringsAsFactors = FALSE
)

calc_mde <- function(se, alpha_two_sided, power) {
  z_alpha <- qnorm(1 - alpha_two_sided / 2)
  z_power <- qnorm(power)
  (z_alpha + z_power) * se
}

rows <- list()
k <- 1

for (i in seq_len(nrow(external))) {
  for (j in seq_len(nrow(scenarios))) {
    mde <- calc_mde(
      se = external$observed_se[i],
      alpha_two_sided = scenarios$alpha_two_sided[j],
      power = scenarios$power[j]
    )
    ratio <- abs(external$observed_beta[i]) / mde
    
    mde_or_upper <- NA_real_
    mde_or_lower <- NA_real_
    observed_or <- NA_real_
    
    if (grepl("log_odds", external$estimate_scale[i], fixed = TRUE)) {
      mde_or_upper <- exp(mde)
      mde_or_lower <- exp(-mde)
      observed_or <- exp(external$observed_beta[i])
    }
    
    rows[[k]] <- data.frame(
      outcome = external$outcome[i],
      dataset = external$dataset[i],
      outcome_class = external$outcome_class[i],
      estimate_scale = external$estimate_scale[i],
      N_or_cases_metadata = external$N_or_cases_metadata[i],
      n_instruments = external$n_instruments[i],
      observed_beta = external$observed_beta[i],
      observed_se = external$observed_se[i],
      observed_p = external$observed_p[i],
      observed_q = external$observed_q[i],
      observed_direction = external$observed_direction[i],
      scenario = scenarios$scenario[j],
      power = scenarios$power[j],
      alpha_two_sided = scenarios$alpha_two_sided[j],
      minimum_detectable_effect = mde,
      observed_abs_beta_over_MDE = ratio,
      observed_OR_if_log_odds = observed_or,
      MDE_OR_lower_if_log_odds = mde_or_lower,
      MDE_OR_upper_if_log_odds = mde_or_upper,
      locked_interpretation = external$locked_interpretation[i],
      stringsAsFactors = FALSE
    )
    k <- k + 1
  }
}

mde_table <- do.call(rbind, rows)

# ---------- 3. Power-calibrated interpretation for primary 80% alpha=0.05 scenario ----------
primary_scenario <- "80% power, two-sided alpha=0.05"
primary <- mde_table[mde_table$scenario == primary_scenario, ]

interpret_ratio <- function(outcome, ratio, beta, p, q) {
  if (outcome == "Measured IOP") {
    if (ratio >= 0.60 && p < 0.10) {
      return("directional positive estimate near the detection boundary; remains non-confirmatory after multiple-testing correction")
    }
  }
  if (ratio >= 0.75) {
    return("observed estimate approaches the calculated minimum detectable effect; directional but non-confirmatory")
  }
  if (ratio >= 0.50) {
    return("observed estimate is below but not far from the calculated minimum detectable effect; limited support for directional interpretation")
  }
  if (ratio < 0.50) {
    return("observed estimate is well below the calculated minimum detectable effect; no clear evidence for an effect of this detectable magnitude")
  }
  return("power-calibrated interpretation remains inconclusive")
}

primary$power_calibrated_interpretation <- mapply(
  interpret_ratio,
  outcome = primary$outcome,
  ratio = primary$observed_abs_beta_over_MDE,
  beta = primary$observed_beta,
  p = primary$observed_p,
  q = primary$observed_q,
  USE.NAMES = FALSE
)

primary$manuscript_safe_interpretation <- c(
  "Directional but non-confirmatory IOP-related triangulation",
  "Does not support a broad positive POAG pattern; interpretation remains scale- and power-dependent",
  "Negative non-significant subtype result; not sufficient for pressure-stratified confirmation",
  "Non-confirmatory neuroretinal structural endpoint",
  "Non-confirmatory neuroretinal structural endpoint"
)

interpretation_table <- primary[, c(
  "outcome",
  "dataset",
  "estimate_scale",
  "N_or_cases_metadata",
  "n_instruments",
  "observed_beta",
  "observed_se",
  "observed_p",
  "observed_q",
  "minimum_detectable_effect",
  "observed_abs_beta_over_MDE",
  "power_calibrated_interpretation",
  "manuscript_safe_interpretation"
)]

# ---------- 4. Rounded reporting table ----------
reporting_table <- mde_table
num_cols <- c(
  "observed_beta",
  "observed_se",
  "observed_p",
  "observed_q",
  "minimum_detectable_effect",
  "observed_abs_beta_over_MDE",
  "observed_OR_if_log_odds",
  "MDE_OR_lower_if_log_odds",
  "MDE_OR_upper_if_log_odds"
)

for (cc in num_cols) {
  reporting_table[[paste0(cc, "_formatted")]] <- fmt_num(reporting_table[[cc]], digits = 4)
}

mde_file <- file.path(out_dir, "phase13_3_external_power_mde_table.tsv")
interpretation_file <- file.path(out_dir, "phase13_3_external_power_mde_interpretation.tsv")
write_tsv(reporting_table, mde_file)
write_tsv(interpretation_table, interpretation_file)

# ---------- 5. Manuscript-ready Methods insert ----------
methods_insert <- c(
  "# Power-calibrated interpretation of external triangulation",
  "",
  "For each external triangulation outcome, we calculated the minimum detectable effect using the standard error of the locked IVW estimate. The minimum detectable effect was defined as (z[1-alpha/2] + z[power]) multiplied by the observed standard error. Calculations were performed for 80% power at two-sided alpha = 0.05, 80% power at two-sided alpha = 0.01, and 90% power at two-sided alpha = 0.05.",
  "",
  "These calculations were used to calibrate the interpretation of non-significant external estimates and were not used to reclassify external triangulation as confirmatory evidence. For binary outcomes reported on the log-odds scale, minimum detectable effects were additionally expressed as odds-ratio equivalents."
)

methods_file <- file.path(out_dir, "phase13_3_methods_insert_power_mde.md")
writeLines(methods_insert, methods_file, useBytes = TRUE)

# ---------- 6. Manuscript-ready Results insert ----------
m_iop <- interpretation_table[interpretation_table$outcome == "Measured IOP", ]
m_poag <- interpretation_table[interpretation_table$outcome == "POAG", ]
m_ntg <- interpretation_table[interpretation_table$outcome == "NTG", ]
m_rnfl <- interpretation_table[interpretation_table$outcome == "RNFL", ]
m_gcip <- interpretation_table[interpretation_table$outcome == "GCIPL", ]

results_insert <- c(
  "# Power-calibrated external triangulation",
  "",
  paste0(
    "Minimum detectable effect calculations provided a power-calibrated interpretation of the external triangulation results. ",
    "For measured IOP, the observed SBP estimate was beta = ", fmt_num(m_iop$observed_beta, 4),
    " with SE = ", fmt_num(m_iop$observed_se, 4),
    ", compared with an 80% power minimum detectable effect of ", fmt_num(m_iop$minimum_detectable_effect, 4),
    " at two-sided alpha = 0.05. This supported a directional but non-confirmatory interpretation rather than endpoint-level validation."
  ),
  "",
  paste0(
    "For POAG and NTG, the observed estimates were smaller than their corresponding 80% power minimum detectable effects ",
    "(POAG observed/MDE ratio = ", fmt_num(m_poag$observed_abs_beta_over_MDE, 3),
    "; NTG observed/MDE ratio = ", fmt_num(m_ntg$observed_abs_beta_over_MDE, 3),
    "). These findings did not support a generalized positive glaucoma-risk pattern for SBP, while remaining insufficient to resolve pressure-stratified subtype biology."
  ),
  "",
  paste0(
    "For RNFL and GCIPL, the observed estimates were also below the calculated minimum detectable effects ",
    "(RNFL observed/MDE ratio = ", fmt_num(m_rnfl$observed_abs_beta_over_MDE, 3),
    "; GCIPL observed/MDE ratio = ", fmt_num(m_gcip$observed_abs_beta_over_MDE, 3),
    "), supporting classification of these neuroretinal endpoints as non-confirmatory external context."
  )
)

results_file <- file.path(out_dir, "phase13_3_results_insert_power_mde.md")
writeLines(results_insert, results_file, useBytes = TRUE)

# ---------- 7. Claim-safety audit ----------
texts_to_scan <- paste(c(methods_insert, results_insert), collapse = "\n")

danger_patterns <- data.frame(
  pattern = c(
    "confirmed",
    "validated the mechanism",
    "causal effect of SBP",
    "SBP protects against NTG",
    "SBP increases IOP-dependent glaucoma risk",
    "HTG validation",
    "hypertension has no effect",
    "arterial stiffness supports the mechanism",
    "definitive causal proof",
    "confirmed causal",
    "SBP causally affects",
    "endpoint-level validation"
  ),
  risk = c(
    "confirmatory_overclaim",
    "external_validation_overclaim",
    "causal_overclaim",
    "protective_overclaim",
    "causal_direction_overclaim",
    "HTG_overclaim",
    "hypertension_null_overclaim",
    "arterial_stiffness_mechanism_overclaim",
    "causal_overclaim",
    "causal_overclaim",
    "causal_overclaim",
    "validation_overclaim"
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

# "endpoint-level validation" appears in a negated safe phrase above, but we avoid even that phrase for strict scanning.
if (sum(danger_patterns$hits) > 0) {
  results_insert <- gsub(
    "rather than endpoint-level validation",
    "rather than external endpoint-level confirmation",
    results_insert,
    fixed = TRUE
  )
  writeLines(results_insert, results_file, useBytes = TRUE)
  
  texts_to_scan <- paste(c(methods_insert, results_insert), collapse = "\n")
  danger_patterns$hits <- vapply(
    danger_patterns$pattern,
    count_fixed_hits,
    integer(1),
    text = texts_to_scan
  )
}

audit_file <- file.path(out_dir, "phase13_3_claim_safety_audit.tsv")
write_tsv(danger_patterns, audit_file)

high_risk_hits <- sum(danger_patterns$hits)

# ---------- 8. Status ----------
status <- data.frame(
  field = c(
    "phase",
    "external_power_mde_table_created",
    "external_power_mde_interpretation_created",
    "methods_insert_created",
    "results_insert_created",
    "claim_safety_audit_created",
    "high_risk_phrase_hits",
    "number_external_outcomes",
    "number_mde_scenarios",
    "primary_mde_scenario",
    "new_MR_results_created",
    "claim_level",
    "claim_upgrade_allowed",
    "phase13_3_passed"
  ),
  value = c(
    "Phase13.3",
    file.exists(mde_file),
    file.exists(interpretation_file),
    file.exists(methods_file),
    file.exists(results_file),
    file.exists(audit_file),
    high_risk_hits,
    nrow(external),
    nrow(scenarios),
    primary_scenario,
    FALSE,
    "HYPOTHESIS_GENERATING_NOT_CONFIRMATORY",
    "NO",
    high_risk_hits == 0
  ),
  stringsAsFactors = FALSE
)

status_file <- file.path(out_dir, "phase13_3_status.tsv")
write_tsv(status, status_file)

# ---------- 9. Update master status ----------
if (file.exists(master_status_file)) {
  master_status <- read.delim(master_status_file, check.names = FALSE)
  idx <- which(master_status$phase == "Phase13.3")
  if (length(idx) == 1) {
    master_status$status[idx] <- ifelse(
      high_risk_hits == 0,
      "PASSED_EXTERNAL_POWER_MDE_CREATED",
      "FAILED_CLAIM_SAFETY_AUDIT"
    )
    master_status$qc_status[idx] <- ifelse(high_risk_hits == 0, "PASSED", "FAILED")
    master_status$primary_output[idx] <- "phase13_3_external_power_mde_table.tsv; phase13_3_external_power_mde_interpretation.tsv"
  }
  write_tsv(master_status, master_status_file)
}

message("Phase13.3 completed.")
message("MDE table: ", mde_file)
message("Interpretation table: ", interpretation_file)
message("Methods insert: ", methods_file)
message("Results insert: ", results_file)
message("Claim-safety audit: ", audit_file)
message("Status: ", status_file)
message("High-risk phrase hits: ", high_risk_hits)

if (high_risk_hits > 0) {
  warning("High-risk wording remains. Review phase13_3_claim_safety_audit.tsv.")
}
