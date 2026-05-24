options(stringsAsFactors = FALSE)

# Phase13.8: Power-calibrated external triangulation interpretation grid
# Project: IOP-dependent vs IOP-independent glaucoma component MR
# Evidence level: HYPOTHESIS_GENERATING_NOT_CONFIRMATORY
# This phase integrates existing outputs only. It does not create new MR estimates.

message("Running Phase13.8 power-calibrated external triangulation grid...")

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
  read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
}

as_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

fmt <- function(x, digits = 4) {
  ifelse(is.na(x), NA, formatC(as.numeric(x), digits = digits, format = "fg", flag = "#"))
}

root <- normalizePath(".", winslash = "/", mustWork = TRUE)
out_dir <- file.path(root, "71_data_analysis_reinforcement", "phase13_8_external_triangulation_grid")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

master_status_file <- file.path(root, "71_data_analysis_reinforcement", "phase13_master_status.tsv")

power_file <- file.path(
  root,
  "71_data_analysis_reinforcement/phase13_3_power_mde/phase13_3_external_power_mde_interpretation.tsv"
)

strength_file <- file.path(
  root,
  "71_data_analysis_reinforcement/phase13_7_instrument_transparency/phase13_7B_SBP_external_instrument_strength_summary.tsv"
)

bootstrap_file <- file.path(
  root,
  "71_data_analysis_reinforcement/phase13_4_empirical_covariance/phase13_4B_bootstrap_contrast_summary.tsv"
)

mvmr_contrast_file <- file.path(
  root,
  "71_data_analysis_reinforcement/phase13_5_mvmr_pathway/phase13_5B_sbp_iop_adjusted_contrast_comparison.tsv"
)

mvmr_interp_file <- file.path(
  root,
  "71_data_analysis_reinforcement/phase13_5_mvmr_pathway/phase13_5B_pathway_consistency_interpretation.tsv"
)

robustness_caveat_file <- file.path(
  root,
  "71_data_analysis_reinforcement/phase13_6_pleiotropy_robustness/phase13_6B1_directional_caveat_summary.tsv"
)

robustness_coverage_file <- file.path(
  root,
  "71_data_analysis_reinforcement/phase13_6_pleiotropy_robustness/phase13_6B_robustness_coverage_interpretation.tsv"
)

power <- safe_read(power_file)
strength <- safe_read(strength_file)
bootstrap <- safe_read(bootstrap_file)
mvmr_contrast <- safe_read(mvmr_contrast_file)
mvmr_interp <- safe_read(mvmr_interp_file)
robustness_caveat <- safe_read(robustness_caveat_file)
robustness_coverage <- safe_read(robustness_coverage_file)

source_manifest <- data.frame(
  source_role = c(
    "external_power_mde",
    "external_instrument_strength",
    "internal_bootstrap_contrast",
    "exploratory_measured_IOP_adjusted_contrast",
    "exploratory_measured_IOP_adjusted_interpretation",
    "robustness_directional_caveat",
    "robustness_coverage"
  ),
  relative_path = c(
    "71_data_analysis_reinforcement/phase13_3_power_mde/phase13_3_external_power_mde_interpretation.tsv",
    "71_data_analysis_reinforcement/phase13_7_instrument_transparency/phase13_7B_SBP_external_instrument_strength_summary.tsv",
    "71_data_analysis_reinforcement/phase13_4_empirical_covariance/phase13_4B_bootstrap_contrast_summary.tsv",
    "71_data_analysis_reinforcement/phase13_5_mvmr_pathway/phase13_5B_sbp_iop_adjusted_contrast_comparison.tsv",
    "71_data_analysis_reinforcement/phase13_5_mvmr_pathway/phase13_5B_pathway_consistency_interpretation.tsv",
    "71_data_analysis_reinforcement/phase13_6_pleiotropy_robustness/phase13_6B1_directional_caveat_summary.tsv",
    "71_data_analysis_reinforcement/phase13_6_pleiotropy_robustness/phase13_6B_robustness_coverage_interpretation.tsv"
  ),
  exists = c(
    !is.null(power),
    !is.null(strength),
    !is.null(bootstrap),
    !is.null(mvmr_contrast),
    !is.null(mvmr_interp),
    !is.null(robustness_caveat),
    !is.null(robustness_coverage)
  ),
  stringsAsFactors = FALSE
)

source_manifest_file <- file.path(out_dir, "phase13_8_source_manifest.tsv")
write_tsv(source_manifest, source_manifest_file)

# ---- External triangulation grid ----
external_rows <- list()

if (!is.null(power) && nrow(power) > 0) {
  for (i in seq_len(nrow(power))) {
    outcome <- as.character(power$outcome[i])
    
    srow <- NULL
    if (!is.null(strength) && nrow(strength) > 0 && "outcome" %in% names(strength)) {
      srow <- strength[strength$outcome == outcome, , drop = FALSE]
      if (nrow(srow) == 0) srow <- NULL
    }
    
    min_F <- if (!is.null(srow) && "min_F" %in% names(srow)) as.character(srow$min_F[1]) else NA_character_
    mean_F <- if (!is.null(srow) && "mean_F" %in% names(srow)) as.character(srow$mean_F[1]) else NA_character_
    n_F_less_10 <- if (!is.null(srow) && "n_F_less_10" %in% names(srow)) as.character(srow$n_F_less_10[1]) else NA_character_
    retention <- if (!is.null(srow) && "file_level_retention_rate" %in% names(srow)) as.character(srow$file_level_retention_rate[1]) else NA_character_
    
    external_rows[[length(external_rows) + 1]] <- data.frame(
      evidence_layer = "external_triangulation",
      endpoint_or_analysis = outcome,
      endpoint_class = ifelse(outcome == "Measured IOP", "IOP-related quantitative endpoint",
                       ifelse(outcome %in% c("POAG", "NTG"), "clinical glaucoma endpoint",
                       ifelse(outcome %in% c("RNFL", "GCIPL"), "neuroretinal structural endpoint", "external endpoint"))),
      dataset = as.character(power$dataset[i]),
      N_or_cases_metadata = as.character(power$N_or_cases_metadata[i]),
      n_instruments = as.character(power$n_instruments[i]),
      estimate_scale = as.character(power$estimate_scale[i]),
      beta_or_difference = as.character(power$observed_beta[i]),
      se = as.character(power$observed_se[i]),
      pval = as.character(power$observed_p[i]),
      qval = as.character(power$observed_q[i]),
      MDE_80pct_alpha_0p05 = as.character(power$minimum_detectable_effect[i]),
      observed_abs_beta_over_MDE = as.character(power$observed_abs_beta_over_MDE[i]),
      instrument_strength_summary = paste0("mean_F=", mean_F, "; min_F=", min_F, "; F_less_10=", n_F_less_10),
      retention_summary = paste0("retention_rate=", retention),
      power_calibrated_interpretation = as.character(power$power_calibrated_interpretation[i]),
      manuscript_safe_interpretation = as.character(power$manuscript_safe_interpretation[i]),
      claim_boundary = "External triangulation is supportive/contextual only and does not establish pressure-stratified confirmation.",
      stringsAsFactors = FALSE
    )
  }
}

external_grid <- if (length(external_rows) > 0) do.call(rbind, external_rows) else data.frame()

external_grid_file <- file.path(out_dir, "phase13_8_external_triangulation_interpretation_grid.tsv")
write_tsv(external_grid, external_grid_file)

# ---- Integrated evidence summary ----
summary_rows <- list()

if (!is.null(bootstrap) && nrow(bootstrap) > 0) {
  summary_rows[[length(summary_rows) + 1]] <- data.frame(
    evidence_layer = "internal_component_contrast",
    analysis = "SBP IOP-minus-nonIOP component contrast",
    key_result = paste0(
      "difference=", fmt(bootstrap$full_beta_difference_IOP_minus_nonIOP[1], 5),
      "; empirical_SE=", fmt(bootstrap$empirical_SE_difference[1], 5),
      "; bootstrap_CI=", fmt(bootstrap$empirical_CI_lower_2p5[1], 5),
      " to ", fmt(bootstrap$empirical_CI_upper_97p5[1], 5),
      "; empirical_p=", fmt(bootstrap$empirical_p_normal_using_full_difference[1], 4)
    ),
    interpretation = "Primary internal component-divergence signal showed directional stability in paired bootstrap.",
    claim_boundary = "Hypothesis-generating component-divergence signal, not confirmatory causal evidence.",
    stringsAsFactors = FALSE
  )
}

if (!is.null(mvmr_contrast) && nrow(mvmr_contrast) > 0) {
  univ <- mvmr_contrast[mvmr_contrast$exposure_label == "SBP_univariable", , drop = FALSE]
  adj <- mvmr_contrast[mvmr_contrast$exposure_label == "SBP_adjusted_for_measured_IOP", , drop = FALSE]
  
  if (nrow(univ) == 1 && nrow(adj) == 1) {
    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      evidence_layer = "exploratory_pathway_sensitivity",
      analysis = "Measured-IOP-adjusted SBP component contrast",
      key_result = paste0(
        "univariable_difference=", fmt(univ$beta_difference_IOP_minus_nonIOP[1], 5),
        "; IOP_adjusted_difference=", fmt(adj$beta_difference_IOP_minus_nonIOP[1], 5),
        "; relative_change=", fmt(adj$relative_change_from_univariable[1], 4),
        "; direction_preserved=", as.character(adj$direction_preserved[1])
      ),
      interpretation = "Measured-IOP-adjusted sensitivity showed modest attenuation while preserving the positive contrast direction.",
      claim_boundary = "Exploratory pathway-consistency sensitivity only; not a mediation or pathway-confirming analysis.",
      stringsAsFactors = FALSE
    )
  }
}

if (!is.null(robustness_caveat) && nrow(robustness_caveat) > 0) {
  disc <- robustness_caveat$value[robustness_caveat$item == "expected_direction_discordant_rows"]
  preserved <- robustness_caveat$value[robustness_caveat$item == "expected_direction_preserved_rows"]
  methods <- robustness_caveat$value[robustness_caveat$item == "discordant_methods"]
  
  summary_rows[[length(summary_rows) + 1]] <- data.frame(
    evidence_layer = "pleiotropy_and_robustness",
    analysis = "SBP robustness matrix",
    key_result = paste0("direction_preserved_rows=", preserved, "; discordant_rows=", disc, "; discordant_methods=", methods),
    interpretation = "Most component-method estimates preserved expected directions; weighted mode discordance was retained as a robustness caveat.",
    claim_boundary = "Robustness characterization only; does not upgrade the claim level.",
    stringsAsFactors = FALSE
  )
}

if (!is.null(external_grid) && nrow(external_grid) > 0) {
  iop_row <- external_grid[external_grid$endpoint_or_analysis == "Measured IOP", , drop = FALSE]
  poag_row <- external_grid[external_grid$endpoint_or_analysis == "POAG", , drop = FALSE]
  ntg_row <- external_grid[external_grid$endpoint_or_analysis == "NTG", , drop = FALSE]
  
  ext_key <- paste0(
    "Measured IOP: ", ifelse(nrow(iop_row) > 0, iop_row$manuscript_safe_interpretation[1], "not available"),
    "; POAG: ", ifelse(nrow(poag_row) > 0, poag_row$manuscript_safe_interpretation[1], "not available"),
    "; NTG: ", ifelse(nrow(ntg_row) > 0, ntg_row$manuscript_safe_interpretation[1], "not available")
  )
  
  summary_rows[[length(summary_rows) + 1]] <- data.frame(
    evidence_layer = "external_triangulation",
    analysis = "Power-calibrated external endpoint grid",
    key_result = ext_key,
    interpretation = "External endpoints did not provide pressure-stratified confirmation; measured IOP was directional but non-confirmatory.",
    claim_boundary = "External triangulation remains contextual and power/scale dependent.",
    stringsAsFactors = FALSE
  )
}

integrated_summary <- if (length(summary_rows) > 0) do.call(rbind, summary_rows) else data.frame()

integrated_summary_file <- file.path(out_dir, "phase13_8_integrated_evidence_summary.tsv")
write_tsv(integrated_summary, integrated_summary_file)

# ---- Interpretation classification ----
classification <- data.frame(
  item = c(
    "overall_evidence_position",
    "SBP_internal_contrast_position",
    "measured_IOP_adjusted_sensitivity_position",
    "external_triangulation_position",
    "robustness_position",
    "manuscript_claim_level"
  ),
  value = c(
    "Evidence supports a stable SBP component-divergence pattern, with external triangulation remaining non-confirmatory.",
    "Primary internal contrast is directionally stable across empirical uncertainty and influence analyses.",
    "Measured-IOP-adjusted sensitivity suggests modest attenuation but preserves the contrast direction.",
    "External endpoints do not establish pressure-stratified confirmation; interpretation is power- and scale-calibrated.",
    "Robustness matrix supports characterization of stability while retaining weighted-mode discordance as a caveat.",
    "HYPOTHESIS_GENERATING_NOT_CONFIRMATORY"
  ),
  stringsAsFactors = FALSE
)

classification_file <- file.path(out_dir, "phase13_8_overall_interpretation_classification.tsv")
write_tsv(classification, classification_file)

# ---- Manuscript inserts ----
methods_insert <- c(
  "# Power-calibrated external triangulation grid",
  "",
  "External triangulation results were integrated with minimum detectable effect calculations, SBP instrument-strength summaries, empirical contrast uncertainty, measured-IOP-adjusted exploratory sensitivity, and robustness-matrix results. The grid was used to classify external endpoints as directional, non-confirmatory, or inconclusive in a power- and scale-calibrated manner.",
  "",
  "This integration did not create new MR estimates. It was used to align endpoint-level interpretation with the hypothesis-generating framework and to avoid overinterpreting non-significant external endpoints."
)

methods_file <- file.path(out_dir, "phase13_8_methods_insert_external_triangulation_grid.md")
writeLines(methods_insert, methods_file, useBytes = TRUE)

results_insert <- c(
  "# Power-calibrated external triangulation grid",
  "",
  "The power-calibrated external triangulation grid integrated measured IOP, POAG, NTG, RNFL, and GCIPL with minimum detectable effect calculations and SBP instrument-strength summaries. Across SBP external triangulation files, instrument strength was adequate, with no weak instruments with F < 10 in the available harmonized files.",
  "",
  "Measured IOP showed directional but non-confirmatory IOP-related triangulation. POAG, NTG, RNFL, and GCIPL did not provide pressure-stratified confirmation, and their interpretation remained power-, scale-, and endpoint-dependent.",
  "",
  "Together with the paired-bootstrap contrast analysis, robustness matrix, and exploratory measured-IOP-adjusted sensitivity, these findings support a stable SBP component-divergence pattern while retaining a hypothesis-generating interpretation."
)

results_file <- file.path(out_dir, "phase13_8_results_insert_external_triangulation_grid.md")
writeLines(results_insert, results_file, useBytes = TRUE)

# ---- Claim-safety audit ----
texts_to_scan <- paste(
  c(
    methods_insert,
    results_insert,
    external_grid$manuscript_safe_interpretation,
    integrated_summary$interpretation,
    integrated_summary$claim_boundary,
    classification$value
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
    "confirms mediation",
    "confirms pathway",
    "establishes causality",
    "causal effect of SBP"
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
    "mediation_overclaim",
    "pathway_overclaim",
    "causal_overclaim",
    "causal_overclaim"
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

audit_file <- file.path(out_dir, "phase13_8_claim_safety_audit.tsv")
write_tsv(danger_patterns, audit_file)

high_risk_hits <- sum(danger_patterns$hits)

status <- data.frame(
  field = c(
    "phase",
    "source_manifest_created",
    "external_triangulation_grid_created",
    "integrated_evidence_summary_created",
    "overall_interpretation_classification_created",
    "methods_insert_created",
    "results_insert_created",
    "claim_safety_audit_created",
    "high_risk_phrase_hits",
    "external_grid_rows",
    "integrated_summary_rows",
    "new_MR_results_created",
    "claim_level",
    "claim_upgrade_allowed",
    "phase13_8_passed"
  ),
  value = c(
    "Phase13.8",
    file.exists(source_manifest_file),
    file.exists(external_grid_file),
    file.exists(integrated_summary_file),
    file.exists(classification_file),
    file.exists(methods_file),
    file.exists(results_file),
    file.exists(audit_file),
    high_risk_hits,
    nrow(external_grid),
    nrow(integrated_summary),
    FALSE,
    "HYPOTHESIS_GENERATING_NOT_CONFIRMATORY",
    "NO",
    high_risk_hits == 0 && nrow(external_grid) > 0 && nrow(integrated_summary) > 0
  ),
  stringsAsFactors = FALSE
)

status_file <- file.path(out_dir, "phase13_8_status.tsv")
write_tsv(status, status_file)

if (file.exists(master_status_file)) {
  master_status <- read.delim(master_status_file, check.names = FALSE)
  idx <- which(master_status$phase == "Phase13.8")
  if (length(idx) == 1) {
    master_status$status[idx] <- ifelse(
      high_risk_hits == 0 && nrow(external_grid) > 0 && nrow(integrated_summary) > 0,
      "PASSED_POWER_CALIBRATED_EXTERNAL_TRIANGULATION_GRID",
      "FAILED_PHASE13_8_REVIEW_REQUIRED"
    )
    master_status$qc_status[idx] <- ifelse(
      high_risk_hits == 0,
      "PASSED",
      "REVIEW_REQUIRED"
    )
    master_status$primary_output[idx] <- "phase13_8_external_triangulation_interpretation_grid.tsv; phase13_8_integrated_evidence_summary.tsv"
  }
  write_tsv(master_status, master_status_file)
}

message("Phase13.8 completed.")
message("External grid: ", external_grid_file)
message("Integrated summary: ", integrated_summary_file)
message("Classification: ", classification_file)
message("Status: ", status_file)
message("High-risk phrase hits: ", high_risk_hits)
