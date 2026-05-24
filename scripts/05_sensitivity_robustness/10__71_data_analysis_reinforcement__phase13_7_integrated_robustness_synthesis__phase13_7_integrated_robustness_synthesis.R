options(stringsAsFactors = FALSE)

message("Running Phase13.7 integrated robustness evidence table and manuscript-safe synthesis...")

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
  tryCatch(
    read.delim(path, check.names = FALSE, stringsAsFactors = FALSE),
    error = function(e) NULL
  )
}

safe_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

fmt <- function(x, digits = 4) {
  x <- safe_num(x)
  ifelse(is.na(x), NA, formatC(x, digits = digits, format = "fg", flag = "#"))
}

status_value <- function(path, field_name) {
  x <- safe_read(path)
  if (is.null(x) || nrow(x) == 0) return(NA_character_)
  if (!all(c("field", "value") %in% names(x))) return(NA_character_)
  y <- x$value[x$field == field_name]
  if (length(y) == 0) return(NA_character_)
  as.character(y[1])
}

metric_value <- function(path, metric_name) {
  x <- safe_read(path)
  if (is.null(x) || nrow(x) == 0) return(NA_character_)
  if (!all(c("metric", "value") %in% names(x))) return(NA_character_)
  y <- x$value[x$metric == metric_name]
  if (length(y) == 0) return(NA_character_)
  as.character(y[1])
}

first_existing <- function(paths) {
  for (p in paths) {
    if (file.exists(p)) return(p)
  }
  return(paths[1])
}

add_row <- function(rows, domain, analysis, status, key_result, manuscript_use,
                    claim_impact, table_destination, source_file) {
  rows[[length(rows) + 1]] <- data.frame(
    evidence_domain = domain,
    analysis = analysis,
    status = status,
    key_result = key_result,
    manuscript_use = manuscript_use,
    claim_impact = claim_impact,
    table_destination = table_destination,
    source_file = source_file,
    stringsAsFactors = FALSE
  )
  rows
}

root <- normalizePath(".", winslash = "/", mustWork = TRUE)

out_dir <- file.path(
  root,
  "71_data_analysis_reinforcement",
  "phase13_7_integrated_robustness_synthesis"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------
# Source files
# ---------------------------

paths <- list(
  primary_family_status = "71_data_analysis_reinforcement/phase13_1_primary_family_lock/phase13_1_status.tsv",
  primary_family_roles = "71_data_analysis_reinforcement/phase13_1_primary_family_lock/phase13_1_primary_secondary_trait_roles.tsv",
  htg_status = first_existing(c(
    "71_data_analysis_reinforcement/phase13_2_htg_rescue/phase13_2B_status.tsv",
    "71_data_analysis_reinforcement/phase13_2_htg_rescue/phase13_2A_status.tsv"
  )),
  htg_decision = first_existing(c(
    "71_data_analysis_reinforcement/phase13_2_htg_rescue/phase13_2B_htg_accession_mapping.tsv",
    "71_data_analysis_reinforcement/phase13_2_htg_rescue/phase13_2A_htg_rescue_decision_table.tsv"
  )),
  power_interpretation = "71_data_analysis_reinforcement/phase13_3_power_mde/phase13_3_external_power_mde_interpretation.tsv",
  bootstrap_summary = "71_data_analysis_reinforcement/phase13_4_empirical_covariance/phase13_4B_bootstrap_contrast_summary.tsv",
  chr_jackknife_summary = "71_data_analysis_reinforcement/phase13_4_empirical_covariance/phase13_4B_leave_one_chromosome_summary.tsv",
  locus_jackknife_summary = "71_data_analysis_reinforcement/phase13_4_empirical_covariance/phase13_4B_locus_proxy_jackknife_summary.tsv",
  mvmr_status = "71_data_analysis_reinforcement/phase13_5_mvmr_pathway/phase13_5B_status.tsv",
  mvmr_contrast = "71_data_analysis_reinforcement/phase13_5_mvmr_pathway/phase13_5B_sbp_iop_adjusted_contrast_comparison.tsv",
  sbp_dbp_status = "71_data_analysis_reinforcement/phase13_5_mvmr_pathway/phase13_5C_A1_status.tsv",
  sbp_dbp_decision = "71_data_analysis_reinforcement/phase13_5_mvmr_pathway/phase13_5C_A1_sbp_dbp_decision_lock.md",
  robustness_status = "71_data_analysis_reinforcement/phase13_6_pleiotropy_robustness/phase13_6B_status.tsv",
  robustness_matrix = "71_data_analysis_reinforcement/phase13_6_pleiotropy_robustness/phase13_6B_SBP_robustness_matrix.tsv",
  robustness_caveat = "71_data_analysis_reinforcement/phase13_6_pleiotropy_robustness/phase13_6B1_directional_caveat_summary.tsv",
  exclusion_status = "71_data_analysis_reinforcement/phase13_6_pleiotropy_robustness/phase13_6C_status.tsv",
  exclusion_summary = "71_data_analysis_reinforcement/phase13_6_pleiotropy_robustness/phase13_6C_exclusion_sensitivity_summary.tsv",
  steiger_status = "71_data_analysis_reinforcement/phase13_6_pleiotropy_robustness/phase13_6D_status.tsv",
  steiger_decision = "71_data_analysis_reinforcement/phase13_6_pleiotropy_robustness/phase13_6D_directionality_decision_table.tsv",
  steiger_audit = "71_data_analysis_reinforcement/phase13_6_pleiotropy_robustness/phase13_6D_steiger_feasibility_file_audit.tsv"
)

source_manifest <- data.frame(
  source_name = names(paths),
  relative_path = unlist(paths, use.names = FALSE),
  file_exists = file.exists(file.path(root, unlist(paths, use.names = FALSE))),
  stringsAsFactors = FALSE
)

write_tsv(
  source_manifest,
  file.path(out_dir, "phase13_7_source_manifest.tsv")
)

# ---------------------------
# Integrated robustness matrix
# ---------------------------

rows <- list()

# Phase13.1 primary family
pf_passed <- status_value(file.path(root, paths$primary_family_status), "phase13_1_passed")
rows <- add_row(
  rows,
  "analysis_family",
  "Biologically prioritized SBP primary exposure with contextual vascular traits",
  ifelse(identical(pf_passed, "TRUE"), "COMPLETED", "CHECK_SOURCE"),
  "SBP retained as the biologically prioritized primary vascular exposure; DBP, hypertension liability, and arterial stiffness retained as contextual or exploratory traits.",
  "Methods: analysis-family definition and multiplicity control.",
  "Clarifies statistical family; does not upgrade causal inference.",
  "Supplementary Methods / Supplementary Table",
  paths$primary_family_roles
)

# Phase13.2 HTG source audit
htg_passed <- status_value(file.path(root, paths$htg_status), "phase13_2B_passed")
if (is.na(htg_passed)) htg_passed <- status_value(file.path(root, paths$htg_status), "phase13_2A_passed")
rows <- add_row(
  rows,
  "external_endpoint_availability",
  "HTG source rescue and accession/source audit",
  ifelse(identical(htg_passed, "TRUE"), "COMPLETED", "CHECK_SOURCE_OR_PENDING"),
  "Structured HTG source review performed; no claim that HTG MR validation was completed unless a verified full HTG summary-statistics file is later recovered.",
  "Supplementary source audit; brief limitation in Discussion.",
  "Defines a validation gap; does not create an external confirmation claim.",
  "Supplementary Table",
  paths$htg_decision
)

# Phase13.3 power/MDE
power <- safe_read(file.path(root, paths$power_interpretation))
power_summary <- if (!is.null(power) && nrow(power) > 0) {
  paste0("Power-calibrated MDE interpretation available for ", nrow(power), " external endpoints.")
} else {
  "Power-calibrated MDE table not found."
}
rows <- add_row(
  rows,
  "external_triangulation_power",
  "Minimum detectable effect calibration for external outcomes",
  ifelse(!is.null(power), "COMPLETED", "MISSING"),
  power_summary,
  "Supplementary Table: Minimum detectable effect analysis for external triangulation outcomes.",
  "Separates non-confirmatory results into power-calibrated interpretation categories.",
  "Supplementary Table",
  paths$power_interpretation
)

# Phase13.4 bootstrap and jackknife
boot <- safe_read(file.path(root, paths$bootstrap_summary))
if (!is.null(boot) && nrow(boot) > 0) {
  boot_result <- paste0(
    "SBP contrast difference=", fmt(boot$full_beta_difference_IOP_minus_nonIOP[1], 5),
    "; empirical SE=", fmt(boot$empirical_SE_difference[1], 5),
    "; 95% bootstrap CI=", fmt(boot$empirical_CI_lower_2p5[1], 5),
    " to ", fmt(boot$empirical_CI_upper_97p5[1], 5),
    "; empirical p=", fmt(boot$empirical_p_normal_using_full_difference[1], 4),
    "; empirical r=", fmt(boot$empirical_correlation_IOP_nonIOP[1], 4)
  )
} else {
  boot_result <- "Bootstrap contrast summary not found."
}
rows <- add_row(
  rows,
  "internal_contrast_uncertainty",
  "Paired SNP bootstrap for SBP component contrast",
  ifelse(!is.null(boot), "COMPLETED", "MISSING"),
  boot_result,
  "Results and Supplementary robustness table.",
  "Supports stability of the component contrast but remains hypothesis-generating.",
  "Supplementary Table",
  paths$bootstrap_summary
)

chr_preserved <- metric_value(file.path(root, paths$chr_jackknife_summary), "all_direction_preserved")
chr_delta <- metric_value(file.path(root, paths$chr_jackknife_summary), "max_abs_delta")
locus_preserved <- metric_value(file.path(root, paths$locus_jackknife_summary), "all_direction_preserved")
locus_delta <- metric_value(file.path(root, paths$locus_jackknife_summary), "max_abs_delta")

rows <- add_row(
  rows,
  "internal_contrast_influence",
  "Leave-one-chromosome and locus-proxy jackknife",
  ifelse(file.exists(file.path(root, paths$chr_jackknife_summary)) && file.exists(file.path(root, paths$locus_jackknife_summary)), "COMPLETED", "CHECK_SOURCE"),
  paste0(
    "Chromosome direction preserved=", chr_preserved,
    "; chromosome max absolute delta=", chr_delta,
    "; locus direction preserved=", locus_preserved,
    "; locus max absolute delta=", locus_delta, "."
  ),
  "Supplementary robustness table.",
  "Addresses dominance by chromosome or local locus-proxy block; no claim upgrade.",
  "Supplementary Table",
  paste(paths$chr_jackknife_summary, paths$locus_jackknife_summary, sep = "; ")
)

# Phase13.5B measured-IOP-adjusted sensitivity
mvmr_passed <- status_value(file.path(root, paths$mvmr_status), "phase13_5B_passed")
mvmr <- safe_read(file.path(root, paths$mvmr_contrast))
mvmr_result <- "Measured-IOP-adjusted sensitivity table not found."
if (!is.null(mvmr) && nrow(mvmr) > 0 && "exposure_label" %in% names(mvmr)) {
  univ <- mvmr[mvmr$exposure_label == "SBP_univariable", , drop = FALSE]
  adj <- mvmr[mvmr$exposure_label == "SBP_adjusted_for_measured_IOP", , drop = FALSE]
  if (nrow(univ) == 1 && nrow(adj) == 1) {
    mvmr_result <- paste0(
      "Univariable contrast=", fmt(univ$beta_difference_IOP_minus_nonIOP[1], 5),
      "; measured-IOP-adjusted contrast=", fmt(adj$beta_difference_IOP_minus_nonIOP[1], 5),
      "; relative change=", fmt(adj$relative_change_from_univariable[1], 4),
      "; direction preserved=", as.character(adj$direction_preserved[1]), "."
    )
  }
}
rows <- add_row(
  rows,
  "exploratory_pathway_sensitivity",
  "Measured-IOP-adjusted SBP component contrast",
  ifelse(identical(mvmr_passed, "TRUE"), "COMPLETED_EXPLORATORY", "CHECK_SOURCE"),
  mvmr_result,
  "Brief Results sentence; Supplementary Methods/Table.",
  "Exploratory pathway-consistency sensitivity only; not mediation proof.",
  "Supplementary Table",
  paths$mvmr_contrast
)

# Phase13.5C_A1 SBP/DBP feasibility lock
sbp_dbp_passed <- status_value(file.path(root, paths$sbp_dbp_status), "phase13_5C_A1_passed")
rows <- add_row(
  rows,
  "bivariable_blood_pressure_feasibility",
  "SBP plus DBP bivariate sensitivity feasibility decision",
  ifelse(identical(sbp_dbp_passed, "TRUE"), "COMPLETED_DECISION_LOCK", "CHECK_SOURCE"),
  "SBP+DBP bivariate sensitivity was not pursued unless additional full cross-exposure association data are recovered, because overlap/collinearity made the model unsuitable for manuscript-level inference.",
  "Supplementary feasibility audit or limitation.",
  "Prevents unstable SBP-specific overinterpretation.",
  "Supplementary Methods",
  paths$sbp_dbp_decision
)

# Phase13.6B robustness matrix
robust_passed <- status_value(file.path(root, paths$robustness_status), "phase13_6B_passed")
caveat <- safe_read(file.path(root, paths$robustness_caveat))
caveat_summary <- "Robustness caveat summary not found."
if (!is.null(caveat) && nrow(caveat) > 0 && all(c("item", "value") %in% names(caveat))) {
  disc <- caveat$value[caveat$item == "expected_direction_discordant_rows"]
  preserved <- caveat$value[caveat$item == "expected_direction_preserved_rows"]
  methods <- caveat$value[caveat$item == "discordant_methods"]
  caveat_summary <- paste0(
    "Direction-preserved rows=", ifelse(length(preserved) > 0, preserved[1], NA),
    "; discordant rows=", ifelse(length(disc) > 0, disc[1], NA),
    "; discordant methods=", ifelse(length(methods) > 0, methods[1], NA), "."
  )
}
rows <- add_row(
  rows,
  "pleiotropy_and_method_robustness",
  "Integrated SBP robustness matrix",
  ifelse(identical(robust_passed, "TRUE"), "COMPLETED", "CHECK_SOURCE"),
  caveat_summary,
  "Supplementary robustness matrix and concise Results synthesis.",
  "Characterizes robustness and caveats; does not establish confirmatory causality.",
  "Supplementary Table",
  paths$robustness_matrix
)

# Phase13.6C instrument exclusion sensitivity
excl_passed <- status_value(file.path(root, paths$exclusion_status), "phase13_6C_passed")
excl_result <- paste0(
  "All exclusion scenarios direction preserved=",
  metric_value(file.path(root, paths$exclusion_summary), "all_exclusion_scenarios_direction_preserved"),
  "; minimum contrast=",
  fmt(metric_value(file.path(root, paths$exclusion_summary), "minimum_exclusion_contrast"), 5),
  "; maximum contrast=",
  fmt(metric_value(file.path(root, paths$exclusion_summary), "maximum_exclusion_contrast"), 5),
  "; maximum absolute delta=",
  fmt(metric_value(file.path(root, paths$exclusion_summary), "maximum_absolute_delta_from_full"), 5),
  "."
)
rows <- add_row(
  rows,
  "instrument_exclusion_sensitivity",
  "SBP instrument-exclusion sensitivity",
  ifelse(identical(excl_passed, "TRUE"), "COMPLETED_EXPLORATORY", "CHECK_SOURCE"),
  excl_result,
  "Supplementary robustness table; brief Results sentence.",
  "Supports instrument-level stability; no claim upgrade.",
  "Supplementary Table",
  paths$exclusion_summary
)

# Phase13.6D Steiger feasibility
steiger_passed <- status_value(file.path(root, paths$steiger_status), "phase13_6D_passed")
steiger_result <- paste0(
  "Formal Steiger feasible files=",
  status_value(file.path(root, paths$steiger_status), "formal_steiger_feasible_files"),
  "; proxy-only possible files=",
  status_value(file.path(root, paths$steiger_status), "proxy_only_possible_files"),
  "; claim upgrade allowed=",
  status_value(file.path(root, paths$steiger_status), "claim_upgrade_allowed"),
  "."
)
rows <- add_row(
  rows,
  "directionality_feasibility",
  "Steiger-style directionality feasibility audit",
  ifelse(identical(steiger_passed, "TRUE"), "COMPLETED_FEASIBILITY_AUDIT", "CHECK_SOURCE"),
  steiger_result,
  "Supplementary feasibility audit and limitation sentence.",
  "Formal Steiger not feasible with verified fields; raw variance proxies diagnostic only.",
  "Supplementary Table",
  paths$steiger_audit
)

integrated_matrix <- do.call(rbind, rows)

matrix_file <- file.path(out_dir, "phase13_7_integrated_robustness_matrix.tsv")
write_tsv(integrated_matrix, matrix_file)

# ---------------------------
# Supplementary table plan
# ---------------------------

supp_plan <- data.frame(
  proposed_table = c(
    "Supplementary Table X1. Analysis-family definition and trait roles",
    "Supplementary Table X2. HTG source rescue and accession audit",
    "Supplementary Table X3. Minimum detectable effect analysis for external triangulation outcomes",
    "Supplementary Table X4. Empirical covariance, bootstrap, and jackknife sensitivity for the SBP component contrast",
    "Supplementary Table X5. Exploratory measured-IOP-adjusted SBP component-contrast sensitivity",
    "Supplementary Table X6. SBP pleiotropy and method-robustness matrix",
    "Supplementary Table X7. SBP instrument-exclusion sensitivity",
    "Supplementary Table X8. Steiger-style directionality feasibility audit"
  ),
  source_phase = c(
    "Phase13.1",
    "Phase13.2",
    "Phase13.3",
    "Phase13.4",
    "Phase13.5B",
    "Phase13.6B",
    "Phase13.6C",
    "Phase13.6D"
  ),
  manuscript_role = c(
    "Clarifies primary/contextual exposure hierarchy.",
    "Documents pressure-stratified glaucoma validation gap.",
    "Power-calibrates non-confirmatory external endpoint results.",
    "Addresses covariance and influence concerns for the internal contrast.",
    "Exploratory pathway-consistency sensitivity; not mediation.",
    "Summarizes method-level and pleiotropy robustness.",
    "Evaluates sensitivity to outlier and IOP-associated instrument exclusions.",
    "Explains why formal Steiger directionality was not used as primary robustness evidence."
  ),
  claim_boundary = c(
    "Does not convert SBP into confirmatory causal exposure.",
    "Does not claim HTG validation was performed.",
    "Does not treat external nulls as definitive absence of effect.",
    "Does not establish causal direction.",
    "Does not confirm mediation through IOP.",
    "Does not upgrade inference beyond hypothesis-generating.",
    "Does not create new primary MR estimates.",
    "Does not establish causal directionality."
  ),
  stringsAsFactors = FALSE
)

supp_plan_file <- file.path(out_dir, "phase13_7_supplementary_table_plan.tsv")
write_tsv(supp_plan, supp_plan_file)

# ---------------------------
# Manuscript-safe Methods insert
# ---------------------------

methods_lines <- c(
  "# Integrated robustness and feasibility synthesis",
  "",
  "We integrated the robustness analyses into a prespecified synthesis framework for the SBP component-contrast signal. This synthesis did not generate new primary MR estimates. Instead, it summarized the analysis-family decision, empirical uncertainty of the SBP IOP-minus-nonIOP component contrast, chromosome- and locus-level influence checks, instrument-exclusion sensitivity, exploratory measured-IOP-adjusted sensitivity, external power-calibrated triangulation, HTG source availability, and Steiger-style directionality feasibility.",
  "",
  "SBP was retained as the biologically prioritized primary vascular exposure for the component-contrast analysis, whereas DBP, hypertension liability, and arterial stiffness were treated as contextual or exploratory vascular traits. External endpoint results were interpreted using minimum detectable effect calibration where available. Directionality screening was considered but not implemented as a formal Steiger robustness analysis because comparable trait-scale variance parameters were not verified across the component and external endpoint files. Raw variance-proxy diagnostics were retained for diagnostic context only.",
  "",
  "All integrated outputs were interpreted within a hypothesis-generating framework. Robustness analyses were used to evaluate stability, influence, feasibility, and interpretation boundaries rather than to establish confirmatory causal evidence."
)

methods_file <- file.path(out_dir, "phase13_7_manuscript_methods_insert.md")
writeLines(methods_lines, methods_file, useBytes = TRUE)

# ---------------------------
# Manuscript-safe Results insert
# ---------------------------

results_lines <- c(
  "# Integrated robustness and feasibility synthesis",
  "",
  paste0(
    "The integrated robustness synthesis supported the stability of the SBP component-divergence pattern while retaining a hypothesis-generating interpretation. The paired bootstrap analysis yielded an IOP-minus-nonIOP contrast of ",
    if (!is.null(boot) && nrow(boot) > 0) fmt(boot$full_beta_difference_IOP_minus_nonIOP[1], 5) else "[not available]",
    " with an empirical SE of ",
    if (!is.null(boot) && nrow(boot) > 0) fmt(boot$empirical_SE_difference[1], 5) else "[not available]",
    "."
  ),
  "",
  paste0(
    "Influence analyses did not indicate that the positive contrast direction depended on a single chromosome or locus-proxy block. Instrument-exclusion sensitivity analyses preserved the positive IOP-minus-nonIOP contrast direction across exclusion scenarios, with the contrast ranging from ",
    fmt(metric_value(file.path(root, paths$exclusion_summary), "minimum_exclusion_contrast"), 5),
    " to ",
    fmt(metric_value(file.path(root, paths$exclusion_summary), "maximum_exclusion_contrast"), 5),
    "."
  ),
  "",
  paste0(
    "In the exploratory measured-IOP-adjusted sensitivity analysis, the SBP component contrast was modestly attenuated but the positive contrast direction was preserved. This analysis was interpreted as pathway-consistency context only and not as evidence confirming mediation through IOP."
  ),
  "",
  paste0(
    "A Steiger-style directionality screen was considered but was not feasible as a formal robustness analysis using the currently verified trait-scale fields. The directionality audit was therefore retained as a feasibility and limitation assessment, and raw variance-proxy diagnostics were not used to upgrade causal directionality claims."
  ),
  "",
  "Overall, these analyses support a stable, hypothesis-generating SBP component-divergence signal while preserving the non-confirmatory interpretation of the external and directionality evidence."
)

results_file <- file.path(out_dir, "phase13_7_manuscript_results_insert.md")
writeLines(results_lines, results_file, useBytes = TRUE)

# ---------------------------
# Claim-safety audit
# ---------------------------

created_text <- paste(
  paste(readLines(methods_file, warn = FALSE), collapse = "\n"),
  paste(readLines(results_file, warn = FALSE), collapse = "\n"),
  paste(apply(integrated_matrix, 1, paste, collapse = " "), collapse = "\n"),
  sep = "\n"
)

high_risk <- data.frame(
  phrase = c(
    "SBP causally affects glaucoma",
    "SBP increases IOP-dependent glaucoma risk",
    "SBP protects against NTG",
    "external validation confirmed",
    "validated the mechanism",
    "confirms mediation",
    "confirms pathway",
    "HTG validation was performed",
    "hypertension has no effect",
    "definitive proof"
  ),
  category = c(
    "causal_overclaim",
    "component_causal_overclaim",
    "protective_overclaim",
    "external_confirmation_overclaim",
    "mechanism_validation_overclaim",
    "mediation_overclaim",
    "pathway_overclaim",
    "HTG_overclaim",
    "hypertension_null_overclaim",
    "proof_overclaim"
  ),
  hits = NA_integer_,
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(high_risk))) {
  high_risk$hits[i] <- length(gregexpr(high_risk$phrase[i], created_text, fixed = TRUE)[[1]])
  if (identical(gregexpr(high_risk$phrase[i], created_text, fixed = TRUE)[[1]], -1L)) {
    high_risk$hits[i] <- 0L
  }
}

audit_file <- file.path(out_dir, "phase13_7_claim_safety_audit.tsv")
write_tsv(high_risk, audit_file)

high_risk_hits <- sum(high_risk$hits, na.rm = TRUE)

# ---------------------------
# Status
# ---------------------------

critical_sources <- c(
  paths$bootstrap_summary,
  paths$exclusion_summary,
  paths$steiger_status
)

critical_found <- file.exists(file.path(root, critical_sources))

status <- data.frame(
  field = c(
    "phase",
    "source_manifest_created",
    "integrated_matrix_created",
    "supplementary_table_plan_created",
    "methods_insert_created",
    "results_insert_created",
    "claim_safety_audit_created",
    "source_files_audited",
    "source_files_found",
    "critical_sources_found",
    "high_risk_phrase_hits",
    "new_primary_MR_estimates_created",
    "claim_level",
    "claim_upgrade_allowed",
    "phase13_7_passed"
  ),
  value = c(
    "Phase13.7",
    file.exists(file.path(out_dir, "phase13_7_source_manifest.tsv")),
    file.exists(matrix_file),
    file.exists(supp_plan_file),
    file.exists(methods_file),
    file.exists(results_file),
    file.exists(audit_file),
    length(paths),
    sum(source_manifest$file_exists),
    paste(critical_found, collapse = ";"),
    high_risk_hits,
    "FALSE",
    "HYPOTHESIS_GENERATING_NOT_CONFIRMATORY",
    "NO",
    all(critical_found) && high_risk_hits == 0
  ),
  stringsAsFactors = FALSE
)

status_file <- file.path(out_dir, "phase13_7_status.tsv")
write_tsv(status, status_file)

# Update master status if available
master_status_file <- file.path(root, "71_data_analysis_reinforcement", "phase13_master_status.tsv")
if (file.exists(master_status_file)) {
  master_status <- safe_read(master_status_file)
  if (!is.null(master_status) && all(c("phase", "status", "primary_output", "claim_upgrade_allowed") %in% names(master_status))) {
    if ("Phase13.7" %in% master_status$phase) {
      idx <- which(master_status$phase == "Phase13.7")
      master_status$status[idx] <- ifelse(all(critical_found) && high_risk_hits == 0, "PASSED", "CHECK_REQUIRED")
      master_status$primary_output[idx] <- "phase13_7_integrated_robustness_matrix.tsv; phase13_7_manuscript_methods_insert.md; phase13_7_manuscript_results_insert.md"
      master_status$claim_upgrade_allowed[idx] <- "NO"
    } else {
      master_status <- rbind(
        master_status,
        data.frame(
          phase = "Phase13.7",
          objective = "Integrated robustness evidence table and manuscript-safe synthesis",
          status = ifelse(all(critical_found) && high_risk_hits == 0, "PASSED", "CHECK_REQUIRED"),
          primary_output = "phase13_7_integrated_robustness_matrix.tsv; phase13_7_manuscript_methods_insert.md; phase13_7_manuscript_results_insert.md",
          manuscript_impact = "Methods/Results robustness synthesis and supplementary table plan",
          claim_upgrade_allowed = "NO",
          stringsAsFactors = FALSE
        )[names(master_status)]
      )
    }
    write_tsv(master_status, master_status_file)
  }
}

message("Phase13.7 completed.")
message("Status: ", status_file)
message("Integrated matrix: ", matrix_file)
message("Methods insert: ", methods_file)
message("Results insert: ", results_file)
