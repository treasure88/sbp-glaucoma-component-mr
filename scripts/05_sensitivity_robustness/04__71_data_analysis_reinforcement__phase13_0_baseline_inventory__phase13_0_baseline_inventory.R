options(stringsAsFactors = FALSE)

# Phase13.0: Baseline inventory and locked-result snapshot
# Project: IOP-dependent vs IOP-independent glaucoma component MR
# Evidence level: HYPOTHESIS_GENERATING_NOT_CONFIRMATORY
# This script does not create new statistical results.

message("Running Phase13.0 baseline inventory...")

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

safe_exists <- function(path) {
  file.exists(path) | dir.exists(path)
}

safe_size <- function(path) {
  if (!safe_exists(path)) return(NA_real_)
  if (dir.exists(path)) return(NA_real_)
  as.numeric(file.info(path)$size)
}

safe_md5 <- function(path) {
  if (!safe_exists(path)) return(NA_character_)
  if (dir.exists(path)) return(NA_character_)
  as.character(tools::md5sum(path))
}

# You are currently in the project root.
root <- normalizePath(".", winslash = "/", mustWork = TRUE)

out_dir <- file.path(root, "71_data_analysis_reinforcement", "phase13_0_baseline_inventory")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("Project root: ", root)
message("Output directory: ", out_dir)

items <- data.frame(
  role = c(
    "final_submission_package_dir",
    "final_semantic_audit_dir",
    "final_s1_patch_dir",
    "final_post_patch_lock_dir",
    "final_delivery_archive_dir",
    "main_manuscript_phase8_11",
    "submission_ready_manuscript_phase10_4",
    "figure1_png",
    "figure1_pdf",
    "figure2A_png",
    "figure2A_pdf",
    "figure2B_png",
    "figure2B_pdf",
    "figure3_png",
    "figure3_pdf",
    "final_S1_GWAS_data_source_transparency",
    "final_S2_external_specificity_triangulation",
    "final_submission_package_archive",
    "final_submission_audit_trail_archive"
  ),
  relative_path = c(
    "66_final_submission_package",
    "67_final_submission_semantic_audit",
    "68_final_s1_semantic_patch",
    "69_final_post_patch_lock",
    "70_final_delivery_archive",
    "54_topjournal_manuscript_revision/phase8_11_topjournal_manuscript_with_figures.md",
    "60_submission_text_revision_control/phase10_4_submission_ready_manuscript.md",
    "48_figures/Figure1_conceptual_model.png",
    "48_figures/Figure1_conceptual_model.pdf",
    "48_figures/Figure2A_internal_component_estimates.png",
    "48_figures/Figure2A_internal_component_estimates.pdf",
    "48_figures/Figure2B_internal_component_contrasts.png",
    "48_figures/Figure2B_internal_component_contrasts.pdf",
    "48_figures/Figure3_external_SBP_triangulation.png",
    "48_figures/Figure3_external_SBP_triangulation.pdf",
    "66_final_submission_package/Supplementary_Table_S1_GWAS_data_source_transparency.tsv",
    "66_final_submission_package/Supplementary_Table_S2_external_specificity_triangulation.tsv",
    "70_final_delivery_archive/final_submission_package_20260523_171756.tar.gz",
    "70_final_delivery_archive/final_submission_audit_trail_20260523_171756.tar.gz"
  ),
  stringsAsFactors = FALSE
)

items$resolved_path <- file.path(root, items$relative_path)
items$exists <- vapply(items$resolved_path, safe_exists, logical(1))
items$is_dir <- dir.exists(items$resolved_path)
items$size_bytes <- vapply(items$resolved_path, safe_size, numeric(1))
items$md5 <- vapply(items$resolved_path, safe_md5, character(1))

inventory_file <- file.path(out_dir, "phase13_0_baseline_file_inventory.tsv")
write_tsv(items, inventory_file)

locked_key_results <- data.frame(
  result_label = c(
    "SBP -> GBS_nonIOPcomponent",
    "SBP -> GBS_IOPcomponent",
    "SBP component contrast: IOP minus nonIOP",
    "DBP component contrast: IOP minus nonIOP",
    "Arterial stiffness component contrast: IOP minus nonIOP",
    "Hypertension liability",
    "SBP -> measured IOP external triangulation",
    "SBP -> POAG external triangulation",
    "SBP -> NTG external triangulation",
    "SBP -> RNFL external triangulation",
    "SBP -> GCIPL external triangulation",
    "HTG external triangulation"
  ),
  beta_or_beta_difference = c(
    -0.014988432,
    0.0077208916,
    0.022709324,
    0.017886294,
    0.5425936144,
    NA,
    0.0109703204093603,
    -0.00395132035025,
    -0.00288639234477424,
    0.00679328223717,
    0.0051631420971,
    NA
  ),
  se = c(
    0.0058051821,
    0.0030961462,
    0.0065792295,
    0.011463401,
    0.7233267392,
    NA,
    0.00573866502418868,
    0.00331537642249,
    0.00278002265855878,
    0.00779146629124,
    0.0102889553768,
    NA
  ),
  p = c(
    0.0098256472,
    0.012641568,
    0.00055713036,
    0.11869003,
    0.4531726434,
    NA,
    0.0559210698671773,
    0.233333201468,
    0.299148017948638,
    0.383269741568,
    0.615798347565,
    NA
  ),
  q = c(
    0.06320784,
    0.06320784,
    0.0027856518,
    0.29672507,
    NA,
    NA,
    0.111842139734355,
    0.466666402936,
    0.598296035897277,
    0.581038533929,
    0.615798347565,
    NA
  ),
  n_instruments = c(
    391,
    391,
    391,
    392,
    3,
    0,
    456,
    456,
    311,
    456,
    456,
    NA
  ),
  interpretation = c(
    "Nominal negative internal component estimate; not FDR-significant individually",
    "Nominal positive internal component estimate; not FDR-significant individually",
    "Primary hypothesis-generating directional component-divergence signal",
    "Contextual comparator; no clear component-divergence evidence",
    "Exploratory low-power directional context only",
    "No strict genome-wide significant instruments; not analyzable, not null evidence",
    "Directional positive borderline external triangulation; non-confirmatory",
    "Negative non-significant broad POAG result; non-confirmatory",
    "Negative non-significant NTG result; non-confirmatory",
    "Non-confirmatory neuroretinal endpoint",
    "Non-confirmatory neuroretinal endpoint",
    "No verified full downloadable HTG summary-statistics dataset locked"
  ),
  claim_level = c(
    "contextual",
    "contextual",
    "main_hypothesis_generating_result",
    "comparator",
    "exploratory_low_power",
    "not_analyzable",
    "directional_context_only",
    "directional_context_only",
    "directional_context_only",
    "contextual_only",
    "contextual_only",
    "limitation"
  ),
  stringsAsFactors = FALSE
)

locked_file <- file.path(out_dir, "phase13_0_locked_key_results.tsv")
write_tsv(locked_key_results, locked_file)

scope <- data.frame(
  phase = c(
    "Phase13.0",
    "Phase13.1",
    "Phase13.2",
    "Phase13.3",
    "Phase13.4",
    "Phase13.5",
    "Phase13.6",
    "Phase13.7",
    "Phase13.8",
    "Phase13.9",
    "Phase13.10"
  ),
  analysis_name = c(
    "Baseline inventory",
    "Primary analysis-family lock",
    "HTG rescue package",
    "Power and minimum detectable effect analysis",
    "Empirical covariance and block-jackknife contrast uncertainty",
    "MVMR pathway-consistency analysis",
    "Pleiotropy and robustness matrix",
    "Instrument strength and harmonization transparency",
    "Power-calibrated external triangulation grid",
    "Figure and table statistical redesign",
    "Final claim-safety audit"
  ),
  planned_output = c(
    "baseline file inventory and locked key results",
    "primary/contextual trait-role register",
    "HTG source inventory and final availability status",
    "external outcome MDE table",
    "bootstrap and block-jackknife contrast summaries",
    "MVMR results and conditional instrument strength",
    "method/outlier/annotation-filter robustness matrix",
    "F-statistic and SNP-retention tables",
    "integrated external interpretation table",
    "updated figure/table manifest",
    "unsafe phrase and claim-upgrade scan"
  ),
  claim_upgrade_allowed = rep("NO", 11),
  evidence_level = rep("HYPOTHESIS_GENERATING_NOT_CONFIRMATORY", 11),
  stringsAsFactors = FALSE
)

scope_file <- file.path(out_dir, "phase13_0_analysis_upgrade_scope.tsv")
write_tsv(scope, scope_file)

master_status <- data.frame(
  phase = scope$phase,
  analysis_name = scope$analysis_name,
  status = c("PASSED_BASELINE_INVENTORY_CREATED", rep("PENDING", 10)),
  primary_output = c(
    "phase13_0_baseline_file_inventory.tsv; phase13_0_locked_key_results.tsv",
    "phase13_1_primary_family_decision_register.tsv",
    "phase13_2D_htg_rescue_final_status.tsv",
    "phase13_3_external_power_mde_table.tsv",
    "phase13_4_bootstrap_and_jackknife_summary.tsv",
    "phase13_5_mvmr_results.tsv",
    "phase13_6_robustness_matrix.tsv",
    "phase13_7_instrument_strength_summary.tsv",
    "phase13_8_external_triangulation_interpretation_grid.tsv",
    "phase13_9_figure_table_manifest.tsv",
    "phase13_10_claim_safety_summary.tsv"
  ),
  qc_status = c("PENDING_MANUAL_REVIEW", rep("PENDING", 10)),
  manuscript_action = c(
    "None; baseline lock only",
    "Methods revision: analysis-family definition",
    "Supplement source audit; optional Results limitation refinement",
    "Supplementary MDE table and external interpretation refinement",
    "Methods/results sensitivity update",
    "Exploratory pathway-consistency section if feasible",
    "Supplementary robustness matrix",
    "Supplementary instrument transparency table",
    "Results/discussion external triangulation refinement",
    "Main figures/tables redesign",
    "Final text lock"
  ),
  claim_upgrade_allowed = rep("NO", 11),
  stringsAsFactors = FALSE
)

master_file <- file.path(root, "71_data_analysis_reinforcement", "phase13_master_status.tsv")
write_tsv(master_status, master_file)

required_roles <- c(
  "main_manuscript_phase8_11",
  "final_S1_GWAS_data_source_transparency",
  "final_S2_external_specificity_triangulation"
)

required_exists <- all(items$exists[items$role %in% required_roles])

recommended_roles <- c(
  "final_submission_package_dir",
  "final_semantic_audit_dir",
  "final_s1_patch_dir",
  "final_post_patch_lock_dir",
  "final_delivery_archive_dir"
)

recommended_exists <- all(items$exists[items$role %in% recommended_roles])

status <- data.frame(
  field = c(
    "phase",
    "project_root",
    "baseline_inventory_created",
    "locked_key_results_created",
    "analysis_upgrade_scope_created",
    "phase13_master_status_created",
    "required_core_files_found",
    "recommended_lock_dirs_found",
    "SBP_primary_contrast_locked",
    "claim_level",
    "new_statistical_results_created",
    "safe_to_proceed_to_phase13_1"
  ),
  value = c(
    "Phase13.0",
    root,
    file.exists(inventory_file),
    file.exists(locked_file),
    file.exists(scope_file),
    file.exists(master_file),
    required_exists,
    recommended_exists,
    TRUE,
    "HYPOTHESIS_GENERATING_NOT_CONFIRMATORY",
    FALSE,
    required_exists
  ),
  stringsAsFactors = FALSE
)

status_file <- file.path(out_dir, "phase13_0_status.tsv")
write_tsv(status, status_file)

message("Phase13.0 completed.")
message("Inventory: ", inventory_file)
message("Locked key results: ", locked_file)
message("Scope: ", scope_file)
message("Master status: ", master_file)
message("Status: ", status_file)

if (!required_exists) {
  warning("Some required core files were not found. Review phase13_0_baseline_file_inventory.tsv before continuing.")
}
