# Phase 10.0: file inventory and master checklist
# Purpose: lock core file paths before Phase 10.3 SNP/locus influence analysis.

options(stringsAsFactors = FALSE)

out_dir <- "../../58_component_divergence_reinforcement"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

inventory_out <- file.path(out_dir, "phase10_0_file_inventory.tsv")
checklist_out <- file.path(out_dir, "phase10_0_master_checklist.tsv")
status_out <- file.path(out_dir, "phase10_0_status.tsv")

files <- data.frame(
  file_role = character(),
  phase = character(),
  path = character(),
  required = logical(),
  notes = character(),
  stringsAsFactors = FALSE
)

add_file <- function(file_role, phase, path, required = TRUE, notes = "") {
  files <<- rbind(
    files,
    data.frame(
      file_role = file_role,
      phase = phase,
      path = path,
      required = required,
      notes = notes,
      stringsAsFactors = FALSE
    )
  )
}

# -----------------------------
# Core manuscript / prior evidence files
# -----------------------------
add_file(
  "current_main_manuscript_markdown",
  "Phase8.11",
  "../../54_topjournal_manuscript_revision/phase8_11_topjournal_manuscript_with_figures.md",
  TRUE,
  "Current manuscript source before Phase 10 rewrite."
)

add_file(
  "clean_manuscript_before_figures",
  "Phase8.7C",
  "../../54_topjournal_manuscript_revision/phase8_7C_topjournal_manuscript_draft_heading_polished.md",
  FALSE,
  "Earlier clean manuscript draft; useful fallback."
)

add_file(
  "final_number_crosscheck",
  "Phase8.8",
  "../../54_topjournal_manuscript_revision/phase8_8_final_number_crosscheck.tsv",
  TRUE,
  "Locked number cross-check."
)

add_file(
  "final_claim_crosscheck",
  "Phase8.8",
  "../../54_topjournal_manuscript_revision/phase8_8_final_claim_crosscheck.tsv",
  TRUE,
  "Locked claim-safety cross-check."
)

add_file(
  "final_figure_captions",
  "Phase8.9",
  "../../55_final_figure_caption_alignment/phase8_9_final_figure_captions.tsv",
  TRUE,
  "Caption source for current figures."
)

add_file(
  "figure_file_inventory",
  "Phase8.9",
  "../../55_final_figure_caption_alignment/phase8_9_figure_file_inventory.tsv",
  TRUE,
  "Existing figure file inventory."
)

add_file(
  "final_package_manifest",
  "Phase8.12",
  "../../56_final_submission_package/phase8_12_final_package_manifest.tsv",
  FALSE,
  "Earlier final-package manifest; supplementary reference."
)

add_file(
  "final_claim_safety_summary",
  "Phase8.12",
  "../../56_final_submission_package/phase8_12_final_claim_safety_summary.tsv",
  FALSE,
  "Earlier claim-safety summary; supplementary reference."
)

# -----------------------------
# Internal SBP / vascular evidence
# -----------------------------
add_file(
  "routeA_master_summary_matrix",
  "Phase6.1",
  "../../44_master_evidence_table/phase6_1_routeA_master_summary_matrix.tsv",
  TRUE,
  "Master evidence matrix."
)

add_file(
  "publication_readiness_flags",
  "Phase6.1",
  "../../44_master_evidence_table/phase6_1_publication_readiness_flags.tsv",
  TRUE,
  "Publication-readiness flags."
)

add_file(
  "permitted_claims",
  "Phase6.2",
  "../../45_final_claims_wording/phase6_2_final_permitted_claims.tsv",
  TRUE,
  "Permitted claim boundaries."
)

add_file(
  "unsafe_claims",
  "Phase6.2",
  "../../45_final_claims_wording/phase6_2_prohibited_or_unsafe_claims.tsv",
  TRUE,
  "Unsafe/prohibited claims."
)

add_file(
  "primary_vascular_component_table",
  "Phase6.3",
  "../../46_figure_table_preparation/table1_primary_vascular_component_evidence.tsv",
  TRUE,
  "Primary vascular component evidence table."
)

add_file(
  "external_routeA_triangulation_table",
  "Phase6.3",
  "../../46_figure_table_preparation/table2_external_routeA_triangulation.tsv",
  TRUE,
  "External triangulation table."
)

add_file(
  "figure2_internal_component_contrast_input",
  "Phase6.3",
  "../../46_figure_table_preparation/figure2_internal_component_contrast_plot_input.tsv",
  TRUE,
  "Input for internal component contrast figure."
)

add_file(
  "figure3_external_validation_input",
  "Phase6.3",
  "../../46_figure_table_preparation/figure3_external_validation_plot_input.tsv",
  TRUE,
  "Input for external triangulation figure."
)

add_file(
  "contrast_audit_summary",
  "Phase8.2",
  "../../49_contrast_audit/phase8_2_overall_summary.tsv",
  TRUE,
  "Contrast-definition and FDR-family audit summary."
)

add_file(
  "sbp_stability_summary",
  "Phase8.3",
  "../../50_sbp_sensitivity/phase8_3_overall_sbp_stability_summary.tsv",
  TRUE,
  "SBP stability summary."
)

# -----------------------------
# Phase 10.2 reinforcement outputs
# -----------------------------
add_file(
  "phase10_2A_2B_status",
  "Phase10.2",
  file.path(out_dir, "phase10_2A_2B_status.tsv"),
  TRUE,
  "Status file for extended covariance and empirical SE analyses."
)

add_file(
  "phase10_2A_extended_covariance_sensitivity",
  "Phase10.2A",
  file.path(out_dir, "phase10_2A_extended_covariance_sensitivity.tsv"),
  TRUE,
  "Extended covariance sensitivity across r values."
)

add_file(
  "phase10_2B_paired_full_data_reproduction",
  "Phase10.2B",
  file.path(out_dir, "phase10_2B_paired_full_data_reproduction.tsv"),
  TRUE,
  "Paired SNP-level data needed for Phase 10.3 influence analyses."
)

add_file(
  "phase10_2B_empirical_contrast_se_summary",
  "Phase10.2B",
  file.path(out_dir, "phase10_2B_empirical_contrast_se_summary.tsv"),
  TRUE,
  "Empirical contrast SE summary."
)

add_file(
  "phase10_2_final_reinforcement_evidence_matrix",
  "Phase10.2",
  file.path(out_dir, "phase10_2_final_reinforcement_evidence_matrix.tsv"),
  TRUE,
  "Final reinforcement evidence matrix."
)

add_file(
  "phase10_2_claim_boundary_matrix",
  "Phase10.2",
  file.path(out_dir, "phase10_2_claim_boundary_matrix.tsv"),
  TRUE,
  "Claim-boundary matrix after Phase 10.2."
)

add_file(
  "phase10_2C_status",
  "Phase10.2C",
  file.path(out_dir, "phase10_2C_status.tsv"),
  TRUE,
  "Manuscript-safe integration text status."
)

add_file(
  "phase10_2C_required_text_audit",
  "Phase10.2C",
  file.path(out_dir, "phase10_2C_required_text_audit.tsv"),
  TRUE,
  "Required text audit."
)

add_file(
  "phase10_2C_unsafe_phrase_audit",
  "Phase10.2C",
  file.path(out_dir, "phase10_2C_unsafe_phrase_audit.tsv"),
  TRUE,
  "Unsafe phrase audit."
)

add_file(
  "phase10_2C_results_insert",
  "Phase10.2C",
  file.path(out_dir, "phase10_2C_results_insert.txt"),
  FALSE,
  "Safe Results insert text."
)

add_file(
  "phase10_2C_discussion_insert",
  "Phase10.2C",
  file.path(out_dir, "phase10_2C_discussion_insert.txt"),
  FALSE,
  "Safe Discussion insert text."
)

add_file(
  "phase10_2C_methods_insert",
  "Phase10.2C",
  file.path(out_dir, "phase10_2C_methods_insert.txt"),
  FALSE,
  "Safe Methods/Supplement insert text."
)

# -----------------------------
# External triangulation / HTG / orthogonal evidence summaries
# -----------------------------
add_file(
  "htg_final_status",
  "Phase8.4",
  "../../51_htg_final_validation_attempt/phase8_4_overall_htg_final_status.tsv",
  FALSE,
  "HTG final availability status."
)

add_file(
  "ldsc_orthogonal_evidence_summary",
  "Phase8.5",
  "../../52_orthogonal_genetic_evidence/phase8_5_overall_orthogonal_evidence_summary.tsv",
  FALSE,
  "Orthogonal genetic evidence feasibility summary."
)

add_file(
  "ntg_htg_external_validation_grid",
  "Phase5.10F",
  "../../42_ntg_htg_validation_integration/phase5_10F_routeA_external_validation_grid.tsv",
  FALSE,
  "NTG/HTG external validation grid."
)

# -----------------------------
# Planned Phase 10 outputs
# -----------------------------
add_file(
  "phase10_3A_leave_one_snp_out_contrast",
  "Phase10.3 planned",
  file.path(out_dir, "phase10_3A_leave_one_snp_out_contrast.tsv"),
  FALSE,
  "Planned Phase 10.3 output."
)

add_file(
  "phase10_3A_leave_one_snp_out_summary",
  "Phase10.3 planned",
  file.path(out_dir, "phase10_3A_leave_one_snp_out_summary.tsv"),
  FALSE,
  "Planned Phase 10.3 output."
)

add_file(
  "phase10_3B_snp_contribution_to_contrast",
  "Phase10.3 planned",
  file.path(out_dir, "phase10_3B_snp_contribution_to_contrast.tsv"),
  FALSE,
  "Planned Phase 10.3 output."
)

add_file(
  "phase10_3B_top_contributor_summary",
  "Phase10.3 planned",
  file.path(out_dir, "phase10_3B_top_contributor_summary.tsv"),
  FALSE,
  "Planned Phase 10.3 output."
)

add_file(
  "phase10_3C_top_contribution_removal_results",
  "Phase10.3 planned",
  file.path(out_dir, "phase10_3C_top_contribution_removal_results.tsv"),
  FALSE,
  "Planned Phase 10.3 output."
)

add_file(
  "phase10_3D_radial_baujat_influence_input",
  "Phase10.3 planned",
  file.path(out_dir, "phase10_3D_radial_baujat_influence_input.tsv"),
  FALSE,
  "Planned Phase 10.3 output."
)

add_file(
  "phase10_4A_sbp_instrument_set_manifest",
  "Phase10.4 planned",
  file.path(out_dir, "phase10_4A_sbp_instrument_set_manifest.tsv"),
  FALSE,
  "Planned Phase 10.4 output."
)

add_file(
  "phase10_4B_alternative_sbp_instrument_results",
  "Phase10.4 planned",
  file.path(out_dir, "phase10_4B_alternative_sbp_instrument_results.tsv"),
  FALSE,
  "Planned Phase 10.4 output."
)

add_file(
  "phase10_4C_alternative_instrument_summary",
  "Phase10.4 planned",
  file.path(out_dir, "phase10_4C_alternative_instrument_summary.tsv"),
  FALSE,
  "Planned Phase 10.4 output."
)

add_file(
  "phase10_1_gwas_data_source_transparency_table",
  "Phase10.1 planned",
  file.path(out_dir, "phase10_1_gwas_data_source_transparency_table.tsv"),
  FALSE,
  "Planned data-source transparency table."
)

add_file(
  "phase10_5_external_specificity_triangulation_table",
  "Phase10.5 planned",
  file.path(out_dir, "phase10_5_external_specificity_triangulation_table.tsv"),
  FALSE,
  "Planned external specificity triangulation table."
)

add_file(
  "phase10_6_HTG_summary_statistics_availability_audit",
  "Phase10.6 planned",
  file.path(out_dir, "phase10_6_HTG_summary_statistics_availability_audit.tsv"),
  FALSE,
  "Planned structured HTG availability audit."
)

add_file(
  "phase10_10_submission_package_manifest",
  "Phase10.10 planned",
  file.path(out_dir, "phase10_10_submission_package_manifest.tsv"),
  FALSE,
  "Planned final submission package manifest."
)

# -----------------------------
# Utility functions
# -----------------------------
is_tabular <- function(path) {
  grepl("\\.(tsv|tab|csv)(\\.gz)?$", path, ignore.case = TRUE)
}

detect_sep <- function(path) {
  if (grepl("\\.csv(\\.gz)?$", path, ignore.case = TRUE)) "," else "\t"
}

safe_first_line <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  con <- NULL
  out <- NA_character_
  tryCatch({
    con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "rt") else file(path, "rt")
    x <- readLines(con, n = 1, warn = FALSE)
    if (length(x) > 0) out <- x[1]
  }, error = function(e) {
    out <<- NA_character_
  }, finally = {
    if (!is.null(con)) close(con)
  })
  out
}

safe_n_columns <- function(path) {
  if (!file.exists(path) || !is_tabular(path)) return(NA_integer_)
  first <- safe_first_line(path)
  if (is.na(first) || !nzchar(first)) return(NA_integer_)
  length(strsplit(first, detect_sep(path), fixed = TRUE)[[1]])
}

safe_n_rows <- function(path, max_size_bytes = 50 * 1024^2) {
  if (!file.exists(path) || !is_tabular(path)) return(NA_integer_)
  info <- file.info(path)
  if (is.na(info$size) || info$size > max_size_bytes) return(NA_integer_)
  if (grepl("\\.gz$", path, ignore.case = TRUE)) return(NA_integer_)
  con <- NULL
  n <- 0L
  tryCatch({
    con <- file(path, "rt")
    repeat {
      x <- readLines(con, n = 10000, warn = FALSE)
      if (length(x) == 0) break
      n <- n + length(x)
    }
  }, error = function(e) {
    n <<- NA_integer_
  }, finally = {
    if (!is.null(con)) close(con)
  })
  if (is.na(n) || n == 0) return(NA_integer_)
  max(n - 1L, 0L)
}

safe_md5 <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  out <- tryCatch(as.character(tools::md5sum(path)), error = function(e) NA_character_)
  unname(out)
}

read_tsv_safe <- function(path) {
  if (!file.exists(path)) return(data.frame())
  tryCatch(
    read.table(
      path,
      sep = "\t",
      header = TRUE,
      quote = "",
      comment.char = "",
      fill = TRUE,
      check.names = FALSE,
      stringsAsFactors = FALSE
    ),
    error = function(e) data.frame()
  )
}

first_status <- function(path) {
  x <- read_tsv_safe(path)
  if (nrow(x) == 0) return(NA_character_)
  if (!"status" %in% names(x)) return("PRESENT_NO_STATUS_COLUMN")
  as.character(x$status[1])
}

status_passed <- function(x) {
  !is.na(x) && grepl("PASS|PASSED|COMPLETED", x, ignore.case = TRUE)
}

# -----------------------------
# Build inventory
# -----------------------------
files$file_exists <- file.exists(files$path)
files$file_size_bytes <- ifelse(files$file_exists, file.info(files$path)$size, NA)
files$n_rows <- vapply(files$path, safe_n_rows, integer(1))
files$n_columns <- vapply(files$path, safe_n_columns, integer(1))
files$md5 <- vapply(files$path, safe_md5, character(1))

files$status <- ifelse(
  files$file_exists & files$required,
  "PASS_REQUIRED_PRESENT",
  ifelse(
    !files$file_exists & files$required,
    "FAIL_REQUIRED_MISSING",
    ifelse(
      files$file_exists & !files$required,
      "PASS_OPTIONAL_PRESENT",
      ifelse(
        grepl("planned", files$phase, ignore.case = TRUE),
        "PLANNED_NOT_YET_CREATED",
        "OPTIONAL_MISSING"
      )
    )
  )
)

inventory <- files[, c(
  "file_role",
  "phase",
  "path",
  "required",
  "file_exists",
  "file_size_bytes",
  "n_rows",
  "n_columns",
  "md5",
  "status",
  "notes"
)]

write.table(
  inventory,
  inventory_out,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = "NA"
)

# -----------------------------
# Build master checklist
# -----------------------------
required_missing <- inventory[inventory$required & !inventory$file_exists, , drop = FALSE]

phase10_2ab_status <- first_status(file.path(out_dir, "phase10_2A_2B_status.tsv"))
phase10_2c_status <- first_status(file.path(out_dir, "phase10_2C_status.tsv"))

phase10_2ab_files <- c(
  file.path(out_dir, "phase10_2A_2B_status.tsv"),
  file.path(out_dir, "phase10_2A_extended_covariance_sensitivity.tsv"),
  file.path(out_dir, "phase10_2B_paired_full_data_reproduction.tsv"),
  file.path(out_dir, "phase10_2B_empirical_contrast_se_summary.tsv")
)

phase10_2c_files <- c(
  file.path(out_dir, "phase10_2C_status.tsv"),
  file.path(out_dir, "phase10_2C_required_text_audit.tsv"),
  file.path(out_dir, "phase10_2C_unsafe_phrase_audit.tsv")
)

phase10_2ab_ready <- all(file.exists(phase10_2ab_files)) && status_passed(phase10_2ab_status)
phase10_2c_ready <- all(file.exists(phase10_2c_files)) && status_passed(phase10_2c_status)

phase10_3_ready <- phase10_2ab_ready &&
  file.exists(file.path(out_dir, "phase10_2B_paired_full_data_reproduction.tsv"))

checklist <- data.frame(
  phase = character(),
  priority = character(),
  checkpoint = character(),
  required_files_or_inputs = character(),
  pass_rule = character(),
  status = character(),
  blocking = logical(),
  next_action = character(),
  planned_outputs = character(),
  notes = character(),
  stringsAsFactors = FALSE
)

add_check <- function(
  phase,
  priority,
  checkpoint,
  required_files_or_inputs,
  pass_rule,
  status,
  blocking,
  next_action,
  planned_outputs,
  notes = ""
) {
  checklist <<- rbind(
    checklist,
    data.frame(
      phase = phase,
      priority = priority,
      checkpoint = checkpoint,
      required_files_or_inputs = required_files_or_inputs,
      pass_rule = pass_rule,
      status = status,
      blocking = blocking,
      next_action = next_action,
      planned_outputs = planned_outputs,
      notes = notes,
      stringsAsFactors = FALSE
    )
  )
}

add_check(
  "Phase10.0",
  "P0",
  "File inventory generated",
  inventory_out,
  "phase10_0_file_inventory.tsv exists and records required/planned files.",
  "PASS",
  FALSE,
  "Review missing required files, if any.",
  inventory_out,
  "Generated by this script."
)

add_check(
  "Phase10.0",
  "P0",
  "Required file presence audit",
  "All rows with required=TRUE in phase10_0_file_inventory.tsv",
  "No required file is missing.",
  ifelse(nrow(required_missing) == 0, "PASS", "FAIL"),
  nrow(required_missing) > 0,
  ifelse(
    nrow(required_missing) == 0,
    "Proceed to Phase 10.3 readiness check.",
    paste("Resolve missing required files:", paste(required_missing$file_role, collapse = "; "))
  ),
  checklist_out,
  paste0("n_required_missing=", nrow(required_missing))
)

add_check(
  "Phase10.2",
  "P0",
  "Extended covariance and empirical SE outputs available",
  paste(phase10_2ab_files, collapse = "; "),
  "Phase10.2A/10.2B status contains PASS/PASSED/COMPLETED and required files exist.",
  ifelse(phase10_2ab_ready, "PASS", "FAIL"),
  !phase10_2ab_ready,
  ifelse(
    phase10_2ab_ready,
    "Use paired SNP-level file as Phase 10.3 input.",
    "Re-run or locate Phase 10.2A/10.2B outputs before Phase 10.3."
  ),
  "phase10_3A_leave_one_snp_out_contrast.tsv; phase10_3B_snp_contribution_to_contrast.tsv",
  paste0("phase10_2A_2B_status=", phase10_2ab_status)
)

add_check(
  "Phase10.2C",
  "P1",
  "Manuscript-safe integration text available",
  paste(phase10_2c_files, collapse = "; "),
  "Phase10.2C status contains PASS/PASSED/COMPLETED and safety audits exist.",
  ifelse(phase10_2c_ready, "PASS", "WARN"),
  FALSE,
  ifelse(
    phase10_2c_ready,
    "Keep text for later manuscript rewrite; do not rewrite manuscript yet.",
    "Regenerate Phase 10.2C text before final manuscript rewrite."
  ),
  "phase10_2C_results_insert.txt; phase10_2C_discussion_insert.txt; phase10_2C_methods_insert.txt",
  paste0("phase10_2C_status=", phase10_2c_status)
)

add_check(
  "Phase10.3",
  "P0",
  "Ready for SNP/locus influence analysis",
  file.path(out_dir, "phase10_2B_paired_full_data_reproduction.tsv"),
  "Paired SNP-level reproduction file exists and Phase10.2A/10.2B passed.",
  ifelse(phase10_3_ready, "READY", "BLOCKED"),
  !phase10_3_ready,
  ifelse(
    phase10_3_ready,
    "Start Phase 10.3 SNP/locus influence analysis.",
    "Do not start Phase 10.3 until Phase 10.2B paired file and status are available."
  ),
  "phase10_3A_leave_one_snp_out_contrast.tsv; phase10_3A_leave_one_snp_out_summary.tsv; phase10_3B_snp_contribution_to_contrast.tsv; phase10_3B_top_contributor_summary.tsv; phase10_3C_top_contribution_removal_results.tsv",
  "Main question: whether a small number of SNPs/loci drive the SBP component contrast."
)

add_check(
  "Phase10.4",
  "P1",
  "Alternative SBP instrument definitions planned",
  "Phase10.3 outputs",
  "Start only after Phase 10.3 shows whether the contrast is not dominated by single SNPs/loci.",
  "PENDING",
  FALSE,
  "Wait for Phase 10.3 results.",
  "phase10_4A_sbp_instrument_set_manifest.tsv; phase10_4B_alternative_sbp_instrument_results.tsv; phase10_4C_alternative_instrument_summary.tsv",
  "Tests whether the signal depends on one SBP instrument definition."
)

add_check(
  "Phase10.1",
  "P1",
  "GWAS data-source transparency table planned",
  "Locked internal and external GWAS source files",
  "Prepare transparent source table after core statistical reinforcements are stable.",
  "PENDING",
  FALSE,
  "Defer until after Phase 10.3/10.4 unless needed earlier for drafting.",
  "phase10_1_gwas_data_source_transparency_table.tsv; phase10_1_data_source_missingness_audit.tsv",
  "This is transparency, not new causal evidence."
)

add_check(
  "Phase10.5_10.6",
  "P2",
  "External specificity and HTG audit planned",
  "Existing external triangulation files and HTG source review",
  "Summarize external evidence as directional but non-confirmatory; HTG unavailable unless verified full summary statistics are locked.",
  "PENDING",
  FALSE,
  "Defer until after Phase 10.3/10.4.",
  "phase10_5_external_specificity_triangulation_table.tsv; phase10_6_HTG_summary_statistics_availability_audit.tsv",
  "Avoid claiming external confirmation."
)

add_check(
  "Phase10.10",
  "P3",
  "Final submission lock planned",
  "Final manuscript, figures, tables, supplement, claim safety audit",
  "Only start after Phase 10.8/10.9 rewrite and figure/table rebuild.",
  "PENDING",
  FALSE,
  "Do not start yet.",
  "phase10_10_submission_package_manifest.tsv; phase10_10_final_claim_safety_audit.tsv",
  "Final lock comes last."
)

write.table(
  checklist,
  checklist_out,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = "NA"
)

overall_status <- ifelse(
  any(checklist$blocking),
  "BLOCKED_REQUIRED_FILES_MISSING_OR_PHASE10_2_NOT_READY",
  "PASSED_READY_FOR_PHASE10_3"
)

status <- data.frame(
  phase = "Phase10.0",
  status = overall_status,
  n_inventory_rows = nrow(inventory),
  n_required_files = sum(inventory$required),
  n_required_missing = nrow(required_missing),
  phase10_2A_2B_status = phase10_2ab_status,
  phase10_2C_status = phase10_2c_status,
  ready_for_phase10_3 = phase10_3_ready,
  file_inventory = inventory_out,
  master_checklist = checklist_out,
  output_directory = out_dir,
  timestamp = as.character(Sys.time()),
  stringsAsFactors = FALSE
)

write.table(
  status,
  status_out,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = "NA"
)

cat("Phase10.0 overall status:", overall_status, "\n")
cat("Inventory:", inventory_out, "\n")
cat("Checklist:", checklist_out, "\n")
cat("Status:", status_out, "\n")
