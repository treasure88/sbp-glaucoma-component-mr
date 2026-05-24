options(stringsAsFactors = FALSE)

message("Running Phase13.6B1 robustness matrix semantic patch...")

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

root <- normalizePath(".", winslash = "/", mustWork = TRUE)
out_dir <- file.path(root, "71_data_analysis_reinforcement", "phase13_6_pleiotropy_robustness")

matrix_file <- file.path(out_dir, "phase13_6B_SBP_robustness_matrix.tsv")
coverage_file <- file.path(out_dir, "phase13_6B_robustness_coverage_interpretation.tsv")
status_file <- file.path(out_dir, "phase13_6B_status.tsv")
audit_file <- file.path(out_dir, "phase13_6B_claim_safety_audit.tsv")
master_status_file <- file.path(root, "71_data_analysis_reinforcement", "phase13_master_status.tsv")

stopifnot(file.exists(matrix_file))
stopifnot(file.exists(coverage_file))

matrix <- read.delim(matrix_file, check.names = FALSE)
coverage <- read.delim(coverage_file, check.names = FALSE)

# Backup before patch.
backup_matrix <- file.path(out_dir, "phase13_6B_SBP_robustness_matrix.before_semantic_patch.tsv")
backup_coverage <- file.path(out_dir, "phase13_6B_robustness_coverage_interpretation.before_semantic_patch.tsv")
write_tsv(matrix, backup_matrix)
write_tsv(coverage, backup_coverage)

# ---- 1. Fix expected-direction flags using exact outcome labels ----
fix_expected_flag <- function(outcome_suffix, beta_or_difference, old_flag) {
  b <- as_num(beta_or_difference)
  if (is.na(b)) return(old_flag)
  
  if (identical(outcome_suffix, "GBS_nonIOPcomponent")) {
    return(paste0("expected_direction_preserved=", b < 0))
  }
  
  if (identical(outcome_suffix, "GBS_IOPcomponent")) {
    return(paste0("expected_direction_preserved=", b > 0))
  }
  
  old_flag
}

idx_component <- matrix$robustness_domain == "component_specific_method_consistency"

matrix$stability_flag[idx_component] <- mapply(
  fix_expected_flag,
  outcome_suffix = matrix$outcome_suffix[idx_component],
  beta_or_difference = matrix$beta_or_difference[idx_component],
  old_flag = matrix$stability_flag[idx_component],
  USE.NAMES = FALSE
)

# ---- 2. Add stricter safe interpretation for discordant robust methods ----
idx_discordant <- idx_component & grepl("expected_direction_preserved=FALSE", matrix$stability_flag, fixed = TRUE)

matrix$manuscript_safe_interpretation[idx_discordant] <- paste0(
  matrix$manuscript_safe_interpretation[idx_discordant],
  " Directionally discordant method-specific estimates are treated as robustness caveats rather than as primary evidence."
)

# ---- 3. Remove inherited risky wording from source-file interpretations ----
safe_replace <- function(x) {
  x <- gsub("Causal estimate from MR-PRESSO should be treated cautiously because it differs from the independent IVW estimate\\.",
            "The MR-PRESSO estimate should be treated cautiously because it differs from the independent IVW estimate.",
            x)
  x <- gsub("causal estimate", "estimate", x, ignore.case = TRUE)
  x <- gsub("confirmatory causal evidence", "confirmatory interpretation", x, ignore.case = TRUE)
  x <- gsub("confirmatory evidence", "confirmatory interpretation", x, ignore.case = TRUE)
  x
}

matrix$manuscript_safe_interpretation <- safe_replace(matrix$manuscript_safe_interpretation)
coverage$claim_safety_boundary <- safe_replace(coverage$claim_safety_boundary)

# ---- 4. Create a compact caveat summary ----
component_methods <- matrix[matrix$robustness_domain == "component_specific_method_consistency", ]

caveat_summary <- data.frame(
  item = c(
    "component_method_rows",
    "expected_direction_preserved_rows",
    "expected_direction_discordant_rows",
    "discordant_methods",
    "interpretation"
  ),
  value = c(
    nrow(component_methods),
    sum(grepl("expected_direction_preserved=TRUE", component_methods$stability_flag, fixed = TRUE)),
    sum(grepl("expected_direction_preserved=FALSE", component_methods$stability_flag, fixed = TRUE)),
    paste(unique(component_methods$method_or_check[
      grepl("expected_direction_preserved=FALSE", component_methods$stability_flag, fixed = TRUE)
    ]), collapse = "; "),
    "Discordant robust-method estimates are retained as caveats; they do not overturn the primary paired-bootstrap and IVW contrast sensitivity interpretation."
  ),
  stringsAsFactors = FALSE
)

caveat_file <- file.path(out_dir, "phase13_6B1_directional_caveat_summary.tsv")
write_tsv(caveat_summary, caveat_file)

# ---- 5. Re-run claim-safety audit ----
texts_to_scan <- paste(
  c(matrix$manuscript_safe_interpretation, coverage$claim_safety_boundary),
  collapse = "\n"
)

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
    "causal estimate",
    "confirmatory causal evidence"
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
    "causal_language",
    "causal_language"
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

write_tsv(danger_patterns, audit_file)

high_risk_hits <- sum(danger_patterns$hits)

# ---- 6. Write patched files ----
write_tsv(matrix, matrix_file)
write_tsv(coverage, coverage_file)

# ---- 7. Patch status ----
status <- data.frame(
  field = c(
    "phase",
    "patch",
    "matrix_backup_created",
    "coverage_backup_created",
    "robustness_matrix_patched",
    "coverage_interpretation_patched",
    "directional_caveat_summary_created",
    "high_risk_phrase_hits",
    "component_method_rows",
    "expected_direction_preserved_rows",
    "expected_direction_discordant_rows",
    "new_primary_MR_estimates_created",
    "claim_level",
    "claim_upgrade_allowed",
    "phase13_6B_passed"
  ),
  value = c(
    "Phase13.6B",
    "Phase13.6B1_semantic_patch",
    file.exists(backup_matrix),
    file.exists(backup_coverage),
    file.exists(matrix_file),
    file.exists(coverage_file),
    file.exists(caveat_file),
    high_risk_hits,
    nrow(component_methods),
    sum(grepl("expected_direction_preserved=TRUE", component_methods$stability_flag, fixed = TRUE)),
    sum(grepl("expected_direction_preserved=FALSE", component_methods$stability_flag, fixed = TRUE)),
    FALSE,
    "HYPOTHESIS_GENERATING_NOT_CONFIRMATORY",
    "NO",
    high_risk_hits == 0
  ),
  stringsAsFactors = FALSE
)

write_tsv(status, status_file)

# ---- 8. Update master ----
if (file.exists(master_status_file)) {
  master_status <- read.delim(master_status_file, check.names = FALSE)
  idx <- which(master_status$phase == "Phase13.6")
  if (length(idx) == 1) {
    master_status$status[idx] <- ifelse(
      high_risk_hits == 0,
      "PASSED_ROBUSTNESS_MATRIX_INTEGRATED_AND_SEMANTIC_PATCHED",
      "FAILED_PHASE13_6B1_REVIEW_REQUIRED"
    )
    master_status$qc_status[idx] <- ifelse(high_risk_hits == 0, "PASSED", "REVIEW_REQUIRED")
    master_status$primary_output[idx] <- "phase13_6B_SBP_robustness_matrix.tsv; phase13_6B_robustness_coverage_interpretation.tsv; phase13_6B1_directional_caveat_summary.tsv"
  }
  write_tsv(master_status, master_status_file)
}

message("Phase13.6B1 patch completed.")
message("High-risk phrase hits after patch: ", high_risk_hits)
message("Status: ", status_file)
message("Caveat summary: ", caveat_file)
