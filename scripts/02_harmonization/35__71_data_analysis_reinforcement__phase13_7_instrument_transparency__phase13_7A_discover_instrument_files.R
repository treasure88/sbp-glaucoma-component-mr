options(stringsAsFactors = FALSE)

# Phase13.7A: Discover SNP-level instrument / harmonization files
# Project: IOP-dependent vs IOP-independent glaucoma component MR
# Evidence level: HYPOTHESIS_GENERATING_NOT_CONFIRMATORY
# This script only inventories files and columns. It does not create new MR estimates.

message("Running Phase13.7A instrument/harmonization file discovery...")

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

safe_exists <- function(path) file.exists(path) | dir.exists(path)

root <- normalizePath(".", winslash = "/", mustWork = TRUE)
out_dir <- file.path(root, "71_data_analysis_reinforcement", "phase13_7_instrument_transparency")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

master_status_file <- file.path(root, "71_data_analysis_reinforcement", "phase13_master_status.tsv")

# Candidate project directories likely to contain analysis-level tables.
candidate_dirs <- c(
  "03_file_inspection",
  "05_harmonization_planning",
  "22_sbp_robustness",
  "24_vascular_panel_integration",
  "25_external_outcome_triangulation",
  "27_external_outcome_standardized",
  "28_external_outcome_triangulation_inputs",
  "29_external_neuroretinal_mr_results",
  "30_external_triangulation_integration",
  "32_poag_external_validation_inputs",
  "33_poag_external_mr_results",
  "36_iop_external_validation_inputs",
  "37_iop_external_mr_results",
  "39_ntg_htg_external_validation",
  "40_ntg_external_validation_inputs",
  "41_ntg_external_mr_results",
  "42_ntg_htg_validation_integration",
  "44_master_evidence_table",
  "46_figure_table_preparation",
  "66_final_submission_package"
)

candidate_dirs <- candidate_dirs[dir.exists(file.path(root, candidate_dirs))]

message("Candidate directories found: ", paste(candidate_dirs, collapse = ", "))

# File extensions to scan. Header only, so gz files are OK.
patterns <- "\\.(tsv|txt|csv|tsv\\.gz|txt\\.gz|csv\\.gz)$"

all_files <- unique(unlist(lapply(candidate_dirs, function(d) {
  list.files(
    file.path(root, d),
    pattern = patterns,
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )
})))

# Avoid massive raw GWAS files if accidentally present in scanned folders.
# Keep files under 500 MB for header/sample inspection.
file_size <- suppressWarnings(as.numeric(file.info(all_files)$size))
keep <- !is.na(file_size) & file_size <= 500 * 1024^2
all_files <- all_files[keep]
file_size <- file_size[keep]

message("Files selected for header scan: ", length(all_files))

open_connection <- function(path) {
  if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    gzfile(path, open = "rt")
  } else {
    file(path, open = "rt")
  }
}

get_header <- function(path) {
  con <- NULL
  out <- NA_character_
  tryCatch({
    con <- open_connection(path)
    tmp <- readLines(con, n = 1, warn = FALSE)
    if (length(tmp) == 0 || is.na(tmp[1]) || !nzchar(tmp[1])) {
      out <- NA_character_
    } else {
      out <- tmp[1]
    }
  }, error = function(e) {
    out <<- NA_character_
  }, finally = {
    if (!is.null(con)) try(close(con), silent = TRUE)
  })
  out
}

guess_sep <- function(header_line, path) {
  if (length(header_line) != 1 || is.na(header_line) || !nzchar(header_line)) return("\t")
  count_fixed <- function(text, pattern) {
    hit <- gregexpr(pattern, text, fixed = TRUE)[[1]]
    if (length(hit) == 1 && hit[1] == -1) return(0L)
    length(hit)
  }
  n_tab <- count_fixed(header_line, "\t")
  n_comma <- count_fixed(header_line, ",")
  if (grepl("\\.csv(\\.gz)?$", path, ignore.case = TRUE)) return(",")
  if (n_comma > n_tab) return(",")
  "\t"
}

split_header <- function(header_line, sep) {
  if (length(header_line) != 1 || is.na(header_line) || !nzchar(header_line)) return(character(0))
  out <- strsplit(header_line, split = sep, fixed = TRUE)[[1]]
  out <- trimws(out)
  out[nzchar(out)]
}

has_any <- function(cols_lower, patterns) {
  any(cols_lower %in% tolower(patterns))
}

has_regex <- function(cols_lower, pattern) {
  any(grepl(pattern, cols_lower, ignore.case = TRUE))
}

classify_file <- function(cols) {
  cl <- tolower(cols)
  
  has_snp <- has_any(cl, c(
    "snp", "rsid", "rs_id", "rs_number", "variant", "variant_id", "markername", "marker_name"
  ))
  
  has_chr <- has_any(cl, c("chr", "chromosome", "chrom"))
  has_pos <- has_any(cl, c("pos", "position", "bp", "base_pair_location"))
  
  has_beta_exp <- has_any(cl, c(
    "beta.exposure", "beta_exposure", "beta_exp", "bx", "b_exp", "beta.x", "beta_x"
  )) || has_regex(cl, "^beta.*exposure$")
  
  has_se_exp <- has_any(cl, c(
    "se.exposure", "se_exposure", "se_exp", "sex", "se_exp", "se.x", "se_x"
  )) || has_regex(cl, "^se.*exposure$")
  
  has_p_exp <- has_any(cl, c(
    "pval.exposure", "p.exposure", "p_exposure", "pval_exposure", "p_exp", "pval_exp"
  )) || has_regex(cl, "p.*exposure")
  
  has_beta_out <- has_any(cl, c(
    "beta.outcome", "beta_outcome", "beta_out", "by", "b_out", "beta.y", "beta_y"
  )) || has_regex(cl, "^beta.*outcome$")
  
  has_se_out <- has_any(cl, c(
    "se.outcome", "se_outcome", "se_out", "sey", "se.y", "se_y"
  )) || has_regex(cl, "^se.*outcome$")
  
  has_p_out <- has_any(cl, c(
    "pval.outcome", "p.outcome", "p_outcome", "pval_outcome", "p_out", "pval_out"
  )) || has_regex(cl, "p.*outcome")
  
  has_ea <- has_any(cl, c(
    "effect_allele", "effect_allele.exposure", "effect_allele_exposure", "ea", "a1", "allele1"
  ))
  
  has_oa <- has_any(cl, c(
    "other_allele", "other_allele.exposure", "other_allele_exposure", "oa", "a2", "allele2"
  ))
  
  has_locus <- has_any(cl, c("locus", "ld_block", "block", "locus_id", "clump_id"))
  
  has_exposure_name <- has_any(cl, c("exposure", "exposure_id", "trait", "trait_name"))
  has_outcome_name <- has_any(cl, c("outcome", "outcome_id"))
  
  can_compute_F <- has_beta_exp && has_se_exp
  looks_harmonized <- has_snp && has_beta_exp && has_se_exp && has_beta_out && has_se_out
  looks_exposure_instrument <- has_snp && has_beta_exp && has_se_exp && !has_beta_out
  looks_result_summary <- has_regex(cl, "beta") && has_regex(cl, "se") && has_regex(cl, "p")
  
  file_class <- "other_table"
  if (looks_harmonized) file_class <- "harmonized_snp_level_candidate"
  if (looks_exposure_instrument) file_class <- "exposure_instrument_candidate"
  if (!looks_harmonized && !looks_exposure_instrument && looks_result_summary) file_class <- "summary_result_candidate"
  
  list(
    has_snp = has_snp,
    has_chr = has_chr,
    has_pos = has_pos,
    has_locus = has_locus,
    has_beta_exposure = has_beta_exp,
    has_se_exposure = has_se_exp,
    has_p_exposure = has_p_exp,
    has_beta_outcome = has_beta_out,
    has_se_outcome = has_se_out,
    has_p_outcome = has_p_out,
    has_effect_allele = has_ea,
    has_other_allele = has_oa,
    has_exposure_name = has_exposure_name,
    has_outcome_name = has_outcome_name,
    can_compute_F = can_compute_F,
    looks_harmonized = looks_harmonized,
    looks_exposure_instrument = looks_exposure_instrument,
    file_class = file_class
  )
}

scan_rows <- list()

for (i in seq_along(all_files)) {
  f <- all_files[i]
  header <- get_header(f)
  sep <- guess_sep(header, f)
  cols <- split_header(header, sep)
  cls <- classify_file(cols)
  
  rel <- sub(paste0("^", gsub("\\\\", "/", root), "/?"), "", normalizePath(f, winslash = "/", mustWork = FALSE))
  
  scan_rows[[i]] <- data.frame(
    relative_path = rel,
    size_bytes = file_size[i],
    n_columns = length(cols),
    separator_guess = ifelse(sep == "\t", "tab", "comma"),
    file_class = cls$file_class,
    has_snp = cls$has_snp,
    has_chr = cls$has_chr,
    has_pos = cls$has_pos,
    has_locus = cls$has_locus,
    has_beta_exposure = cls$has_beta_exposure,
    has_se_exposure = cls$has_se_exposure,
    has_p_exposure = cls$has_p_exposure,
    has_beta_outcome = cls$has_beta_outcome,
    has_se_outcome = cls$has_se_outcome,
    has_p_outcome = cls$has_p_outcome,
    has_effect_allele = cls$has_effect_allele,
    has_other_allele = cls$has_other_allele,
    has_exposure_name = cls$has_exposure_name,
    has_outcome_name = cls$has_outcome_name,
    can_compute_F = cls$can_compute_F,
    looks_harmonized = cls$looks_harmonized,
    looks_exposure_instrument = cls$looks_exposure_instrument,
    columns_preview = paste(head(cols, 80), collapse = "|"),
    stringsAsFactors = FALSE
  )
}

scan_table <- if (length(scan_rows) > 0) do.call(rbind, scan_rows) else data.frame()

scan_file <- file.path(out_dir, "phase13_7A_candidate_file_column_scan.tsv")
write_tsv(scan_table, scan_file)

# Prioritize likely usable files.
priority <- scan_table
if (nrow(priority) > 0) {
  priority$priority_score <- 
    4 * as.integer(priority$looks_harmonized) +
    3 * as.integer(priority$looks_exposure_instrument) +
    2 * as.integer(priority$can_compute_F) +
    1 * as.integer(priority$has_locus) +
    1 * as.integer(priority$has_chr & priority$has_pos) +
    1 * as.integer(priority$has_effect_allele & priority$has_other_allele)
  
  priority <- priority[order(-priority$priority_score, priority$relative_path), ]
  priority <- priority[priority$priority_score > 0, ]
} else {
  priority$priority_score <- numeric(0)
}

priority_file <- file.path(out_dir, "phase13_7A_prioritized_instrument_file_candidates.tsv")
write_tsv(priority, priority_file)

# Create recommendation table for next phase.
n_harmonized <- sum(scan_table$looks_harmonized, na.rm = TRUE)
n_exposure_instr <- sum(scan_table$looks_exposure_instrument, na.rm = TRUE)
n_can_compute_F <- sum(scan_table$can_compute_F, na.rm = TRUE)
n_locus <- sum(scan_table$has_locus | (scan_table$has_chr & scan_table$has_pos), na.rm = TRUE)

recommendations <- data.frame(
  item = c(
    "harmonized_snp_level_candidates",
    "exposure_instrument_candidates",
    "files_with_beta_exposure_and_se_exposure",
    "files_with_locus_or_chr_pos",
    "recommended_next_step"
  ),
  value = c(
    n_harmonized,
    n_exposure_instr,
    n_can_compute_F,
    n_locus,
    ifelse(
      n_can_compute_F > 0,
      "Proceed to Phase13.7B to compute F-statistic and harmonization retention from prioritized candidates.",
      "Manually identify the SBP instrument/harmonization input files before Phase13.7B."
    )
  ),
  stringsAsFactors = FALSE
)

recommendation_file <- file.path(out_dir, "phase13_7A_recommendations.tsv")
write_tsv(recommendations, recommendation_file)

# Status
status <- data.frame(
  field = c(
    "phase",
    "candidate_directories_scanned",
    "files_header_scanned",
    "candidate_file_column_scan_created",
    "prioritized_candidates_created",
    "recommendations_created",
    "harmonized_snp_level_candidates",
    "exposure_instrument_candidates",
    "files_with_beta_exposure_and_se_exposure",
    "files_with_locus_or_chr_pos",
    "new_MR_results_created",
    "claim_level",
    "claim_upgrade_allowed",
    "safe_to_proceed_to_phase13_7B"
  ),
  value = c(
    "Phase13.7A",
    length(candidate_dirs),
    length(all_files),
    file.exists(scan_file),
    file.exists(priority_file),
    file.exists(recommendation_file),
    n_harmonized,
    n_exposure_instr,
    n_can_compute_F,
    n_locus,
    FALSE,
    "HYPOTHESIS_GENERATING_NOT_CONFIRMATORY",
    "NO",
    n_can_compute_F > 0
  ),
  stringsAsFactors = FALSE
)

status_file <- file.path(out_dir, "phase13_7A_status.tsv")
write_tsv(status, status_file)

# Update master status as in-progress/discovery complete rather than full Phase13.7 complete.
if (file.exists(master_status_file)) {
  master_status <- read.delim(master_status_file, check.names = FALSE)
  idx <- which(master_status$phase == "Phase13.7")
  if (length(idx) == 1) {
    master_status$status[idx] <- "PHASE13_7A_FILE_DISCOVERY_COMPLETED"
    master_status$qc_status[idx] <- ifelse(n_can_compute_F > 0, "READY_FOR_PHASE13_7B", "NEEDS_MANUAL_FILE_SELECTION")
    master_status$primary_output[idx] <- "phase13_7A_candidate_file_column_scan.tsv; phase13_7A_prioritized_instrument_file_candidates.tsv"
  }
  write_tsv(master_status, master_status_file)
}

message("Phase13.7A completed.")
message("Column scan: ", scan_file)
message("Prioritized candidates: ", priority_file)
message("Recommendations: ", recommendation_file)
message("Status: ", status_file)
