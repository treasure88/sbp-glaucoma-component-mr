options(stringsAsFactors = FALSE)

# Phase13.5A: MVMR feasibility scan
# Project: IOP-dependent vs IOP-independent glaucoma component MR
# Evidence level: HYPOTHESIS_GENERATING_NOT_CONFIRMATORY
# This phase only inventories local data availability for possible SBP+IOP MVMR.
# It does not generate new MR estimates.

message("Running Phase13.5A MVMR feasibility scan...")

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

root <- normalizePath(".", winslash = "/", mustWork = TRUE)
out_dir <- file.path(root, "71_data_analysis_reinforcement", "phase13_5_mvmr_pathway")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

master_status_file <- file.path(root, "71_data_analysis_reinforcement", "phase13_master_status.tsv")

candidate_dirs <- c(
  "05_harmonization_planning",
  "16_mr_input_datasets",
  "28_external_outcome_triangulation_inputs",
  "36_iop_external_validation_inputs",
  "37_iop_external_mr_results",
  "44_master_evidence_table",
  "71_data_analysis_reinforcement"
)

candidate_dirs <- candidate_dirs[dir.exists(file.path(root, candidate_dirs))]

patterns <- "\\.(tsv|txt|csv|tsv\\.gz|txt\\.gz|csv\\.gz|R|r)$"

all_files <- unique(unlist(lapply(candidate_dirs, function(d) {
  list.files(
    file.path(root, d),
    pattern = patterns,
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )
})))

file_size <- suppressWarnings(as.numeric(file.info(all_files)$size))
keep <- !is.na(file_size) & file_size <= 500 * 1024^2
all_files <- all_files[keep]
file_size <- file_size[keep]

open_connection <- function(path) {
  if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    gzfile(path, open = "rt")
  } else {
    file(path, open = "rt")
  }
}

get_header <- function(path) {
  if (grepl("\\.(R|r)$", path, ignore.case = TRUE)) return(NA_character_)
  con <- NULL
  out <- NA_character_
  tryCatch({
    con <- open_connection(path)
    tmp <- readLines(con, n = 1, warn = FALSE)
    if (length(tmp) > 0 && nzchar(tmp[1])) out <- tmp[1]
  }, error = function(e) {
    out <<- NA_character_
  }, finally = {
    if (!is.null(con)) try(close(con), silent = TRUE)
  })
  out
}

split_header <- function(header_line) {
  if (length(header_line) != 1 || is.na(header_line) || !nzchar(header_line)) return(character(0))
  h <- trimws(header_line)
  if (grepl("\t", h, fixed = TRUE)) {
    out <- strsplit(h, "\t", fixed = TRUE)[[1]]
  } else if (grepl("|", h, fixed = TRUE)) {
    out <- strsplit(h, "|", fixed = TRUE)[[1]]
  } else if (grepl(",", h, fixed = TRUE)) {
    out <- strsplit(h, ",", fixed = TRUE)[[1]]
  } else {
    out <- strsplit(h, "\\s+")[[1]]
  }
  trimws(out)
}

has_any <- function(cols_lower, targets) {
  any(cols_lower %in% tolower(targets))
}

has_regex <- function(cols_lower, pat) {
  any(grepl(pat, cols_lower, ignore.case = TRUE))
}

classify_mvmr_file <- function(rel, cols) {
  rel_low <- tolower(rel)
  cl <- tolower(cols)
  
  has_snp <- has_any(cl, c("snp", "rsid", "variant", "variant_id", "match_key"))
  has_beta <- has_any(cl, c("beta", "b", "beta_exposure", "beta_outcome", "beta_outcome_aligned")) ||
    has_regex(cl, "^beta")
  has_se <- has_any(cl, c("se", "stderr", "se_exposure", "se_outcome"))
  has_p <- has_any(cl, c("p", "pval", "pvalue", "pval_exposure", "pval_outcome"))
  has_chr_pos <- has_any(cl, c("chr", "chromosome", "exposure_chr")) &&
    has_any(cl, c("pos", "position", "bp", "exposure_pos"))
  
  has_sbp_path <- grepl("sbp", rel_low)
  has_iop_path <- grepl("iopcc|iop", rel_low)
  has_component_path <- grepl("gbs|noniop|iopcomponent|component", rel_low)
  has_index_path <- grepl("index|by_rsid|by_chrpos|outcome_indexes|exposure", rel_low)
  has_pairwise_path <- grepl("pairwise|mr_input|harmonized", rel_low)
  
  file_class <- "other"
  priority_score <- 0
  
  if (has_sbp_path && has_component_path && has_pairwise_path && has_snp && has_beta && has_se) {
    file_class <- "SBP_component_harmonized_input"
    priority_score <- priority_score + 10
  }
  
  if (has_sbp_path && has_iop_path && has_pairwise_path && has_snp && has_beta && has_se) {
    file_class <- "SBP_instruments_with_IOP_association"
    priority_score <- priority_score + 10
  }
  
  if (has_iop_path && has_index_path && has_snp && has_beta && has_se) {
    file_class <- "IOP_association_index_candidate"
    priority_score <- priority_score + 9
  }
  
  if (has_sbp_path && has_index_path && has_snp && has_beta && has_se) {
    file_class <- "SBP_association_index_candidate"
    priority_score <- priority_score + 8
  }
  
  if (has_component_path && has_index_path && has_snp && has_beta && has_se) {
    file_class <- "component_outcome_index_candidate"
    priority_score <- priority_score + 8
  }
  
  if (has_iop_path && has_pairwise_path && has_snp && has_beta && has_se && !has_sbp_path) {
    file_class <- "IOP_pairwise_or_instrument_candidate"
    priority_score <- priority_score + 7
  }
  
  if (has_sbp_path) priority_score <- priority_score + 2
  if (has_iop_path) priority_score <- priority_score + 2
  if (has_component_path) priority_score <- priority_score + 2
  if (has_snp) priority_score <- priority_score + 1
  if (has_beta && has_se) priority_score <- priority_score + 1
  if (has_chr_pos) priority_score <- priority_score + 1
  
  list(
    file_class = file_class,
    priority_score = priority_score,
    has_snp = has_snp,
    has_beta = has_beta,
    has_se = has_se,
    has_p = has_p,
    has_chr_pos = has_chr_pos,
    has_sbp_path = has_sbp_path,
    has_iop_path = has_iop_path,
    has_component_path = has_component_path,
    has_index_path = has_index_path,
    has_pairwise_path = has_pairwise_path
  )
}

rows <- list()

for (i in seq_along(all_files)) {
  f <- all_files[i]
  rel <- sub(paste0("^", gsub("\\\\", "/", root), "/?"), "", normalizePath(f, winslash = "/", mustWork = FALSE))
  header <- get_header(f)
  cols <- split_header(header)
  cls <- classify_mvmr_file(rel, cols)
  
  rows[[i]] <- data.frame(
    relative_path = rel,
    size_bytes = file_size[i],
    n_columns = length(cols),
    file_class = cls$file_class,
    priority_score = cls$priority_score,
    has_snp = cls$has_snp,
    has_beta = cls$has_beta,
    has_se = cls$has_se,
    has_p = cls$has_p,
    has_chr_pos = cls$has_chr_pos,
    has_sbp_path = cls$has_sbp_path,
    has_iop_path = cls$has_iop_path,
    has_component_path = cls$has_component_path,
    has_index_path = cls$has_index_path,
    has_pairwise_path = cls$has_pairwise_path,
    columns_preview = paste(head(cols, 80), collapse = "|"),
    stringsAsFactors = FALSE
  )
}

scan <- do.call(rbind, rows)
scan <- scan[order(-scan$priority_score, scan$relative_path), ]

scan_file <- file.path(out_dir, "phase13_5A_mvmr_feasibility_file_scan.tsv")
write_tsv(scan, scan_file)

priority <- scan[scan$priority_score > 0, ]
priority_file <- file.path(out_dir, "phase13_5A_mvmr_priority_candidates.tsv")
write_tsv(priority, priority_file)

n_sbp_component_inputs <- sum(priority$file_class == "SBP_component_harmonized_input", na.rm = TRUE)
n_sbp_iop_assoc <- sum(priority$file_class == "SBP_instruments_with_IOP_association", na.rm = TRUE)
n_iop_index <- sum(priority$file_class == "IOP_association_index_candidate", na.rm = TRUE)
n_sbp_index <- sum(priority$file_class == "SBP_association_index_candidate", na.rm = TRUE)
n_component_index <- sum(priority$file_class == "component_outcome_index_candidate", na.rm = TRUE)
n_iop_pairwise <- sum(priority$file_class == "IOP_pairwise_or_instrument_candidate", na.rm = TRUE)

feasibility <- data.frame(
  analysis_option = c(
    "SBP-instrument IOP-adjusted sensitivity",
    "Full union-instrument SBP+IOP MVMR",
    "No MVMR; narrative limitation only"
  ),
  required_inputs = c(
    "SBP-component harmonized inputs plus IOP association for the same SBP instruments",
    "SBP and IOP exposure association indexes plus component outcome indexes for the union of instruments",
    "Triggered only if IOP association data cannot be paired to SBP instruments"
  ),
  local_evidence_found = c(
    paste0("SBP_component_inputs=", n_sbp_component_inputs, "; SBP_instruments_with_IOP_association=", n_sbp_iop_assoc, "; IOP_index=", n_iop_index),
    paste0("SBP_index=", n_sbp_index, "; IOP_index=", n_iop_index, "; component_index=", n_component_index, "; IOP_pairwise=", n_iop_pairwise),
    "Use only if both above options are infeasible"
  ),
  feasibility_status = c(
    ifelse(n_sbp_component_inputs >= 2 && (n_sbp_iop_assoc >= 1 || n_iop_index >= 1), "POTENTIALLY_FEASIBLE", "NOT_YET_FEASIBLE"),
    ifelse(n_sbp_index >= 1 && n_iop_index >= 1 && n_component_index >= 2, "POTENTIALLY_FEASIBLE", "NOT_YET_FEASIBLE"),
    ifelse(n_sbp_component_inputs >= 2 && (n_sbp_iop_assoc >= 1 || n_iop_index >= 1), "NOT_NEEDED_NOW", "POSSIBLE_FALLBACK")
  ),
  claim_boundary = c(
    "If performed, label as exploratory pathway-consistency sensitivity, not mediation proof.",
    "If performed, label as exploratory MVMR sensitivity, not pathway confirmation.",
    "If not feasible, report as data-availability limitation."
  ),
  stringsAsFactors = FALSE
)

feasibility_file <- file.path(out_dir, "phase13_5A_mvmr_feasibility_decision_table.tsv")
write_tsv(feasibility, feasibility_file)

recommendation <- if (feasibility$feasibility_status[1] == "POTENTIALLY_FEASIBLE") {
  "Proceed to Phase13.5B: build SBP-instrument IOP-adjusted sensitivity dataset."
} else if (feasibility$feasibility_status[2] == "POTENTIALLY_FEASIBLE") {
  "Proceed to Phase13.5B-full: build union-instrument SBP+IOP MVMR dataset."
} else {
  "Do not run MVMR yet; manually recover IOP association data or report as limitation."
}

recommendation_table <- data.frame(
  item = c(
    "SBP_component_harmonized_inputs",
    "SBP_instruments_with_IOP_association",
    "IOP_association_index_candidates",
    "SBP_association_index_candidates",
    "component_outcome_index_candidates",
    "IOP_pairwise_or_instrument_candidates",
    "recommendation"
  ),
  value = c(
    n_sbp_component_inputs,
    n_sbp_iop_assoc,
    n_iop_index,
    n_sbp_index,
    n_component_index,
    n_iop_pairwise,
    recommendation
  ),
  stringsAsFactors = FALSE
)

recommendation_file <- file.path(out_dir, "phase13_5A_recommendation.tsv")
write_tsv(recommendation_table, recommendation_file)

status <- data.frame(
  field = c(
    "phase",
    "files_scanned",
    "priority_candidates_created",
    "feasibility_decision_table_created",
    "SBP_component_harmonized_inputs",
    "SBP_instruments_with_IOP_association",
    "IOP_association_index_candidates",
    "SBP_association_index_candidates",
    "component_outcome_index_candidates",
    "IOP_pairwise_or_instrument_candidates",
    "recommended_next_step",
    "new_MR_results_created",
    "claim_level",
    "claim_upgrade_allowed",
    "safe_to_proceed_to_phase13_5B"
  ),
  value = c(
    "Phase13.5A",
    nrow(scan),
    file.exists(priority_file),
    file.exists(feasibility_file),
    n_sbp_component_inputs,
    n_sbp_iop_assoc,
    n_iop_index,
    n_sbp_index,
    n_component_index,
    n_iop_pairwise,
    recommendation,
    FALSE,
    "HYPOTHESIS_GENERATING_NOT_CONFIRMATORY",
    "NO",
    feasibility$feasibility_status[1] == "POTENTIALLY_FEASIBLE" ||
      feasibility$feasibility_status[2] == "POTENTIALLY_FEASIBLE"
  ),
  stringsAsFactors = FALSE
)

status_file <- file.path(out_dir, "phase13_5A_status.tsv")
write_tsv(status, status_file)

if (file.exists(master_status_file)) {
  master_status <- read.delim(master_status_file, check.names = FALSE)
  idx <- which(master_status$phase == "Phase13.5")
  if (length(idx) == 1) {
    master_status$status[idx] <- "PHASE13_5A_MVMR_FEASIBILITY_SCAN_COMPLETED"
    master_status$qc_status[idx] <- ifelse(
      feasibility$feasibility_status[1] == "POTENTIALLY_FEASIBLE" ||
        feasibility$feasibility_status[2] == "POTENTIALLY_FEASIBLE",
      "READY_FOR_PHASE13_5B",
      "NEEDS_DATA_RECOVERY_OR_LIMITATION"
    )
    master_status$primary_output[idx] <- "phase13_5A_mvmr_feasibility_decision_table.tsv; phase13_5A_mvmr_priority_candidates.tsv"
  }
  write_tsv(master_status, master_status_file)
}

message("Phase13.5A completed.")
message("Priority candidates: ", priority_file)
message("Feasibility decision: ", feasibility_file)
message("Recommendation: ", recommendation_file)
message("Status: ", status_file)
message("Recommendation: ", recommendation)
