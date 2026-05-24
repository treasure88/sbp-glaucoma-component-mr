options(stringsAsFactors = FALSE)

# Phase13.6A: Robustness and pleiotropy source inventory
# Project: IOP-dependent vs IOP-independent glaucoma component MR
# Evidence level: HYPOTHESIS_GENERATING_NOT_CONFIRMATORY
# This phase discovers existing robustness/pleiotropy outputs.
# It does not create new MR estimates.

message("Running Phase13.6A robustness source inventory...")

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
out_dir <- file.path(root, "71_data_analysis_reinforcement", "phase13_6_pleiotropy_robustness")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

master_status_file <- file.path(root, "71_data_analysis_reinforcement", "phase13_master_status.tsv")

candidate_dirs <- c(
  "05_harmonization_planning",
  "16_mr_input_datasets",
  "22_sbp_robustness",
  "24_vascular_panel_integration",
  "29_external_neuroretinal_mr_results",
  "33_poag_external_mr_results",
  "37_iop_external_mr_results",
  "41_ntg_external_mr_results",
  "44_master_evidence_table",
  "46_figure_table_preparation",
  "66_final_submission_package",
  "71_data_analysis_reinforcement"
)

candidate_dirs <- candidate_dirs[dir.exists(file.path(root, candidate_dirs))]

patterns <- "\\.(tsv|txt|csv|tsv\\.gz|txt\\.gz|csv\\.gz|R|r|Rmd|rmd)$"

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
  if (grepl("\\.(R|r|Rmd|rmd)$", path)) return(NA_character_)
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

classify_robustness <- function(rel, cols) {
  rel_low <- tolower(rel)
  cols_low <- tolower(cols)
  
  has_sbp <- grepl("sbp", rel_low) || any(grepl("sbp", cols_low))
  has_component <- grepl("gbs|noniop|iopcomponent|component", rel_low) || any(grepl("gbs|noniop|iopcomponent|component", cols_low))
  
  category <- "other"
  method_family <- NA_character_
  priority <- 0
  
  if (grepl("mre?gger|egger", rel_low) || any(grepl("egger", cols_low))) {
    category <- "pleiotropy_method_result"
    method_family <- "MR-Egger"
    priority <- priority + 8
  }
  
  if (grepl("weighted.*median|median", rel_low) || any(grepl("weighted.*median|median", cols_low))) {
    category <- "robust_method_result"
    method_family <- "weighted_median"
    priority <- priority + 7
  }
  
  if (grepl("weighted.*mode|mode", rel_low) || any(grepl("weighted.*mode|mode", cols_low))) {
    category <- "robust_method_result"
    method_family <- "weighted_mode"
    priority <- priority + 7
  }
  
  if (grepl("mrpresso|mr-presso|presso", rel_low) || any(grepl("mrpresso|mr-presso|presso", cols_low))) {
    category <- "pleiotropy_outlier_result"
    method_family <- "MR-PRESSO"
    priority <- priority + 9
  }
  
  if (grepl("mrraps|mr-raps|raps", rel_low) || any(grepl("mrraps|mr-raps|raps", cols_low))) {
    category <- "robust_method_result"
    method_family <- "MR-RAPS"
    priority <- priority + 8
  }
  
  if (grepl("radial", rel_low) || any(grepl("radial", cols_low))) {
    category <- "radial_outlier_result"
    method_family <- "radial_MR"
    priority <- priority + 8
  }
  
  if (grepl("leave.*one|single_snp|influence", rel_low) ||
      any(grepl("left_out|leave_one|single_snp|delta_from_full", cols_low))) {
    category <- "influence_result"
    method_family <- "leave_one_or_single_snp"
    priority <- priority + 7
  }
  
  if (grepl("steiger", rel_low) || any(grepl("steiger", cols_low))) {
    category <- "directionality_result"
    method_family <- "Steiger"
    priority <- priority + 7
  }
  
  if (grepl("heterogeneity|cochran|\\bQ\\b", rel_low) ||
      any(grepl("heterogeneity|cochran|q_pval|q_df|phi", cols_low))) {
    if (category == "other") category <- "heterogeneity_result"
    if (is.na(method_family)) method_family <- "heterogeneity"
    priority <- priority + 5
  }
  
  if (grepl("method|ivw|weighted|egger|median|mode", rel_low) ||
      any(cols_low %in% c("method", "beta", "se", "pval", "p"))) {
    if (category == "other") category <- "method_summary_candidate"
    if (is.na(method_family)) method_family <- "method_summary"
    priority <- priority + 3
  }
  
  if (has_sbp) priority <- priority + 3
  if (has_component) priority <- priority + 3
  
  list(
    category = category,
    method_family = method_family,
    priority_score = priority,
    has_sbp = has_sbp,
    has_component = has_component
  )
}

rows <- list()

for (i in seq_along(all_files)) {
  f <- all_files[i]
  rel <- sub(paste0("^", gsub("\\\\", "/", root), "/?"), "", normalizePath(f, winslash = "/", mustWork = FALSE))
  header <- get_header(f)
  cols <- split_header(header)
  cls <- classify_robustness(rel, cols)
  
  rows[[i]] <- data.frame(
    relative_path = rel,
    size_bytes = file_size[i],
    extension = sub("^.*\\.", ".", basename(rel)),
    n_columns = length(cols),
    category = cls$category,
    method_family = cls$method_family,
    priority_score = cls$priority_score,
    has_sbp = cls$has_sbp,
    has_component = cls$has_component,
    columns_preview = paste(head(cols, 80), collapse = "|"),
    stringsAsFactors = FALSE
  )
}

inventory <- do.call(rbind, rows)
inventory <- inventory[order(-inventory$priority_score, inventory$relative_path), ]

inventory_file <- file.path(out_dir, "phase13_6A_robustness_source_inventory.tsv")
write_tsv(inventory, inventory_file)

priority <- inventory[inventory$priority_score > 0, ]
priority_file <- file.path(out_dir, "phase13_6A_prioritized_robustness_sources.tsv")
write_tsv(priority, priority_file)

# Summarize method coverage.
coverage <- data.frame(
  method_family = c(
    "MR-Egger",
    "weighted_median",
    "weighted_mode",
    "MR-PRESSO",
    "MR-RAPS",
    "radial_MR",
    "leave_one_or_single_snp",
    "Steiger",
    "heterogeneity",
    "method_summary"
  ),
  n_candidate_files = NA_integer_,
  n_SBP_component_candidate_files = NA_integer_,
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(coverage))) {
  mf <- coverage$method_family[i]
  coverage$n_candidate_files[i] <- sum(priority$method_family == mf, na.rm = TRUE)
  coverage$n_SBP_component_candidate_files[i] <- sum(
    priority$method_family == mf &
      priority$has_sbp &
      priority$has_component,
    na.rm = TRUE
  )
}

coverage_file <- file.path(out_dir, "phase13_6A_method_coverage_summary.tsv")
write_tsv(coverage, coverage_file)

# Script trace for robustness generation.
script_files <- all_files[grepl("\\.(R|r|Rmd|rmd)$", all_files)]
terms <- c(
  "mr_egger",
  "egger",
  "weighted_median",
  "weighted_mode",
  "mrpresso",
  "MRPRESSO",
  "mrraps",
  "radial",
  "leave_one",
  "single_snp",
  "steiger",
  "heterogeneity",
  "phase4_2",
  "SBP__GBS_nonIOPcomponent",
  "SBP__GBS_IOPcomponent"
)

script_hits <- list()
k <- 1

for (s in script_files) {
  rel <- sub(paste0("^", gsub("\\\\", "/", root), "/?"), "", normalizePath(s, winslash = "/", mustWork = FALSE))
  txt <- tryCatch(readLines(s, warn = FALSE, encoding = "UTF-8"), error = function(e) character(0))
  if (length(txt) == 0) next
  
  for (term in terms) {
    idx <- grep(term, txt, ignore.case = TRUE, fixed = TRUE)
    if (length(idx) > 0) {
      for (ii in head(idx, 30)) {
        script_hits[[k]] <- data.frame(
          script_path = rel,
          search_term = term,
          line_number = ii,
          line_text = trimws(txt[ii]),
          stringsAsFactors = FALSE
        )
        k <- k + 1
      }
    }
  }
}

script_hits_df <- if (length(script_hits) > 0) do.call(rbind, script_hits) else data.frame(
  script_path = character(0),
  search_term = character(0),
  line_number = integer(0),
  line_text = character(0)
)

script_hits_file <- file.path(out_dir, "phase13_6A_robustness_script_trace_hits.tsv")
write_tsv(script_hits_df, script_hits_file)

status <- data.frame(
  field = c(
    "phase",
    "files_scanned",
    "prioritized_sources_created",
    "method_coverage_created",
    "script_trace_hits_created",
    "MR_PRESSO_SBP_component_candidate_files",
    "MR_RAPS_SBP_component_candidate_files",
    "radial_SBP_component_candidate_files",
    "leave_one_SBP_component_candidate_files",
    "egger_SBP_component_candidate_files",
    "new_MR_results_created",
    "claim_level",
    "claim_upgrade_allowed",
    "safe_to_proceed_to_phase13_6B_matrix"
  ),
  value = c(
    "Phase13.6A",
    nrow(inventory),
    file.exists(priority_file),
    file.exists(coverage_file),
    file.exists(script_hits_file),
    coverage$n_SBP_component_candidate_files[coverage$method_family == "MR-PRESSO"],
    coverage$n_SBP_component_candidate_files[coverage$method_family == "MR-RAPS"],
    coverage$n_SBP_component_candidate_files[coverage$method_family == "radial_MR"],
    coverage$n_SBP_component_candidate_files[coverage$method_family == "leave_one_or_single_snp"],
    coverage$n_SBP_component_candidate_files[coverage$method_family == "MR-Egger"],
    FALSE,
    "HYPOTHESIS_GENERATING_NOT_CONFIRMATORY",
    "NO",
    sum(coverage$n_SBP_component_candidate_files, na.rm = TRUE) > 0
  ),
  stringsAsFactors = FALSE
)

status_file <- file.path(out_dir, "phase13_6A_status.tsv")
write_tsv(status, status_file)

if (file.exists(master_status_file)) {
  master_status <- read.delim(master_status_file, check.names = FALSE)
  idx <- which(master_status$phase == "Phase13.6")
  if (length(idx) == 1) {
    master_status$status[idx] <- "PHASE13_6A_ROBUSTNESS_SOURCE_INVENTORY_COMPLETED"
    master_status$qc_status[idx] <- ifelse(
      sum(coverage$n_SBP_component_candidate_files, na.rm = TRUE) > 0,
      "READY_FOR_PHASE13_6B",
      "NEEDS_MANUAL_REVIEW"
    )
    master_status$primary_output[idx] <- "phase13_6A_prioritized_robustness_sources.tsv; phase13_6A_method_coverage_summary.tsv"
  }
  write_tsv(master_status, master_status_file)
}

message("Phase13.6A completed.")
message("Priority sources: ", priority_file)
message("Coverage summary: ", coverage_file)
message("Script trace hits: ", script_hits_file)
message("Status: ", status_file)
