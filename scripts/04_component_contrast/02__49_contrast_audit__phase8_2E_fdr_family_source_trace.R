# Phase 8.2E: Trace source of component-contrast FDR family
# Purpose: identify TSV files containing contrast p/q columns and infer the locked FDR family

options(stringsAsFactors = FALSE)

project_root <- "../.."
out_dir <- "../../49_contrast_audit"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

out_file <- file.path(out_dir, "phase8_2E_fdr_family_source_trace.tsv")
candidate_rows_file <- file.path(out_dir, "phase8_2E_fdr_candidate_rows.tsv")
status_file <- file.path(out_dir, "phase8_2E_fdr_family_source_trace_status.tsv")

read_tsv_safe <- function(path) {
  tryCatch(
    read.delim(
      path,
      header = TRUE,
      sep = "\t",
      quote = "",
      comment.char = "",
      check.names = FALSE,
      fill = TRUE
    ),
    error = function(e) NULL
  )
}

to_num <- function(x) suppressWarnings(as.numeric(x))

write_status <- function(status, note) {
  x <- data.frame(
    phase = "Phase8.2E",
    status = status,
    note = note,
    timestamp = as.character(Sys.time()),
    stringsAsFactors = FALSE
  )
  write.table(x, status_file, sep = "\t", quote = FALSE, row.names = FALSE)
}

tsv_files <- list.files(project_root, pattern = "\\.tsv$", recursive = TRUE, full.names = TRUE)

trace_rows <- list()
candidate_rows <- list()

for (f in tsv_files) {
  d <- read_tsv_safe(f)
  if (is.null(d) || nrow(d) == 0) next

  cols <- names(d)

  has_contrast_p <- any(cols %in% c("contrast_p_r0", "contrast_p", "p_contrast", "contrast_pval"))
  has_contrast_q <- any(cols %in% c("contrast_q_r0", "contrast_q", "q_contrast", "contrast_qval"))

  if (has_contrast_p || has_contrast_q) {
    p_col <- intersect(c("contrast_p_r0", "contrast_p", "p_contrast", "contrast_pval"), cols)
    q_col <- intersect(c("contrast_q_r0", "contrast_q", "q_contrast", "contrast_qval"), cols)

    p_col <- ifelse(length(p_col) > 0, p_col[1], NA)
    q_col <- ifelse(length(q_col) > 0, q_col[1], NA)

    trace_rows[[length(trace_rows) + 1]] <- data.frame(
      file_path = f,
      n_rows = nrow(d),
      p_column = p_col,
      q_column = q_col,
      columns = paste(cols, collapse = ";"),
      stringsAsFactors = FALSE
    )

    id_col <- intersect(c("exposure_id", "trait", "trait_id", "exposure", "outcome_id"), cols)
    id_col <- ifelse(length(id_col) > 0, id_col[1], NA)

    p_vals <- if (!is.na(p_col)) to_num(d[[p_col]]) else rep(NA_real_, nrow(d))
    q_vals <- if (!is.na(q_col)) to_num(d[[q_col]]) else rep(NA_real_, nrow(d))
    id_vals <- if (!is.na(id_col)) as.character(d[[id_col]]) else paste0("row_", seq_len(nrow(d)))

    tmp <- data.frame(
      file_path = f,
      row_index = seq_len(nrow(d)),
      id = id_vals,
      p_column = p_col,
      q_column = q_col,
      p_value = p_vals,
      q_value = q_vals,
      stringsAsFactors = FALSE
    )

    tmp <- tmp[!is.na(tmp$p_value) | !is.na(tmp$q_value), , drop = FALSE]

    if (nrow(tmp) > 0) {
      candidate_rows[[length(candidate_rows) + 1]] <- tmp
    }
  }
}

trace <- if (length(trace_rows) > 0) do.call(rbind, trace_rows) else data.frame()
candidates <- if (length(candidate_rows) > 0) do.call(rbind, candidate_rows) else data.frame()

if (nrow(trace) > 0) {
  write.table(trace, out_file, sep = "\t", quote = FALSE, row.names = FALSE)
}

if (nrow(candidates) > 0) {
  candidates <- candidates[order(candidates$file_path, candidates$p_value), , drop = FALSE]
  write.table(candidates, candidate_rows_file, sep = "\t", quote = FALSE, row.names = FALSE)
}

if (nrow(trace) == 0) {
  write_status("NO_CONTRAST_FDR_FILES_FOUND", "No TSV files with recognized contrast p/q columns were found.")
} else {
  note <- paste0(
    "Found ", nrow(trace),
    " TSV files with recognized contrast p/q columns. Review phase8_2E_fdr_candidate_rows.tsv to identify the locked 5-test FDR family."
  )
  write_status("PASSED_TRACE_CREATED", note)
}
