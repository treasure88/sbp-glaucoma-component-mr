# Phase 8.3A: Find candidate SBP SNP-level harmonized files
# Purpose: locate files suitable for leave-one-chromosome-out or alternative instrument sensitivity

options(stringsAsFactors = FALSE)

project_root <- "../.."
out_dir <- "../../50_sbp_sensitivity"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

out_file <- file.path(out_dir, "phase8_3A_sbp_snp_level_input_inventory.tsv")
status_file <- file.path(out_dir, "phase8_3A_sbp_snp_level_input_inventory_status.tsv")

write_status <- function(status, note) {
  x <- data.frame(
    phase = "Phase8.3A",
    status = status,
    note = note,
    timestamp = as.character(Sys.time()),
    stringsAsFactors = FALSE
  )
  write.table(x, status_file, sep = "\t", quote = FALSE, row.names = FALSE)
}

read_tsv_safe <- function(path, n_max = 50) {
  tryCatch(
    read.delim(
      path,
      header = TRUE,
      sep = "\t",
      quote = "",
      comment.char = "",
      check.names = FALSE,
      fill = TRUE,
      nrows = n_max
    ),
    error = function(e) NULL
  )
}

files <- list.files(
  project_root,
  pattern = "\\.(tsv|txt|csv)(\\.gz)?$",
  recursive = TRUE,
  full.names = TRUE
)

candidate_rows <- list()

for (f in files) {
  # Skip huge raw external outcomes unless filename suggests harmonized or SBP/component
  fname <- basename(f)
  lower_path <- tolower(f)

  likely_relevant <- grepl("sbp|harmon|component|gbs|noniop|iop|mr_input|instrument|allele", lower_path)

  if (!likely_relevant) next

  d <- read_tsv_safe(f, n_max = 30)
  if (is.null(d) || ncol(d) == 0) next

  cols <- names(d)
  cols_lower <- tolower(cols)

  has_snp <- any(cols_lower %in% c("snp", "rsid", "rs_id", "variant", "markername"))
  has_beta_exp <- any(grepl("beta.*exposure|beta_exposure|beta\\.exposure|exposure.*beta|bx|beta_sbp", cols_lower))
  has_se_exp <- any(grepl("se.*exposure|se_exposure|se\\.exposure|exposure.*se|se_sbp|sebx", cols_lower))
  has_beta_out <- any(grepl("beta.*outcome|beta_outcome|beta\\.outcome|outcome.*beta|by", cols_lower))
  has_se_out <- any(grepl("se.*outcome|se_outcome|se\\.outcome|outcome.*se|seby", cols_lower))
  has_chr <- any(cols_lower %in% c("chr", "chrom", "chromosome", "#chrom"))
  has_pos <- any(cols_lower %in% c("pos", "position", "bp", "base_pair_location"))

  contains_sbp_text <- any(grepl("SBP", as.matrix(d), ignore.case = FALSE), na.rm = TRUE)
  contains_iop_text <- any(grepl("IOP|nonIOP|GBS", as.matrix(d), ignore.case = FALSE), na.rm = TRUE)

  score <- sum(c(
    has_snp,
    has_beta_exp,
    has_se_exp,
    has_beta_out,
    has_se_out,
    has_chr,
    has_pos,
    contains_sbp_text,
    contains_iop_text
  ))

  if (score >= 3) {
    candidate_rows[[length(candidate_rows) + 1]] <- data.frame(
      file_path = f,
      file_name = fname,
      n_preview_rows = nrow(d),
      n_columns = ncol(d),
      has_snp = has_snp,
      has_beta_exposure_like = has_beta_exp,
      has_se_exposure_like = has_se_exp,
      has_beta_outcome_like = has_beta_out,
      has_se_outcome_like = has_se_out,
      has_chr = has_chr,
      has_pos = has_pos,
      contains_sbp_text = contains_sbp_text,
      contains_iop_or_component_text = contains_iop_text,
      suitability_score = score,
      columns = paste(cols, collapse = ";"),
      stringsAsFactors = FALSE
    )
  }
}

inventory <- if (length(candidate_rows) > 0) do.call(rbind, candidate_rows) else data.frame()

if (nrow(inventory) > 0) {
  inventory <- inventory[order(-inventory$suitability_score, inventory$file_path), , drop = FALSE]
  write.table(inventory, out_file, sep = "\t", quote = FALSE, row.names = FALSE)
  write_status(
    "PASSED_INVENTORY_CREATED",
    paste0("Found ", nrow(inventory), " candidate SNP-level or harmonized input files. Review inventory before Phase 8.3B.")
  )
} else {
  write.table(inventory, out_file, sep = "\t", quote = FALSE, row.names = FALSE)
  write_status(
    "NO_CANDIDATE_FILES_FOUND",
    "No likely SBP SNP-level harmonized input files found by heuristic scan."
  )
}
