options(stringsAsFactors = FALSE)

message("Running Phase13.8A measured IOP instrument-strength label patch...")

write_tsv <- function(x, file) {
  write.table(x, file = file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE, na = "")
}

root <- normalizePath(".", winslash = "/", mustWork = TRUE)
out_dir <- file.path(root, "71_data_analysis_reinforcement", "phase13_8_external_triangulation_grid")

grid_file <- file.path(out_dir, "phase13_8_external_triangulation_interpretation_grid.tsv")
status_file <- file.path(out_dir, "phase13_8_status.tsv")
strength_file <- file.path(root, "71_data_analysis_reinforcement/phase13_7_instrument_transparency/phase13_7B_SBP_external_instrument_strength_summary.tsv")

stopifnot(file.exists(grid_file))
stopifnot(file.exists(strength_file))

grid <- read.delim(grid_file, check.names = FALSE)
strength <- read.delim(strength_file, check.names = FALSE)

backup_file <- paste0(grid_file, ".before_phase13_8A_patch")
write_tsv(grid, backup_file)

# Map external grid endpoint labels to strength-table outcome labels.
endpoint_to_strength_outcome <- c(
  "Measured IOP" = "IOP",
  "POAG" = "POAG",
  "NTG" = "NTG",
  "RNFL" = "RNFL",
  "GCIPL" = "GCIPL"
)

for (i in seq_len(nrow(grid))) {
  endpoint <- as.character(grid$endpoint_or_analysis[i])
  strength_outcome <- endpoint_to_strength_outcome[[endpoint]]
  if (is.null(strength_outcome)) next
  
  srow <- strength[strength$outcome == strength_outcome & strength$exposure == "SBP", , drop = FALSE]
  if (nrow(srow) < 1) next
  
  grid$instrument_strength_summary[i] <- paste0(
    "mean_F=", srow$mean_F[1],
    "; min_F=", srow$min_F[1],
    "; F_less_10=", srow$n_F_less_10[1]
  )
  grid$retention_summary[i] <- paste0("retention_rate=", srow$file_level_retention_rate[1])
}

write_tsv(grid, grid_file)

measured_iop_row <- grid[grid$endpoint_or_analysis == "Measured IOP", , drop = FALSE]
measured_iop_strength_patched <- nrow(measured_iop_row) == 1 &&
  !grepl("mean_F=NA", measured_iop_row$instrument_strength_summary[1], fixed = TRUE)

patch_summary <- data.frame(
  item = c(
    "patch",
    "grid_backup_created",
    "measured_iop_strength_patched",
    "measured_iop_instrument_strength_summary",
    "measured_iop_retention_summary",
    "new_MR_results_created",
    "claim_level",
    "claim_upgrade_allowed"
  ),
  value = c(
    "Phase13.8A_measured_IOP_strength_label_patch",
    file.exists(backup_file),
    measured_iop_strength_patched,
    ifelse(nrow(measured_iop_row) == 1, measured_iop_row$instrument_strength_summary[1], NA),
    ifelse(nrow(measured_iop_row) == 1, measured_iop_row$retention_summary[1], NA),
    FALSE,
    "HYPOTHESIS_GENERATING_NOT_CONFIRMATORY",
    "NO"
  ),
  stringsAsFactors = FALSE
)

patch_summary_file <- file.path(out_dir, "phase13_8A_patch_measured_iop_strength_summary.tsv")
write_tsv(patch_summary, patch_summary_file)

# Patch status in place.
if (file.exists(status_file)) {
  status <- read.delim(status_file, check.names = FALSE)
  add_or_set <- function(df, field, value) {
    idx <- which(df$field == field)
    if (length(idx) == 0) {
      df <- rbind(df, data.frame(field = field, value = as.character(value), stringsAsFactors = FALSE))
    } else {
      df$value[idx[1]] <- as.character(value)
    }
    df
  }
  status <- add_or_set(status, "patch", "Phase13.8A_measured_IOP_strength_label_patch")
  status <- add_or_set(status, "measured_iop_strength_patched", measured_iop_strength_patched)
  status <- add_or_set(status, "phase13_8_passed", measured_iop_strength_patched)
  write_tsv(status, status_file)
}

message("Phase13.8A patch completed.")
message("Measured IOP strength patched: ", measured_iop_strength_patched)
message("Patch summary: ", patch_summary_file)
