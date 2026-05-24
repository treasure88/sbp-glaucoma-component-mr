options(stringsAsFactors = FALSE)

script_path <- "../../65_metadata_completion_final_lock/phase10_9_metadata_completion_final_lock.R"
txt <- paste(readLines(script_path, warn = FALSE), collapse = "\n")

start_pat <- "# -----------------------------\n# Final external triangulation table lock\n# -----------------------------"
end_pat <- "# -----------------------------\n# Supplementary numbering final candidate\n# -----------------------------"

start <- regexpr(start_pat, txt, fixed = TRUE)[1]
end <- regexpr(end_pat, txt, fixed = TRUE)[1]

if (start < 0 || end < 0 || end <= start) {
  stop("Could not locate external-lock block boundaries in Phase 10.9 script.")
}

new_block <- paste(c(
"# -----------------------------",
"# Final external triangulation table lock",
"# -----------------------------",
"",
"ext_path <- if (file.exists(paths$external_10_8)) paths$external_10_8 else paths$external_10_5",
"external <- read_tsv(ext_path)",
"",
"locked_n <- data.frame(",
"  outcome = c(\"Measured IOP\", \"HTG\", \"NTG\", \"POAG\", \"RNFL\", \"GCIPL\"),",
"  locked_n_instruments = c(\"456\", \"NA\", \"311\", \"456\", \"456\", \"456\"),",
"  n_instrument_source = c(",
"    \"Phase5.9E IOP external MR result\",",
"    \"Not available; HTG MR not performed\",",
"    \"Phase5.10E NTG external MR result\",",
"    \"Phase5.8E POAG external MR result\",",
"    \"Phase5.6 RNFL external MR result\",",
"    \"Phase5.6 GCIPL external MR result\"",
"  ),",
"  stringsAsFactors = FALSE",
")",
"",
"if (!\"outcome\" %in% names(external)) stop(\"External triangulation table lacks outcome column.\")",
"",
"# Phase 10.8 may already contain lock-helper columns. Drop them before merge",
"# so merge() cannot create .x/.y columns and leave locked_n_instruments undefined.",
"drop_lock_cols <- c(",
"  \"locked_n_instruments\", \"locked_n_instruments.x\", \"locked_n_instruments.y\",",
"  \"n_instrument_source\", \"n_instrument_source.x\", \"n_instrument_source.y\",",
"  \"instrument_count_lock_status\"",
")",
"external <- external[, setdiff(names(external), drop_lock_cols), drop = FALSE]",
"",
"external_final <- merge(locked_n, external, by = \"outcome\", all.x = TRUE, sort = FALSE)",
"external_final <- external_final[match(locked_n$outcome, external_final$outcome), , drop = FALSE]",
"",
"if (!\"n_instruments_original\" %in% names(external_final)) {",
"  if (\"n_instruments\" %in% names(external_final)) {",
"    external_final$n_instruments_original <- as.character(external_final$n_instruments)",
"  } else {",
"    external_final$n_instruments_original <- NA_character_",
"  }",
"}",
"",
"n_orig <- as.character(external_final$n_instruments_original)",
"need_n_original <- is.na(n_orig) |",
"  trimws(n_orig) == \"\" |",
"  trimws(n_orig) == \"NA\" |",
"  grepl(unresolved_regex, n_orig, ignore.case = TRUE)",
"",
"if (any(need_n_original)) {",
"  n_orig[need_n_original] <- as.character(external_final$locked_n_instruments[need_n_original])",
"}",
"",
"external_final$n_instruments_original <- n_orig",
"external_final$n_instruments <- as.character(external_final$locked_n_instruments)",
"",
"external_final$instrument_count_lock_status <- ifelse(",
"  as.character(external_final$n_instruments) == as.character(external_final$locked_n_instruments),",
"  \"PASS_LOCKED\",",
"  \"FAIL_CHECK\"",
")",
"",
"write_tsv(",
"  external_final,",
"  file.path(out_dir, \"phase10_9_external_specificity_triangulation_table_final_locked.tsv\")",
")",
"",
"external_qc <- data.frame(",
"  check = c(\"all_expected_outcomes_present\", \"instrument_counts_locked\", \"no_unresolved_n_instrument_placeholders\"),",
"  status = c(",
"    ifelse(all(locked_n$outcome %in% external_final$outcome), \"PASS\", \"FAIL\"),",
"    ifelse(all(external_final$instrument_count_lock_status == \"PASS_LOCKED\"), \"PASS\", \"FAIL\"),",
"    ifelse(!any(grepl(unresolved_regex, external_final$n_instruments, ignore.case = TRUE)), \"PASS\", \"FAIL\")",
"  ),",
"  stringsAsFactors = FALSE",
")",
"",
"write_tsv(external_qc, file.path(out_dir, \"phase10_9_external_table_qc.tsv\"))",
""
), collapse = "\n")

txt2 <- paste0(substr(txt, 1, start - 1), new_block, substr(txt, end, nchar(txt)))
writeLines(txt2, script_path, useBytes = TRUE)

cat("Patched Phase 10.9 external-lock block:\n")
cat(script_path, "\n")
