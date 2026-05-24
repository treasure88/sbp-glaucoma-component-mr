options(stringsAsFactors = FALSE)

message("Running Phase13.7B patch: claim-safety audit and jackknife summary extraction...")

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

safe_read <- function(path) {
  if (!file.exists(path)) return(NULL)
  tryCatch(
    read.delim(path, check.names = FALSE, stringsAsFactors = FALSE),
    error = function(e) NULL
  )
}

safe_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

fmt <- function(x, digits = 5) {
  x <- safe_num(x)
  ifelse(is.na(x), NA, formatC(x, digits = digits, format = "fg", flag = "#"))
}

backup_file <- function(path) {
  if (file.exists(path)) {
    bak <- paste0(path, ".pre_phase13_7B_patch.bak")
    file.copy(path, bak, overwrite = TRUE)
  }
}

root <- normalizePath(".", winslash = "/", mustWork = TRUE)

out_dir <- file.path(
  root,
  "71_data_analysis_reinforcement",
  "phase13_7_integrated_robustness_synthesis"
)

matrix_file <- file.path(out_dir, "phase13_7_integrated_robustness_matrix.tsv")
status_file <- file.path(out_dir, "phase13_7_status.tsv")
audit_file <- file.path(out_dir, "phase13_7_claim_safety_audit.tsv")
methods_file <- file.path(out_dir, "phase13_7_manuscript_methods_insert.md")
results_file <- file.path(out_dir, "phase13_7_manuscript_results_insert.md")

backup_file(matrix_file)
backup_file(status_file)
backup_file(audit_file)

# ---------------------------
# Robust value extraction
# ---------------------------

extract_kv_value <- function(df, keys) {
  if (is.null(df) || nrow(df) == 0) return(NA_character_)

  possible_key_cols <- c("metric", "field", "item", "key", "statistic", "parameter", "name")
  possible_value_cols <- c("value", "result", "estimate", "summary")

  key_col <- intersect(possible_key_cols, names(df))
  value_col <- intersect(possible_value_cols, names(df))

  if (length(key_col) > 0 && length(value_col) > 0) {
    kc <- key_col[1]
    vc <- value_col[1]
    key_lower <- tolower(trimws(as.character(df[[kc]])))
    for (k in keys) {
      hit <- which(key_lower == tolower(k))
      if (length(hit) > 0) return(as.character(df[[vc]][hit[1]]))
    }
    for (k in keys) {
      hit <- grep(tolower(k), key_lower, fixed = TRUE)
      if (length(hit) > 0) return(as.character(df[[vc]][hit[1]]))
    }
  }

  for (k in keys) {
    if (k %in% names(df)) return(as.character(df[[k]][1]))
  }

  return(NA_character_)
}

extract_jackknife_summary <- function(path, full_difference = NA_real_) {
  df <- safe_read(path)
  if (is.null(df) || nrow(df) == 0) {
    return(list(preserved = "SOURCE_NOT_READABLE", max_delta = NA_character_))
  }

  preserved <- extract_kv_value(
    df,
    c(
      "all_direction_preserved",
      "all_contrast_direction_preserved",
      "direction_preserved_all",
      "all_positive_direction_preserved",
      "all_direction_stable",
      "direction_preserved"
    )
  )

  max_delta <- extract_kv_value(
    df,
    c(
      "max_abs_delta",
      "maximum_absolute_delta",
      "max_abs_delta_from_full",
      "maximum_absolute_delta_from_full",
      "maximum_delta_from_full"
    )
  )

  # If not key-value style, infer from wide/detail table.
  if (is.na(preserved)) {
    direction_cols <- grep("direction.*preserved|preserved.*direction", names(df), ignore.case = TRUE, value = TRUE)
    if (length(direction_cols) > 0) {
      vals <- toupper(trimws(as.character(df[[direction_cols[1]]])))
      vals <- vals[!is.na(vals) & vals != ""]
      if (length(vals) > 0) {
        preserved <- as.character(all(vals %in% c("TRUE", "YES", "1", "PRESERVED")))
      }
    }
  }

  if (is.na(max_delta)) {
    delta_cols <- grep("delta", names(df), ignore.case = TRUE, value = TRUE)
    if (length(delta_cols) > 0) {
      nums <- safe_num(df[[delta_cols[1]]])
      if (sum(!is.na(nums)) > 0) {
        max_delta <- fmt(max(abs(nums), na.rm = TRUE), 5)
      }
    }
  }

  if (is.na(preserved) || is.na(max_delta)) {
    beta_cols <- grep("beta.*difference|contrast|IOP_minus_nonIOP|difference", names(df), ignore.case = TRUE, value = TRUE)
    beta_cols <- beta_cols[!grepl("se|p$|p_value|q|ci|lower|upper", beta_cols, ignore.case = TRUE)]

    if (length(beta_cols) > 0) {
      nums <- safe_num(df[[beta_cols[1]]])
      nums <- nums[!is.na(nums)]

      if (length(nums) > 0 && is.na(preserved)) {
        preserved <- as.character(all(nums > 0))
      }

      if (length(nums) > 0 && is.na(max_delta) && !is.na(full_difference)) {
        max_delta <- fmt(max(abs(nums - full_difference), na.rm = TRUE), 5)
      }
    }
  }

  if (is.na(preserved)) preserved <- "AVAILABLE_SEE_SOURCE"
  if (is.na(max_delta)) max_delta <- "AVAILABLE_SEE_SOURCE"

  list(preserved = preserved, max_delta = max_delta)
}

# Get full difference from bootstrap if possible.
bootstrap_path <- file.path(
  root,
  "71_data_analysis_reinforcement",
  "phase13_4_empirical_covariance",
  "phase13_4B_bootstrap_contrast_summary.tsv"
)

boot <- safe_read(bootstrap_path)
full_difference <- NA_real_
if (!is.null(boot) && "full_beta_difference_IOP_minus_nonIOP" %in% names(boot)) {
  full_difference <- safe_num(boot$full_beta_difference_IOP_minus_nonIOP[1])
}

chr_path <- file.path(
  root,
  "71_data_analysis_reinforcement",
  "phase13_4_empirical_covariance",
  "phase13_4B_leave_one_chromosome_summary.tsv"
)

locus_path <- file.path(
  root,
  "71_data_analysis_reinforcement",
  "phase13_4_empirical_covariance",
  "phase13_4B_locus_proxy_jackknife_summary.tsv"
)

chr_sum <- extract_jackknife_summary(chr_path, full_difference)
locus_sum <- extract_jackknife_summary(locus_path, full_difference)

# ---------------------------
# Patch integrated matrix
# ---------------------------

matrix <- safe_read(matrix_file)

if (!is.null(matrix) && nrow(matrix) > 0) {
  idx <- which(matrix$evidence_domain == "internal_contrast_influence")

  if (length(idx) == 1) {
    matrix$key_result[idx] <- paste0(
      "Chromosome direction preserved=", chr_sum$preserved,
      "; chromosome max absolute delta=", chr_sum$max_delta,
      "; locus direction preserved=", locus_sum$preserved,
      "; locus max absolute delta=", locus_sum$max_delta,
      "."
    )
  }

  write_tsv(matrix, matrix_file)
}

# ---------------------------
# Correct claim-safety audit
# ---------------------------

created_text_parts <- c()

if (file.exists(methods_file)) {
  created_text_parts <- c(created_text_parts, readLines(methods_file, warn = FALSE))
}
if (file.exists(results_file)) {
  created_text_parts <- c(created_text_parts, readLines(results_file, warn = FALSE))
}
if (!is.null(matrix)) {
  created_text_parts <- c(
    created_text_parts,
    apply(matrix, 1, paste, collapse = " ")
  )
}

created_text <- paste(created_text_parts, collapse = "\n")

count_fixed_hits <- function(phrase, text) {
  hit <- gregexpr(phrase, text, fixed = TRUE)[[1]]
  if (length(hit) == 1 && hit[1] == -1) return(0L)
  length(hit)
}

high_risk <- data.frame(
  phrase = c(
    "SBP causally affects glaucoma",
    "SBP increases IOP-dependent glaucoma risk",
    "SBP protects against NTG",
    "external validation confirmed",
    "validated the mechanism",
    "confirms mediation",
    "confirms pathway",
    "HTG validation was performed",
    "hypertension has no effect",
    "definitive proof"
  ),
  category = c(
    "causal_overclaim",
    "component_causal_overclaim",
    "protective_overclaim",
    "external_confirmation_overclaim",
    "mechanism_validation_overclaim",
    "mediation_overclaim",
    "pathway_overclaim",
    "HTG_overclaim",
    "hypertension_null_overclaim",
    "proof_overclaim"
  ),
  hits = NA_integer_,
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(high_risk))) {
  high_risk$hits[i] <- count_fixed_hits(high_risk$phrase[i], created_text)
}

write_tsv(high_risk, audit_file)

high_risk_hits <- sum(high_risk$hits, na.rm = TRUE)

# ---------------------------
# Patch status
# ---------------------------

status <- safe_read(status_file)

if (!is.null(status) && all(c("field", "value") %in% names(status))) {
  set_status <- function(field, value) {
    idx <- which(status$field == field)
    if (length(idx) == 0) {
      status <<- rbind(
        status,
        data.frame(field = field, value = as.character(value), stringsAsFactors = FALSE)
      )
    } else {
      status$value[idx[1]] <<- as.character(value)
    }
  }

  critical_raw <- status$value[status$field == "critical_sources_found"]
  critical_ok <- FALSE
  if (length(critical_raw) > 0) {
    critical_vals <- trimws(unlist(strsplit(critical_raw[1], ";", fixed = TRUE)))
    critical_ok <- all(critical_vals == "TRUE")
  }

  set_status("jackknife_summary_patch_created", "TRUE")
  set_status("high_risk_phrase_hits", high_risk_hits)
  set_status("new_primary_MR_estimates_created", "FALSE")
  set_status("claim_level", "HYPOTHESIS_GENERATING_NOT_CONFIRMATORY")
  set_status("claim_upgrade_allowed", "NO")
  set_status("phase13_7_passed", critical_ok && high_risk_hits == 0)

  write_tsv(status, status_file)
}

patch_report <- data.frame(
  item = c(
    "false_positive_claim_audit_fixed",
    "corrected_high_risk_phrase_hits",
    "chromosome_jackknife_direction_preserved",
    "chromosome_jackknife_max_abs_delta",
    "locus_jackknife_direction_preserved",
    "locus_jackknife_max_abs_delta",
    "new_primary_MR_estimates_created",
    "claim_upgrade_allowed"
  ),
  value = c(
    "TRUE",
    high_risk_hits,
    chr_sum$preserved,
    chr_sum$max_delta,
    locus_sum$preserved,
    locus_sum$max_delta,
    "FALSE",
    "NO"
  ),
  stringsAsFactors = FALSE
)

write_tsv(
  patch_report,
  file.path(out_dir, "phase13_7B_patch_report.tsv")
)

message("Phase13.7B patch completed.")
message("Corrected high-risk phrase hits: ", high_risk_hits)
message("Chromosome jackknife: direction preserved=", chr_sum$preserved, "; max delta=", chr_sum$max_delta)
message("Locus jackknife: direction preserved=", locus_sum$preserved, "; max delta=", locus_sum$max_delta)
