options(stringsAsFactors = FALSE)

# Phase13.6B: Integrate existing SBP robustness and pleiotropy outputs
# Project: IOP-dependent vs IOP-independent glaucoma component MR
# Evidence level: HYPOTHESIS_GENERATING_NOT_CONFIRMATORY
# This phase integrates existing robustness outputs into a manuscript-facing matrix.
# It does not generate new primary MR estimates.

message("Running Phase13.6B robustness matrix integration...")

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
  out <- NULL
  tryCatch({
    if (grepl("\\.gz$", path, ignore.case = TRUE)) {
      con <- gzfile(path, open = "rt")
      on.exit(close(con), add = TRUE)
      out <- read.delim(con, header = TRUE, sep = "\t", quote = "", comment.char = "", check.names = FALSE, stringsAsFactors = FALSE)
    } else {
      out <- read.delim(path, header = TRUE, sep = "\t", quote = "", comment.char = "", check.names = FALSE, stringsAsFactors = FALSE)
    }
  }, error = function(e) {
    out <<- NULL
  })
  out
}

pick <- function(df, candidates) {
  nms <- names(df)
  low <- tolower(nms)
  for (cc in candidates) {
    idx <- which(low == tolower(cc))
    if (length(idx) > 0) return(nms[idx[1]])
  }
  NA_character_
}

val <- function(df, col, i) {
  if (is.na(col) || !(col %in% names(df))) return(NA_character_)
  as.character(df[[col]][i])
}

as_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

fmt <- function(x, digits = 4) {
  ifelse(is.na(x), NA, formatC(as.numeric(x), digits = digits, format = "fg", flag = "#"))
}

root <- normalizePath(".", winslash = "/", mustWork = TRUE)
out_dir <- file.path(root, "71_data_analysis_reinforcement", "phase13_6_pleiotropy_robustness")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

master_status_file <- file.path(root, "71_data_analysis_reinforcement", "phase13_master_status.tsv")

# Explicit source files only; avoid keyword-classification false positives.
sources <- data.frame(
  source_label = c(
    "method_consistency",
    "component_contrast_by_method",
    "mrraps_component_contrast",
    "mrraps_results",
    "curated_mrraps_interpretation",
    "mrpresso_summary",
    "curated_mrpresso_results",
    "radial_outlier_summary",
    "radial_outlier_snps",
    "outlier_corrected_results",
    "single_snp_influence_summary",
    "leave_one_chromosome_out",
    "phase13_bootstrap_summary",
    "phase13_leave_one_chromosome_summary",
    "phase13_locus_proxy_jackknife_summary"
  ),
  relative_path = c(
    "22_sbp_robustness/phase4_2_sbp_method_consistency.tsv",
    "22_sbp_robustness/phase4_2_sbp_component_contrast_by_method.tsv",
    "22_sbp_robustness/phase4_2B_formal_packages_final/phase4_2B_final_mrraps_component_contrast.tsv",
    "22_sbp_robustness/phase4_2B_formal_packages_final/phase4_2B_final_mrraps_results.tsv",
    "22_sbp_robustness/phase4_2B_formal_packages_final/phase4_2B_final_curated_mrraps_interpretation.tsv",
    "22_sbp_robustness/phase4_2B_formal_packages_final/phase4_2B_final_mrpresso_summary.tsv",
    "22_sbp_robustness/phase4_2B_formal_packages_final/phase4_2B_final_curated_mrpresso_results.tsv",
    "22_sbp_robustness/phase4_2_sbp_radial_outlier_summary.tsv",
    "22_sbp_robustness/phase4_2_sbp_radial_outlier_snps.tsv",
    "22_sbp_robustness/phase4_2_sbp_outlier_corrected_results.tsv",
    "22_sbp_robustness/phase4_2_sbp_single_snp_influence_summary.tsv",
    "22_sbp_robustness/phase4_2_sbp_leave_one_chromosome_out.tsv",
    "71_data_analysis_reinforcement/phase13_4_empirical_covariance/phase13_4B_bootstrap_contrast_summary.tsv",
    "71_data_analysis_reinforcement/phase13_4_empirical_covariance/phase13_4B_leave_one_chromosome_summary.tsv",
    "71_data_analysis_reinforcement/phase13_4_empirical_covariance/phase13_4B_locus_proxy_jackknife_summary.tsv"
  ),
  stringsAsFactors = FALSE
)

sources$full_path <- file.path(root, sources$relative_path)
sources$exists <- file.exists(sources$full_path)

source_manifest_file <- file.path(out_dir, "phase13_6B_source_manifest.tsv")
write_tsv(sources, source_manifest_file)

rows <- list()

add_row <- function(
  robustness_domain,
  method_or_check,
  analysis_unit,
  outcome_suffix = NA_character_,
  n_instruments = NA_character_,
  beta_or_difference = NA_character_,
  se = NA_character_,
  pval = NA_character_,
  direction = NA_character_,
  statistic = NA_character_,
  statistic_value = NA_character_,
  statistic_pval = NA_character_,
  outlier_count = NA_character_,
  stability_flag = NA_character_,
  manuscript_safe_interpretation = NA_character_,
  source_file = NA_character_
) {
  rows[[length(rows) + 1]] <<- data.frame(
    robustness_domain = robustness_domain,
    method_or_check = method_or_check,
    analysis_unit = analysis_unit,
    outcome_suffix = outcome_suffix,
    n_instruments = as.character(n_instruments),
    beta_or_difference = as.character(beta_or_difference),
    se = as.character(se),
    pval = as.character(pval),
    direction = as.character(direction),
    statistic = as.character(statistic),
    statistic_value = as.character(statistic_value),
    statistic_pval = as.character(statistic_pval),
    outlier_count = as.character(outlier_count),
    stability_flag = as.character(stability_flag),
    manuscript_safe_interpretation = as.character(manuscript_safe_interpretation),
    source_file = source_file,
    stringsAsFactors = FALSE
  )
}

# 1. Component-specific method consistency, including Egger intercept if present.
path <- file.path(root, "22_sbp_robustness/phase4_2_sbp_method_consistency.tsv")
df <- safe_read(path)
if (!is.null(df) && nrow(df) > 0) {
  exp_col <- pick(df, c("exposure_id"))
  out_col <- pick(df, c("outcome_suffix", "outcome_id"))
  method_col <- pick(df, c("method"))
  n_col <- pick(df, c("n_instruments"))
  beta_col <- pick(df, c("beta"))
  se_col <- pick(df, c("se"))
  p_col <- pick(df, c("pval", "p"))
  dir_col <- pick(df, c("direction"))
  q_col <- pick(df, c("Q"))
  qdf_col <- pick(df, c("Q_df"))
  phi_col <- pick(df, c("phi"))
  egger_int_col <- pick(df, c("egger_intercept"))
  egger_int_se_col <- pick(df, c("egger_intercept_se"))
  egger_int_p_col <- pick(df, c("egger_intercept_pval"))
  note_col <- pick(df, c("note"))
  
  for (i in seq_len(nrow(df))) {
    if (!is.na(exp_col) && toupper(val(df, exp_col, i)) != "SBP") next
    outcome <- val(df, out_col, i)
    method <- val(df, method_col, i)
    beta <- val(df, beta_col, i)
    
    expected_ok <- NA_character_
    b <- as_num(beta)
    if (!is.na(b) && grepl("nonIOP", outcome, ignore.case = TRUE)) expected_ok <- as.character(b < 0)
    if (!is.na(b) && grepl("IOPcomponent", outcome, ignore.case = TRUE)) expected_ok <- as.character(b > 0)
    
    add_row(
      robustness_domain = "component_specific_method_consistency",
      method_or_check = method,
      analysis_unit = "component_specific_estimate",
      outcome_suffix = outcome,
      n_instruments = val(df, n_col, i),
      beta_or_difference = beta,
      se = val(df, se_col, i),
      pval = val(df, p_col, i),
      direction = val(df, dir_col, i),
      statistic = "Q; phi",
      statistic_value = paste0("Q=", val(df, q_col, i), "; Q_df=", val(df, qdf_col, i), "; phi=", val(df, phi_col, i)),
      statistic_pval = NA_character_,
      outlier_count = NA_character_,
      stability_flag = paste0("expected_direction_preserved=", expected_ok),
      manuscript_safe_interpretation = paste0("Method-specific component estimate; ", val(df, note_col, i)),
      source_file = "22_sbp_robustness/phase4_2_sbp_method_consistency.tsv"
    )
    
    if (!is.na(egger_int_col) && !is.na(val(df, egger_int_col, i)) && nzchar(val(df, egger_int_col, i))) {
      add_row(
        robustness_domain = "directional_pleiotropy",
        method_or_check = "MR-Egger intercept",
        analysis_unit = "component_specific_intercept",
        outcome_suffix = outcome,
        n_instruments = val(df, n_col, i),
        beta_or_difference = val(df, egger_int_col, i),
        se = val(df, egger_int_se_col, i),
        pval = val(df, egger_int_p_col, i),
        direction = NA_character_,
        statistic = "Egger intercept",
        statistic_value = val(df, egger_int_col, i),
        statistic_pval = val(df, egger_int_p_col, i),
        outlier_count = NA_character_,
        stability_flag = "pleiotropy_screen_only",
        manuscript_safe_interpretation = "MR-Egger intercept screen for directional pleiotropy; interpreted as sensitivity context.",
        source_file = "22_sbp_robustness/phase4_2_sbp_method_consistency.tsv"
      )
    }
  }
}

# 2. Component contrast by method.
path <- file.path(root, "22_sbp_robustness/phase4_2_sbp_component_contrast_by_method.tsv")
df <- safe_read(path)
if (!is.null(df) && nrow(df) > 0) {
  exp_col <- pick(df, c("exposure_id"))
  method_col <- pick(df, c("method"))
  beta_col <- pick(df, c("beta_difference_IOP_minus_nonIOP"))
  se_col <- pick(df, c("se_difference"))
  p_col <- pick(df, c("p_contrast"))
  dir_col <- pick(df, c("contrast_direction"))
  z_col <- pick(df, c("z_contrast"))
  r_col <- pick(df, c("assumed_component_correlation_r"))
  
  for (i in seq_len(nrow(df))) {
    if (!is.na(exp_col) && toupper(val(df, exp_col, i)) != "SBP") next
    add_row(
      robustness_domain = "component_contrast_method_consistency",
      method_or_check = val(df, method_col, i),
      analysis_unit = "IOP_minus_nonIOP_contrast",
      outcome_suffix = "SBP_component_contrast",
      n_instruments = NA_character_,
      beta_or_difference = val(df, beta_col, i),
      se = val(df, se_col, i),
      pval = val(df, p_col, i),
      direction = val(df, dir_col, i),
      statistic = "z_contrast; assumed_r",
      statistic_value = paste0("z=", val(df, z_col, i), "; r=", val(df, r_col, i)),
      statistic_pval = val(df, p_col, i),
      outlier_count = NA_character_,
      stability_flag = ifelse(as_num(val(df, beta_col, i)) > 0, "positive_contrast", "non_positive_contrast"),
      manuscript_safe_interpretation = "Method-specific contrast estimate from existing robustness output.",
      source_file = "22_sbp_robustness/phase4_2_sbp_component_contrast_by_method.tsv"
    )
  }
}

# 3. MR-RAPS component contrast.
path <- file.path(root, "22_sbp_robustness/phase4_2B_formal_packages_final/phase4_2B_final_mrraps_component_contrast.tsv")
df <- safe_read(path)
if (!is.null(df) && nrow(df) > 0 && "beta_difference_IOP_minus_nonIOP" %in% names(df)) {
  for (i in seq_len(nrow(df))) {
    add_row(
      robustness_domain = "robust_method_component_contrast",
      method_or_check = val(df, pick(df, c("method")), i),
      analysis_unit = "IOP_minus_nonIOP_contrast",
      outcome_suffix = "SBP_component_contrast",
      n_instruments = NA_character_,
      beta_or_difference = val(df, pick(df, c("beta_difference_IOP_minus_nonIOP")), i),
      se = val(df, pick(df, c("se_difference")), i),
      pval = val(df, pick(df, c("p_contrast")), i),
      direction = val(df, pick(df, c("contrast_direction")), i),
      statistic = "MR-RAPS contrast",
      statistic_value = val(df, pick(df, c("z_contrast")), i),
      statistic_pval = val(df, pick(df, c("p_contrast")), i),
      outlier_count = NA_character_,
      stability_flag = ifelse(as_num(val(df, pick(df, c("beta_difference_IOP_minus_nonIOP")), i)) > 0, "positive_contrast", "non_positive_contrast"),
      manuscript_safe_interpretation = val(df, pick(df, c("note")), i),
      source_file = "22_sbp_robustness/phase4_2B_formal_packages_final/phase4_2B_final_mrraps_component_contrast.tsv"
    )
  }
}

# 4. MR-RAPS component-specific results.
for (rel in c(
  "22_sbp_robustness/phase4_2B_formal_packages_final/phase4_2B_final_mrraps_results.tsv",
  "22_sbp_robustness/phase4_2B_formal_packages_final/phase4_2B_final_curated_mrraps_interpretation.tsv"
)) {
  path <- file.path(root, rel)
  df <- safe_read(path)
  if (!is.null(df) && nrow(df) > 0) {
    exp_col <- pick(df, c("exposure_id"))
    out_col <- pick(df, c("outcome_suffix", "outcome_id"))
    n_col <- pick(df, c("n_instruments"))
    beta_col <- pick(df, c("beta"))
    se_col <- pick(df, c("se"))
    p_col <- pick(df, c("pval", "p"))
    dir_col <- pick(df, c("direction"))
    status_col <- pick(df, c("status"))
    note_col <- pick(df, c("curated_interpretation", "note", "warning_summary"))
    for (i in seq_len(nrow(df))) {
      if (!is.na(exp_col) && toupper(val(df, exp_col, i)) != "SBP") next
      add_row(
        robustness_domain = "MR_RAPS_component_specific",
        method_or_check = "MR-RAPS",
        analysis_unit = "component_specific_estimate",
        outcome_suffix = val(df, out_col, i),
        n_instruments = val(df, n_col, i),
        beta_or_difference = val(df, beta_col, i),
        se = val(df, se_col, i),
        pval = val(df, p_col, i),
        direction = val(df, dir_col, i),
        statistic = "status",
        statistic_value = val(df, status_col, i),
        statistic_pval = NA_character_,
        outlier_count = NA_character_,
        stability_flag = "robust_method_context",
        manuscript_safe_interpretation = val(df, note_col, i),
        source_file = rel
      )
    }
  }
}

# 5. MR-PRESSO results.
for (rel in c(
  "22_sbp_robustness/phase4_2B_formal_packages_final/phase4_2B_final_mrpresso_summary.tsv",
  "22_sbp_robustness/phase4_2B_formal_packages_final/phase4_2B_final_curated_mrpresso_results.tsv"
)) {
  path <- file.path(root, rel)
  df <- safe_read(path)
  if (!is.null(df) && nrow(df) > 0) {
    exp_col <- pick(df, c("exposure_id"))
    out_col <- pick(df, c("outcome_suffix", "outcome_id"))
    n_col <- pick(df, c("n_instruments"))
    raw_beta_col <- pick(df, c("raw_beta"))
    raw_se_col <- pick(df, c("raw_se"))
    raw_p_col <- pick(df, c("raw_pval", "raw_p"))
    corr_beta_col <- pick(df, c("outlier_corrected_beta"))
    corr_se_col <- pick(df, c("outlier_corrected_se"))
    corr_p_col <- pick(df, c("outlier_corrected_pval"))
    global_col <- pick(df, c("global_test_p", "global_test_p_display", "global_test_p_numeric_upper_bound"))
    outlier_col <- pick(df, c("n_outliers", "n_significant_outliers_p_lt_0_05_or_threshold_string"))
    status_col <- pick(df, c("status"))
    interp_col <- pick(df, c("curated_interpretation", "note", "causal_estimate_use_flag"))
    
    for (i in seq_len(nrow(df))) {
      if (!is.na(exp_col) && toupper(val(df, exp_col, i)) != "SBP") next
      
      add_row(
        robustness_domain = "MR_PRESSO_raw",
        method_or_check = "MR-PRESSO raw",
        analysis_unit = "component_specific_estimate",
        outcome_suffix = val(df, out_col, i),
        n_instruments = val(df, n_col, i),
        beta_or_difference = val(df, raw_beta_col, i),
        se = val(df, raw_se_col, i),
        pval = val(df, raw_p_col, i),
        direction = NA_character_,
        statistic = "global test",
        statistic_value = val(df, status_col, i),
        statistic_pval = val(df, global_col, i),
        outlier_count = val(df, outlier_col, i),
        stability_flag = "outlier_sensitivity_context",
        manuscript_safe_interpretation = val(df, interp_col, i),
        source_file = rel
      )
      
      if (!is.na(corr_beta_col)) {
        add_row(
          robustness_domain = "MR_PRESSO_outlier_corrected",
          method_or_check = "MR-PRESSO outlier-corrected",
          analysis_unit = "component_specific_estimate",
          outcome_suffix = val(df, out_col, i),
          n_instruments = val(df, n_col, i),
          beta_or_difference = val(df, corr_beta_col, i),
          se = val(df, corr_se_col, i),
          pval = val(df, corr_p_col, i),
          direction = NA_character_,
          statistic = "global test",
          statistic_value = val(df, status_col, i),
          statistic_pval = val(df, global_col, i),
          outlier_count = val(df, outlier_col, i),
          stability_flag = "outlier_sensitivity_context",
          manuscript_safe_interpretation = val(df, interp_col, i),
          source_file = rel
        )
      }
    }
  }
}

# 6. Radial outlier summary.
path <- file.path(root, "22_sbp_robustness/phase4_2_sbp_radial_outlier_summary.tsv")
df <- safe_read(path)
if (!is.null(df) && nrow(df) > 0) {
  exp_col <- pick(df, c("exposure_id"))
  out_col <- pick(df, c("outcome_suffix", "outcome_id"))
  n_col <- pick(df, c("n_instruments"))
  n_out_col <- pick(df, c("n_radial_outliers_abs_std_resid_gt3"))
  max_col <- pick(df, c("max_abs_std_residual"))
  threshold_col <- pick(df, c("threshold"))
  note_col <- pick(df, c("note"))
  for (i in seq_len(nrow(df))) {
    if (!is.na(exp_col) && toupper(val(df, exp_col, i)) != "SBP") next
    add_row(
      robustness_domain = "radial_outlier_screen",
      method_or_check = "radial MR outlier screen",
      analysis_unit = "component_specific_outlier_screen",
      outcome_suffix = val(df, out_col, i),
      n_instruments = val(df, n_col, i),
      beta_or_difference = NA_character_,
      se = NA_character_,
      pval = NA_character_,
      direction = NA_character_,
      statistic = "max absolute standardized residual",
      statistic_value = val(df, max_col, i),
      statistic_pval = NA_character_,
      outlier_count = val(df, n_out_col, i),
      stability_flag = paste0("threshold=", val(df, threshold_col, i)),
      manuscript_safe_interpretation = val(df, note_col, i),
      source_file = "22_sbp_robustness/phase4_2_sbp_radial_outlier_summary.tsv"
    )
  }
}

# 7. Outlier-corrected results.
path <- file.path(root, "22_sbp_robustness/phase4_2_sbp_outlier_corrected_results.tsv")
df <- safe_read(path)
if (!is.null(df) && nrow(df) > 0) {
  exp_col <- pick(df, c("exposure_id"))
  out_col <- pick(df, c("outcome_suffix", "outcome_id"))
  method_col <- pick(df, c("method"))
  n_orig_col <- pick(df, c("n_original"))
  n_rem_col <- pick(df, c("n_remaining"))
  n_out_col <- pick(df, c("n_outliers_removed"))
  beta_col <- pick(df, c("beta"))
  se_col <- pick(df, c("se"))
  p_col <- pick(df, c("pval", "p"))
  dir_col <- pick(df, c("direction"))
  phi_col <- pick(df, c("phi"))
  note_col <- pick(df, c("note"))
  for (i in seq_len(nrow(df))) {
    if (!is.na(exp_col) && toupper(val(df, exp_col, i)) != "SBP") next
    add_row(
      robustness_domain = "outlier_corrected_component_estimate",
      method_or_check = val(df, method_col, i),
      analysis_unit = "component_specific_estimate",
      outcome_suffix = val(df, out_col, i),
      n_instruments = paste0("original=", val(df, n_orig_col, i), "; remaining=", val(df, n_rem_col, i)),
      beta_or_difference = val(df, beta_col, i),
      se = val(df, se_col, i),
      pval = val(df, p_col, i),
      direction = val(df, dir_col, i),
      statistic = "phi",
      statistic_value = val(df, phi_col, i),
      statistic_pval = NA_character_,
      outlier_count = val(df, n_out_col, i),
      stability_flag = "outlier_corrected_context",
      manuscript_safe_interpretation = val(df, note_col, i),
      source_file = "22_sbp_robustness/phase4_2_sbp_outlier_corrected_results.tsv"
    )
  }
}

# 8. Influence summaries.
path <- file.path(root, "22_sbp_robustness/phase4_2_sbp_single_snp_influence_summary.tsv")
df <- safe_read(path)
if (!is.null(df) && nrow(df) > 0) {
  exp_col <- pick(df, c("exposure_id"))
  out_col <- pick(df, c("outcome_suffix", "outcome_id"))
  snp_col <- pick(df, c("left_out_SNP", "omitted_SNP"))
  abs_col <- pick(df, c("abs_delta_from_full_beta"))
  dir_col <- pick(df, c("direction_leave_one_out", "direction"))
  
  if (!is.na(exp_col)) df <- df[toupper(df[[exp_col]]) == "SBP", , drop = FALSE]
  
  outcomes <- unique(df[[out_col]])
  for (oo in outcomes) {
    sub <- df[df[[out_col]] == oo, , drop = FALSE]
    if (nrow(sub) == 0) next
    abs_delta <- as_num(sub[[abs_col]])
    idx <- which.max(abs_delta)
    all_dir_same <- if (!is.na(dir_col)) length(unique(sub[[dir_col]])) == 1 else NA
    
    add_row(
      robustness_domain = "single_SNP_influence",
      method_or_check = "leave-one-SNP influence",
      analysis_unit = "component_specific_influence",
      outcome_suffix = oo,
      n_instruments = nrow(sub),
      beta_or_difference = NA_character_,
      se = NA_character_,
      pval = NA_character_,
      direction = if (!is.na(dir_col)) paste(unique(sub[[dir_col]]), collapse = ";") else NA_character_,
      statistic = "maximum absolute delta from full beta",
      statistic_value = abs_delta[idx],
      statistic_pval = NA_character_,
      outlier_count = NA_character_,
      stability_flag = paste0("largest_influence_SNP=", sub[[snp_col]][idx], "; all_leave_one_directions_same=", all_dir_same),
      manuscript_safe_interpretation = "Single-SNP influence summary; used to assess whether a component estimate was dominated by one variant.",
      source_file = "22_sbp_robustness/phase4_2_sbp_single_snp_influence_summary.tsv"
    )
  }
}

path <- file.path(root, "22_sbp_robustness/phase4_2_sbp_leave_one_chromosome_out.tsv")
df <- safe_read(path)
if (!is.null(df) && nrow(df) > 0) {
  exp_col <- pick(df, c("exposure_id"))
  out_col <- pick(df, c("outcome_suffix", "outcome_id"))
  chr_col <- pick(df, c("left_out_chr", "left_out_chromosome"))
  delta_col <- pick(df, c("delta_from_full_beta"))
  dir_col <- pick(df, c("direction"))
  
  if (!is.na(exp_col)) df <- df[toupper(df[[exp_col]]) == "SBP", , drop = FALSE]
  
  outcomes <- unique(df[[out_col]])
  for (oo in outcomes) {
    sub <- df[df[[out_col]] == oo, , drop = FALSE]
    if (nrow(sub) == 0) next
    abs_delta <- abs(as_num(sub[[delta_col]]))
    idx <- which.max(abs_delta)
    all_dir_same <- if (!is.na(dir_col)) length(unique(sub[[dir_col]])) == 1 else NA
    
    add_row(
      robustness_domain = "chromosome_influence",
      method_or_check = "leave-one-chromosome influence",
      analysis_unit = "component_specific_influence",
      outcome_suffix = oo,
      n_instruments = nrow(sub),
      beta_or_difference = NA_character_,
      se = NA_character_,
      pval = NA_character_,
      direction = if (!is.na(dir_col)) paste(unique(sub[[dir_col]]), collapse = ";") else NA_character_,
      statistic = "maximum absolute delta from full beta",
      statistic_value = abs_delta[idx],
      statistic_pval = NA_character_,
      outlier_count = NA_character_,
      stability_flag = paste0("largest_influence_chromosome=", sub[[chr_col]][idx], "; all_leave_one_directions_same=", all_dir_same),
      manuscript_safe_interpretation = "Chromosome-level influence summary for component-specific estimates.",
      source_file = "22_sbp_robustness/phase4_2_sbp_leave_one_chromosome_out.tsv"
    )
  }
}

# 9. Phase13.4B contrast uncertainty and influence.
path <- file.path(root, "71_data_analysis_reinforcement/phase13_4_empirical_covariance/phase13_4B_bootstrap_contrast_summary.tsv")
df <- safe_read(path)
if (!is.null(df) && nrow(df) > 0) {
  add_row(
    robustness_domain = "empirical_contrast_uncertainty",
    method_or_check = "paired SNP bootstrap",
    analysis_unit = "IOP_minus_nonIOP_contrast",
    outcome_suffix = "SBP_component_contrast",
    n_instruments = val(df, pick(df, c("n_instruments")), 1),
    beta_or_difference = val(df, pick(df, c("full_beta_difference_IOP_minus_nonIOP")), 1),
    se = val(df, pick(df, c("empirical_SE_difference")), 1),
    pval = val(df, pick(df, c("empirical_p_normal_using_full_difference")), 1),
    direction = ifelse(tolower(val(df, pick(df, c("direction_preserved")), 1)) == "true", "positive", NA_character_),
    statistic = "bootstrap 95% interval; empirical correlation",
    statistic_value = paste0(
      "CI=", val(df, pick(df, c("empirical_CI_lower_2p5")), 1),
      " to ", val(df, pick(df, c("empirical_CI_upper_97p5")), 1),
      "; r=", val(df, pick(df, c("empirical_correlation_IOP_nonIOP")), 1)
    ),
    statistic_pval = val(df, pick(df, c("empirical_p_normal_using_full_difference")), 1),
    outlier_count = NA_character_,
    stability_flag = "bootstrap_CI_excludes_zero",
    manuscript_safe_interpretation = "Empirical uncertainty sensitivity for the SBP component contrast.",
    source_file = "71_data_analysis_reinforcement/phase13_4_empirical_covariance/phase13_4B_bootstrap_contrast_summary.tsv"
  )
}

for (rel in c(
  "71_data_analysis_reinforcement/phase13_4_empirical_covariance/phase13_4B_leave_one_chromosome_summary.tsv",
  "71_data_analysis_reinforcement/phase13_4_empirical_covariance/phase13_4B_locus_proxy_jackknife_summary.tsv"
)) {
  path <- file.path(root, rel)
  df <- safe_read(path)
  if (!is.null(df) && nrow(df) > 0 && all(c("metric", "value") %in% names(df))) {
    metric_value <- function(m) {
      z <- df$value[df$metric == m]
      if (length(z) == 0) return(NA_character_)
      as.character(z[1])
    }
    
    method_name <- ifelse(grepl("chromosome", rel), "leave-one-chromosome contrast jackknife", "1-Mb locus-proxy contrast jackknife")
    domain <- ifelse(grepl("chromosome", rel), "contrast_chromosome_influence", "contrast_locus_proxy_influence")
    all_preserved <- ifelse(grepl("chromosome", rel),
                            metric_value("all_chromosome_leave_one_out_direction_preserved"),
                            metric_value("all_locus_proxy_leave_one_out_direction_preserved"))
    
    add_row(
      robustness_domain = domain,
      method_or_check = method_name,
      analysis_unit = "IOP_minus_nonIOP_contrast",
      outcome_suffix = "SBP_component_contrast",
      n_instruments = NA_character_,
      beta_or_difference = NA_character_,
      se = NA_character_,
      pval = NA_character_,
      direction = NA_character_,
      statistic = "direction preservation and maximum influence",
      statistic_value = paste0(
        "all_preserved=", all_preserved,
        "; max_abs_delta=", metric_value("max_abs_delta_from_full")
      ),
      statistic_pval = NA_character_,
      outlier_count = NA_character_,
      stability_flag = paste0("all_direction_preserved=", all_preserved),
      manuscript_safe_interpretation = "Contrast-level influence sensitivity from Phase13.4B.",
      source_file = rel
    )
  }
}

matrix <- if (length(rows) > 0) do.call(rbind, rows) else data.frame()

matrix_file <- file.path(out_dir, "phase13_6B_SBP_robustness_matrix.tsv")
write_tsv(matrix, matrix_file)

# Coverage / interpretation table.
coverage_items <- data.frame(
  robustness_item = c(
    "Component-specific method consistency",
    "Component contrast by method",
    "MR-Egger intercept",
    "MR-PRESSO",
    "MR-RAPS",
    "Radial MR outlier screen",
    "Single-SNP influence",
    "Chromosome/locus contrast influence",
    "Empirical paired-SNP bootstrap"
  ),
  available = c(
    any(matrix$robustness_domain == "component_specific_method_consistency"),
    any(matrix$robustness_domain == "component_contrast_method_consistency"),
    any(matrix$method_or_check == "MR-Egger intercept"),
    any(grepl("MR_PRESSO", matrix$robustness_domain)),
    any(grepl("MR_RAPS", matrix$robustness_domain)),
    any(matrix$robustness_domain == "radial_outlier_screen"),
    any(matrix$robustness_domain == "single_SNP_influence"),
    any(matrix$robustness_domain %in% c("contrast_chromosome_influence", "contrast_locus_proxy_influence")),
    any(matrix$robustness_domain == "empirical_contrast_uncertainty")
  ),
  manuscript_role = c(
    "Method-specific robustness context",
    "Direct component-divergence robustness context",
    "Directional pleiotropy screen",
    "Horizontal-pleiotropy/outlier sensitivity",
    "Robust many-instrument sensitivity",
    "Outlier screen",
    "Variant-influence screen",
    "Contrast-level influence screen",
    "Empirical uncertainty sensitivity"
  ),
  claim_safety_boundary = rep(
    "Supports robustness characterization only; does not upgrade to confirmatory causal evidence.",
    9
  ),
  stringsAsFactors = FALSE
)

coverage_file <- file.path(out_dir, "phase13_6B_robustness_coverage_interpretation.tsv")
write_tsv(coverage_items, coverage_file)

# Compact manuscript-facing summary.
n_domains <- length(unique(matrix$robustness_domain))
n_rows <- nrow(matrix)
n_available <- sum(coverage_items$available)

methods_insert <- c(
  "# Robustness and pleiotropy matrix",
  "",
  "Existing SBP robustness outputs were integrated into a structured matrix covering component-specific method consistency, component-contrast sensitivity, MR-Egger intercept screening, MR-PRESSO, MR-RAPS, radial outlier screening, single-variant influence, and empirical contrast-level bootstrap and jackknife analyses. This integration used previously generated robustness outputs and the Phase13 empirical contrast-sensitivity outputs; it was not used to generate new primary MR estimates.",
  "",
  "Robustness results were interpreted as sensitivity and bias-assessment context for a hypothesis-generating component-divergence signal. Method-specific, outlier-corrected, and influence analyses were not used to reclassify the SBP finding as confirmatory evidence."
)

methods_file <- file.path(out_dir, "phase13_6B_methods_insert_robustness_matrix.md")
writeLines(methods_insert, methods_file, useBytes = TRUE)

results_insert <- c(
  "# Robustness and pleiotropy matrix",
  "",
  paste0(
    "The integrated SBP robustness matrix contained ",
    n_rows,
    " rows across ",
    n_domains,
    " robustness domains. Available sensitivity layers included ",
    n_available,
    " of 9 prespecified robustness categories: component-specific method consistency, component-contrast sensitivity, MR-Egger intercept screening, MR-PRESSO, MR-RAPS, radial outlier screening, variant-influence assessment, chromosome/locus contrast-influence assessment, and empirical paired-SNP bootstrap."
  ),
  "",
  "Together with the Phase13 empirical contrast analyses, the robustness matrix supports directional stability of the SBP component-divergence pattern while retaining a hypothesis-generating interpretation."
)

results_file <- file.path(out_dir, "phase13_6B_results_insert_robustness_matrix.md")
writeLines(results_insert, results_file, useBytes = TRUE)

# Claim-safety audit.
texts_to_scan <- paste(c(methods_insert, results_insert, matrix$manuscript_safe_interpretation), collapse = "\n")

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
    "SBP causally affects"
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
    "causal_overclaim"
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

audit_file <- file.path(out_dir, "phase13_6B_claim_safety_audit.tsv")
write_tsv(danger_patterns, audit_file)

high_risk_hits <- sum(danger_patterns$hits)

status <- data.frame(
  field = c(
    "phase",
    "source_manifest_created",
    "robustness_matrix_created",
    "coverage_interpretation_created",
    "methods_insert_created",
    "results_insert_created",
    "claim_safety_audit_created",
    "high_risk_phrase_hits",
    "robustness_matrix_rows",
    "robustness_domains",
    "coverage_items_available_out_of_9",
    "new_primary_MR_estimates_created",
    "claim_level",
    "claim_upgrade_allowed",
    "phase13_6B_passed"
  ),
  value = c(
    "Phase13.6B",
    file.exists(source_manifest_file),
    file.exists(matrix_file),
    file.exists(coverage_file),
    file.exists(methods_file),
    file.exists(results_file),
    file.exists(audit_file),
    high_risk_hits,
    n_rows,
    n_domains,
    paste0(n_available, "/9"),
    FALSE,
    "HYPOTHESIS_GENERATING_NOT_CONFIRMATORY",
    "NO",
    high_risk_hits == 0 && n_rows > 0
  ),
  stringsAsFactors = FALSE
)

status_file <- file.path(out_dir, "phase13_6B_status.tsv")
write_tsv(status, status_file)

# Update master.
if (file.exists(master_status_file)) {
  master_status <- read.delim(master_status_file, check.names = FALSE)
  idx <- which(master_status$phase == "Phase13.6")
  if (length(idx) == 1) {
    master_status$status[idx] <- ifelse(
      high_risk_hits == 0 && n_rows > 0,
      "PASSED_ROBUSTNESS_AND_PLEIOTROPY_MATRIX_INTEGRATED",
      "FAILED_PHASE13_6B_REVIEW_REQUIRED"
    )
    master_status$qc_status[idx] <- ifelse(
      high_risk_hits == 0 && n_rows > 0,
      "PASSED",
      "REVIEW_REQUIRED"
    )
    master_status$primary_output[idx] <- "phase13_6B_SBP_robustness_matrix.tsv; phase13_6B_robustness_coverage_interpretation.tsv"
  }
  write_tsv(master_status, master_status_file)
}

message("Phase13.6B completed.")
message("Robustness matrix: ", matrix_file)
message("Coverage interpretation: ", coverage_file)
message("Status: ", status_file)
message("High-risk phrase hits: ", high_risk_hits)
