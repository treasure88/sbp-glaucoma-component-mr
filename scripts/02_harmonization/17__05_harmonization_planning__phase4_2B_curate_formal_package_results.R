OUTDIR <- "../../22_sbp_robustness/phase4_2B_formal_packages_final"

CURATED_OUT <- file.path(OUTDIR, "phase4_2B_final_curated_formal_package_interpretation.tsv")
MRPRESSO_CURATED_OUT <- file.path(OUTDIR, "phase4_2B_final_curated_mrpresso_results.tsv")
MRPRESSO_OUTLIER_CURATED_OUT <- file.path(OUTDIR, "phase4_2B_final_curated_mrpresso_significant_outliers.tsv")
MRRAPS_CURATED_OUT <- file.path(OUTDIR, "phase4_2B_final_curated_mrraps_interpretation.tsv")
STATUS_OUT <- file.path(OUTDIR, "phase4_2B_final_curated_status.tsv")

p_to_numeric <- function(x) {
  if (length(x) == 0 || is.null(x) || is.na(x)) return(NA_real_)
  x <- as.character(x)[1]
  x <- trimws(x)
  x <- gsub("<", "", x, fixed = TRUE)
  suppressWarnings(as.numeric(x))
}

p_display <- function(x) {
  if (length(x) == 0 || is.null(x) || is.na(x)) return("NA")
  as.character(x)[1]
}

dir_label <- function(beta) {
  beta <- suppressWarnings(as.numeric(beta))
  if (length(beta) == 0 || is.na(beta)) return("NA")
  if (beta > 0) return("positive")
  if (beta < 0) return("negative")
  "zero"
}

extract_mrpresso <- function(outcome_suffix) {
  path <- file.path(OUTDIR, paste0("MRPRESSO_raw_", outcome_suffix, ".rds"))
  obj <- readRDS(path)
  res <- obj$result

  main <- res[["Main MR results"]]
  gp <- res[["MR-PRESSO results"]][["Global Test"]][["Pvalue"]]
  rss <- res[["MR-PRESSO results"]][["Global Test"]][["RSSobs"]]

  raw <- main[main[["MR Analysis"]] == "Raw", , drop = FALSE]
  corr <- main[main[["MR Analysis"]] == "Outlier-corrected", , drop = FALSE]

  outlier_tab <- res[["MR-PRESSO results"]][["Outlier Test"]]
  distortion <- res[["MR-PRESSO results"]][["Distortion Test"]]

  n_outlier_tests <- 0
  n_significant_outliers <- 0
  outlier_indices <- "NA"

  sig_outliers <- data.frame()

  if (is.data.frame(outlier_tab)) {
    n_outlier_tests <- nrow(outlier_tab)

    p_raw <- as.character(outlier_tab$Pvalue)
    p_num <- sapply(p_raw, p_to_numeric)
    is_sig <- grepl("^<", p_raw) | (!is.na(p_num) & p_num < 0.05)

    n_significant_outliers <- sum(is_sig, na.rm = TRUE)

    if (n_significant_outliers > 0) {
      sig_outliers <- outlier_tab[is_sig, , drop = FALSE]
      sig_outliers$outlier_index <- which(is_sig)
      sig_outliers$exposure_id <- "SBP"
      sig_outliers$outcome_suffix <- outcome_suffix
      sig_outliers <- sig_outliers[, c("exposure_id", "outcome_suffix", "outlier_index", "RSSobs", "Pvalue")]
      outlier_indices <- paste(which(is_sig), collapse = ";")
    }
  }

  distortion_p <- NA
  distortion_coeff <- NA

  if (!is.null(distortion)) {
    distortion_p <- tryCatch(distortion[["Pvalue"]], error = function(e) NA)
    distortion_coeff <- tryCatch(as.numeric(distortion[["Distortion Coefficient"]][1]), error = function(e) NA)
  }

  corrected_beta <- if (nrow(corr) > 0) corr[["Causal Estimate"]][1] else NA
  corrected_p <- if (nrow(corr) > 0) corr[["P-value"]][1] else NA

  causal_estimate_use_flag <- "USABLE"
  if (outcome_suffix == "GBS_IOPcomponent") {
    causal_estimate_use_flag <- "USE_FOR_GLOBAL_OUTLIER_STATUS_ONLY_REVIEW_CAUSAL_ESTIMATE"
  }

  interpretation <- "NA"
  if (outcome_suffix == "GBS_nonIOPcomponent") {
    interpretation <- "MR-PRESSO global test indicated horizontal pleiotropy/outlier structure; one significant outlier was identified; outlier-corrected estimate remained negative and nominally significant."
  }
  if (outcome_suffix == "GBS_IOPcomponent") {
    interpretation <- "MR-PRESSO global test did not indicate outlier pleiotropy; no outliers were identified. Causal estimate from MR-PRESSO should be treated cautiously because it differs from the independent IVW estimate."
  }

  row <- data.frame(
    exposure_id = "SBP",
    outcome_suffix = outcome_suffix,
    package = "MRPRESSO",
    status = "OK_CURATED",
    n_instruments = obj$n_instruments,
    mode = obj$mode,
    nb_distribution = obj$nb_distribution,
    global_test_RSSobs = rss,
    global_test_p_display = p_display(gp),
    global_test_p_numeric_upper_bound = p_to_numeric(gp),
    raw_beta = raw[["Causal Estimate"]][1],
    raw_se = raw[["Sd"]][1],
    raw_pval = raw[["P-value"]][1],
    outlier_corrected_beta = corrected_beta,
    outlier_corrected_pval = corrected_p,
    n_outlier_tests = n_outlier_tests,
    n_significant_outliers_p_lt_0_05_or_threshold_string = n_significant_outliers,
    significant_outlier_indices = outlier_indices,
    distortion_pval = distortion_p,
    distortion_coefficient = distortion_coeff,
    causal_estimate_use_flag = causal_estimate_use_flag,
    curated_interpretation = interpretation,
    stringsAsFactors = FALSE
  )

  list(row = row, sig_outliers = sig_outliers)
}

extract_mrraps <- function(outcome_suffix) {
  path <- file.path(OUTDIR, paste0("MRRAPS_raw_", outcome_suffix, ".rds"))
  obj <- readRDS(path)
  res <- obj$result

  beta <- res$beta.hat
  se <- res$beta.se
  pval <- 2 * pnorm(abs(beta / se), lower.tail = FALSE)

  expected_direction <- ifelse(outcome_suffix == "GBS_nonIOPcomponent", "negative", "positive")
  observed_direction <- dir_label(beta)

  direction_consistency <- ifelse(observed_direction == expected_direction, "CONSISTENT_WITH_PRIMARY_DIRECTION", "DISCORDANT_WITH_PRIMARY_DIRECTION")

  interpretation <- ifelse(
    outcome_suffix == "GBS_nonIOPcomponent",
    "MR-RAPS yielded a positive nonIOP estimate, discordant with IVW and MR-PRESSO outlier-corrected negative estimates.",
    "MR-RAPS yielded a positive IOP estimate, consistent with the primary positive IOP direction."
  )

  data.frame(
    exposure_id = "SBP",
    outcome_suffix = outcome_suffix,
    package = "mr.raps",
    status = "OK_CURATED_WITH_WARNINGS",
    n_instruments = obj$n_instruments,
    beta = beta,
    se = se,
    pval = pval,
    direction = observed_direction,
    expected_primary_direction = expected_direction,
    direction_consistency = direction_consistency,
    tau2_hat = res$tau2.hat,
    n_warnings = length(obj$warnings),
    warning_summary = paste(unique(obj$warnings), collapse = " | "),
    curated_interpretation = interpretation,
    stringsAsFactors = FALSE
  )
}

mrp_non <- extract_mrpresso("GBS_nonIOPcomponent")
mrp_iop <- extract_mrpresso("GBS_IOPcomponent")
mrp_rows <- rbind(mrp_non$row, mrp_iop$row)

outlier_rows <- rbind(mrp_non$sig_outliers, mrp_iop$sig_outliers)
if (is.null(outlier_rows) || nrow(outlier_rows) == 0) {
  outlier_rows <- data.frame(
    exposure_id = character(),
    outcome_suffix = character(),
    outlier_index = character(),
    RSSobs = character(),
    Pvalue = character(),
    stringsAsFactors = FALSE
  )
}

raps_rows <- rbind(
  extract_mrraps("GBS_nonIOPcomponent"),
  extract_mrraps("GBS_IOPcomponent")
)

contrast <- read.delim(
  file.path(OUTDIR, "phase4_2B_final_mrraps_component_contrast.tsv"),
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

contrast_r0 <- contrast[contrast$assumed_component_correlation_r == 0, , drop = FALSE]

curated_rows <- data.frame(
  evidence_item = c(
    "MR-PRESSO nonIOP outlier-corrected result",
    "MR-PRESSO IOP global/outlier result",
    "MR-RAPS nonIOP result",
    "MR-RAPS IOP result",
    "MR-RAPS component contrast r=0",
    "Overall Phase 4.2B interpretation"
  ),
  evidence_status = c(
    "SUPPORTIVE_FOR_NONIOP_NEGATIVE_AFTER_OUTLIER_REMOVAL",
    "NO_MRPRESSO_OUTLIERS_FOR_IOP",
    "DISCORDANT_DIRECTION",
    "SUPPORTIVE_FOR_IOP_POSITIVE",
    "SUPPORTIVE_FOR_IOP_MORE_POSITIVE_THAN_NONIOP",
    "MIXED_BUT_COMPONENT_CONTRAST_SUPPORTED"
  ),
  key_result = c(
    paste0("beta=", mrp_rows$outlier_corrected_beta[mrp_rows$outcome_suffix == "GBS_nonIOPcomponent"],
           "; p=", mrp_rows$outlier_corrected_pval[mrp_rows$outcome_suffix == "GBS_nonIOPcomponent"],
           "; outliers=", mrp_rows$n_significant_outliers_p_lt_0_05_or_threshold_string[mrp_rows$outcome_suffix == "GBS_nonIOPcomponent"]),
    paste0("global_p=", mrp_rows$global_test_p_display[mrp_rows$outcome_suffix == "GBS_IOPcomponent"],
           "; outliers=", mrp_rows$n_significant_outliers_p_lt_0_05_or_threshold_string[mrp_rows$outcome_suffix == "GBS_IOPcomponent"]),
    paste0("beta=", raps_rows$beta[raps_rows$outcome_suffix == "GBS_nonIOPcomponent"],
           "; p=", raps_rows$pval[raps_rows$outcome_suffix == "GBS_nonIOPcomponent"],
           "; direction=", raps_rows$direction[raps_rows$outcome_suffix == "GBS_nonIOPcomponent"]),
    paste0("beta=", raps_rows$beta[raps_rows$outcome_suffix == "GBS_IOPcomponent"],
           "; p=", raps_rows$pval[raps_rows$outcome_suffix == "GBS_IOPcomponent"],
           "; direction=", raps_rows$direction[raps_rows$outcome_suffix == "GBS_IOPcomponent"]),
    paste0("beta_diff=", contrast_r0$beta_difference_IOP_minus_nonIOP[1],
           "; p=", contrast_r0$p_contrast[1],
           "; direction=", contrast_r0$contrast_direction[1]),
    "Formal packages support component-level divergence/contrast, but MR-RAPS does not reproduce the IVW nonIOP-negative direction."
  ),
  manuscript_use = c(
    "Can report as MR-PRESSO outlier-corrected sensitivity supporting nonIOP negative direction.",
    "Can report as no MR-PRESSO outlier evidence for IOP; avoid overinterpreting MR-PRESSO causal estimate.",
    "Report as discordant robust estimator; weakens claim of universal direction consistency.",
    "Report as supportive for positive IOP direction.",
    "Report as supportive formal-package contrast sensitivity.",
    "Use cautious wording: sensitivity evidence is mixed, but contrast remains supported."
  ),
  stringsAsFactors = FALSE
)

write.table(mrp_rows, MRPRESSO_CURATED_OUT, sep = "\t", row.names = FALSE, quote = FALSE)
write.table(outlier_rows, MRPRESSO_OUTLIER_CURATED_OUT, sep = "\t", row.names = FALSE, quote = FALSE)
write.table(raps_rows, MRRAPS_CURATED_OUT, sep = "\t", row.names = FALSE, quote = FALSE)
write.table(curated_rows, CURATED_OUT, sep = "\t", row.names = FALSE, quote = FALSE)

con <- file(STATUS_OUT, "w")
writeLines("phase\tstatus\tkey_result\tnote", con)
writeLines(paste(
  "Phase 4.2B-final curated interpretation",
  "PASSED_WITH_MIXED_FORMAL_PACKAGE_EVIDENCE",
  "MRPRESSO_nonIOP_supportive;MRRAPS_contrast_supportive;MRRAPS_nonIOP_direction_discordant",
  "Formal package outputs curated after raw-object inspection",
  sep = "\t"
), con)
writeLines(paste(
  "MRPRESSO_nonIOP",
  "SUPPORTIVE",
  "outlier_corrected_negative",
  "One significant outlier; corrected beta remains negative",
  sep = "\t"
), con)
writeLines(paste(
  "MRPRESSO_IOP",
  "NO_OUTLIER",
  "global_p_1",
  "No outlier identified; causal estimate from MR-PRESSO should not replace independent IVW",
  sep = "\t"
), con)
writeLines(paste(
  "MR_RAPS_nonIOP",
  "DISCORDANT",
  "positive",
  "Does not reproduce primary nonIOP negative direction",
  sep = "\t"
), con)
writeLines(paste(
  "MR_RAPS_contrast",
  "SUPPORTIVE",
  "IOP_MORE_POSITIVE_THAN_NONIOP",
  "MR-RAPS component contrast remains significant",
  sep = "\t"
), con)
close(con)

cat("===== Phase 4.2B curated interpretation completed =====\n")
cat("Wrote:", STATUS_OUT, "\n")
cat("Wrote:", MRPRESSO_CURATED_OUT, "\n")
cat("Wrote:", MRPRESSO_OUTLIER_CURATED_OUT, "\n")
cat("Wrote:", MRRAPS_CURATED_OUT, "\n")
cat("Wrote:", CURATED_OUT, "\n")
