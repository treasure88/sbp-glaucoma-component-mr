# Phase 8.3C: SBP instrument-quality subset sensitivity
# Purpose: evaluate whether the SBP component-divergence contrast depends on weak instruments or palindromic SNPs
# Method: base R IVW with multiplicative random-effects SE
# Language policy: English labels for statistical outputs

options(stringsAsFactors = FALSE)

out_dir <- "../../50_sbp_sensitivity"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

non_file <- "../../16_mr_input_datasets/pairwise/SBP__GBS_nonIOPcomponent.mr_input.tsv.gz"
iop_file <- "../../16_mr_input_datasets/pairwise/SBP__GBS_IOPcomponent.mr_input.tsv.gz"

subset_counts_file <- file.path(out_dir, "phase8_3C_instrument_quality_subset_counts.tsv")
results_file <- file.path(out_dir, "phase8_3C_instrument_quality_sensitivity_results.tsv")
summary_file <- file.path(out_dir, "phase8_3C_instrument_quality_sensitivity_summary.tsv")
status_file <- file.path(out_dir, "phase8_3C_instrument_quality_sensitivity_status.tsv")

read_tsv <- function(path) {
  read.delim(
    path,
    header = TRUE,
    sep = "\t",
    quote = "",
    comment.char = "",
    check.names = FALSE,
    fill = TRUE
  )
}

to_num <- function(x) suppressWarnings(as.numeric(x))

write_status <- function(status, note) {
  x <- data.frame(
    phase = "Phase8.3C",
    status = status,
    note = note,
    timestamp = as.character(Sys.time()),
    stringsAsFactors = FALSE
  )
  write.table(x, status_file, sep = "\t", quote = FALSE, row.names = FALSE)
}

ivw_random <- function(d) {
  bx <- to_num(d$beta_exposure)
  by <- to_num(d$beta_outcome)
  se_by <- to_num(d$se_outcome)

  keep <- is.finite(bx) & is.finite(by) & is.finite(se_by) & se_by > 0 & bx != 0
  d <- d[keep, , drop = FALSE]
  bx <- bx[keep]
  by <- by[keep]
  se_by <- se_by[keep]

  n <- length(bx)

  if (n < 3) {
    return(data.frame(
      n_instruments = n,
      beta = NA_real_,
      se_fixed = NA_real_,
      se_random = NA_real_,
      p_random = NA_real_,
      Q = NA_real_,
      Q_df = NA_integer_,
      Q_pval = NA_real_,
      phi = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  w <- 1 / (se_by^2)
  beta <- sum(w * bx * by) / sum(w * bx^2)
  se_fixed <- sqrt(1 / sum(w * bx^2))
  residual <- by - beta * bx
  Q <- sum(w * residual^2)
  Q_df <- n - 1
  phi <- max(1, Q / Q_df)
  se_random <- se_fixed * sqrt(phi)
  p_random <- 2 * pnorm(-abs(beta / se_random))
  Q_pval <- pchisq(Q, df = Q_df, lower.tail = FALSE)

  data.frame(
    n_instruments = n,
    beta = beta,
    se_fixed = se_fixed,
    se_random = se_random,
    p_random = p_random,
    Q = Q,
    Q_df = Q_df,
    Q_pval = Q_pval,
    phi = phi,
    stringsAsFactors = FALSE
  )
}

estimate_component <- function(d, component_id, subset_name) {
  est <- ivw_random(d)
  data.frame(
    component = component_id,
    subset_name = subset_name,
    n_instruments = est$n_instruments,
    beta = est$beta,
    se = est$se_random,
    p = est$p_random,
    Q = est$Q,
    Q_df = est$Q_df,
    Q_pval = est$Q_pval,
    phi = est$phi,
    stringsAsFactors = FALSE
  )
}

is_palindromic_pair <- function(a1, a2) {
  a1 <- toupper(as.character(a1))
  a2 <- toupper(as.character(a2))
  pair <- paste0(a1, a2)
  pair %in% c("AT", "TA", "CG", "GC")
}

get_f_stat <- function(d) {
  if ("instrument_F_stat" %in% names(d)) {
    f <- to_num(d$instrument_F_stat)
  } else if ("F_stat" %in% names(d)) {
    f <- to_num(d$F_stat)
  } else {
    f <- to_num(d$beta_exposure)^2 / to_num(d$se_exposure)^2
  }
  f
}

tryCatch({

  if (!file.exists(non_file)) {
    write_status("FAILED_INPUT_MISSING", paste("Missing file:", non_file))
    stop("Missing nonIOP file")
  }

  if (!file.exists(iop_file)) {
    write_status("FAILED_INPUT_MISSING", paste("Missing file:", iop_file))
    stop("Missing IOP file")
  }

  non <- read_tsv(non_file)
  iop <- read_tsv(iop_file)

  required <- c("SNP", "beta_exposure", "se_exposure", "beta_outcome", "se_outcome", "effect_allele", "other_allele")

  miss_non <- setdiff(required, names(non))
  miss_iop <- setdiff(required, names(iop))

  if (length(miss_non) > 0 || length(miss_iop) > 0) {
    write_status(
      "FAILED_COLUMN_QC",
      paste(
        "Missing nonIOP columns:", paste(miss_non, collapse = ";"),
        "| Missing IOP columns:", paste(miss_iop, collapse = ";")
      )
    )
    stop("Missing required columns")
  }

  shared_snps <- intersect(non$SNP, iop$SNP)

  non <- non[non$SNP %in% shared_snps, , drop = FALSE]
  iop <- iop[iop$SNP %in% shared_snps, , drop = FALSE]

  non <- non[order(non$SNP), , drop = FALSE]
  iop <- iop[order(iop$SNP), , drop = FALSE]

  non$F_for_subset <- get_f_stat(non)
  iop$F_for_subset <- get_f_stat(iop)

  # Use exposure alleles from nonIOP file; exposure should be the same across both component outcomes.
  pal <- is_palindromic_pair(non$effect_allele, non$other_allele)

  subset_definitions <- list(
    FULL_SHARED = rep(TRUE, nrow(non)),
    F_GE_30 = non$F_for_subset >= 30 & iop$F_for_subset >= 30,
    F_GE_50 = non$F_for_subset >= 50 & iop$F_for_subset >= 50,
    NONPALINDROMIC = !pal,
    NONPALINDROMIC_AND_F_GE_30 = !pal & non$F_for_subset >= 30 & iop$F_for_subset >= 30,
    NONPALINDROMIC_AND_F_GE_50 = !pal & non$F_for_subset >= 50 & iop$F_for_subset >= 50
  )

  count_rows <- list()
  result_rows <- list()

  full_beta_diff <- NA_real_

  for (subset_name in names(subset_definitions)) {
    keep <- subset_definitions[[subset_name]]
    keep[is.na(keep)] <- FALSE

    n_subset <- sum(keep)
    run_status <- ifelse(n_subset >= 50, "RUN", ifelse(n_subset >= 3, "LOW_N_RUN_WITH_CAUTION", "SKIP_TOO_FEW"))

    count_rows[[length(count_rows) + 1]] <- data.frame(
      subset_name = subset_name,
      n_shared_snps = length(shared_snps),
      n_subset = n_subset,
      n_excluded = length(shared_snps) - n_subset,
      run_status = run_status,
      min_F_nonIOP = suppressWarnings(min(non$F_for_subset[keep], na.rm = TRUE)),
      min_F_IOP = suppressWarnings(min(iop$F_for_subset[keep], na.rm = TRUE)),
      n_palindromic_excluded = ifelse(grepl("NONPALINDROMIC", subset_name), sum(pal, na.rm = TRUE), 0),
      stringsAsFactors = FALSE
    )

    if (run_status == "SKIP_TOO_FEW") next

    non_sub <- non[keep, , drop = FALSE]
    iop_sub <- iop[keep, , drop = FALSE]

    a <- estimate_component(non_sub, "GBS_nonIOPcomponent", subset_name)
    b <- estimate_component(iop_sub, "GBS_IOPcomponent", subset_name)

    beta_diff <- b$beta - a$beta
    se_diff <- sqrt(b$se^2 + a$se^2)
    z <- beta_diff / se_diff
    p <- 2 * pnorm(-abs(z))

    if (subset_name == "FULL_SHARED") {
      full_beta_diff <- beta_diff
    }

    result_rows[[length(result_rows) + 1]] <- data.frame(
      exposure_id = "SBP",
      subset_name = subset_name,
      n_nonIOP = a$n_instruments,
      beta_nonIOP = a$beta,
      se_nonIOP = a$se,
      p_nonIOP = a$p,
      n_IOP = b$n_instruments,
      beta_IOP = b$beta,
      se_IOP = b$se,
      p_IOP = b$p,
      beta_difference_IOP_minus_nonIOP = beta_diff,
      se_difference_r0 = se_diff,
      z_contrast = z,
      p_contrast_r0 = p,
      contrast_direction = ifelse(beta_diff > 0, "IOP_MORE_POSITIVE_THAN_NONIOP",
                                  ifelse(beta_diff < 0, "IOP_MORE_NEGATIVE_THAN_NONIOP", "NO_DIFFERENCE")),
      nonIOP_direction = ifelse(a$beta < 0, "negative", ifelse(a$beta > 0, "positive", "zero")),
      IOP_direction = ifelse(b$beta < 0, "negative", ifelse(b$beta > 0, "positive", "zero")),
      subset_run_status = run_status,
      stringsAsFactors = FALSE
    )
  }

  counts <- do.call(rbind, count_rows)
  write.table(counts, subset_counts_file, sep = "\t", quote = FALSE, row.names = FALSE)

  if (length(result_rows) == 0) {
    write_status("FAILED_NO_RUNNABLE_SUBSETS", "No instrument-quality subset had enough SNPs for IVW analysis.")
    stop("No runnable subsets")
  }

  res <- do.call(rbind, result_rows)

  full_row <- res[res$subset_name == "FULL_SHARED", , drop = FALSE]
  full_beta <- full_row$beta_difference_IOP_minus_nonIOP

  res$contrast_direction_preserved_vs_full <- ifelse(
    sign(res$beta_difference_IOP_minus_nonIOP) == sign(full_beta),
    "YES",
    "NO"
  )

  res$contrast_beta_percent_of_full <- 100 * res$beta_difference_IOP_minus_nonIOP / full_beta
  res$absolute_percent_attenuation_vs_full <- abs(100 - res$contrast_beta_percent_of_full)

  res$nominal_contrast_p_lt_0_05 <- ifelse(
    is.finite(res$p_contrast_r0) & res$p_contrast_r0 < 0.05,
    "YES",
    "NO"
  )

  res$component_sign_pattern_preserved <- ifelse(
    res$nonIOP_direction == "negative" & res$IOP_direction == "positive",
    "YES",
    "NO"
  )

  write.table(res, results_file, sep = "\t", quote = FALSE, row.names = FALSE)

  sensitivity_rows <- res[res$subset_name != "FULL_SHARED", , drop = FALSE]

  direction_all <- all(sensitivity_rows$contrast_direction_preserved_vs_full == "YES")
  sign_all <- all(sensitivity_rows$component_sign_pattern_preserved == "YES")
  nominal_all <- all(sensitivity_rows$nominal_contrast_p_lt_0_05 == "YES")
  max_atten <- max(sensitivity_rows$absolute_percent_attenuation_vs_full, na.rm = TRUE)

  summary <- data.frame(
    metric = c(
      "n_shared_snps_full",
      "n_runnable_sensitivity_subsets",
      "contrast_direction_preserved_all_sensitivity_subsets",
      "component_sign_pattern_preserved_all_sensitivity_subsets",
      "nominal_contrast_p_lt_0_05_all_sensitivity_subsets",
      "max_absolute_percent_attenuation_vs_full",
      "full_shared_contrast_beta",
      "full_shared_contrast_se",
      "full_shared_contrast_p"
    ),
    value = c(
      length(shared_snps),
      nrow(sensitivity_rows),
      direction_all,
      sign_all,
      nominal_all,
      max_atten,
      full_row$beta_difference_IOP_minus_nonIOP,
      full_row$se_difference_r0,
      full_row$p_contrast_r0
    ),
    stringsAsFactors = FALSE
  )

  write.table(summary, summary_file, sep = "\t", quote = FALSE, row.names = FALSE)

  if (direction_all && sign_all) {
    write_status(
      "PASSED_INSTRUMENT_QUALITY_DIRECTION_STABLE",
      paste0(
        "SBP contrast direction and component sign pattern preserved across ",
        nrow(sensitivity_rows),
        " runnable instrument-quality subsets. Nominal p<0.05 all=",
        nominal_all,
        "."
      )
    )
  } else {
    write_status(
      "CHECK_REQUIRED",
      "SBP contrast direction or component sign pattern was not preserved in at least one instrument-quality subset."
    )
  }

}, error = function(e) {
  write_status("FAILED", paste("Error:", conditionMessage(e)))
  stop(e)
})
