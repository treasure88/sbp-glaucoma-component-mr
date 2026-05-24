suppressPackageStartupMessages({
  has_mrpresso <- requireNamespace("MRPRESSO", quietly = TRUE)
  has_mrraps <- requireNamespace("mr.raps", quietly = TRUE)
})

OUTDIR <- "../../22_sbp_robustness/phase4_2B_formal_packages"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

STATUS_OUT <- file.path(OUTDIR, "phase4_2B_status.tsv")
RUNTIME_OUT <- file.path(OUTDIR, "phase4_2B_runtime_log.tsv")
MRPRESSO_SUMMARY_OUT <- file.path(OUTDIR, "phase4_2B_mrpresso_summary.tsv")
MRPRESSO_OUTLIER_OUT <- file.path(OUTDIR, "phase4_2B_mrpresso_outlier_table.tsv")
MRRAPS_OUT <- file.path(OUTDIR, "phase4_2B_mrraps_results.tsv")
CONTRAST_OUT <- file.path(OUTDIR, "phase4_2B_mrraps_component_contrast.tsv")

INPUTS <- list(
  GBS_nonIOPcomponent = "../../16_mr_input_datasets/pairwise/SBP__GBS_nonIOPcomponent.mr_input.tsv.gz",
  GBS_IOPcomponent = "../../16_mr_input_datasets/pairwise/SBP__GBS_IOPcomponent.mr_input.tsv.gz"
)

R_VALUES <- c(-0.25, 0, 0.25, 0.5)

now_str <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")

pick_col <- function(nms, candidates) {
  hit <- candidates[candidates %in% nms]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

normal_p <- function(z) {
  2 * pnorm(-abs(z))
}

fmt <- function(x) {
  if (length(x) == 0 || is.null(x) || is.na(x)) return("NA")
  if (is.numeric(x)) return(format(x, digits = 8, scientific = TRUE, trim = TRUE))
  as.character(x)
}

get_first <- function(x, candidates) {
  for (nm in candidates) {
    if (!is.null(x[[nm]])) return(x[[nm]])
  }
  NA
}

read_mr_input <- function(path, outcome_suffix) {
  d <- read.delim(gzfile(path), sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

  nms <- names(d)

  col_snp <- pick_col(nms, c("SNP", "rsid", "variant"))
  col_chr <- pick_col(nms, c("chr", "CHR", "chromosome"))
  col_bx <- pick_col(nms, c("beta_exposure", "beta.exposure", "exposure_beta", "beta_exp"))
  col_sx <- pick_col(nms, c("se_exposure", "se.exposure", "exposure_se", "se_exp"))
  col_by <- pick_col(nms, c("beta_outcome_harmonized", "beta_outcome", "beta.outcome", "outcome_beta", "beta_out"))
  col_sy <- pick_col(nms, c("se_outcome", "se.outcome", "outcome_se", "se_out"))

  required <- c(col_snp, col_bx, col_sx, col_by, col_sy)
  if (any(is.na(required))) {
    stop(paste("Missing required columns in", path, "header:", paste(nms, collapse = "|")))
  }

  out <- data.frame(
    outcome_suffix = outcome_suffix,
    SNP = d[[col_snp]],
    chr = if (!is.na(col_chr)) d[[col_chr]] else NA,
    BetaExposure = suppressWarnings(as.numeric(d[[col_bx]])),
    SdExposure = suppressWarnings(as.numeric(d[[col_sx]])),
    BetaOutcome = suppressWarnings(as.numeric(d[[col_by]])),
    SdOutcome = suppressWarnings(as.numeric(d[[col_sy]])),
    stringsAsFactors = FALSE
  )

  out <- out[is.finite(out$BetaExposure) &
             is.finite(out$SdExposure) &
             is.finite(out$BetaOutcome) &
             is.finite(out$SdOutcome) &
             out$BetaExposure != 0 &
             out$SdExposure > 0 &
             out$SdOutcome > 0, ]

  out
}

safe_write <- function(df, path) {
  write.table(df, path, sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
}

run_mrpresso_one <- function(dat, outcome_suffix) {
  if (!has_mrpresso) {
    return(list(
      summary = data.frame(
        exposure_id = "SBP",
        outcome_suffix = outcome_suffix,
        package = "MRPRESSO",
        status = "PACKAGE_NOT_AVAILABLE",
        n_instruments = nrow(dat),
        global_test_p = NA,
        raw_beta = NA,
        raw_pval = NA,
        outlier_corrected_beta = NA,
        outlier_corrected_pval = NA,
        n_outliers = NA,
        distortion_pval = NA,
        note = "MRPRESSO package not available",
        stringsAsFactors = FALSE
      ),
      outliers = data.frame()
    ))
  }

  mr_presso_fun <- getFromNamespace("mr_presso", "MRPRESSO")

  set.seed(20260517)

  cat("Running MR-PRESSO for", outcome_suffix, "with NbDistribution=200...\n"); flush.console()
  res <- tryCatch({
    mr_presso_fun(
      BetaOutcome = "BetaOutcome",
      BetaExposure = "BetaExposure",
      SdOutcome = "SdOutcome",
      SdExposure = "SdExposure",
      OUTLIERtest = TRUE,
      DISTORTIONtest = TRUE,
      data = dat,
      NbDistribution = 200,
      SignifThreshold = 0.05
    )
  }, error = function(e) e)

  saveRDS(res, file.path(OUTDIR, paste0("MRPRESSO_raw_", outcome_suffix, ".rds")))
  capture.output(str(res), file = file.path(OUTDIR, paste0("MRPRESSO_structure_", outcome_suffix, ".txt")))

  if (inherits(res, "error")) {
    return(list(
      summary = data.frame(
        exposure_id = "SBP",
        outcome_suffix = outcome_suffix,
        package = "MRPRESSO",
        status = "ERROR",
        n_instruments = nrow(dat),
        global_test_p = NA,
        raw_beta = NA,
        raw_pval = NA,
        outlier_corrected_beta = NA,
        outlier_corrected_pval = NA,
        n_outliers = NA,
        distortion_pval = NA,
        note = conditionMessage(res),
        stringsAsFactors = FALSE
      ),
      outliers = data.frame()
    ))
  }

  global_p <- NA
  raw_beta <- NA
  raw_p <- NA
  corrected_beta <- NA
  corrected_p <- NA
  distortion_p <- NA
  n_outliers <- NA

  rp <- res[["MR-PRESSO results"]]
  main <- res[["Main MR results"]]

  if (!is.null(rp[["Global Test"]])) {
    gt <- rp[["Global Test"]]
    global_p <- suppressWarnings(as.numeric(gt[grep("P", names(gt), ignore.case = TRUE)[1]]))
  }

  if (!is.null(main) && is.data.frame(main)) {
    beta_col <- grep("Causal|Estimate|Beta", names(main), ignore.case = TRUE, value = TRUE)[1]
    p_col <- grep("P", names(main), ignore.case = TRUE, value = TRUE)[1]

    raw_row <- grep("Raw", main[[1]], ignore.case = TRUE)
    corr_row <- grep("Outlier", main[[1]], ignore.case = TRUE)

    if (!is.na(beta_col) && length(raw_row) > 0) raw_beta <- suppressWarnings(as.numeric(main[raw_row[1], beta_col]))
    if (!is.na(p_col) && length(raw_row) > 0) raw_p <- suppressWarnings(as.numeric(main[raw_row[1], p_col]))
    if (!is.na(beta_col) && length(corr_row) > 0) corrected_beta <- suppressWarnings(as.numeric(main[corr_row[1], beta_col]))
    if (!is.na(p_col) && length(corr_row) > 0) corrected_p <- suppressWarnings(as.numeric(main[corr_row[1], p_col]))
  }

  outlier_df <- data.frame()
  if (!is.null(rp[["Outlier Test"]]) && is.data.frame(rp[["Outlier Test"]])) {
    outlier_df <- rp[["Outlier Test"]]
    outlier_df$exposure_id <- "SBP"
    outlier_df$outcome_suffix <- outcome_suffix
    if (nrow(outlier_df) == nrow(dat)) {
      outlier_df$SNP <- dat$SNP
      outlier_df$chr <- dat$chr
    }
    pcols <- grep("P", names(outlier_df), ignore.case = TRUE, value = TRUE)
    if (length(pcols) > 0) {
      pp <- suppressWarnings(as.numeric(outlier_df[[pcols[1]]]))
      n_outliers <- sum(is.finite(pp) & pp < 0.05, na.rm = TRUE)
    }
  }

  if (!is.null(rp[["Distortion Test"]])) {
    dt <- rp[["Distortion Test"]]
    if (is.list(dt)) {
      pnm <- grep("P", names(dt), ignore.case = TRUE, value = TRUE)[1]
      if (!is.na(pnm)) distortion_p <- suppressWarnings(as.numeric(dt[[pnm]]))
    }
  }

  list(
    summary = data.frame(
      exposure_id = "SBP",
      outcome_suffix = outcome_suffix,
      package = "MRPRESSO",
      status = "COMPLETED",
      n_instruments = nrow(dat),
      global_test_p = global_p,
      raw_beta = raw_beta,
      raw_pval = raw_p,
      outlier_corrected_beta = corrected_beta,
      outlier_corrected_pval = corrected_p,
      n_outliers = n_outliers,
      distortion_pval = distortion_p,
      note = "MR-PRESSO run with NbDistribution=1000; raw object saved as RDS",
      stringsAsFactors = FALSE
    ),
    outliers = outlier_df
  )
}

run_mrraps_one <- function(dat, outcome_suffix) {
  if (!has_mrraps) {
    return(data.frame(
      exposure_id = "SBP",
      outcome_suffix = outcome_suffix,
      package = "mr.raps",
      status = "PACKAGE_NOT_AVAILABLE",
      n_instruments = nrow(dat),
      beta = NA,
      se = NA,
      pval = NA,
      direction = NA,
      note = "mr.raps package not available",
      stringsAsFactors = FALSE
    ))
  }

  raps_fun <- getFromNamespace("mr.raps", "mr.raps")

  cat("Running MR-RAPS for", outcome_suffix, "...\n"); flush.console()
  res <- tryCatch({
    raps_fun(
      b_exp = dat$BetaExposure,
      b_out = dat$BetaOutcome,
      se_exp = dat$SdExposure,
      se_out = dat$SdOutcome,
      over.dispersion = TRUE,
      loss.function = "huber"
    )
  }, error = function(e) e)

  saveRDS(res, file.path(OUTDIR, paste0("mr_raps_raw_", outcome_suffix, ".rds")))
  capture.output(str(res), file = file.path(OUTDIR, paste0("mr_raps_structure_", outcome_suffix, ".txt")))

  if (inherits(res, "error")) {
    return(data.frame(
      exposure_id = "SBP",
      outcome_suffix = outcome_suffix,
      package = "mr.raps",
      status = "ERROR",
      n_instruments = nrow(dat),
      beta = NA,
      se = NA,
      pval = NA,
      direction = NA,
      note = conditionMessage(res),
      stringsAsFactors = FALSE
    ))
  }

  beta <- suppressWarnings(as.numeric(get_first(res, c("beta.hat", "beta", "estimate"))))
  se <- suppressWarnings(as.numeric(get_first(res, c("beta.se", "se", "std.error"))))
  pval <- suppressWarnings(as.numeric(get_first(res, c("beta.p.value", "p.value", "pval"))))

  if (!is.finite(pval) && is.finite(beta) && is.finite(se) && se > 0) {
    pval <- normal_p(beta / se)
  }

  data.frame(
    exposure_id = "SBP",
    outcome_suffix = outcome_suffix,
    package = "mr.raps",
    status = "COMPLETED",
    n_instruments = nrow(dat),
    beta = beta,
    se = se,
    pval = pval,
    direction = ifelse(beta > 0, "positive", ifelse(beta < 0, "negative", "zero")),
    note = "mr.raps with over.dispersion=TRUE and loss.function=huber",
    stringsAsFactors = FALSE
  )
}

contrast_rows_from_raps <- function(raps_results) {
  non <- raps_results[raps_results$outcome_suffix == "GBS_nonIOPcomponent", ]
  iop <- raps_results[raps_results$outcome_suffix == "GBS_IOPcomponent", ]

  if (nrow(non) == 0 || nrow(iop) == 0) return(data.frame())
  if (!is.finite(non$beta[1]) || !is.finite(non$se[1]) || !is.finite(iop$beta[1]) || !is.finite(iop$se[1])) return(data.frame())

  out <- lapply(R_VALUES, function(r) {
    beta_diff <- iop$beta[1] - non$beta[1]
    var_diff <- iop$se[1]^2 + non$se[1]^2 - 2 * r * iop$se[1] * non$se[1]
    se_diff <- sqrt(var_diff)
    z <- beta_diff / se_diff
    p <- normal_p(z)

    data.frame(
      exposure_id = "SBP",
      method = "mr.raps",
      assumed_component_correlation_r = r,
      beta_nonIOP = non$beta[1],
      se_nonIOP = non$se[1],
      beta_IOP = iop$beta[1],
      se_IOP = iop$se[1],
      beta_difference_IOP_minus_nonIOP = beta_diff,
      se_difference = se_diff,
      z_contrast = z,
      p_contrast = p,
      contrast_direction = ifelse(beta_diff > 0, "IOP_MORE_POSITIVE_THAN_NONIOP", "IOP_MORE_NEGATIVE_THAN_NONIOP"),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, out)
}

cat("===== Phase 4.2B formal MR-PRESSO and MR-RAPS =====\n")
cat("Start:", now_str(), "\n")
t0 <- Sys.time()

datasets <- lapply(names(INPUTS), function(outcome) {
  read_mr_input(INPUTS[[outcome]], outcome)
})
names(datasets) <- names(INPUTS)

mrpresso_results <- lapply(names(datasets), function(outcome) {
  run_mrpresso_one(datasets[[outcome]], outcome)
})

mrpresso_summary <- do.call(rbind, lapply(mrpresso_results, `[[`, "summary"))
mrpresso_outliers <- do.call(rbind, lapply(mrpresso_results, `[[`, "outliers"))

if (is.null(mrpresso_outliers) || nrow(mrpresso_outliers) == 0) {
  mrpresso_outliers <- data.frame(note = "No parseable MR-PRESSO outlier table or no outliers reported")
}

mrraps_results <- do.call(rbind, lapply(names(datasets), function(outcome) {
  run_mrraps_one(datasets[[outcome]], outcome)
}))

mrraps_contrast <- contrast_rows_from_raps(mrraps_results)
if (nrow(mrraps_contrast) == 0) {
  mrraps_contrast <- data.frame(note = "MR-RAPS contrast not available")
}

safe_write(mrpresso_summary, MRPRESSO_SUMMARY_OUT)
safe_write(mrpresso_outliers, MRPRESSO_OUTLIER_OUT)
safe_write(mrraps_results, MRRAPS_OUT)
safe_write(mrraps_contrast, CONTRAST_OUT)

elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

mrraps_non <- mrraps_results[mrraps_results$outcome_suffix == "GBS_nonIOPcomponent", ]
mrraps_iop <- mrraps_results[mrraps_results$outcome_suffix == "GBS_IOPcomponent", ]

mrraps_direction_call <- "NA"
if (nrow(mrraps_non) > 0 && nrow(mrraps_iop) > 0) {
  mrraps_direction_call <- paste0(
    "nonIOP=", mrraps_non$direction[1],
    ";IOP=", mrraps_iop$direction[1]
  )
}

con <- file(STATUS_OUT, "w")
writeLines("phase\tstatus\tkey_result\tnote", con)
writeLines(paste(
  "Phase 4.2B formal MR-PRESSO and MR-RAPS robustness",
  ifelse(has_mrpresso || has_mrraps, "COMPLETED_WITH_AVAILABLE_PACKAGES", "NO_PACKAGES_AVAILABLE"),
  "Formal package-based robustness attempted",
  "MR-PRESSO and mr.raps results written when packages are available",
  sep = "\t"
), con)
writeLines(paste("MRPRESSO_available", "INFO", has_mrpresso, "Package availability at runtime", sep = "\t"), con)
writeLines(paste("mr_raps_available", "INFO", has_mrraps, "Package availability at runtime", sep = "\t"), con)
writeLines(paste("MRPRESSO_summary_rows", "INFO", nrow(mrpresso_summary), "Rows in MR-PRESSO summary table", sep = "\t"), con)
writeLines(paste("MR_RAPS_direction_call", "INFO", mrraps_direction_call, "Expected supportive pattern is nonIOP=negative;IOP=positive", sep = "\t"), con)
writeLines(paste("runtime", "INFO", paste0(round(elapsed, 3), "s"), "Phase 4.2B runtime", sep = "\t"), con)
close(con)

runtime_con <- file(RUNTIME_OUT, "w")
writeLines("phase\tstart_time\tend_time\telapsed_seconds\telapsed_human", runtime_con)
writeLines(paste("Phase 4.2B", "NA", now_str(), round(elapsed, 6), paste0(round(elapsed, 1), "s"), sep = "\t"), runtime_con)
close(runtime_con)

cat("===== Phase 4.2B completed =====\n")
cat("End:", now_str(), "\n")
cat("Elapsed:", round(elapsed, 1), "s\n")
cat("Wrote:", STATUS_OUT, "\n")
cat("Wrote:", MRPRESSO_SUMMARY_OUT, "\n")
cat("Wrote:", MRPRESSO_OUTLIER_OUT, "\n")
cat("Wrote:", MRRAPS_OUT, "\n")
cat("Wrote:", CONTRAST_OUT, "\n")
