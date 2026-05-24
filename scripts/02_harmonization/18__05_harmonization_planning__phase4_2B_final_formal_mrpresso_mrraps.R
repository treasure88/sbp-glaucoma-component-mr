suppressPackageStartupMessages({
  has_mrpresso <- requireNamespace("MRPRESSO", quietly = TRUE)
  has_mrraps <- requireNamespace("mr.raps", quietly = TRUE)
})

OUTDIR <- "../../22_sbp_robustness/phase4_2B_formal_packages_final"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

STATUS_OUT <- file.path(OUTDIR, "phase4_2B_final_status.tsv")
MRPRESSO_SUMMARY_OUT <- file.path(OUTDIR, "phase4_2B_final_mrpresso_summary.tsv")
MRPRESSO_OUTLIER_OUT <- file.path(OUTDIR, "phase4_2B_final_mrpresso_outlier_table.tsv")
MRRAPS_OUT <- file.path(OUTDIR, "phase4_2B_final_mrraps_results.tsv")
CONTRAST_OUT <- file.path(OUTDIR, "phase4_2B_final_mrraps_component_contrast.tsv")
RUNTIME_OUT <- file.path(OUTDIR, "phase4_2B_final_runtime_log.tsv")

INPUTS <- list(
  GBS_nonIOPcomponent = "../../16_mr_input_datasets/pairwise/SBP__GBS_nonIOPcomponent.mr_input.tsv.gz",
  GBS_IOPcomponent = "../../16_mr_input_datasets/pairwise/SBP__GBS_IOPcomponent.mr_input.tsv.gz"
)

R_VALUES <- c(-0.25, 0, 0.25, 0.5)

NB_DISTRIBUTION <- as.integer(Sys.getenv("MRPRESSO_NB", "10000"))
MRPRESSO_MODE <- Sys.getenv("MRPRESSO_MODE", "full")
# MRPRESSO_MODE:
#   full   = global + outlier + distortion test
#   global = global test only, faster

now_str <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")

safe_num <- function(x) {
  y <- suppressWarnings(as.numeric(x))
  if (length(y) == 0) return(NA_real_)
  if (all(is.na(y))) return(NA_real_)
  y[1]
}

p_from_z <- function(beta, se) {
  beta <- safe_num(beta)
  se <- safe_num(se)
  if (is.na(beta) || is.na(se) || se <= 0) return(NA_real_)
  2 * pnorm(abs(beta / se), lower.tail = FALSE)
}

dir_label <- function(beta) {
  beta <- safe_num(beta)
  if (is.na(beta)) return("NA")
  if (beta > 0) return("positive")
  if (beta < 0) return("negative")
  "zero"
}

read_input <- function(path) {
  dat <- read.delim(
    gzfile(path),
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  needed <- c("SNP", "beta_exposure", "se_exposure", "beta_outcome", "se_outcome")
  miss <- setdiff(needed, names(dat))
  if (length(miss) > 0) {
    stop("Missing required columns: ", paste(miss, collapse = ", "))
  }

  dat$beta_exposure <- safe_num(dat$beta_exposure)
  dat$se_exposure <- safe_num(dat$se_exposure)
  dat$beta_outcome <- safe_num(dat$beta_outcome)
  dat$se_outcome <- safe_num(dat$se_outcome)

  dat <- dat[
    is.finite(dat$beta_exposure) &
      is.finite(dat$se_exposure) &
      is.finite(dat$beta_outcome) &
      is.finite(dat$se_outcome) &
      dat$se_exposure > 0 &
      dat$se_outcome > 0,
  ]

  dat
}

extract_unlisted_num <- function(obj, patterns) {
  u <- tryCatch(unlist(obj, recursive = TRUE), error = function(e) NULL)
  if (is.null(u)) return(NA_real_)

  nms <- names(u)
  if (is.null(nms)) return(NA_real_)

  for (pat in patterns) {
    hit <- grep(pat, nms, ignore.case = TRUE)
    if (length(hit) > 0) {
      val <- safe_num(u[hit[1]])
      if (!is.na(val)) return(val)
    }
  }

  NA_real_
}

extract_from_df <- function(df, patterns) {
  if (!is.data.frame(df)) return(NA_real_)
  nms <- names(df)
  for (pat in patterns) {
    hit <- grep(pat, nms, ignore.case = TRUE)
    if (length(hit) > 0) {
      val <- safe_num(df[[hit[1]]][1])
      if (!is.na(val)) return(val)
    }
  }
  NA_real_
}

extract_mrpresso_global_p <- function(res) {
  p <- tryCatch(res[["MR-PRESSO results"]][["Global Test"]][["Pvalue"]], error = function(e) NA)
  p <- safe_num(p)
  if (!is.na(p)) return(p)

  p <- tryCatch(res[["MR-PRESSO results"]][["Global Test"]][["P-value"]], error = function(e) NA)
  p <- safe_num(p)
  if (!is.na(p)) return(p)

  extract_unlisted_num(res, c("Global.*Pvalue", "Global.*P.value", "Global.*P-value"))
}

extract_mrpresso_distortion_p <- function(res) {
  p <- tryCatch(res[["MR-PRESSO results"]][["Distortion Test"]][["Pvalue"]], error = function(e) NA)
  p <- safe_num(p)
  if (!is.na(p)) return(p)

  p <- tryCatch(res[["MR-PRESSO results"]][["Distortion Test"]][["P-value"]], error = function(e) NA)
  p <- safe_num(p)
  if (!is.na(p)) return(p)

  extract_unlisted_num(res, c("Distortion.*Pvalue", "Distortion.*P.value", "Distortion.*P-value"))
}

extract_mrpresso_main <- function(res) {
  main <- tryCatch(res[["Main MR results"]], error = function(e) NULL)

  out <- list(
    raw_beta = NA_real_,
    raw_se = NA_real_,
    raw_pval = NA_real_,
    corrected_beta = NA_real_,
    corrected_se = NA_real_,
    corrected_pval = NA_real_
  )

  if (!is.data.frame(main)) {
    return(out)
  }

  nms <- names(main)
  method_col <- NA_character_

  possible_method_cols <- c("MR Analysis", "MR.Analysis", "Method", "method")
  for (cc in possible_method_cols) {
    if (cc %in% nms) method_col <- cc
  }

  beta_col <- nms[grep("Causal|Estimate|Beta|b", nms, ignore.case = TRUE)][1]
  se_col <- nms[grep("^Sd$|Std|SE|se", nms, ignore.case = TRUE)][1]
  p_col <- nms[grep("Pvalue|P.value|P-value|P", nms, ignore.case = TRUE)][1]

  if (is.na(beta_col)) beta_col <- NA_character_
  if (is.na(se_col)) se_col <- NA_character_
  if (is.na(p_col)) p_col <- NA_character_

  raw_idx <- 1
  corrected_idx <- NA_integer_

  if (!is.na(method_col)) {
    method_values <- as.character(main[[method_col]])
    raw_hits <- grep("Raw", method_values, ignore.case = TRUE)
    corr_hits <- grep("Outlier|Correct", method_values, ignore.case = TRUE)
    if (length(raw_hits) > 0) raw_idx <- raw_hits[1]
    if (length(corr_hits) > 0) corrected_idx <- corr_hits[1]
  } else {
    if (nrow(main) >= 2) corrected_idx <- 2
  }

  if (!is.na(beta_col)) out$raw_beta <- safe_num(main[[beta_col]][raw_idx])
  if (!is.na(se_col)) out$raw_se <- safe_num(main[[se_col]][raw_idx])
  if (!is.na(p_col)) out$raw_pval <- safe_num(main[[p_col]][raw_idx])

  if (!is.na(corrected_idx)) {
    if (!is.na(beta_col)) out$corrected_beta <- safe_num(main[[beta_col]][corrected_idx])
    if (!is.na(se_col)) out$corrected_se <- safe_num(main[[se_col]][corrected_idx])
    if (!is.na(p_col)) out$corrected_pval <- safe_num(main[[p_col]][corrected_idx])
  }

  out
}

extract_mrpresso_outliers <- function(res, outcome_suffix) {
  tab <- tryCatch(res[["MR-PRESSO results"]][["Outlier Test"]], error = function(e) NULL)
  if (!is.data.frame(tab)) {
    return(data.frame(
      exposure_id = character(),
      outcome_suffix = character(),
      outlier_row = character(),
      stringsAsFactors = FALSE
    ))
  }

  tab2 <- tab
  tab2$exposure_id <- "SBP"
  tab2$outcome_suffix <- outcome_suffix
  tab2$outlier_row <- rownames(tab2)
  rownames(tab2) <- NULL

  # Put identifiers first.
  tab2 <- tab2[, c("exposure_id", "outcome_suffix", "outlier_row", setdiff(names(tab2), c("exposure_id", "outcome_suffix", "outlier_row"))), drop = FALSE]
  tab2
}

run_mrpresso_one <- function(outcome_suffix, input_dat) {
  raw_file <- file.path(OUTDIR, paste0("MRPRESSO_raw_", outcome_suffix, ".rds"))
  struct_file <- file.path(OUTDIR, paste0("MRPRESSO_structure_", outcome_suffix, ".txt"))

  if (!has_mrpresso) {
    return(list(
      summary = data.frame(
        exposure_id = "SBP",
        outcome_suffix = outcome_suffix,
        package = "MRPRESSO",
        status = "PACKAGE_NOT_AVAILABLE",
        n_instruments = nrow(input_dat),
        mode = MRPRESSO_MODE,
        nb_distribution = NB_DISTRIBUTION,
        global_test_p = NA,
        raw_beta = NA,
        raw_se = NA,
        raw_pval = NA,
        outlier_corrected_beta = NA,
        outlier_corrected_se = NA,
        outlier_corrected_pval = NA,
        n_outliers = NA,
        distortion_pval = NA,
        n_warnings = NA,
        note = "MRPRESSO package not available",
        stringsAsFactors = FALSE
      ),
      outliers = data.frame()
    ))
  }

  mr_presso_fun <- getFromNamespace("mr_presso", "MRPRESSO")

  dat <- data.frame(
    BetaExposure = input_dat$beta_exposure,
    BetaOutcome = input_dat$beta_outcome,
    SdExposure = input_dat$se_exposure,
    SdOutcome = input_dat$se_outcome,
    SNP = input_dat$SNP,
    stringsAsFactors = FALSE
  )

  dat <- dat[
    is.finite(dat$BetaExposure) &
      is.finite(dat$BetaOutcome) &
      is.finite(dat$SdExposure) &
      is.finite(dat$SdOutcome) &
      dat$SdExposure > 0 &
      dat$SdOutcome > 0,
  ]

  outlier_test <- identical(tolower(MRPRESSO_MODE), "full")
  distortion_test <- identical(tolower(MRPRESSO_MODE), "full")

  if (!file.exists(raw_file)) {
    cat("Running MR-PRESSO for", outcome_suffix,
        "mode=", MRPRESSO_MODE,
        "NbDistribution=", NB_DISTRIBUTION, "\n")
    flush.console()

    warn_vec <- character()

    res <- withCallingHandlers(
      tryCatch(
        mr_presso_fun(
          BetaOutcome = "BetaOutcome",
          BetaExposure = "BetaExposure",
          SdOutcome = "SdOutcome",
          SdExposure = "SdExposure",
          data = dat,
          OUTLIERtest = outlier_test,
          DISTORTIONtest = distortion_test,
          NbDistribution = NB_DISTRIBUTION,
          SignifThreshold = 0.05,
          seed = 12345
        ),
        error = function(e) e
      ),
      warning = function(w) {
        warn_vec <<- c(warn_vec, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )

    saveRDS(
      list(
        result = res,
        warnings = warn_vec,
        mode = MRPRESSO_MODE,
        nb_distribution = NB_DISTRIBUTION,
        n_instruments = nrow(dat),
        time = now_str()
      ),
      raw_file
    )

    sink(struct_file)
    cat("MR-PRESSO structure for", outcome_suffix, "\n")
    cat("mode:", MRPRESSO_MODE, "\n")
    cat("NbDistribution:", NB_DISTRIBUTION, "\n")
    cat("warnings:\n")
    print(warn_vec)
    cat("\nresult structure:\n")
    print(str(res))
    sink()
  } else {
    cat("Skipping existing MR-PRESSO raw result for", outcome_suffix, "\n")
    flush.console()
  }

  obj <- readRDS(raw_file)
  res <- obj$result
  warn_vec <- obj$warnings

  if (inherits(res, "error")) {
    return(list(
      summary = data.frame(
        exposure_id = "SBP",
        outcome_suffix = outcome_suffix,
        package = "MRPRESSO",
        status = "ERROR",
        n_instruments = nrow(dat),
        mode = obj$mode,
        nb_distribution = obj$nb_distribution,
        global_test_p = NA,
        raw_beta = NA,
        raw_se = NA,
        raw_pval = NA,
        outlier_corrected_beta = NA,
        outlier_corrected_se = NA,
        outlier_corrected_pval = NA,
        n_outliers = NA,
        distortion_pval = NA,
        n_warnings = length(warn_vec),
        note = conditionMessage(res),
        stringsAsFactors = FALSE
      ),
      outliers = data.frame()
    ))
  }

  main <- extract_mrpresso_main(res)
  global_p <- extract_mrpresso_global_p(res)
  distortion_p <- extract_mrpresso_distortion_p(res)
  outlier_tab <- extract_mrpresso_outliers(res, outcome_suffix)
  n_out <- nrow(outlier_tab)

  list(
    summary = data.frame(
      exposure_id = "SBP",
      outcome_suffix = outcome_suffix,
      package = "MRPRESSO",
      status = "OK",
      n_instruments = nrow(dat),
      mode = obj$mode,
      nb_distribution = obj$nb_distribution,
      global_test_p = global_p,
      raw_beta = main$raw_beta,
      raw_se = main$raw_se,
      raw_pval = main$raw_pval,
      outlier_corrected_beta = main$corrected_beta,
      outlier_corrected_se = main$corrected_se,
      outlier_corrected_pval = main$corrected_pval,
      n_outliers = n_out,
      distortion_pval = distortion_p,
      n_warnings = length(warn_vec),
      note = ifelse(length(warn_vec) > 0, paste(unique(warn_vec), collapse = " | "), "MR-PRESSO completed"),
      stringsAsFactors = FALSE
    ),
    outliers = outlier_tab
  )
}

extract_raps_field <- function(res, patterns) {
  if (is.data.frame(res)) {
    val <- extract_from_df(res, patterns)
    if (!is.na(val)) return(val)
  }

  u <- tryCatch(unlist(res, recursive = TRUE), error = function(e) NULL)
  if (is.null(u)) return(NA_real_)

  nms <- names(u)
  if (is.null(nms)) return(NA_real_)

  for (pat in patterns) {
    hit <- grep(pat, nms, ignore.case = TRUE)
    if (length(hit) > 0) {
      val <- safe_num(u[hit[1]])
      if (!is.na(val)) return(val)
    }
  }

  NA_real_
}

run_mrraps_one <- function(outcome_suffix, input_dat) {
  raw_file <- file.path(OUTDIR, paste0("MRRAPS_raw_", outcome_suffix, ".rds"))
  struct_file <- file.path(OUTDIR, paste0("MRRAPS_structure_", outcome_suffix, ".txt"))

  if (!has_mrraps) {
    return(data.frame(
      exposure_id = "SBP",
      outcome_suffix = outcome_suffix,
      package = "mr.raps",
      status = "PACKAGE_NOT_AVAILABLE",
      n_instruments = nrow(input_dat),
      beta = NA,
      se = NA,
      pval = NA,
      direction = "NA",
      over_dispersion = "TRUE",
      loss_function = "huber",
      n_warnings = NA,
      note = "mr.raps package not available",
      stringsAsFactors = FALSE
    ))
  }

  raps_fun <- getFromNamespace("mr.raps", "mr.raps")

  dat <- data.frame(
    beta.exposure = input_dat$beta_exposure,
    beta.outcome = input_dat$beta_outcome,
    se.exposure = input_dat$se_exposure,
    se.outcome = input_dat$se_outcome,
    stringsAsFactors = FALSE
  )

  dat <- dat[
    is.finite(dat$beta.exposure) &
      is.finite(dat$beta.outcome) &
      is.finite(dat$se.exposure) &
      is.finite(dat$se.outcome) &
      dat$se.exposure > 0 &
      dat$se.outcome > 0,
  ]

  if (!file.exists(raw_file)) {
    cat("Running MR-RAPS for", outcome_suffix, "\n")
    flush.console()

    warn_vec <- character()

    res <- withCallingHandlers(
      tryCatch(
        raps_fun(
          dat,
          diagnostics = FALSE,
          over.dispersion = TRUE,
          loss.function = "huber"
        ),
        error = function(e) e
      ),
      warning = function(w) {
        warn_vec <<- c(warn_vec, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )

    saveRDS(
      list(
        result = res,
        warnings = warn_vec,
        n_instruments = nrow(dat),
        time = now_str()
      ),
      raw_file
    )

    sink(struct_file)
    cat("MR-RAPS structure for", outcome_suffix, "\n")
    cat("warnings:\n")
    print(warn_vec)
    cat("\nresult structure:\n")
    print(str(res))
    sink()
  } else {
    cat("Skipping existing MR-RAPS raw result for", outcome_suffix, "\n")
    flush.console()
  }

  obj <- readRDS(raw_file)
  res <- obj$result
  warn_vec <- obj$warnings

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
      direction = "NA",
      over_dispersion = "TRUE",
      loss_function = "huber",
      n_warnings = length(warn_vec),
      note = conditionMessage(res),
      stringsAsFactors = FALSE
    ))
  }

  beta <- extract_raps_field(res, c("^beta\\.hat$", "beta.hat", "beta_hat", "estimate", "Estimate", "^beta$"))
  se <- extract_raps_field(res, c("beta\\.se", "se\\.beta", "se_beta", "std", "Std", "^se$"))
  pval <- extract_raps_field(res, c("p\\.value", "pvalue", "pval", "P\\.value"))

  if (is.na(pval) && !is.na(beta) && !is.na(se) && se > 0) {
    pval <- p_from_z(beta, se)
  }

  data.frame(
    exposure_id = "SBP",
    outcome_suffix = outcome_suffix,
    package = "mr.raps",
    status = "OK",
    n_instruments = nrow(dat),
    beta = beta,
    se = se,
    pval = pval,
    direction = dir_label(beta),
    over_dispersion = "TRUE",
    loss_function = "huber",
    n_warnings = length(warn_vec),
    note = ifelse(length(warn_vec) > 0, paste(unique(warn_vec), collapse = " | "), "MR-RAPS completed"),
    stringsAsFactors = FALSE
  )
}

component_contrast_from_mrraps <- function(mrraps_rows) {
  non <- mrraps_rows[mrraps_rows$outcome_suffix == "GBS_nonIOPcomponent" & mrraps_rows$status == "OK", ]
  iop <- mrraps_rows[mrraps_rows$outcome_suffix == "GBS_IOPcomponent" & mrraps_rows$status == "OK", ]

  if (nrow(non) == 0 || nrow(iop) == 0) {
    return(data.frame(
      exposure_id = character(),
      method = character(),
      assumed_component_correlation_r = numeric(),
      beta_nonIOP = numeric(),
      se_nonIOP = numeric(),
      beta_IOP = numeric(),
      se_IOP = numeric(),
      beta_difference_IOP_minus_nonIOP = numeric(),
      se_difference = numeric(),
      z_contrast = numeric(),
      p_contrast = numeric(),
      contrast_direction = character(),
      note = character(),
      stringsAsFactors = FALSE
    ))
  }

  out <- list()

  for (r in R_VALUES) {
    beta_non <- safe_num(non$beta[1])
    se_non <- safe_num(non$se[1])
    beta_iop <- safe_num(iop$beta[1])
    se_iop <- safe_num(iop$se[1])

    beta_diff <- beta_iop - beta_non
    se_diff <- sqrt(se_iop^2 + se_non^2 - 2 * r * se_iop * se_non)
    z <- beta_diff / se_diff
    p <- 2 * pnorm(abs(z), lower.tail = FALSE)

    out[[length(out) + 1]] <- data.frame(
      exposure_id = "SBP",
      method = "MR_RAPS",
      assumed_component_correlation_r = r,
      beta_nonIOP = beta_non,
      se_nonIOP = se_non,
      beta_IOP = beta_iop,
      se_IOP = se_iop,
      beta_difference_IOP_minus_nonIOP = beta_diff,
      se_difference = se_diff,
      z_contrast = z,
      p_contrast = p,
      contrast_direction = ifelse(beta_diff > 0, "IOP_MORE_POSITIVE_THAN_NONIOP", "IOP_MORE_NEGATIVE_THAN_NONIOP"),
      note = "MR-RAPS component contrast; covariance unknown; sensitivity over assumed r",
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, out)
}

write_empty_if_needed <- function(path, header) {
  if (!file.exists(path)) {
    con <- file(path, "w")
    writeLines(paste(header, collapse = "\t"), con)
    close(con)
  }
}

main <- function() {
  start <- Sys.time()

  cat("===== Phase 4.2B-final formal MR-PRESSO and MR-RAPS =====\n")
  cat("Start:", now_str(), "\n")
  cat("MRPRESSO available:", has_mrpresso, "\n")
  cat("mr.raps available:", has_mrraps, "\n")
  cat("MRPRESSO_MODE:", MRPRESSO_MODE, "\n")
  cat("MRPRESSO_NB:", NB_DISTRIBUTION, "\n")
  flush.console()

  input_data <- lapply(INPUTS, read_input)

  mrpresso_summaries <- list()
  mrpresso_outliers <- list()
  mrraps_rows <- list()

  for (outcome_suffix in names(input_data)) {
    dat <- input_data[[outcome_suffix]]

    mp <- run_mrpresso_one(outcome_suffix, dat)
    mrpresso_summaries[[outcome_suffix]] <- mp$summary
    mrpresso_outliers[[outcome_suffix]] <- mp$outliers

    rr <- run_mrraps_one(outcome_suffix, dat)
    mrraps_rows[[outcome_suffix]] <- rr
  }

  mp_sum <- do.call(rbind, mrpresso_summaries)
  rr_sum <- do.call(rbind, mrraps_rows)
  contrast <- component_contrast_from_mrraps(rr_sum)

  write.table(mp_sum, MRPRESSO_SUMMARY_OUT, sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(rr_sum, MRRAPS_OUT, sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(contrast, CONTRAST_OUT, sep = "\t", row.names = FALSE, quote = FALSE)

  outlier_all <- do.call(rbind, mrpresso_outliers)
  if (is.null(outlier_all) || nrow(outlier_all) == 0) {
    write.table(
      data.frame(
        exposure_id = character(),
        outcome_suffix = character(),
        outlier_row = character(),
        stringsAsFactors = FALSE
      ),
      MRPRESSO_OUTLIER_OUT,
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )
  } else {
    write.table(outlier_all, MRPRESSO_OUTLIER_OUT, sep = "\t", row.names = FALSE, quote = FALSE)
  }

  runtime <- as.numeric(difftime(Sys.time(), start, units = "secs"))

  mpr_ok <- sum(mp_sum$status == "OK")
  raps_ok <- sum(rr_sum$status == "OK")
  raps_non <- rr_sum[rr_sum$outcome_suffix == "GBS_nonIOPcomponent", ]
  raps_iop <- rr_sum[rr_sum$outcome_suffix == "GBS_IOPcomponent", ]

  raps_call <- paste0(
    "nonIOP=", ifelse(nrow(raps_non) > 0, raps_non$direction[1], "NA"),
    ";IOP=", ifelse(nrow(raps_iop) > 0, raps_iop$direction[1], "NA")
  )

  con <- file(STATUS_OUT, "w")
  writeLines("phase\tstatus\tkey_result\tnote", con)
  writeLines(paste(
    "Phase 4.2B-final formal MR-PRESSO and MR-RAPS",
    ifelse(raps_ok == 2 || mpr_ok > 0, "COMPLETED_WITH_RESULTS", "COMPLETED_WITH_ERRORS"),
    paste0("MRPRESSO_OK=", mpr_ok, "/2;MR_RAPS_OK=", raps_ok, "/2"),
    "Formal package-based robustness rerun with corrected MR-PRESSO NbDistribution and MR-RAPS input format",
    sep = "\t"
  ), con)
  writeLines(paste("MRPRESSO_available", "INFO", has_mrpresso, "Package availability at runtime", sep = "\t"), con)
  writeLines(paste("mr_raps_available", "INFO", has_mrraps, "Package availability at runtime", sep = "\t"), con)
  writeLines(paste("MRPRESSO_mode", "INFO", MRPRESSO_MODE, "full=global+outlier+distortion; global=global test only", sep = "\t"), con)
  writeLines(paste("MRPRESSO_NbDistribution", "INFO", NB_DISTRIBUTION, "Use >=10000 for 391 instruments", sep = "\t"), con)
  writeLines(paste("MR_RAPS_direction_call", "INFO", raps_call, "Expected supportive pattern is nonIOP=negative;IOP=positive", sep = "\t"), con)
  writeLines(paste("MR_RAPS_contrast_rows", "INFO", nrow(contrast), "Rows in MR-RAPS contrast table", sep = "\t"), con)
  writeLines(paste("runtime", "INFO", paste0(round(runtime, 3), "s"), "Phase 4.2B-final runtime", sep = "\t"), con)
  close(con)

  write.table(
    data.frame(
      phase = "Phase 4.2B-final",
      start_time = format(start, "%Y-%m-%d %H:%M:%S"),
      end_time = now_str(),
      elapsed_seconds = round(runtime, 3),
      elapsed_human = paste0(round(runtime, 1), "s"),
      stringsAsFactors = FALSE
    ),
    RUNTIME_OUT,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  cat("===== Phase 4.2B-final completed =====\n")
  cat("End:", now_str(), "\n")
  cat("Wrote:", STATUS_OUT, "\n")
  cat("Wrote:", MRPRESSO_SUMMARY_OUT, "\n")
  cat("Wrote:", MRPRESSO_OUTLIER_OUT, "\n")
  cat("Wrote:", MRRAPS_OUT, "\n")
  cat("Wrote:", CONTRAST_OUT, "\n")
}

main()
