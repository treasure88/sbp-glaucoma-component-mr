OUTDIR <- "../../37_iop_external_mr_results"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

INPUTS <- list(
  SBP = "../../36_iop_external_validation_inputs/mr_input/SBP__IOP.external_mr_input.tsv.gz",
  ART_STIFFNESS = "../../36_iop_external_validation_inputs/mr_input/ART_STIFFNESS__IOP.external_mr_input.tsv.gz"
)

RESULTS_OUT <- file.path(OUTDIR, "phase5_9E_IOP_external_mr_results.tsv")
HET_OUT <- file.path(OUTDIR, "phase5_9E_IOP_external_heterogeneity.tsv")
EGGER_OUT <- file.path(OUTDIR, "phase5_9E_IOP_external_egger_intercept.tsv")
SINGLE_OUT <- file.path(OUTDIR, "phase5_9E_IOP_external_single_snp_ratios.tsv.gz")
LOO_OUT <- file.path(OUTDIR, "phase5_9E_IOP_external_leave_one_out.tsv.gz")
STATUS_OUT <- file.path(OUTDIR, "phase5_9E_status.tsv")
RUNTIME_OUT <- file.path(OUTDIR, "phase5_9E_runtime_log.tsv")

safe_num <- function(x) suppressWarnings(as.numeric(x))

direction_label <- function(beta) {
  if (is.na(beta)) return("NA")
  if (beta > 0) return("positive")
  if (beta < 0) return("negative")
  "zero"
}

p_from_z <- function(beta, se) {
  if (is.na(beta) || is.na(se) || se <= 0) return(NA_real_)
  2 * pnorm(abs(beta / se), lower.tail = FALSE)
}

read_mr_input <- function(path) {
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
  dat$pval_exposure <- safe_num(dat$pval_exposure)
  dat$pval_outcome <- safe_num(dat$pval_outcome)

  dat <- dat[
    is.finite(dat$beta_exposure) &
      is.finite(dat$se_exposure) &
      is.finite(dat$beta_outcome) &
      is.finite(dat$se_outcome) &
      dat$se_exposure > 0 &
      dat$se_outcome > 0 &
      dat$beta_exposure != 0,
  ]

  dat
}

ivw_result <- function(dat, exposure_id, random = FALSE) {
  bx <- dat$beta_exposure
  by <- dat$beta_outcome
  sy <- dat$se_outcome
  w <- 1 / sy^2

  n <- length(bx)

  beta <- sum(w * bx * by) / sum(w * bx^2)
  se_fixed <- sqrt(1 / sum(w * bx^2))

  residual <- by - beta * bx
  Q <- sum(w * residual^2)
  Q_df <- n - 1
  Q_pval <- ifelse(Q_df > 0, pchisq(Q, df = Q_df, lower.tail = FALSE), NA_real_)
  phi <- ifelse(Q_df > 0, max(1, Q / Q_df), 1)

  se <- if (random) se_fixed * sqrt(phi) else se_fixed
  pval <- p_from_z(beta, se)

  data.frame(
    exposure_id = exposure_id,
    outcome_id = "IOP",
    outcome_scale = "continuous_IOP",
    method = ifelse(random, "IVW_multiplicative_random_effects", "IVW_fixed_effects"),
    n_instruments = n,
    beta = beta,
    se = se,
    pval = pval,
    qval_bh_ivw_random = NA_real_,
    direction = direction_label(beta),
    Q = Q,
    Q_df = Q_df,
    Q_pval = Q_pval,
    phi = ifelse(random, phi, 1),
    note = ifelse(random, "Primary external IOP MR method", "External IOP MR; fixed effects"),
    stringsAsFactors = FALSE
  )
}

weighted_median <- function(dat, exposure_id) {
  bx <- dat$beta_exposure
  by <- dat$beta_outcome
  sy <- dat$se_outcome

  ratio <- by / bx
  ratio_se <- abs(sy / bx)
  w <- 1 / ratio_se^2

  ok <- is.finite(ratio) & is.finite(w) & w > 0
  ratio <- ratio[ok]
  w <- w[ok]

  if (length(ratio) == 0) {
    beta <- NA_real_
  } else {
    ord <- order(ratio)
    ratio <- ratio[ord]
    w <- w[ord]
    cw <- cumsum(w) / sum(w)
    beta <- ratio[which(cw >= 0.5)[1]]
  }

  data.frame(
    exposure_id = exposure_id,
    outcome_id = "IOP",
    outcome_scale = "continuous_IOP",
    method = "weighted_median_ratio_descriptive",
    n_instruments = nrow(dat),
    beta = beta,
    se = NA_real_,
    pval = NA_real_,
    qval_bh_ivw_random = NA_real_,
    direction = direction_label(beta),
    Q = NA_real_,
    Q_df = NA_real_,
    Q_pval = NA_real_,
    phi = NA_real_,
    note = "Descriptive weighted median of SNP ratio estimates",
    stringsAsFactors = FALSE
  )
}

egger_result <- function(dat, exposure_id) {
  bx <- dat$beta_exposure
  by <- dat$beta_outcome
  sy <- dat$se_outcome
  w <- 1 / sy^2
  n <- length(bx)

  if (n < 3) {
    return(list(
      slope = data.frame(
        exposure_id = exposure_id,
        outcome_id = "IOP",
        outcome_scale = "continuous_IOP",
        method = "MR_Egger_slope",
        n_instruments = n,
        beta = NA_real_,
        se = NA_real_,
        pval = NA_real_,
        qval_bh_ivw_random = NA_real_,
        direction = "NA",
        Q = NA_real_,
        Q_df = NA_real_,
        Q_pval = NA_real_,
        phi = NA_real_,
        note = "MR-Egger not run; fewer than 3 instruments",
        stringsAsFactors = FALSE
      ),
      intercept = data.frame(
        exposure_id = exposure_id,
        outcome_id = "IOP",
        outcome_scale = "continuous_IOP",
        n_instruments = n,
        egger_intercept = NA_real_,
        egger_intercept_se = NA_real_,
        egger_intercept_pval = NA_real_,
        pleiotropy_flag = "NOT_RUN_N_LT_3",
        egger_slope = NA_real_,
        egger_slope_se = NA_real_,
        egger_slope_pval = NA_real_,
        note = "MR-Egger not run; fewer than 3 instruments",
        stringsAsFactors = FALSE
      )
    ))
  }

  fit <- lm(by ~ bx, weights = w)
  sm <- summary(fit)

  intercept <- coef(sm)[1, 1]
  intercept_se <- coef(sm)[1, 2]
  intercept_p <- coef(sm)[1, 4]

  slope <- coef(sm)[2, 1]
  slope_se <- coef(sm)[2, 2]
  slope_p <- coef(sm)[2, 4]

  residual <- residuals(fit)
  Q <- sum(w * residual^2)
  Q_df <- n - 2
  Q_pval <- ifelse(Q_df > 0, pchisq(Q, df = Q_df, lower.tail = FALSE), NA_real_)
  phi <- ifelse(Q_df > 0, max(1, Q / Q_df), 1)

  list(
    slope = data.frame(
      exposure_id = exposure_id,
      outcome_id = "IOP",
      outcome_scale = "continuous_IOP",
      method = "MR_Egger_slope",
      n_instruments = n,
      beta = slope,
      se = slope_se,
      pval = slope_p,
      qval_bh_ivw_random = NA_real_,
      direction = direction_label(slope),
      Q = Q,
      Q_df = Q_df,
      Q_pval = Q_pval,
      phi = phi,
      note = "MR-Egger slope; low power when instrument count is small",
      stringsAsFactors = FALSE
    ),
    intercept = data.frame(
      exposure_id = exposure_id,
      outcome_id = "IOP",
      outcome_scale = "continuous_IOP",
      n_instruments = n,
      egger_intercept = intercept,
      egger_intercept_se = intercept_se,
      egger_intercept_pval = intercept_p,
      pleiotropy_flag = ifelse(!is.na(intercept_p) && intercept_p < 0.05, "EGGER_INTERCEPT_P_LT_0.05", "NO_EGGER_INTERCEPT_EVIDENCE"),
      egger_slope = slope,
      egger_slope_se = slope_se,
      egger_slope_pval = slope_p,
      note = "Egger intercept test for directional pleiotropy",
      stringsAsFactors = FALSE
    )
  )
}

single_snp_rows <- function(dat, exposure_id) {
  bx <- dat$beta_exposure
  by <- dat$beta_outcome
  sy <- dat$se_outcome

  beta <- by / bx
  se <- abs(sy / bx)
  pval <- mapply(p_from_z, beta, se)

  data.frame(
    exposure_id = exposure_id,
    outcome_id = "IOP",
    outcome_scale = "continuous_IOP",
    SNP = dat$SNP,
    beta_ratio = beta,
    se_ratio = se,
    pval_ratio = pval,
    direction = sapply(beta, direction_label),
    beta_exposure = dat$beta_exposure,
    se_exposure = dat$se_exposure,
    beta_outcome = dat$beta_outcome,
    se_outcome = dat$se_outcome,
    harmonization_action = dat$harmonization_action,
    palindromic_possible = dat$palindromic_possible,
    outcome_n = dat$outcome_n,
    note = "Single-SNP Wald ratio",
    stringsAsFactors = FALSE
  )
}

leave_one_out_rows <- function(dat, exposure_id) {
  out <- list()

  if (nrow(dat) < 3) {
    return(data.frame(
      exposure_id = character(),
      outcome_id = character(),
      outcome_scale = character(),
      omitted_SNP = character(),
      n_instruments = integer(),
      beta = numeric(),
      se = numeric(),
      pval = numeric(),
      direction = character(),
      Q = numeric(),
      Q_df = numeric(),
      Q_pval = numeric(),
      phi = numeric(),
      note = character(),
      stringsAsFactors = FALSE
    ))
  }

  for (i in seq_len(nrow(dat))) {
    sub <- dat[-i, ]
    r <- ivw_result(sub, exposure_id, random = TRUE)
    out[[length(out) + 1]] <- data.frame(
      exposure_id = exposure_id,
      outcome_id = "IOP",
      outcome_scale = "continuous_IOP",
      omitted_SNP = dat$SNP[i],
      n_instruments = nrow(sub),
      beta = r$beta,
      se = r$se,
      pval = r$pval,
      direction = r$direction,
      Q = r$Q,
      Q_df = r$Q_df,
      Q_pval = r$Q_pval,
      phi = r$phi,
      note = "Leave-one-out IVW multiplicative random-effects",
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, out)
}

start <- Sys.time()

results <- list()
heterogeneity <- list()
egger_intercepts <- list()
single_all <- list()
loo_all <- list()
runtime_rows <- list()

for (exposure_id in names(INPUTS)) {
  t0 <- Sys.time()
  dat <- read_mr_input(INPUTS[[exposure_id]])

  ivw_fixed <- ivw_result(dat, exposure_id, random = FALSE)
  ivw_random <- ivw_result(dat, exposure_id, random = TRUE)
  wm <- weighted_median(dat, exposure_id)
  eg <- egger_result(dat, exposure_id)

  results[[exposure_id]] <- rbind(ivw_fixed, ivw_random, wm, eg$slope)

  heterogeneity[[exposure_id]] <- data.frame(
    exposure_id = exposure_id,
    outcome_id = "IOP",
    outcome_scale = "continuous_IOP",
    method = "IVW",
    n_instruments = nrow(dat),
    Q = ivw_random$Q,
    Q_df = ivw_random$Q_df,
    Q_pval = ivw_random$Q_pval,
    phi = ivw_random$phi,
    heterogeneity_flag = ifelse(!is.na(ivw_random$Q_pval) && ivw_random$Q_pval < 0.05, "YES", "NO"),
    note = "Cochran Q heterogeneity for IVW",
    stringsAsFactors = FALSE
  )

  egger_intercepts[[exposure_id]] <- eg$intercept
  single_all[[exposure_id]] <- single_snp_rows(dat, exposure_id)
  loo_all[[exposure_id]] <- leave_one_out_rows(dat, exposure_id)

  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  runtime_rows[[exposure_id]] <- data.frame(
    exposure_id = exposure_id,
    outcome_id = "IOP",
    elapsed_seconds = round(elapsed, 3),
    elapsed_human = paste0(round(elapsed, 1), "s"),
    stringsAsFactors = FALSE
  )
}

res <- do.call(rbind, results)

idx <- res$method == "IVW_multiplicative_random_effects"
res$qval_bh_ivw_random[idx] <- p.adjust(res$pval[idx], method = "BH")

het <- do.call(rbind, heterogeneity)
egger <- do.call(rbind, egger_intercepts)
single <- do.call(rbind, single_all)
loo <- do.call(rbind, loo_all)
runtime <- do.call(rbind, runtime_rows)

write.table(res, RESULTS_OUT, sep = "\t", row.names = FALSE, quote = FALSE)
write.table(het, HET_OUT, sep = "\t", row.names = FALSE, quote = FALSE)
write.table(egger, EGGER_OUT, sep = "\t", row.names = FALSE, quote = FALSE)

gz1 <- gzfile(SINGLE_OUT, "wt")
write.table(single, gz1, sep = "\t", row.names = FALSE, quote = FALSE)
close(gz1)

gz2 <- gzfile(LOO_OUT, "wt")
write.table(loo, gz2, sep = "\t", row.names = FALSE, quote = FALSE)
close(gz2)

ivw_random <- res[res$method == "IVW_multiplicative_random_effects", ]
nominal <- sum(ivw_random$pval < 0.05, na.rm = TRUE)
fdr <- sum(ivw_random$qval_bh_ivw_random < 0.05, na.rm = TRUE)

sbp <- ivw_random[ivw_random$exposure_id == "SBP", ]
art <- ivw_random[ivw_random$exposure_id == "ART_STIFFNESS", ]

elapsed_all <- as.numeric(difftime(Sys.time(), start, units = "secs"))

con <- file(STATUS_OUT, "w")
writeLines("phase\tstatus\tkey_result\tnote", con)
writeLines(paste(
  "Phase 5.9E IOP external MR",
  "PASSED",
  paste0("ivw_random_tests=", nrow(ivw_random), ";nominal=", nominal, ";fdr=", fdr),
  "SBP and ART_STIFFNESS tested against measured IOP",
  sep = "\t"
), con)
writeLines(paste("outcome", "INFO", "IOP", "External measured intraocular pressure validation outcome", sep = "\t"), con)
writeLines(paste("outcome_scale", "DOCUMENTED", "continuous_IOP", "MR beta is effect on measured IOP scale", sep = "\t"), con)
writeLines(paste("SBP_IOP_direction", "INFO", ifelse(nrow(sbp) > 0, sbp$direction[1], "NA"), "Key route-A validation direction", sep = "\t"), con)
writeLines(paste("ART_STIFFNESS_IOP_direction", "INFO", ifelse(nrow(art) > 0, art$direction[1], "NA"), "Exploratory vascular-stiffness validation direction", sep = "\t"), con)
writeLines(paste("runtime", "INFO", paste0(round(elapsed_all, 3), "s"), "Phase 5.9E runtime", sep = "\t"), con)
close(con)

write.table(runtime, RUNTIME_OUT, sep = "\t", row.names = FALSE, quote = FALSE)

cat("===== Phase 5.9E completed =====\n")
cat("Wrote:", STATUS_OUT, "\n")
cat("Wrote:", RESULTS_OUT, "\n")
cat("Wrote:", HET_OUT, "\n")
cat("Wrote:", EGGER_OUT, "\n")
cat("Wrote:", SINGLE_OUT, "\n")
cat("Wrote:", LOO_OUT, "\n")
