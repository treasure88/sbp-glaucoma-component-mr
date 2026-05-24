# Phase 8.3B: SBP leave-one-chromosome-out component-contrast sensitivity
# Purpose: assess whether SBP IOP-minus-nonIOP contrast is driven by a single chromosome
# Method: base R IVW with multiplicative random-effects SE
# Language policy: English labels for statistical outputs

options(stringsAsFactors = FALSE)

out_dir <- "../../50_sbp_sensitivity"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

non_file <- "../../16_mr_input_datasets/pairwise/SBP__GBS_nonIOPcomponent.mr_input.tsv.gz"
iop_file <- "../../16_mr_input_datasets/pairwise/SBP__GBS_IOPcomponent.mr_input.tsv.gz"

results_file <- file.path(out_dir, "phase8_3B_leave_one_chromosome_out_results.tsv")
component_file <- file.path(out_dir, "phase8_3B_leave_one_chromosome_out_component_estimates.tsv")
reproduction_file <- file.path(out_dir, "phase8_3B_full_data_reproduction_check.tsv")
summary_file <- file.path(out_dir, "phase8_3B_leave_one_chromosome_out_summary.tsv")
status_file <- file.path(out_dir, "phase8_3B_leave_one_chromosome_out_status.tsv")

target <- list(
  beta_nonIOP = -0.014988432,
  se_nonIOP = 0.0058051821,
  beta_IOP = 0.0077208916,
  se_IOP = 0.0030961462,
  contrast_beta = 0.022709324,
  contrast_se = 0.0065792295,
  contrast_p = 0.00055713036
)

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
    phase = "Phase8.3B",
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

estimate_component <- function(d, component_id, excluded_chr) {
  dd <- d

  if (!identical(excluded_chr, "NONE")) {
    dd <- dd[as.character(dd$chr) != as.character(excluded_chr), , drop = FALSE]
  }

  est <- ivw_random(dd)

  data.frame(
    component = component_id,
    excluded_chr = excluded_chr,
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

  required <- c("SNP", "chr", "pos", "beta_exposure", "se_exposure", "beta_outcome", "se_outcome")

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

  non$chr <- as.character(non$chr)
  iop$chr <- as.character(iop$chr)

  chrs <- sort(unique(c(non$chr, iop$chr)))
  chrs <- chrs[!is.na(chrs) & chrs != ""]

  suppressWarnings({
    chr_numeric <- as.numeric(chrs)
  })

  if (all(!is.na(chr_numeric))) {
    chrs <- chrs[order(chr_numeric)]
  }

  exclusion_set <- c("NONE", chrs)

  component_rows <- list()

  for (excl in exclusion_set) {
    component_rows[[length(component_rows) + 1]] <- estimate_component(non, "GBS_nonIOPcomponent", excl)
    component_rows[[length(component_rows) + 1]] <- estimate_component(iop, "GBS_IOPcomponent", excl)
  }

  comp <- do.call(rbind, component_rows)
  write.table(comp, component_file, sep = "\t", quote = FALSE, row.names = FALSE)

  results_rows <- list()

  full_contrast_beta <- NA_real_

  for (excl in exclusion_set) {
    a <- comp[comp$excluded_chr == excl & comp$component == "GBS_nonIOPcomponent", , drop = FALSE]
    b <- comp[comp$excluded_chr == excl & comp$component == "GBS_IOPcomponent", , drop = FALSE]

    beta_diff <- b$beta - a$beta
    se_diff <- sqrt(b$se^2 + a$se^2)
    z <- beta_diff / se_diff
    p <- 2 * pnorm(-abs(z))

    if (excl == "NONE") {
      full_contrast_beta <- beta_diff
    }

    results_rows[[length(results_rows) + 1]] <- data.frame(
      exposure_id = "SBP",
      excluded_chr = excl,
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
      stringsAsFactors = FALSE
    )
  }

  res <- do.call(rbind, results_rows)

  full_row <- res[res$excluded_chr == "NONE", , drop = FALSE]
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

  res$influential_chr_flag <- ifelse(
    res$excluded_chr != "NONE" &
      (
        res$contrast_direction_preserved_vs_full == "NO" |
          res$component_sign_pattern_preserved == "NO" |
          res$absolute_percent_attenuation_vs_full >= 50
      ),
    "YES",
    "NO"
  )

  write.table(res, results_file, sep = "\t", quote = FALSE, row.names = FALSE)

  # Reproduction check against stored primary values
  reproduction <- data.frame(
    quantity = c(
      "beta_nonIOP",
      "se_nonIOP",
      "beta_IOP",
      "se_IOP",
      "contrast_beta",
      "contrast_se",
      "contrast_p"
    ),
    target_value = c(
      target$beta_nonIOP,
      target$se_nonIOP,
      target$beta_IOP,
      target$se_IOP,
      target$contrast_beta,
      target$contrast_se,
      target$contrast_p
    ),
    recomputed_value = c(
      full_row$beta_nonIOP,
      full_row$se_nonIOP,
      full_row$beta_IOP,
      full_row$se_IOP,
      full_row$beta_difference_IOP_minus_nonIOP,
      full_row$se_difference_r0,
      full_row$p_contrast_r0
    ),
    stringsAsFactors = FALSE
  )

  reproduction$absolute_delta <- abs(reproduction$target_value - reproduction$recomputed_value)
  reproduction$relative_delta <- reproduction$absolute_delta / pmax(abs(reproduction$target_value), .Machine$double.eps)
  reproduction$match_status <- ifelse(reproduction$relative_delta < 1e-4 | reproduction$absolute_delta < 1e-8, "PASS", "CHECK")

  write.table(reproduction, reproduction_file, sep = "\t", quote = FALSE, row.names = FALSE)

  loco <- res[res$excluded_chr != "NONE", , drop = FALSE]

  n_exclusions <- nrow(loco)
  direction_all <- all(loco$contrast_direction_preserved_vs_full == "YES")
  sign_all <- all(loco$component_sign_pattern_preserved == "YES")
  nominal_all <- all(loco$nominal_contrast_p_lt_0_05 == "YES")
  n_influential <- sum(loco$influential_chr_flag == "YES", na.rm = TRUE)
  max_atten <- max(loco$absolute_percent_attenuation_vs_full, na.rm = TRUE)

  min_p <- min(loco$p_contrast_r0, na.rm = TRUE)
  max_p <- max(loco$p_contrast_r0, na.rm = TRUE)

  summary <- data.frame(
    metric = c(
      "n_leave_one_chromosome_tests",
      "contrast_direction_preserved_all",
      "component_sign_pattern_preserved_all",
      "nominal_contrast_p_lt_0_05_all",
      "n_influential_chromosome_flags",
      "max_absolute_percent_attenuation_vs_full",
      "min_leave_one_chr_p_contrast",
      "max_leave_one_chr_p_contrast",
      "full_data_contrast_beta",
      "full_data_contrast_se",
      "full_data_contrast_p",
      "full_data_n_nonIOP",
      "full_data_n_IOP"
    ),
    value = c(
      n_exclusions,
      direction_all,
      sign_all,
      nominal_all,
      n_influential,
      max_atten,
      min_p,
      max_p,
      full_row$beta_difference_IOP_minus_nonIOP,
      full_row$se_difference_r0,
      full_row$p_contrast_r0,
      full_row$n_nonIOP,
      full_row$n_IOP
    ),
    stringsAsFactors = FALSE
  )

  write.table(summary, summary_file, sep = "\t", quote = FALSE, row.names = FALSE)

  reproduction_pass <- all(reproduction$match_status == "PASS")

  if (reproduction_pass && direction_all) {
    write_status(
      "PASSED_LOCO_DIRECTION_STABLE",
      paste0(
        "Full-data MR estimates reproduced; SBP contrast direction preserved in all ",
        n_exclusions,
        " leave-one-chromosome analyses. Nominal p<0.05 all=",
        nominal_all,
        "; influential chromosome flags=",
        n_influential,
        "."
      )
    )
  } else if (!reproduction_pass && direction_all) {
    write_status(
      "PASSED_DIRECTION_STABLE_REPRODUCTION_CHECK",
      "SBP contrast direction preserved, but full-data estimates differ from stored target values. Review reproduction check before manuscript use."
    )
  } else {
    write_status(
      "CHECK_REQUIRED",
      "SBP contrast direction was not preserved in at least one leave-one-chromosome analysis or reproduction failed. Review outputs."
    )
  }

}, error = function(e) {
  write_status("FAILED", paste("Error:", conditionMessage(e)))
  stop(e)
})
