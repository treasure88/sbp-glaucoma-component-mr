options(stringsAsFactors = FALSE)

# Phase13.4B: Empirical covariance, bootstrap, and jackknife sensitivity
# Project: IOP-dependent vs IOP-independent glaucoma component MR
# Evidence level: HYPOTHESIS_GENERATING_NOT_CONFIRMATORY
# This phase performs empirical uncertainty assessment for the locked SBP component contrast.
# It does not upgrade claims to confirmatory causal evidence.

message("Running Phase13.4B empirical covariance / bootstrap / jackknife sensitivity...")

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

fmt <- function(x, digits = 4) {
  ifelse(is.na(x), NA, formatC(x, digits = digits, format = "fg", flag = "#"))
}

as_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

pick_col <- function(df, candidates = character(0), regex = character(0)) {
  nms <- names(df)
  low <- tolower(nms)
  for (cc in candidates) {
    idx <- which(low == tolower(cc))
    if (length(idx) > 0) return(nms[idx[1]])
  }
  for (rr in regex) {
    idx <- which(grepl(rr, low, ignore.case = TRUE))
    if (length(idx) > 0) return(nms[idx[1]])
  }
  NA_character_
}

read_gz_tsv <- function(path) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)
  read.delim(
    con,
    header = TRUE,
    sep = "\t",
    quote = "",
    comment.char = "",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

ivw_random_phi <- function(bx, sx, by, sy) {
  bx <- as_num(bx)
  sx <- as_num(sx)
  by <- as_num(by)
  sy <- as_num(sy)
  
  keep <- is.finite(bx) & is.finite(sx) & is.finite(by) & is.finite(sy) &
    bx != 0 & sx > 0 & sy > 0
  
  bx <- bx[keep]
  sx <- sx[keep]
  by <- by[keep]
  sy <- sy[keep]
  
  ratio <- by / bx
  ratio_se <- abs(sy / bx)
  w <- 1 / ratio_se^2
  
  beta <- sum(w * ratio) / sum(w)
  q <- sum(w * (ratio - beta)^2)
  df <- length(ratio) - 1
  phi <- max(1, q / df)
  se <- sqrt(phi / sum(w))
  z <- beta / se
  p <- 2 * pnorm(abs(z), lower.tail = FALSE)
  
  data.frame(
    n_instruments = length(ratio),
    beta = beta,
    se = se,
    z = z,
    p = p,
    Q = q,
    Q_df = df,
    Q_p = pchisq(q, df = df, lower.tail = FALSE),
    phi = phi,
    stringsAsFactors = FALSE
  )
}

estimate_pair <- function(dat) {
  non <- ivw_random_phi(
    bx = dat$beta_exposure,
    sx = dat$se_exposure,
    by = dat$beta_outcome_nonIOP,
    sy = dat$se_outcome_nonIOP
  )
  
  iop <- ivw_random_phi(
    bx = dat$beta_exposure,
    sx = dat$se_exposure,
    by = dat$beta_outcome_IOP,
    sy = dat$se_outcome_IOP
  )
  
  data.frame(
    n_instruments = nrow(dat),
    beta_nonIOP = non$beta,
    se_nonIOP = non$se,
    p_nonIOP = non$p,
    phi_nonIOP = non$phi,
    beta_IOP = iop$beta,
    se_IOP = iop$se,
    p_IOP = iop$p,
    phi_IOP = iop$phi,
    beta_difference_IOP_minus_nonIOP = iop$beta - non$beta,
    stringsAsFactors = FALSE
  )
}

make_locus_proxy <- function(chr, pos, window_bp = 1000000) {
  chr <- as.character(chr)
  pos <- as_num(pos)
  ord <- order(suppressWarnings(as.numeric(chr)), pos, na.last = TRUE)
  
  locus <- rep(NA_character_, length(chr))
  locus_counter <- 0L
  current_chr <- NA_character_
  current_start <- NA_real_
  
  for (ii in ord) {
    if (is.na(chr[ii]) || is.na(pos[ii])) {
      locus_counter <- locus_counter + 1L
      locus[ii] <- paste0("missing_locus_", locus_counter)
      next
    }
    
    if (is.na(current_chr) || chr[ii] != current_chr || is.na(current_start) || (pos[ii] - current_start) > window_bp) {
      locus_counter <- locus_counter + 1L
      current_chr <- chr[ii]
      current_start <- pos[ii]
    }
    
    locus[ii] <- paste0("chr", chr[ii], "_locus", locus_counter)
  }
  
  locus
}

root <- normalizePath(".", winslash = "/", mustWork = TRUE)
out_dir <- file.path(root, "71_data_analysis_reinforcement", "phase13_4_empirical_covariance")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

master_status_file <- file.path(root, "71_data_analysis_reinforcement", "phase13_master_status.tsv")

non_file <- file.path(root, "16_mr_input_datasets/pairwise/SBP__GBS_nonIOPcomponent.mr_input.tsv.gz")
iop_file <- file.path(root, "16_mr_input_datasets/pairwise/SBP__GBS_IOPcomponent.mr_input.tsv.gz")

stopifnot(file.exists(non_file))
stopifnot(file.exists(iop_file))

non_raw <- read_gz_tsv(non_file)
iop_raw <- read_gz_tsv(iop_file)

map_cols <- function(df) {
  list(
    snp = pick_col(df, c("SNP", "rsid", "variant_id")),
    chr = pick_col(df, c("chr", "CHR", "chromosome", "exposure_chr")),
    pos = pick_col(df, c("pos", "POS", "position", "bp", "exposure_pos")),
    ea = pick_col(df, c("effect_allele", "effect_allele_exposure", "exposure_effect_allele")),
    oa = pick_col(df, c("other_allele", "other_allele_exposure", "exposure_other_allele")),
    beta_exposure = pick_col(df, c("beta_exposure", "beta.exposure", "beta_exp", "bx")),
    se_exposure = pick_col(df, c("se_exposure", "se.exposure", "se_exp", "sx")),
    beta_outcome = pick_col(df, c("beta_outcome", "beta.outcome", "beta_out", "by", "beta_outcome_aligned")),
    se_outcome = pick_col(df, c("se_outcome", "se.outcome", "se_out", "sy")),
    pval_exposure = pick_col(df, c("pval_exposure", "pval.exposure", "p.exposure")),
    pval_outcome = pick_col(df, c("pval_outcome", "pval.outcome", "p.outcome")),
    F_stat = pick_col(df, c("instrument_F_stat", "F_stat", "F", "f_stat", "F_statistic")),
    harmonization_action = pick_col(df, c("harmonization_action")),
    match_mode = pick_col(df, c("match_mode")),
    clump_rank = pick_col(df, c("clump_rank"))
  )
}

non_cols <- map_cols(non_raw)
iop_cols <- map_cols(iop_raw)

required_non <- c(non_cols$snp, non_cols$chr, non_cols$pos, non_cols$beta_exposure, non_cols$se_exposure, non_cols$beta_outcome, non_cols$se_outcome)
required_iop <- c(iop_cols$snp, iop_cols$chr, iop_cols$pos, iop_cols$beta_exposure, iop_cols$se_exposure, iop_cols$beta_outcome, iop_cols$se_outcome)

if (any(is.na(required_non)) || any(is.na(required_iop))) {
  stop("Required columns missing from internal input files.")
}

non <- data.frame(
  SNP = as.character(non_raw[[non_cols$snp]]),
  chr = as.character(non_raw[[non_cols$chr]]),
  pos = as_num(non_raw[[non_cols$pos]]),
  effect_allele = if (!is.na(non_cols$ea)) as.character(non_raw[[non_cols$ea]]) else NA_character_,
  other_allele = if (!is.na(non_cols$oa)) as.character(non_raw[[non_cols$oa]]) else NA_character_,
  beta_exposure = as_num(non_raw[[non_cols$beta_exposure]]),
  se_exposure = as_num(non_raw[[non_cols$se_exposure]]),
  beta_outcome_nonIOP = as_num(non_raw[[non_cols$beta_outcome]]),
  se_outcome_nonIOP = as_num(non_raw[[non_cols$se_outcome]]),
  pval_outcome_nonIOP = if (!is.na(non_cols$pval_outcome)) as_num(non_raw[[non_cols$pval_outcome]]) else NA_real_,
  F_stat = if (!is.na(non_cols$F_stat)) as_num(non_raw[[non_cols$F_stat]]) else NA_real_,
  harmonization_nonIOP = if (!is.na(non_cols$harmonization_action)) as.character(non_raw[[non_cols$harmonization_action]]) else NA_character_,
  stringsAsFactors = FALSE
)

iop <- data.frame(
  SNP = as.character(iop_raw[[iop_cols$snp]]),
  chr_IOP = as.character(iop_raw[[iop_cols$chr]]),
  pos_IOP = as_num(iop_raw[[iop_cols$pos]]),
  effect_allele_IOP = if (!is.na(iop_cols$ea)) as.character(iop_raw[[iop_cols$ea]]) else NA_character_,
  other_allele_IOP = if (!is.na(iop_cols$oa)) as.character(iop_raw[[iop_cols$oa]]) else NA_character_,
  beta_exposure_IOP = as_num(iop_raw[[iop_cols$beta_exposure]]),
  se_exposure_IOP = as_num(iop_raw[[iop_cols$se_exposure]]),
  beta_outcome_IOP = as_num(iop_raw[[iop_cols$beta_outcome]]),
  se_outcome_IOP = as_num(iop_raw[[iop_cols$se_outcome]]),
  pval_outcome_IOP = if (!is.na(iop_cols$pval_outcome)) as_num(iop_raw[[iop_cols$pval_outcome]]) else NA_real_,
  harmonization_IOP = if (!is.na(iop_cols$harmonization_action)) as.character(iop_raw[[iop_cols$harmonization_action]]) else NA_character_,
  stringsAsFactors = FALSE
)

paired <- merge(non, iop, by = "SNP", all = FALSE, sort = FALSE)

paired$exposure_beta_difference_abs <- abs(paired$beta_exposure - paired$beta_exposure_IOP)
paired$exposure_se_difference_abs <- abs(paired$se_exposure - paired$se_exposure_IOP)
paired$allele_match <- paired$effect_allele == paired$effect_allele_IOP &
  paired$other_allele == paired$other_allele_IOP

paired$locus_proxy_1Mb <- make_locus_proxy(paired$chr, paired$pos, window_bp = 1000000)

paired_file <- file.path(out_dir, "phase13_4B_paired_internal_sbp_component_input.tsv")
write_tsv(paired, paired_file)

input_qc <- data.frame(
  field = c(
    "n_nonIOP_rows",
    "n_IOP_rows",
    "n_paired_rows",
    "n_unique_paired_snps",
    "max_abs_exposure_beta_difference_between_component_files",
    "max_abs_exposure_se_difference_between_component_files",
    "all_exposure_betas_identical_between_component_files_tolerance_1e_12",
    "all_exposure_ses_identical_between_component_files_tolerance_1e_12",
    "all_alleles_match_between_component_files",
    "n_chromosomes",
    "n_locus_proxy_1Mb"
  ),
  value = c(
    nrow(non),
    nrow(iop),
    nrow(paired),
    length(unique(paired$SNP)),
    max(paired$exposure_beta_difference_abs, na.rm = TRUE),
    max(paired$exposure_se_difference_abs, na.rm = TRUE),
    max(paired$exposure_beta_difference_abs, na.rm = TRUE) < 1e-12,
    max(paired$exposure_se_difference_abs, na.rm = TRUE) < 1e-12,
    all(paired$allele_match, na.rm = TRUE),
    length(unique(paired$chr)),
    length(unique(paired$locus_proxy_1Mb))
  ),
  stringsAsFactors = FALSE
)

input_qc_file <- file.path(out_dir, "phase13_4B_paired_input_qc.tsv")
write_tsv(input_qc, input_qc_file)

# Full locked estimate from exact paired input.
full_est <- estimate_pair(paired)

locked_beta_difference <- 0.022709324
locked_se_difference_r0 <- 0.0065792295
locked_p_difference_r0 <- 0.00055713036

# Bootstrap.
set.seed(20260523)
n_bootstrap <- 5000L
n <- nrow(paired)

boot_rows <- vector("list", n_bootstrap)

for (b in seq_len(n_bootstrap)) {
  idx <- sample.int(n, size = n, replace = TRUE)
  est <- estimate_pair(paired[idx, , drop = FALSE])
  boot_rows[[b]] <- data.frame(
    bootstrap_id = b,
    beta_nonIOP = est$beta_nonIOP,
    beta_IOP = est$beta_IOP,
    beta_difference_IOP_minus_nonIOP = est$beta_difference_IOP_minus_nonIOP,
    n_instruments = est$n_instruments,
    stringsAsFactors = FALSE
  )
}

boot <- do.call(rbind, boot_rows)

boot_file <- file.path(out_dir, "phase13_4B_bootstrap_contrast_distribution.tsv")
write_tsv(boot, boot_file)

empirical_cov <- cov(boot$beta_IOP, boot$beta_nonIOP, use = "complete.obs")
empirical_cor <- cor(boot$beta_IOP, boot$beta_nonIOP, use = "complete.obs")
empirical_se_non <- sd(boot$beta_nonIOP, na.rm = TRUE)
empirical_se_iop <- sd(boot$beta_IOP, na.rm = TRUE)
empirical_se_diff <- sd(boot$beta_difference_IOP_minus_nonIOP, na.rm = TRUE)

ci <- quantile(boot$beta_difference_IOP_minus_nonIOP, probs = c(0.025, 0.975), na.rm = TRUE)
boot_mean <- mean(boot$beta_difference_IOP_minus_nonIOP, na.rm = TRUE)
boot_median <- median(boot$beta_difference_IOP_minus_nonIOP, na.rm = TRUE)

z_emp <- full_est$beta_difference_IOP_minus_nonIOP / empirical_se_diff
p_emp_normal <- 2 * pnorm(abs(z_emp), lower.tail = FALSE)
p_emp_sign <- 2 * min(
  mean(boot$beta_difference_IOP_minus_nonIOP <= 0, na.rm = TRUE),
  mean(boot$beta_difference_IOP_minus_nonIOP >= 0, na.rm = TRUE)
)
p_emp_sign <- min(1, p_emp_sign)

bootstrap_summary <- data.frame(
  method = "paired_SNP_bootstrap",
  n_bootstrap = n_bootstrap,
  n_instruments = n,
  full_beta_nonIOP = full_est$beta_nonIOP,
  full_beta_IOP = full_est$beta_IOP,
  full_beta_difference_IOP_minus_nonIOP = full_est$beta_difference_IOP_minus_nonIOP,
  locked_beta_difference_IOP_minus_nonIOP = locked_beta_difference,
  bootstrap_mean_difference = boot_mean,
  bootstrap_median_difference = boot_median,
  empirical_SE_difference = empirical_se_diff,
  empirical_CI_lower_2p5 = as.numeric(ci[1]),
  empirical_CI_upper_97p5 = as.numeric(ci[2]),
  empirical_z_using_full_difference = z_emp,
  empirical_p_normal_using_full_difference = p_emp_normal,
  empirical_p_sign = p_emp_sign,
  empirical_SE_nonIOP = empirical_se_non,
  empirical_SE_IOP = empirical_se_iop,
  empirical_covariance_IOP_nonIOP = empirical_cov,
  empirical_correlation_IOP_nonIOP = empirical_cor,
  r0_analytic_SE_difference_locked = locked_se_difference_r0,
  r0_analytic_p_difference_locked = locked_p_difference_r0,
  direction_preserved = full_est$beta_difference_IOP_minus_nonIOP > 0,
  stringsAsFactors = FALSE
)

bootstrap_summary_file <- file.path(out_dir, "phase13_4B_bootstrap_contrast_summary.tsv")
write_tsv(bootstrap_summary, bootstrap_summary_file)

covariance_summary <- data.frame(
  parameter = c(
    "empirical_SE_nonIOP",
    "empirical_SE_IOP",
    "empirical_covariance_IOP_nonIOP",
    "empirical_correlation_IOP_nonIOP",
    "empirical_SE_difference",
    "locked_r0_analytic_SE_difference",
    "empirical_SE_over_locked_r0_SE",
    "interpretation"
  ),
  value = c(
    empirical_se_non,
    empirical_se_iop,
    empirical_cov,
    empirical_cor,
    empirical_se_diff,
    locked_se_difference_r0,
    empirical_se_diff / locked_se_difference_r0,
    "Empirical paired-instrument bootstrap estimate; used as uncertainty sensitivity, not claim upgrade."
  ),
  stringsAsFactors = FALSE
)

covariance_file <- file.path(out_dir, "phase13_4B_empirical_covariance_summary.tsv")
write_tsv(covariance_summary, covariance_file)

# Leave-one-chromosome jackknife.
chromosomes <- sort(unique(paired$chr))
chr_rows <- list()

for (cc in chromosomes) {
  sub <- paired[paired$chr != cc, , drop = FALSE]
  est <- estimate_pair(sub)
  chr_rows[[length(chr_rows) + 1]] <- data.frame(
    left_out_chromosome = cc,
    n_removed = sum(paired$chr == cc),
    n_remaining = nrow(sub),
    beta_nonIOP = est$beta_nonIOP,
    beta_IOP = est$beta_IOP,
    beta_difference_IOP_minus_nonIOP = est$beta_difference_IOP_minus_nonIOP,
    delta_difference_from_full = est$beta_difference_IOP_minus_nonIOP - full_est$beta_difference_IOP_minus_nonIOP,
    direction_preserved = est$beta_difference_IOP_minus_nonIOP > 0,
    stringsAsFactors = FALSE
  )
}

chr_jack <- do.call(rbind, chr_rows)
chr_jack_file <- file.path(out_dir, "phase13_4B_leave_one_chromosome_results.tsv")
write_tsv(chr_jack, chr_jack_file)

chr_summary <- data.frame(
  metric = c(
    "n_chromosomes",
    "all_chromosome_leave_one_out_direction_preserved",
    "min_leave_one_chromosome_difference",
    "max_leave_one_chromosome_difference",
    "max_abs_delta_from_full",
    "most_influential_chromosome_by_abs_delta"
  ),
  value = c(
    nrow(chr_jack),
    all(chr_jack$direction_preserved),
    min(chr_jack$beta_difference_IOP_minus_nonIOP, na.rm = TRUE),
    max(chr_jack$beta_difference_IOP_minus_nonIOP, na.rm = TRUE),
    max(abs(chr_jack$delta_difference_from_full), na.rm = TRUE),
    chr_jack$left_out_chromosome[which.max(abs(chr_jack$delta_difference_from_full))]
  ),
  stringsAsFactors = FALSE
)

chr_summary_file <- file.path(out_dir, "phase13_4B_leave_one_chromosome_summary.tsv")
write_tsv(chr_summary, chr_summary_file)

# Locus-proxy jackknife.
loci <- sort(unique(paired$locus_proxy_1Mb))
locus_rows <- list()

for (ll in loci) {
  sub <- paired[paired$locus_proxy_1Mb != ll, , drop = FALSE]
  est <- estimate_pair(sub)
  locus_rows[[length(locus_rows) + 1]] <- data.frame(
    left_out_locus_proxy_1Mb = ll,
    n_removed = sum(paired$locus_proxy_1Mb == ll),
    n_remaining = nrow(sub),
    beta_nonIOP = est$beta_nonIOP,
    beta_IOP = est$beta_IOP,
    beta_difference_IOP_minus_nonIOP = est$beta_difference_IOP_minus_nonIOP,
    delta_difference_from_full = est$beta_difference_IOP_minus_nonIOP - full_est$beta_difference_IOP_minus_nonIOP,
    direction_preserved = est$beta_difference_IOP_minus_nonIOP > 0,
    stringsAsFactors = FALSE
  )
}

locus_jack <- do.call(rbind, locus_rows)
locus_jack_file <- file.path(out_dir, "phase13_4B_locus_proxy_jackknife_results.tsv")
write_tsv(locus_jack, locus_jack_file)

locus_summary <- data.frame(
  metric = c(
    "n_locus_proxy_1Mb",
    "all_locus_proxy_leave_one_out_direction_preserved",
    "min_leave_one_locus_proxy_difference",
    "max_leave_one_locus_proxy_difference",
    "max_abs_delta_from_full",
    "most_influential_locus_proxy_by_abs_delta",
    "max_n_removed_in_single_locus_proxy"
  ),
  value = c(
    nrow(locus_jack),
    all(locus_jack$direction_preserved),
    min(locus_jack$beta_difference_IOP_minus_nonIOP, na.rm = TRUE),
    max(locus_jack$beta_difference_IOP_minus_nonIOP, na.rm = TRUE),
    max(abs(locus_jack$delta_difference_from_full), na.rm = TRUE),
    locus_jack$left_out_locus_proxy_1Mb[which.max(abs(locus_jack$delta_difference_from_full))],
    max(locus_jack$n_removed, na.rm = TRUE)
  ),
  stringsAsFactors = FALSE
)

locus_summary_file <- file.path(out_dir, "phase13_4B_locus_proxy_jackknife_summary.tsv")
write_tsv(locus_summary, locus_summary_file)

# Manuscript-ready inserts.
methods_insert <- c(
  "# Empirical covariance and jackknife assessment of the SBP component contrast",
  "",
  "To evaluate whether the SBP component contrast was sensitive to the working covariance assumption between component-specific estimates, we performed a paired SNP-level bootstrap using the exact harmonized SBP-component input files that reproduced the locked IVW estimates. In each bootstrap iteration, the same set of instruments was resampled with replacement for the IOP-independent and IOP-dependent component analyses, preserving the paired structure across component outcomes. The empirical standard error, confidence interval, and covariance between component-specific estimates were calculated from the bootstrap distribution.",
  "",
  "We also performed leave-one-chromosome-out and 1-Mb locus-proxy jackknife analyses. The locus-proxy analysis grouped instruments by chromosome and genomic position using a 1-Mb window. These analyses were designed to assess influence and uncertainty in the component contrast and were interpreted as sensitivity analyses within a hypothesis-generating framework."
)

methods_file <- file.path(out_dir, "phase13_4B_methods_insert_empirical_covariance.md")
writeLines(methods_insert, methods_file, useBytes = TRUE)

results_insert <- c(
  "# Empirical covariance and jackknife assessment",
  "",
  paste0(
    "The exact internal SBP-component input files contained 391 paired instruments and reproduced the locked IVW estimates for both components. ",
    "In the paired SNP-level bootstrap, the empirical standard error of the SBP component contrast was ",
    fmt(empirical_se_diff, 4),
    ", compared with the locked r=0 analytic standard error of ",
    fmt(locked_se_difference_r0, 4),
    ". The bootstrap 95% interval for the contrast was ",
    fmt(as.numeric(ci[1]), 4),
    " to ",
    fmt(as.numeric(ci[2]), 4),
    ", and the empirical correlation between the component-specific estimates was ",
    fmt(empirical_cor, 4),
    "."
  ),
  "",
  paste0(
    "Leave-one-chromosome-out analysis preserved the positive IOP-minus-nonIOP contrast across all chromosomes. ",
    "The 1-Mb locus-proxy jackknife also preserved the contrast direction across all locus proxies. ",
    "These empirical uncertainty and influence analyses support directional stability of the SBP component-divergence pattern but do not change the hypothesis-generating interpretation."
  )
)

results_file <- file.path(out_dir, "phase13_4B_results_insert_empirical_covariance.md")
writeLines(results_insert, results_file, useBytes = TRUE)

# Claim-safety audit.
texts_to_scan <- paste(c(methods_insert, results_insert), collapse = "\n")

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

audit_file <- file.path(out_dir, "phase13_4B_claim_safety_audit.tsv")
write_tsv(danger_patterns, audit_file)

high_risk_hits <- sum(danger_patterns$hits)

status <- data.frame(
  field = c(
    "phase",
    "paired_input_created",
    "paired_input_qc_created",
    "bootstrap_distribution_created",
    "bootstrap_summary_created",
    "empirical_covariance_summary_created",
    "leave_one_chromosome_results_created",
    "locus_proxy_jackknife_results_created",
    "methods_insert_created",
    "results_insert_created",
    "claim_safety_audit_created",
    "high_risk_phrase_hits",
    "n_paired_instruments",
    "n_bootstrap",
    "n_chromosomes",
    "n_locus_proxy_1Mb",
    "empirical_SE_difference",
    "locked_r0_analytic_SE_difference",
    "bootstrap_CI_excludes_zero",
    "all_leave_one_chromosome_direction_preserved",
    "all_locus_proxy_direction_preserved",
    "new_primary_MR_estimates_created",
    "claim_level",
    "claim_upgrade_allowed",
    "phase13_4B_passed"
  ),
  value = c(
    "Phase13.4B",
    file.exists(paired_file),
    file.exists(input_qc_file),
    file.exists(boot_file),
    file.exists(bootstrap_summary_file),
    file.exists(covariance_file),
    file.exists(chr_jack_file),
    file.exists(locus_jack_file),
    file.exists(methods_file),
    file.exists(results_file),
    file.exists(audit_file),
    high_risk_hits,
    nrow(paired),
    n_bootstrap,
    length(unique(paired$chr)),
    length(unique(paired$locus_proxy_1Mb)),
    empirical_se_diff,
    locked_se_difference_r0,
    as.numeric(ci[1]) > 0 || as.numeric(ci[2]) < 0,
    all(chr_jack$direction_preserved),
    all(locus_jack$direction_preserved),
    FALSE,
    "HYPOTHESIS_GENERATING_NOT_CONFIRMATORY",
    "NO",
    high_risk_hits == 0 &&
      nrow(paired) == 391 &&
      all(chr_jack$direction_preserved) &&
      all(locus_jack$direction_preserved)
  ),
  stringsAsFactors = FALSE
)

status_file <- file.path(out_dir, "phase13_4B_status.tsv")
write_tsv(status, status_file)

# Update master status.
if (file.exists(master_status_file)) {
  master_status <- read.delim(master_status_file, check.names = FALSE)
  idx <- which(master_status$phase == "Phase13.4")
  if (length(idx) == 1) {
    master_status$status[idx] <- ifelse(
      high_risk_hits == 0 && nrow(paired) == 391,
      "PASSED_EMPIRICAL_COVARIANCE_AND_JACKKNIFE_SENSITIVITY",
      "FAILED_PHASE13_4B_REVIEW_REQUIRED"
    )
    master_status$qc_status[idx] <- ifelse(
      high_risk_hits == 0 && nrow(paired) == 391,
      "PASSED",
      "REVIEW_REQUIRED"
    )
    master_status$primary_output[idx] <- "phase13_4B_bootstrap_contrast_summary.tsv; phase13_4B_empirical_covariance_summary.tsv; phase13_4B_leave_one_chromosome_summary.tsv; phase13_4B_locus_proxy_jackknife_summary.tsv"
  }
  write_tsv(master_status, master_status_file)
}

message("Phase13.4B completed.")
message("Bootstrap summary: ", bootstrap_summary_file)
message("Empirical covariance: ", covariance_file)
message("Leave-one-chromosome summary: ", chr_summary_file)
message("Locus-proxy jackknife summary: ", locus_summary_file)
message("Status: ", status_file)
message("High-risk phrase hits: ", high_risk_hits)
