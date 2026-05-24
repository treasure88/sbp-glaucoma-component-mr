# Phase 10.3: SNP/locus influence analysis for SBP component contrast

options(stringsAsFactors = FALSE)

out_dir <- "../../58_component_divergence_reinforcement"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

non_file <- "../../16_mr_input_datasets/pairwise/SBP__GBS_nonIOPcomponent.mr_input.tsv.gz"
iop_file <- "../../16_mr_input_datasets/pairwise/SBP__GBS_IOPcomponent.mr_input.tsv.gz"
phase10_2_status_file <- file.path(out_dir, "phase10_2A_2B_status.tsv")

out_loo_snp <- file.path(out_dir, "phase10_3A_leave_one_snp_out_contrast.tsv")
out_loo_snp_summary <- file.path(out_dir, "phase10_3A_leave_one_snp_out_summary.tsv")
out_contrib <- file.path(out_dir, "phase10_3B_snp_contribution_to_contrast.tsv")
out_top <- file.path(out_dir, "phase10_3B_top_contributor_summary.tsv")
out_removal <- file.path(out_dir, "phase10_3C_top_contribution_removal_results.tsv")
out_locus <- file.path(out_dir, "phase10_3C_leave_one_locus_out_contrast.tsv")
out_radial <- file.path(out_dir, "phase10_3D_radial_baujat_influence_input.tsv")
out_status <- file.path(out_dir, "phase10_3_status.tsv")

read_tsv <- function(path) {
  if (grepl("\\.gz$", path)) {
    con <- gzfile(path, "rt")
    on.exit(close(con), add = TRUE)
    read.delim(con, sep = "\t", header = TRUE, quote = "", comment.char = "",
               check.names = FALSE, fill = TRUE)
  } else {
    read.delim(path, sep = "\t", header = TRUE, quote = "", comment.char = "",
               check.names = FALSE, fill = TRUE)
  }
}

write_tsv <- function(x, path) {
  write.table(x, path, sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
}

to_num <- function(x) suppressWarnings(as.numeric(x))

status_passed <- function(path) {
  if (!file.exists(path)) return(FALSE)
  x <- read_tsv(path)
  if (!"status" %in% names(x) || nrow(x) == 0) return(FALSE)
  grepl("PASS|PASSED|COMPLETED", as.character(x$status[1]), ignore.case = TRUE)
}

ivw_random <- function(bx, by, se_by) {
  bx <- to_num(bx); by <- to_num(by); se_by <- to_num(se_by)
  keep <- is.finite(bx) & is.finite(by) & is.finite(se_by) & se_by > 0 & bx != 0
  bx <- bx[keep]; by <- by[keep]; se_by <- se_by[keep]
  n <- length(bx)

  if (n < 3) {
    return(data.frame(
      n = n, beta = NA_real_, se_fixed = NA_real_, se_random = NA_real_,
      p_random = NA_real_, Q = NA_real_, Q_df = NA_integer_,
      Q_pval = NA_real_, phi = NA_real_
    ))
  }

  w <- 1 / se_by^2
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
    n = n,
    beta = beta,
    se_fixed = se_fixed,
    se_random = se_random,
    p_random = p_random,
    Q = Q,
    Q_df = Q_df,
    Q_pval = Q_pval,
    phi = phi
  )
}

contrast_from_rows <- function(d, idx) {
  dd <- d[idx, , drop = FALSE]

  non <- ivw_random(dd$beta_exposure_non, dd$beta_outcome_non, dd$se_outcome_non)
  iop <- ivw_random(dd$beta_exposure_iop, dd$beta_outcome_iop, dd$se_outcome_iop)

  beta_diff <- iop$beta - non$beta
  se_r0 <- sqrt(iop$se_random^2 + non$se_random^2)
  z <- beta_diff / se_r0
  p <- 2 * pnorm(-abs(z))

  data.frame(
    n_paired = nrow(dd),
    n_nonIOP = non$n,
    beta_nonIOP = non$beta,
    se_nonIOP = non$se_random,
    p_nonIOP = non$p_random,
    n_IOP = iop$n,
    beta_IOP = iop$beta,
    se_IOP = iop$se_random,
    p_IOP = iop$p_random,
    beta_difference_IOP_minus_nonIOP = beta_diff,
    se_difference_r0 = se_r0,
    z_contrast = z,
    p_contrast_r0 = p,
    direction = ifelse(
      beta_diff > 0,
      "IOP_MORE_POSITIVE_THAN_NONIOP",
      ifelse(beta_diff < 0, "IOP_MORE_NEGATIVE_THAN_NONIOP", "NO_DIFFERENCE")
    )
  )
}

make_locus_id <- function(chr, pos, window_bp = 1000000) {
  chr <- as.character(chr)
  pos <- to_num(pos)
  locus_id <- rep(NA_character_, length(chr))

  for (cc in unique(chr)) {
    idx <- which(chr == cc & is.finite(pos))
    if (length(idx) == 0) next

    idx <- idx[order(pos[idx])]
    cluster_start <- pos[idx[1]]
    cluster_last <- pos[idx[1]]
    members <- idx[1]

    for (k in seq_along(idx)) {
      current <- idx[k]
      if (k == 1) {
        members <- current
        cluster_start <- pos[current]
        cluster_last <- pos[current]
      } else if ((pos[current] - cluster_last) <= window_bp) {
        members <- c(members, current)
        cluster_last <- pos[current]
      } else {
        locus_id[members] <- paste0("chr", cc, ":", min(pos[members]), "-", max(pos[members]))
        members <- current
        cluster_start <- pos[current]
        cluster_last <- pos[current]
      }
    }

    locus_id[members] <- paste0("chr", cc, ":", min(pos[members]), "-", max(pos[members]))
  }

  locus_id[is.na(locus_id)] <- paste0("single_snp_locus_", seq_len(sum(is.na(locus_id))))
  locus_id
}

tryCatch({

  if (!status_passed(phase10_2_status_file)) {
    stop("Phase 10.2A/10.2B status is not PASS/PASSED/COMPLETED. Do not start Phase 10.3.")
  }

  if (!file.exists(non_file)) stop(paste("Missing nonIOP MR input:", non_file))
  if (!file.exists(iop_file)) stop(paste("Missing IOP MR input:", iop_file))

  non <- read_tsv(non_file)
  iop <- read_tsv(iop_file)

  required_mr <- c("SNP", "chr", "pos", "beta_exposure", "se_exposure", "beta_outcome", "se_outcome")
  miss_non <- setdiff(required_mr, names(non))
  miss_iop <- setdiff(required_mr, names(iop))

  if (length(miss_non) > 0) stop(paste("Missing nonIOP columns:", paste(miss_non, collapse = "; ")))
  if (length(miss_iop) > 0) stop(paste("Missing IOP columns:", paste(miss_iop, collapse = "; ")))

  if ("include_in_main" %in% names(non)) {
    non <- non[toupper(as.character(non$include_in_main)) %in% c("TRUE", "T", "1", "YES"), , drop = FALSE]
  }

  if ("include_in_main" %in% names(iop)) {
    iop <- iop[toupper(as.character(iop$include_in_main)) %in% c("TRUE", "T", "1", "YES"), , drop = FALSE]
  }

  paired <- merge(non, iop, by = "SNP", suffixes = c("_non", "_iop"))

  if (nrow(paired) < 30) stop(paste("Too few paired SNPs:", nrow(paired)))

  paired$chr <- if ("chr_non" %in% names(paired)) paired$chr_non else paired$chr
  paired$pos <- if ("pos_non" %in% names(paired)) paired$pos_non else paired$pos

  paired$chr <- as.character(paired$chr)
  paired$pos <- to_num(paired$pos)
  paired$locus_id_1mb_greedy <- make_locus_id(paired$chr, paired$pos, window_bp = 1000000)

  full <- contrast_from_rows(paired, seq_len(nrow(paired)))
  full_beta <- full$beta_difference_IOP_minus_nonIOP
  full_p <- full$p_contrast_r0

  loo_rows <- vector("list", nrow(paired))

  for (j in seq_len(nrow(paired))) {
    idx <- setdiff(seq_len(nrow(paired)), j)
    tmp <- contrast_from_rows(paired, idx)
    tmp$left_out_SNP <- paired$SNP[j]
    tmp$left_out_chr <- paired$chr[j]
    tmp$left_out_pos <- paired$pos[j]
    tmp$left_out_locus_id <- paired$locus_id_1mb_greedy[j]
    tmp$delta_vs_full <- tmp$beta_difference_IOP_minus_nonIOP - full_beta
    tmp$abs_delta_vs_full <- abs(tmp$delta_vs_full)
    tmp$relative_abs_delta_vs_full <- abs(tmp$delta_vs_full) / abs(full_beta)
    tmp$direction_preserved_vs_full <- ifelse(sign(tmp$beta_difference_IOP_minus_nonIOP) == sign(full_beta), "YES", "NO")
    tmp$nominal_p_lt_0_05 <- ifelse(!is.na(tmp$p_contrast_r0) & tmp$p_contrast_r0 < 0.05, "YES", "NO")
    loo_rows[[j]] <- tmp
  }

  loo <- do.call(rbind, loo_rows)
  loo <- loo[, c(
    "left_out_SNP", "left_out_chr", "left_out_pos", "left_out_locus_id",
    setdiff(names(loo), c("left_out_SNP", "left_out_chr", "left_out_pos", "left_out_locus_id"))
  )]
  loo <- loo[order(-loo$abs_delta_vs_full), , drop = FALSE]
  write_tsv(loo, out_loo_snp)

  paired$w_non <- 1 / to_num(paired$se_outcome_non)^2
  paired$w_iop <- 1 / to_num(paired$se_outcome_iop)^2

  denom_non <- sum(paired$w_non * to_num(paired$beta_exposure_non)^2, na.rm = TRUE)
  denom_iop <- sum(paired$w_iop * to_num(paired$beta_exposure_iop)^2, na.rm = TRUE)

  paired$wald_nonIOP <- to_num(paired$beta_outcome_non) / to_num(paired$beta_exposure_non)
  paired$wald_IOP <- to_num(paired$beta_outcome_iop) / to_num(paired$beta_exposure_iop)
  paired$wald_difference_IOP_minus_nonIOP <- paired$wald_IOP - paired$wald_nonIOP

  paired$fixed_ivw_contribution_nonIOP <- paired$w_non * to_num(paired$beta_exposure_non) * to_num(paired$beta_outcome_non) / denom_non
  paired$fixed_ivw_contribution_IOP <- paired$w_iop * to_num(paired$beta_exposure_iop) * to_num(paired$beta_outcome_iop) / denom_iop
  paired$fixed_ivw_contribution_contrast <- paired$fixed_ivw_contribution_IOP - paired$fixed_ivw_contribution_nonIOP
  paired$abs_fixed_ivw_contribution_contrast <- abs(paired$fixed_ivw_contribution_contrast)

  loo_key <- loo[, c(
    "left_out_SNP", "delta_vs_full", "abs_delta_vs_full",
    "relative_abs_delta_vs_full", "p_contrast_r0",
    "direction_preserved_vs_full", "nominal_p_lt_0_05"
  )]
  names(loo_key) <- c(
    "SNP", "leave_one_out_delta_vs_full", "leave_one_out_abs_delta_vs_full",
    "leave_one_out_relative_abs_delta_vs_full", "leave_one_out_p_contrast_r0",
    "leave_one_out_direction_preserved", "leave_one_out_nominal_p_lt_0_05"
  )

  contrib <- merge(paired, loo_key, by = "SNP", all.x = TRUE)

  contrib$rank_by_leave_one_out_abs_delta <- rank(-contrib$leave_one_out_abs_delta_vs_full, ties.method = "min")
  contrib$rank_by_abs_fixed_ivw_contribution <- rank(-contrib$abs_fixed_ivw_contribution_contrast, ties.method = "min")
  contrib$share_of_total_abs_fixed_contribution_percent <- 100 * contrib$abs_fixed_ivw_contribution_contrast /
    sum(contrib$abs_fixed_ivw_contribution_contrast, na.rm = TRUE)

  contrib_out <- contrib[, c(
    "SNP", "chr", "pos", "locus_id_1mb_greedy",
    "wald_nonIOP", "wald_IOP", "wald_difference_IOP_minus_nonIOP",
    "fixed_ivw_contribution_nonIOP", "fixed_ivw_contribution_IOP",
    "fixed_ivw_contribution_contrast", "abs_fixed_ivw_contribution_contrast",
    "share_of_total_abs_fixed_contribution_percent",
    "leave_one_out_delta_vs_full", "leave_one_out_abs_delta_vs_full",
    "leave_one_out_relative_abs_delta_vs_full",
    "leave_one_out_p_contrast_r0",
    "leave_one_out_direction_preserved",
    "leave_one_out_nominal_p_lt_0_05",
    "rank_by_leave_one_out_abs_delta",
    "rank_by_abs_fixed_ivw_contribution"
  )]

  contrib_out <- contrib_out[order(contrib_out$rank_by_leave_one_out_abs_delta), , drop = FALSE]
  write_tsv(contrib_out, out_contrib)

  loci <- unique(paired$locus_id_1mb_greedy)
  locus_rows <- vector("list", length(loci))

  for (k in seq_along(loci)) {
    loc <- loci[k]
    idx <- which(paired$locus_id_1mb_greedy != loc)
    removed <- paired[paired$locus_id_1mb_greedy == loc, , drop = FALSE]
    tmp <- contrast_from_rows(paired, idx)
    tmp$left_out_locus_id <- loc
    tmp$n_removed_snps <- nrow(removed)
    tmp$removed_snps <- paste(removed$SNP, collapse = ";")
    tmp$removed_chr <- paste(unique(removed$chr), collapse = ";")
    tmp$removed_pos_min <- min(removed$pos, na.rm = TRUE)
    tmp$removed_pos_max <- max(removed$pos, na.rm = TRUE)
    tmp$delta_vs_full <- tmp$beta_difference_IOP_minus_nonIOP - full_beta
    tmp$abs_delta_vs_full <- abs(tmp$delta_vs_full)
    tmp$relative_abs_delta_vs_full <- abs(tmp$delta_vs_full) / abs(full_beta)
    tmp$direction_preserved_vs_full <- ifelse(sign(tmp$beta_difference_IOP_minus_nonIOP) == sign(full_beta), "YES", "NO")
    tmp$nominal_p_lt_0_05 <- ifelse(!is.na(tmp$p_contrast_r0) & tmp$p_contrast_r0 < 0.05, "YES", "NO")
    locus_rows[[k]] <- tmp
  }

  locus <- do.call(rbind, locus_rows)
  locus <- locus[, c(
    "left_out_locus_id", "n_removed_snps", "removed_snps", "removed_chr",
    "removed_pos_min", "removed_pos_max",
    setdiff(names(locus), c(
      "left_out_locus_id", "n_removed_snps", "removed_snps",
      "removed_chr", "removed_pos_min", "removed_pos_max"
    ))
  )]

  locus <- locus[order(-locus$abs_delta_vs_full), , drop = FALSE]
  write_tsv(locus, out_locus)

  top_snp <- head(contrib_out[order(contrib_out$rank_by_leave_one_out_abs_delta), , drop = FALSE], 20)
  top_snp_summary <- data.frame(
    unit_type = "SNP",
    rank = seq_len(nrow(top_snp)),
    unit_id = top_snp$SNP,
    chr = top_snp$chr,
    pos_or_range = as.character(top_snp$pos),
    locus_id = top_snp$locus_id_1mb_greedy,
    n_snps_in_unit = 1,
    abs_delta_vs_full = top_snp$leave_one_out_abs_delta_vs_full,
    relative_abs_delta_vs_full = top_snp$leave_one_out_relative_abs_delta_vs_full,
    p_after_removal = top_snp$leave_one_out_p_contrast_r0,
    direction_preserved_after_removal = top_snp$leave_one_out_direction_preserved,
    nominal_p_lt_0_05_after_removal = top_snp$leave_one_out_nominal_p_lt_0_05,
    stringsAsFactors = FALSE
  )

  top_locus <- head(locus, 20)
  top_locus_summary <- data.frame(
    unit_type = "LOCUS_1MB_GREEDY",
    rank = seq_len(nrow(top_locus)),
    unit_id = top_locus$left_out_locus_id,
    chr = top_locus$removed_chr,
    pos_or_range = paste0(top_locus$removed_pos_min, "-", top_locus$removed_pos_max),
    locus_id = top_locus$left_out_locus_id,
    n_snps_in_unit = top_locus$n_removed_snps,
    abs_delta_vs_full = top_locus$abs_delta_vs_full,
    relative_abs_delta_vs_full = top_locus$relative_abs_delta_vs_full,
    p_after_removal = top_locus$p_contrast_r0,
    direction_preserved_after_removal = top_locus$direction_preserved_vs_full,
    nominal_p_lt_0_05_after_removal = top_locus$nominal_p_lt_0_05,
    stringsAsFactors = FALSE
  )

  top_summary <- rbind(top_snp_summary, top_locus_summary)
  write_tsv(top_summary, out_top)

  removal_rows <- list()

  add_removal <- function(label, unit_type, ids, snp_ids) {
    idx <- which(!paired$SNP %in% snp_ids)
    tmp <- contrast_from_rows(paired, idx)
    tmp$removal_label <- label
    tmp$unit_type <- unit_type
    tmp$removed_units <- paste(ids, collapse = ";")
    tmp$removed_snps <- paste(snp_ids, collapse = ";")
    tmp$n_removed_snps <- length(unique(snp_ids))
    tmp$delta_vs_full <- tmp$beta_difference_IOP_minus_nonIOP - full_beta
    tmp$abs_delta_vs_full <- abs(tmp$delta_vs_full)
    tmp$relative_abs_delta_vs_full <- abs(tmp$delta_vs_full) / abs(full_beta)
    tmp$direction_preserved_vs_full <- ifelse(sign(tmp$beta_difference_IOP_minus_nonIOP) == sign(full_beta), "YES", "NO")
    tmp$nominal_p_lt_0_05 <- ifelse(!is.na(tmp$p_contrast_r0) & tmp$p_contrast_r0 < 0.05, "YES", "NO")
    removal_rows[[length(removal_rows) + 1]] <<- tmp
  }

  snp_by_loo <- contrib_out[order(contrib_out$rank_by_leave_one_out_abs_delta), "SNP"]
  snp_by_fixed <- contrib_out[order(contrib_out$rank_by_abs_fixed_ivw_contribution), "SNP"]

  for (kk in c(1, 3, 5, 10)) {
    add_removal(
      paste0("remove_top_", kk, "_SNPs_by_leave_one_out_delta"),
      "SNP",
      head(snp_by_loo, kk),
      head(snp_by_loo, kk)
    )

    add_removal(
      paste0("remove_top_", kk, "_SNPs_by_fixed_ivw_contribution"),
      "SNP",
      head(snp_by_fixed, kk),
      head(snp_by_fixed, kk)
    )
  }

  locus_ids_by_delta <- locus$left_out_locus_id
  for (kk in c(1, 3, 5)) {
    ids <- head(locus_ids_by_delta, kk)
    snps <- paired$SNP[paired$locus_id_1mb_greedy %in% ids]
    add_removal(
      paste0("remove_top_", kk, "_loci_by_leave_one_locus_out_delta"),
      "LOCUS_1MB_GREEDY",
      ids,
      snps
    )
  }

  removal <- do.call(rbind, removal_rows)
  removal <- removal[, c(
    "removal_label", "unit_type", "removed_units", "removed_snps", "n_removed_snps",
    setdiff(names(removal), c("removal_label", "unit_type", "removed_units", "removed_snps", "n_removed_snps"))
  )]
  write_tsv(removal, out_removal)

  full_non <- ivw_random(paired$beta_exposure_non, paired$beta_outcome_non, paired$se_outcome_non)
  full_iop <- ivw_random(paired$beta_exposure_iop, paired$beta_outcome_iop, paired$se_outcome_iop)

  radial <- data.frame(
    SNP = paired$SNP,
    chr = paired$chr,
    pos = paired$pos,
    locus_id_1mb_greedy = paired$locus_id_1mb_greedy,
    beta_exposure_non = to_num(paired$beta_exposure_non),
    beta_outcome_non = to_num(paired$beta_outcome_non),
    se_outcome_non = to_num(paired$se_outcome_non),
    beta_exposure_iop = to_num(paired$beta_exposure_iop),
    beta_outcome_iop = to_num(paired$beta_outcome_iop),
    se_outcome_iop = to_num(paired$se_outcome_iop),
    stringsAsFactors = FALSE
  )

  radial$weighted_residual_nonIOP <- (radial$beta_outcome_non - full_non$beta * radial$beta_exposure_non) / radial$se_outcome_non
  radial$weighted_residual_IOP <- (radial$beta_outcome_iop - full_iop$beta * radial$beta_exposure_iop) / radial$se_outcome_iop
  radial$Q_contribution_nonIOP <- radial$weighted_residual_nonIOP^2
  radial$Q_contribution_IOP <- radial$weighted_residual_IOP^2
  radial$Q_contribution_sum <- radial$Q_contribution_nonIOP + radial$Q_contribution_IOP
  radial$leave_one_out_abs_delta_vs_full <- contrib$leave_one_out_abs_delta_vs_full[match(radial$SNP, contrib$SNP)]
  radial$rank_by_Q_contribution_sum <- rank(-radial$Q_contribution_sum, ties.method = "min")
  radial$rank_by_leave_one_out_abs_delta <- rank(-radial$leave_one_out_abs_delta_vs_full, ties.method = "min")

  radial <- radial[order(radial$rank_by_leave_one_out_abs_delta), , drop = FALSE]
  write_tsv(radial, out_radial)

  loo_summary <- data.frame(
    metric = c(
      "full_data_n_paired",
      "full_data_beta_difference",
      "full_data_se_r0",
      "full_data_p_r0",
      "leave_one_snp_n_tests",
      "leave_one_snp_direction_preserved_percent",
      "leave_one_snp_nominal_p_lt_0_05_percent",
      "leave_one_snp_min_beta_difference",
      "leave_one_snp_max_beta_difference",
      "leave_one_snp_max_abs_delta_vs_full",
      "leave_one_snp_max_relative_abs_delta_vs_full",
      "leave_one_locus_n_tests",
      "leave_one_locus_direction_preserved_percent",
      "leave_one_locus_nominal_p_lt_0_05_percent",
      "leave_one_locus_max_abs_delta_vs_full",
      "leave_one_locus_max_relative_abs_delta_vs_full"
    ),
    value = c(
      full$n_paired,
      full$beta_difference_IOP_minus_nonIOP,
      full$se_difference_r0,
      full$p_contrast_r0,
      nrow(loo),
      mean(loo$direction_preserved_vs_full == "YES", na.rm = TRUE) * 100,
      mean(loo$nominal_p_lt_0_05 == "YES", na.rm = TRUE) * 100,
      min(loo$beta_difference_IOP_minus_nonIOP, na.rm = TRUE),
      max(loo$beta_difference_IOP_minus_nonIOP, na.rm = TRUE),
      max(loo$abs_delta_vs_full, na.rm = TRUE),
      max(loo$relative_abs_delta_vs_full, na.rm = TRUE),
      nrow(locus),
      mean(locus$direction_preserved_vs_full == "YES", na.rm = TRUE) * 100,
      mean(locus$nominal_p_lt_0_05 == "YES", na.rm = TRUE) * 100,
      max(locus$abs_delta_vs_full, na.rm = TRUE),
      max(locus$relative_abs_delta_vs_full, na.rm = TRUE)
    ),
    stringsAsFactors = FALSE
  )

  write_tsv(loo_summary, out_loo_snp_summary)

  single_snp_sign_flip <- any(loo$direction_preserved_vs_full != "YES", na.rm = TRUE)
  single_snp_nominal_loss <- any(loo$nominal_p_lt_0_05 != "YES", na.rm = TRUE)
  locus_sign_flip <- any(locus$direction_preserved_vs_full != "YES", na.rm = TRUE)
  locus_nominal_loss <- any(locus$nominal_p_lt_0_05 != "YES", na.rm = TRUE)

  max_snp_relative_change <- max(loo$relative_abs_delta_vs_full, na.rm = TRUE)
  max_locus_relative_change <- max(locus$relative_abs_delta_vs_full, na.rm = TRUE)

  top1_snp_removal <- removal[removal$removal_label == "remove_top_1_SNPs_by_leave_one_out_delta", , drop = FALSE]
  top3_snp_removal <- removal[removal$removal_label == "remove_top_3_SNPs_by_leave_one_out_delta", , drop = FALSE]
  top1_locus_removal <- removal[removal$removal_label == "remove_top_1_loci_by_leave_one_locus_out_delta", , drop = FALSE]

  passed_no_single_snp_or_locus_dominance <-
    !single_snp_sign_flip &&
    !locus_sign_flip &&
    max_snp_relative_change < 0.25 &&
    max_locus_relative_change < 0.35 &&
    top1_snp_removal$direction_preserved_vs_full[1] == "YES" &&
    top3_snp_removal$direction_preserved_vs_full[1] == "YES" &&
    top1_locus_removal$direction_preserved_vs_full[1] == "YES"

  status_label <- if (passed_no_single_snp_or_locus_dominance) {
    "PASSED_NO_SINGLE_SNP_OR_LOCUS_DOMINANCE"
  } else if (!single_snp_sign_flip && !locus_sign_flip) {
    "CHECK_TOP_CONTRIBUTORS_BUT_DIRECTION_STABLE"
  } else {
    "CHECK_POSSIBLE_SINGLE_SNP_OR_LOCUS_DOMINANCE"
  }

  status <- data.frame(
    phase = "Phase10.3",
    status = status_label,
    n_paired_snps = nrow(paired),
    n_loci_1mb_greedy = length(unique(paired$locus_id_1mb_greedy)),
    full_beta_difference_IOP_minus_nonIOP = full_beta,
    full_se_r0 = full$se_difference_r0,
    full_p_r0 = full_p,
    leave_one_snp_any_direction_flip = single_snp_sign_flip,
    leave_one_snp_any_nominal_loss = single_snp_nominal_loss,
    leave_one_snp_max_relative_abs_delta = max_snp_relative_change,
    leave_one_locus_any_direction_flip = locus_sign_flip,
    leave_one_locus_any_nominal_loss = locus_nominal_loss,
    leave_one_locus_max_relative_abs_delta = max_locus_relative_change,
    top1_snp_removal_direction_preserved = top1_snp_removal$direction_preserved_vs_full[1],
    top3_snp_removal_direction_preserved = top3_snp_removal$direction_preserved_vs_full[1],
    top1_locus_removal_direction_preserved = top1_locus_removal$direction_preserved_vs_full[1],
    output_directory = out_dir,
    timestamp = as.character(Sys.time()),
    stringsAsFactors = FALSE
  )

  write_tsv(status, out_status)

  cat("Phase10.3 status:", status_label, "\n")
  cat("Paired SNPs:", nrow(paired), "\n")
  cat("1Mb greedy loci:", length(unique(paired$locus_id_1mb_greedy)), "\n")
  cat("Outputs written to:", out_dir, "\n")

}, error = function(e) {

  status <- data.frame(
    phase = "Phase10.3",
    status = "FAILED",
    error_message = conditionMessage(e),
    output_directory = out_dir,
    timestamp = as.character(Sys.time()),
    stringsAsFactors = FALSE
  )

  write_tsv(status, out_status)
  stop(e)
})
