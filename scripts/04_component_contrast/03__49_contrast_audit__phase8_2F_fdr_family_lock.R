# Phase 8.2F: Lock FDR family for the primary component-contrast analysis
# Purpose: document that SBP/DBP q values derive from the original five-exposure contrast family

options(stringsAsFactors = FALSE)

out_dir <- "../../49_contrast_audit"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

phase4_1_file <- "../../21_sbp_component_contrast/phase4_1_component_contrast_all_exposures.tsv"
integrated_file <- "../../24_vascular_panel_integration/phase4_8_vascular_panel_integrated_evidence_table.tsv"

lock_file <- file.path(out_dir, "phase8_2F_fdr_family_lock.tsv")
comparison_file <- file.path(out_dir, "phase8_2F_integrated_q_comparison.tsv")
status_file <- file.path(out_dir, "phase8_2F_fdr_family_lock_status.tsv")

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
    phase = "Phase8.2F",
    status = status,
    note = note,
    timestamp = as.character(Sys.time()),
    stringsAsFactors = FALSE
  )
  write.table(x, status_file, sep = "\t", quote = FALSE, row.names = FALSE)
}

tryCatch({

  if (!file.exists(phase4_1_file)) {
    write_status("FAILED", paste("Missing Phase 4.1 file:", phase4_1_file))
    stop("Missing Phase 4.1 file")
  }

  d <- read_tsv(phase4_1_file)

  required <- c(
    "exposure_id",
    "assumed_component_correlation_r",
    "p_contrast",
    "q_contrast_within_r",
    "beta_difference_IOP_minus_nonIOP",
    "se_difference",
    "contrast_direction"
  )

  missing <- setdiff(required, names(d))
  if (length(missing) > 0) {
    write_status("FAILED", paste("Missing columns:", paste(missing, collapse = ";")))
    stop("Missing required columns")
  }

  d$assumed_component_correlation_r_num <- to_num(d$assumed_component_correlation_r)
  d$p_contrast_num <- to_num(d$p_contrast)
  d$q_contrast_within_r_num <- to_num(d$q_contrast_within_r)

  # Primary stored contrast in later tables corresponds to r=0.
  d0 <- d[!is.na(d$assumed_component_correlation_r_num) & d$assumed_component_correlation_r_num == 0, , drop = FALSE]

  d0 <- d0[order(d0$p_contrast_num), , drop = FALSE]

  d0$bh_q_recomputed_within_r0_family <- p.adjust(d0$p_contrast_num, method = "BH")
  d0$p_rank_within_r0_family <- rank(d0$p_contrast_num, ties.method = "first")
  d0$fdr_family_id <- "PHASE4_1_ORIGINAL_FIVE_EXPOSURE_COMPONENT_CONTRAST_FAMILY_R0"
  d0$fdr_family_size <- nrow(d0)
  d0$fdr_family_members_ordered_by_p <- paste(d0$exposure_id, collapse = ";")
  d0$q_matches_recomputed_bh <- ifelse(
    abs(d0$q_contrast_within_r_num - d0$bh_q_recomputed_within_r0_family) < 1e-8,
    "YES",
    "NO"
  )

  lock <- data.frame(
    fdr_family_id = d0$fdr_family_id,
    source_file = phase4_1_file,
    assumed_component_correlation_r = d0$assumed_component_correlation_r_num,
    exposure_id = d0$exposure_id,
    p_contrast = d0$p_contrast_num,
    q_contrast_within_r_from_source = d0$q_contrast_within_r_num,
    bh_q_recomputed_within_r0_family = d0$bh_q_recomputed_within_r0_family,
    q_matches_recomputed_bh = d0$q_matches_recomputed_bh,
    p_rank_within_family = d0$p_rank_within_r0_family,
    fdr_family_size = d0$fdr_family_size,
    fdr_family_members_ordered_by_p = d0$fdr_family_members_ordered_by_p,
    beta_difference_IOP_minus_nonIOP = to_num(d0$beta_difference_IOP_minus_nonIOP),
    se_difference = to_num(d0$se_difference),
    contrast_direction = d0$contrast_direction,
    interpretation_note = "Primary r=0 component-contrast FDR family from Phase 4.1; later vascular-focused tables retained these conservative q values where available.",
    stringsAsFactors = FALSE
  )

  write.table(lock, lock_file, sep = "\t", quote = FALSE, row.names = FALSE)

  # Compare with integrated evidence table if available
  if (file.exists(integrated_file)) {
    g <- read_tsv(integrated_file)

    if (all(c("exposure_id", "contrast_p_r0", "contrast_q_r0") %in% names(g))) {
      g$contrast_p_r0_num <- to_num(g$contrast_p_r0)
      g$contrast_q_r0_num <- to_num(g$contrast_q_r0)

      comp <- merge(
        g[, c("exposure_id", "contrast_p_r0_num", "contrast_q_r0_num")],
        lock[, c("exposure_id", "p_contrast", "q_contrast_within_r_from_source", "bh_q_recomputed_within_r0_family")],
        by = "exposure_id",
        all.x = TRUE
      )

      comp$p_matches_phase4_1 <- ifelse(
        !is.na(comp$contrast_p_r0_num) & !is.na(comp$p_contrast) &
          abs(comp$contrast_p_r0_num - comp$p_contrast) < 1e-8,
        "YES",
        ifelse(is.na(comp$p_contrast), "NO_PHASE4_1_MATCH", "NO")
      )

      comp$q_matches_phase4_1 <- ifelse(
        !is.na(comp$contrast_q_r0_num) & !is.na(comp$q_contrast_within_r_from_source) &
          abs(comp$contrast_q_r0_num - comp$q_contrast_within_r_from_source) < 1e-8,
        "YES",
        ifelse(is.na(comp$contrast_q_r0_num), "INTEGRATED_Q_MISSING", "NO")
      )

      comp$note <- ifelse(
        comp$exposure_id %in% c("SBP", "DBP"),
        "Integrated q is inherited from Phase 4.1 five-exposure family.",
        "Not part of Phase 4.1 five-exposure q lock or q unavailable in integrated table."
      )

      write.table(comp, comparison_file, sep = "\t", quote = FALSE, row.names = FALSE)
    }
  }

  all_lock_match <- all(lock$q_matches_recomputed_bh == "YES")
  family_size_ok <- unique(lock$fdr_family_size) == 5
  has_sbp <- "SBP" %in% lock$exposure_id
  has_dbp <- "DBP" %in% lock$exposure_id

  if (all_lock_match && family_size_ok && has_sbp && has_dbp) {
    write_status(
      "PASSED_FDR_FAMILY_LOCKED",
      paste0(
        "Locked r=0 FDR family as original five-exposure component-contrast family: ",
        paste(lock$exposure_id, collapse = ";"),
        ". Recomputed BH q values match source q values."
      )
    )
  } else {
    write_status(
      "CHECK_REQUIRED",
      "FDR family lock created, but at least one expected condition was not met. Review lock file."
    )
  }

}, error = function(e) {
  write_status("FAILED", paste("Error:", conditionMessage(e)))
  stop(e)
})
