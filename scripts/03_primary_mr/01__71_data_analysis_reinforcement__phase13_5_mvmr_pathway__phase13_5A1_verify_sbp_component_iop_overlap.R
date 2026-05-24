options(stringsAsFactors = FALSE)

# Phase13.5A1: Exact overlap verification for exploratory SBP+IOP MVMR sensitivity
# Project: IOP-dependent vs IOP-independent glaucoma component MR
# Evidence level: HYPOTHESIS_GENERATING_NOT_CONFIRMATORY
# This phase builds an overlap dataset only. It does not run MVMR estimates.

message("Running Phase13.5A1 exact SBP-component + measured IOP overlap verification...")

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

as_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

pick_col <- function(df, candidates) {
  nms <- names(df)
  low <- tolower(nms)
  for (cc in candidates) {
    idx <- which(low == tolower(cc))
    if (length(idx) > 0) return(nms[idx[1]])
  }
  NA_character_
}

read_any <- function(path) {
  if (!file.exists(path)) stop("Missing file: ", path)
  if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    con <- gzfile(path, open = "rt")
    on.exit(close(con), add = TRUE)
    read.delim(con, header = TRUE, sep = "\t", quote = "", comment.char = "", check.names = FALSE)
  } else {
    read.delim(path, header = TRUE, sep = "\t", quote = "", comment.char = "", check.names = FALSE)
  }
}

root <- normalizePath(".", winslash = "/", mustWork = TRUE)
out_dir <- file.path(root, "71_data_analysis_reinforcement", "phase13_5_mvmr_pathway")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

master_status_file <- file.path(root, "71_data_analysis_reinforcement", "phase13_master_status.tsv")

internal_file <- file.path(
  root,
  "71_data_analysis_reinforcement/phase13_4_empirical_covariance/phase13_4B_paired_internal_sbp_component_input.tsv"
)

iop_candidates <- c(
  "36_iop_external_validation_inputs/mr_input/SBP__IOP.external_mr_input.tsv.gz",
  "36_iop_external_validation_inputs/pairwise_harmonized/SBP__IOP.targeted_external_harmonized.tsv.gz",
  "16_mr_input_datasets/pairwise/SBP__IOPcc_coordinate_subset.mr_input.tsv.gz"
)

iop_candidates_full <- file.path(root, iop_candidates)
iop_exists <- file.exists(iop_candidates_full)

manifest <- data.frame(
  role = c("paired_internal_component_input", paste0("candidate_measured_IOP_association_", seq_along(iop_candidates))),
  relative_path = c(
    "71_data_analysis_reinforcement/phase13_4_empirical_covariance/phase13_4B_paired_internal_sbp_component_input.tsv",
    iop_candidates
  ),
  exists = c(file.exists(internal_file), iop_exists),
  stringsAsFactors = FALSE
)

manifest_file <- file.path(out_dir, "phase13_5A1_exact_overlap_input_manifest.tsv")
write_tsv(manifest, manifest_file)

if (!file.exists(internal_file)) {
  stop("Missing paired internal component input from Phase13.4B.")
}

if (!any(iop_exists)) {
  stop("No candidate measured IOP association file found.")
}

# Prefer external measured IOP MR input with 456 SBP instruments.
iop_file <- iop_candidates_full[which(iop_exists)[1]]
iop_relative <- iop_candidates[which(iop_exists)[1]]

internal <- read_any(internal_file)
iop <- read_any(iop_file)

# Internal paired input columns from Phase13.4B.
required_internal <- c(
  "SNP",
  "chr",
  "pos",
  "effect_allele",
  "other_allele",
  "beta_exposure",
  "se_exposure",
  "beta_outcome_nonIOP",
  "se_outcome_nonIOP",
  "beta_outcome_IOP",
  "se_outcome_IOP"
)

missing_internal <- setdiff(required_internal, names(internal))
if (length(missing_internal) > 0) {
  stop("Missing required internal columns: ", paste(missing_internal, collapse = ", "))
}

# IOP association columns.
snp_col <- pick_col(iop, c("SNP", "rsid", "variant_id"))
chr_col <- pick_col(iop, c("chr", "CHR", "chromosome"))
pos_col <- pick_col(iop, c("pos", "POS", "position", "bp"))
ea_col <- pick_col(iop, c("effect_allele", "effect_allele_exposure", "exposure_effect_allele"))
oa_col <- pick_col(iop, c("other_allele", "other_allele_exposure", "exposure_other_allele"))
beta_sbp_col <- pick_col(iop, c("beta_exposure", "beta.exposure", "beta_exp", "bx"))
se_sbp_col <- pick_col(iop, c("se_exposure", "se.exposure", "se_exp", "sx"))
beta_iop_col <- pick_col(iop, c("beta_outcome", "beta.outcome", "beta_out", "by", "beta_outcome_aligned"))
se_iop_col <- pick_col(iop, c("se_outcome", "se.outcome", "se_out", "sy"))
p_iop_col <- pick_col(iop, c("pval_outcome", "pval.outcome", "p.outcome", "minus_log10_p_outcome"))
f_col <- pick_col(iop, c("F_stat", "instrument_F_stat", "F", "f_stat"))
harm_col <- pick_col(iop, c("harmonization_action"))
pal_col <- pick_col(iop, c("palindromic_possible"))

required_iop <- c(snp_col, beta_sbp_col, se_sbp_col, beta_iop_col, se_iop_col)
if (any(is.na(required_iop))) {
  stop("Measured IOP file missing required SNP-level columns.")
}

iop_std <- data.frame(
  SNP = as.character(iop[[snp_col]]),
  chr_iop_assoc = if (!is.na(chr_col)) as.character(iop[[chr_col]]) else NA_character_,
  pos_iop_assoc = if (!is.na(pos_col)) as_num(iop[[pos_col]]) else NA_real_,
  effect_allele_iop_assoc = if (!is.na(ea_col)) as.character(iop[[ea_col]]) else NA_character_,
  other_allele_iop_assoc = if (!is.na(oa_col)) as.character(iop[[oa_col]]) else NA_character_,
  beta_SBP_from_iop_file = as_num(iop[[beta_sbp_col]]),
  se_SBP_from_iop_file = as_num(iop[[se_sbp_col]]),
  beta_measured_IOP = as_num(iop[[beta_iop_col]]),
  se_measured_IOP = as_num(iop[[se_iop_col]]),
  pval_measured_IOP = if (!is.na(p_iop_col)) as_num(iop[[p_iop_col]]) else NA_real_,
  F_stat_iop_file = if (!is.na(f_col)) as_num(iop[[f_col]]) else NA_real_,
  harmonization_iop_file = if (!is.na(harm_col)) as.character(iop[[harm_col]]) else NA_character_,
  palindromic_possible_iop_file = if (!is.na(pal_col)) as.character(iop[[pal_col]]) else NA_character_,
  stringsAsFactors = FALSE
)

internal_std <- data.frame(
  SNP = as.character(internal$SNP),
  chr = as.character(internal$chr),
  pos = as_num(internal$pos),
  effect_allele = as.character(internal$effect_allele),
  other_allele = as.character(internal$other_allele),
  beta_SBP = as_num(internal$beta_exposure),
  se_SBP = as_num(internal$se_exposure),
  beta_GBS_nonIOPcomponent = as_num(internal$beta_outcome_nonIOP),
  se_GBS_nonIOPcomponent = as_num(internal$se_outcome_nonIOP),
  beta_GBS_IOPcomponent = as_num(internal$beta_outcome_IOP),
  se_GBS_IOPcomponent = as_num(internal$se_outcome_IOP),
  locus_proxy_1Mb = if ("locus_proxy_1Mb" %in% names(internal)) as.character(internal$locus_proxy_1Mb) else NA_character_,
  stringsAsFactors = FALSE
)

merged <- merge(internal_std, iop_std, by = "SNP", all.x = TRUE, sort = FALSE)

merged$has_measured_IOP_association <- is.finite(merged$beta_measured_IOP) & is.finite(merged$se_measured_IOP)
merged$allele_match_internal_vs_iop_file <- merged$effect_allele == merged$effect_allele_iop_assoc &
  merged$other_allele == merged$other_allele_iop_assoc

merged$SBP_beta_abs_diff_internal_vs_iop_file <- abs(merged$beta_SBP - merged$beta_SBP_from_iop_file)
merged$SBP_se_abs_diff_internal_vs_iop_file <- abs(merged$se_SBP - merged$se_SBP_from_iop_file)

overlap_dataset_file <- file.path(out_dir, "phase13_5A1_SBP_component_measured_IOP_overlap_dataset.tsv")
write_tsv(merged, overlap_dataset_file)

n_internal <- nrow(internal_std)
n_iop <- nrow(iop_std)
n_overlap <- sum(merged$has_measured_IOP_association, na.rm = TRUE)
overlap_rate <- n_overlap / n_internal

max_beta_diff <- max(merged$SBP_beta_abs_diff_internal_vs_iop_file[merged$has_measured_IOP_association], na.rm = TRUE)
max_se_diff <- max(merged$SBP_se_abs_diff_internal_vs_iop_file[merged$has_measured_IOP_association], na.rm = TRUE)

all_alleles_match <- all(merged$allele_match_internal_vs_iop_file[merged$has_measured_IOP_association], na.rm = TRUE)
all_sbp_beta_match <- max_beta_diff < 1e-12
all_sbp_se_match <- max_se_diff < 1e-12

# Correlation between SBP and measured IOP SNP-exposure columns among overlapping instruments.
cor_sbp_iop <- suppressWarnings(cor(
  merged$beta_SBP[merged$has_measured_IOP_association],
  merged$beta_measured_IOP[merged$has_measured_IOP_association],
  use = "complete.obs"
))

overlap_summary <- data.frame(
  field = c(
    "internal_component_instruments",
    "measured_IOP_association_file",
    "measured_IOP_file_rows",
    "overlapping_internal_instruments_with_IOP_association",
    "overlap_rate",
    "all_alleles_match_for_overlapping_instruments",
    "all_SBP_betas_match_between_internal_and_IOP_files",
    "all_SBP_ses_match_between_internal_and_IOP_files",
    "max_abs_SBP_beta_difference",
    "max_abs_SBP_se_difference",
    "correlation_between_SBP_and_measured_IOP_SNP_associations",
    "minimum_required_overlap_for_exploratory_MVMR",
    "safe_to_build_phase13_5B_MVMR_dataset",
    "claim_boundary"
  ),
  value = c(
    n_internal,
    iop_relative,
    n_iop,
    n_overlap,
    overlap_rate,
    all_alleles_match,
    all_sbp_beta_match,
    all_sbp_se_match,
    max_beta_diff,
    max_se_diff,
    cor_sbp_iop,
    ">=300 overlapping instruments and allele/SBP-column consistency",
    n_overlap >= 300 && all_alleles_match && all_sbp_beta_match && all_sbp_se_match,
    "Any downstream MVMR is exploratory pathway-consistency sensitivity only, not mediation proof or causal confirmation."
  ),
  stringsAsFactors = FALSE
)

overlap_summary_file <- file.path(out_dir, "phase13_5A1_overlap_summary.tsv")
write_tsv(overlap_summary, overlap_summary_file)

# Missing list.
missing <- merged[!merged$has_measured_IOP_association, c("SNP", "chr", "pos", "effect_allele", "other_allele")]
missing_file <- file.path(out_dir, "phase13_5A1_internal_instruments_missing_measured_IOP.tsv")
write_tsv(missing, missing_file)

# MVMR readiness table.
readiness <- data.frame(
  component_outcome = c("GBS_nonIOPcomponent", "GBS_IOPcomponent"),
  available_columns = c(
    "beta_SBP; se_SBP; beta_measured_IOP; se_measured_IOP; beta_GBS_nonIOPcomponent; se_GBS_nonIOPcomponent",
    "beta_SBP; se_SBP; beta_measured_IOP; se_measured_IOP; beta_GBS_IOPcomponent; se_GBS_IOPcomponent"
  ),
  n_candidate_instruments = c(n_overlap, n_overlap),
  model_type = c(
    "exploratory MVMR weighted regression among overlapping SBP instruments",
    "exploratory MVMR weighted regression among overlapping SBP instruments"
  ),
  interpretation_boundary = c(
    "Assesses whether SBP coefficient attenuates after including measured IOP SNP association; not mediation proof.",
    "Assesses whether SBP coefficient attenuates after including measured IOP SNP association; not mediation proof."
  ),
  stringsAsFactors = FALSE
)

readiness_file <- file.path(out_dir, "phase13_5A1_mvmr_readiness_table.tsv")
write_tsv(readiness, readiness_file)

safe_to_proceed <- n_overlap >= 300 && all_alleles_match && all_sbp_beta_match && all_sbp_se_match

status <- data.frame(
  field = c(
    "phase",
    "exact_overlap_input_manifest_created",
    "overlap_dataset_created",
    "overlap_summary_created",
    "mvmr_readiness_table_created",
    "missing_instruments_file_created",
    "internal_component_instruments",
    "overlapping_instruments",
    "overlap_rate",
    "all_alleles_match",
    "all_SBP_columns_match",
    "safe_to_proceed_to_phase13_5B",
    "new_MR_results_created",
    "claim_level",
    "claim_upgrade_allowed"
  ),
  value = c(
    "Phase13.5A1",
    file.exists(manifest_file),
    file.exists(overlap_dataset_file),
    file.exists(overlap_summary_file),
    file.exists(readiness_file),
    file.exists(missing_file),
    n_internal,
    n_overlap,
    overlap_rate,
    all_alleles_match,
    all_sbp_beta_match && all_sbp_se_match,
    safe_to_proceed,
    FALSE,
    "HYPOTHESIS_GENERATING_NOT_CONFIRMATORY",
    "NO"
  ),
  stringsAsFactors = FALSE
)

status_file <- file.path(out_dir, "phase13_5A1_status.tsv")
write_tsv(status, status_file)

if (file.exists(master_status_file)) {
  master_status <- read.delim(master_status_file, check.names = FALSE)
  idx <- which(master_status$phase == "Phase13.5")
  if (length(idx) == 1) {
    master_status$status[idx] <- ifelse(
      safe_to_proceed,
      "PHASE13_5A1_EXACT_OVERLAP_VERIFIED_READY_FOR_EXPLORATORY_MVMR",
      "PHASE13_5A1_OVERLAP_VERIFIED_BUT_MVMR_NOT_READY"
    )
    master_status$qc_status[idx] <- ifelse(safe_to_proceed, "READY_FOR_PHASE13_5B", "LIMITATION_OR_REVIEW_REQUIRED")
    master_status$primary_output[idx] <- "phase13_5A1_overlap_summary.tsv; phase13_5A1_SBP_component_measured_IOP_overlap_dataset.tsv"
  }
  write_tsv(master_status, master_status_file)
}

message("Phase13.5A1 completed.")
message("Overlap summary: ", overlap_summary_file)
message("Overlap dataset: ", overlap_dataset_file)
message("Readiness table: ", readiness_file)
message("Status: ", status_file)
message("safe_to_proceed_to_phase13_5B = ", safe_to_proceed)
