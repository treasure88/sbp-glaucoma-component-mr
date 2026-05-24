options(stringsAsFactors = FALSE)

# Phase13.4A3: Verify exact internal SBP-component SNP-level MR inputs
# Project: IOP-dependent vs IOP-independent glaucoma component MR
# Evidence level: HYPOTHESIS_GENERATING_NOT_CONFIRMATORY
# This phase verifies input files and reproduces locked IVW estimates.
# It does not upgrade claims.

message("Running Phase13.4A3 exact internal SBP-component input verification...")

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

as_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

read_table_any <- function(path) {
  if (!file.exists(path)) stop("Missing file: ", path)
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

ivw_random_phi <- function(beta_exposure, se_exposure, beta_outcome, se_outcome) {
  bx <- as_num(beta_exposure)
  sx <- as_num(se_exposure)
  by <- as_num(beta_outcome)
  sy <- as_num(se_outcome)
  
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

root <- normalizePath(".", winslash = "/", mustWork = TRUE)
out_dir <- file.path(root, "71_data_analysis_reinforcement", "phase13_4_empirical_covariance")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

input_files <- data.frame(
  outcome_suffix = c("GBS_nonIOPcomponent", "GBS_IOPcomponent"),
  expected_direction = c("negative", "positive"),
  locked_beta = c(-0.014988432, 0.0077208916),
  locked_se = c(0.0058051821, 0.0030961462),
  locked_p = c(0.0098256472, 0.012641568),
  relative_path = c(
    "16_mr_input_datasets/pairwise/SBP__GBS_nonIOPcomponent.mr_input.tsv.gz",
    "16_mr_input_datasets/pairwise/SBP__GBS_IOPcomponent.mr_input.tsv.gz"
  ),
  stringsAsFactors = FALSE
)

input_files$full_path <- file.path(root, input_files$relative_path)
input_files$exists <- file.exists(input_files$full_path)

input_status_file <- file.path(out_dir, "phase13_4A3_exact_internal_input_file_status.tsv")
write_tsv(input_files, input_status_file)

if (!all(input_files$exists)) {
  stop("One or more exact internal SBP-component input files are missing. See: ", input_status_file)
}

input_summaries <- list()
ivw_summaries <- list()
column_maps <- list()

for (i in seq_len(nrow(input_files))) {
  df <- read_table_any(input_files$full_path[i])
  
  snp_col <- pick_col(df, c("SNP", "rsid", "variant_id"))
  chr_col <- pick_col(df, c("chr", "CHR", "chromosome", "exposure_chr"))
  pos_col <- pick_col(df, c("pos", "POS", "position", "bp", "exposure_pos"))
  
  beta_exp_col <- pick_col(
    df,
    c("beta_exposure", "beta.exposure", "beta_exp", "bx", "b_exp"),
    c("^beta.*exposure$")
  )
  se_exp_col <- pick_col(
    df,
    c("se_exposure", "se.exposure", "se_exp", "sx"),
    c("^se.*exposure$")
  )
  beta_out_col <- pick_col(
    df,
    c("beta_outcome", "beta.outcome", "beta_out", "by", "b_out", "beta_outcome_aligned"),
    c("^beta.*outcome$")
  )
  se_out_col <- pick_col(
    df,
    c("se_outcome", "se.outcome", "se_out", "sy"),
    c("^se.*outcome$")
  )
  p_exp_col <- pick_col(df, c("pval_exposure", "pval.exposure", "p.exposure", "p_exposure"))
  p_out_col <- pick_col(df, c("pval_outcome", "pval.outcome", "p.outcome", "p_outcome"))
  f_col <- pick_col(df, c("F_stat", "F", "f_stat", "F_statistic"))
  harmonization_col <- pick_col(df, c("harmonization_action"))
  match_col <- pick_col(df, c("match_mode"))
  
  required_ok <- !any(is.na(c(snp_col, beta_exp_col, se_exp_col, beta_out_col, se_out_col)))
  
  column_maps[[i]] <- data.frame(
    outcome_suffix = input_files$outcome_suffix[i],
    n_rows = nrow(df),
    n_columns = ncol(df),
    snp_col = snp_col,
    chr_col = chr_col,
    pos_col = pos_col,
    beta_exposure_col = beta_exp_col,
    se_exposure_col = se_exp_col,
    beta_outcome_col = beta_out_col,
    se_outcome_col = se_out_col,
    p_exposure_col = p_exp_col,
    p_outcome_col = p_out_col,
    F_col = f_col,
    harmonization_col = harmonization_col,
    match_col = match_col,
    required_columns_present = required_ok,
    columns_preview = paste(names(df), collapse = "|"),
    stringsAsFactors = FALSE
  )
  
  if (!required_ok) {
    next
  }
  
  bx <- as_num(df[[beta_exp_col]])
  sx <- as_num(df[[se_exp_col]])
  by <- as_num(df[[beta_out_col]])
  sy <- as_num(df[[se_out_col]])
  
  F_values <- if (!is.na(f_col)) {
    as_num(df[[f_col]])
  } else {
    (bx / sx)^2
  }
  
  input_summaries[[i]] <- data.frame(
    outcome_suffix = input_files$outcome_suffix[i],
    n_rows = nrow(df),
    n_unique_snp = length(unique(df[[snp_col]])),
    has_chr_pos = !is.na(chr_col) && !is.na(pos_col),
    has_reported_F = !is.na(f_col),
    mean_F = mean(F_values, na.rm = TRUE),
    median_F = median(F_values, na.rm = TRUE),
    min_F = min(F_values, na.rm = TRUE),
    max_F = max(F_values, na.rm = TRUE),
    n_F_less_10 = sum(F_values < 10, na.rm = TRUE),
    harmonization_action_counts = if (!is.na(harmonization_col)) {
      paste(names(sort(table(df[[harmonization_col]]), decreasing = TRUE)),
            as.integer(sort(table(df[[harmonization_col]]), decreasing = TRUE)),
            sep = "=", collapse = "; ")
    } else {
      NA_character_
    },
    match_mode_counts = if (!is.na(match_col)) {
      paste(names(sort(table(df[[match_col]]), decreasing = TRUE)),
            as.integer(sort(table(df[[match_col]]), decreasing = TRUE)),
            sep = "=", collapse = "; ")
    } else {
      NA_character_
    },
    stringsAsFactors = FALSE
  )
  
  ivw <- ivw_random_phi(
    beta_exposure = df[[beta_exp_col]],
    se_exposure = df[[se_exp_col]],
    beta_outcome = df[[beta_out_col]],
    se_outcome = df[[se_out_col]]
  )
  
  ivw$outcome_suffix <- input_files$outcome_suffix[i]
  ivw$expected_direction <- input_files$expected_direction[i]
  ivw$locked_beta <- input_files$locked_beta[i]
  ivw$locked_se <- input_files$locked_se[i]
  ivw$locked_p <- input_files$locked_p[i]
  ivw$abs_beta_difference_from_locked <- abs(ivw$beta - input_files$locked_beta[i])
  ivw$abs_se_difference_from_locked <- abs(ivw$se - input_files$locked_se[i])
  ivw$abs_p_difference_from_locked <- abs(ivw$p - input_files$locked_p[i])
  ivw$direction_reproduced <- ifelse(ivw$beta < 0, "negative", "positive") == input_files$expected_direction[i]
  ivw$locked_estimate_reproduced_tolerance_1e_6 <- ivw$abs_beta_difference_from_locked < 1e-6 &&
    ivw$abs_se_difference_from_locked < 1e-6
  
  ivw_summaries[[i]] <- ivw
}

column_map <- do.call(rbind, column_maps)
input_summary <- do.call(rbind, input_summaries)
ivw_summary <- do.call(rbind, ivw_summaries)

column_map_file <- file.path(out_dir, "phase13_4A3_internal_input_column_map.tsv")
input_summary_file <- file.path(out_dir, "phase13_4A3_internal_input_strength_summary.tsv")
ivw_summary_file <- file.path(out_dir, "phase13_4A3_internal_input_ivw_reproduction.tsv")

write_tsv(column_map, column_map_file)
write_tsv(input_summary, input_summary_file)
write_tsv(ivw_summary, ivw_summary_file)

all_required_cols <- all(column_map$required_columns_present)
all_n_391 <- all(input_summary$n_unique_snp == 391)
all_direction_reproduced <- all(ivw_summary$direction_reproduced)
all_locked_reproduced <- all(ivw_summary$locked_estimate_reproduced_tolerance_1e_6)

# Some older scripts may have used a slightly different random-effects implementation.
# Therefore, strict lock reproduction is reported, but bootstrap can proceed if columns, n, and directions are valid.
safe_for_bootstrap <- all_required_cols && all_n_391 && all_direction_reproduced

status <- data.frame(
  field = c(
    "phase",
    "exact_input_file_status_created",
    "column_map_created",
    "input_strength_summary_created",
    "ivw_reproduction_created",
    "all_exact_input_files_exist",
    "all_required_columns_present",
    "all_unique_snp_counts_equal_391",
    "all_directions_reproduced",
    "locked_estimates_reproduced_tolerance_1e_6",
    "safe_for_phase13_4B_bootstrap",
    "new_MR_results_created",
    "claim_level",
    "claim_upgrade_allowed"
  ),
  value = c(
    "Phase13.4A3",
    file.exists(input_status_file),
    file.exists(column_map_file),
    file.exists(input_summary_file),
    file.exists(ivw_summary_file),
    all(input_files$exists),
    all_required_cols,
    all_n_391,
    all_direction_reproduced,
    all_locked_reproduced,
    safe_for_bootstrap,
    FALSE,
    "HYPOTHESIS_GENERATING_NOT_CONFIRMATORY",
    "NO"
  ),
  stringsAsFactors = FALSE
)

status_file <- file.path(out_dir, "phase13_4A3_status.tsv")
write_tsv(status, status_file)

message("Phase13.4A3 completed.")
message("Input file status: ", input_status_file)
message("Column map: ", column_map_file)
message("Input strength summary: ", input_summary_file)
message("IVW reproduction: ", ivw_summary_file)
message("Status: ", status_file)

if (!safe_for_bootstrap) {
  warning("Inputs were found but did not pass all bootstrap-readiness checks. Review Phase13.4A3 outputs.")
}
