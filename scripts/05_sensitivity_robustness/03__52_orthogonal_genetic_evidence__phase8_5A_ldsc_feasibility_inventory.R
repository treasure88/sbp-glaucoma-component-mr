# Phase 8.5A: LDSC / genetic correlation feasibility inventory
# Purpose: identify whether local genome-wide summary statistics exist for SBP and glaucoma components
# Note: LDSC requires genome-wide sumstats, not clumped instruments or MR input files.

options(stringsAsFactors = FALSE)

project_root <- "../.."
out_dir <- "../../52_orthogonal_genetic_evidence"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

inventory_file <- file.path(out_dir, "phase8_5A_ldsc_feasibility_inventory.tsv")
decision_file <- file.path(out_dir, "phase8_5A_ldsc_feasibility_decision.tsv")
status_file <- file.path(out_dir, "phase8_5A_ldsc_feasibility_status.tsv")

write_status <- function(status, note) {
  x <- data.frame(
    phase = "Phase8.5A",
    status = status,
    note = note,
    timestamp = as.character(Sys.time()),
    stringsAsFactors = FALSE
  )
  write.table(x, status_file, sep = "\t", quote = FALSE, row.names = FALSE)
}

read_preview <- function(path, n_max = 30) {
  tryCatch(
    read.delim(
      path,
      header = TRUE,
      sep = "\t",
      quote = "",
      comment.char = "",
      check.names = FALSE,
      fill = TRUE,
      nrows = n_max
    ),
    error = function(e) NULL
  )
}

classify_cols <- function(cols) {
  lower <- tolower(cols)

  has_snp <- any(lower %in% c("snp", "rsid", "rs_id", "variant", "markername"))
  has_chr <- any(lower %in% c("chr", "chrom", "chromosome", "#chrom"))
  has_pos <- any(lower %in% c("pos", "position", "bp", "base_pair_location"))
  has_a1 <- any(lower %in% c("a1", "effect_allele", "ea", "allele1", "tested_allele"))
  has_a2 <- any(lower %in% c("a2", "other_allele", "oa", "allele2", "non_tested_allele"))
  has_beta <- any(lower %in% c("beta", "b", "effect", "log_odds", "logor") | grepl("beta", lower))
  has_se <- any(lower %in% c("se", "stderr", "standard_error") | grepl("^se$|standard", lower))
  has_p <- any(lower %in% c("p", "pval", "p_value", "pvalue") | grepl("^p$|pval|p.value|p_value", lower))
  has_z <- any(lower %in% c("z", "zscore", "z_score"))
  has_n <- any(lower %in% c("n", "samplesize", "sample_size", "n_total"))
  has_eaf <- any(lower %in% c("eaf", "af", "freq", "effect_allele_frequency"))

  score <- sum(c(has_snp, has_chr, has_pos, has_a1, has_a2, has_beta, has_se, has_p, has_z, has_n, has_eaf))

  list(
    has_snp = has_snp,
    has_chr = has_chr,
    has_pos = has_pos,
    has_a1 = has_a1,
    has_a2 = has_a2,
    has_beta = has_beta,
    has_se = has_se,
    has_p = has_p,
    has_z = has_z,
    has_n = has_n,
    has_eaf = has_eaf,
    score = score
  )
}

files <- list.files(
  project_root,
  pattern = "\\.(tsv|txt|csv|gz)$",
  recursive = TRUE,
  full.names = TRUE
)

lower_files <- tolower(files)
base_lower <- tolower(basename(files))

trait_hit <- grepl("sbp|systolic|gbs_iopcomponent|gbs_noniopcomponent|noniop|iopcomponent|glaucoma", lower_files)

# Exclude obvious derived MR-only or result/metadata files unless they might be raw/standardized genome-wide files.
exclude_derived <- grepl(
  "mr_input|analysis_ready|harmonized|pilot|clumped|candidate_instruments|plink_clump|results|status|summary|manifest|figure|table|caption|runtime|qc|lock|audit|sensitivity|validation_grid",
  lower_files
)

# Keep raw/standardized large files and possible component files.
potential <- files[trait_hit & (!exclude_derived | grepl("standardized|raw|sumstat|gwas|component", lower_files))]

rows <- list()

for (f in potential) {
  info <- file.info(f)
  d <- read_preview(f, n_max = 30)

  if (is.null(d)) {
    rows[[length(rows) + 1]] <- data.frame(
      file_path = f,
      file_name = basename(f),
      file_size_bytes = info$size,
      readable_as_table = FALSE,
      n_columns = NA_integer_,
      columns = NA_character_,
      trait_guess = NA_character_,
      has_snp = NA,
      has_chr = NA,
      has_pos = NA,
      has_a1 = NA,
      has_a2 = NA,
      has_beta = NA,
      has_se = NA,
      has_p = NA,
      has_z = NA,
      has_n = NA,
      has_eaf = NA,
      sumstats_column_score = NA_integer_,
      is_large_file = info$size > 10000000,
      likely_genomewide_sumstats = FALSE,
      likely_mr_only = grepl("mr_input|clumped|candidate|analysis_ready|harmonized", tolower(f)),
      notes = "Unreadable preview.",
      stringsAsFactors = FALSE
    )
    next
  }

  cc <- classify_cols(names(d))
  lf <- tolower(f)

  trait_guess <- ifelse(
    grepl("sbp|systolic", lf),
    "SBP",
    ifelse(
      grepl("gbs_noniopcomponent|noniop", lf),
      "GBS_nonIOPcomponent",
      ifelse(
        grepl("gbs_iopcomponent|iopcomponent", lf),
        "GBS_IOPcomponent",
        "UNCLASSIFIED_GLAUCOMA_RELATED"
      )
    )
  )

  is_large <- info$size > 10000000
  has_ldsc_core <- cc$has_snp && cc$has_a1 && cc$has_a2 && (cc$has_z || (cc$has_beta && cc$has_se) || cc$has_p)
  likely_mr_only <- grepl("mr_input|clumped|candidate|analysis_ready|harmonized", tolower(f))
  likely_genomewide <- is_large && has_ldsc_core && !likely_mr_only

  rows[[length(rows) + 1]] <- data.frame(
    file_path = f,
    file_name = basename(f),
    file_size_bytes = info$size,
    readable_as_table = TRUE,
    n_columns = ncol(d),
    columns = paste(names(d), collapse = ";"),
    trait_guess = trait_guess,
    has_snp = cc$has_snp,
    has_chr = cc$has_chr,
    has_pos = cc$has_pos,
    has_a1 = cc$has_a1,
    has_a2 = cc$has_a2,
    has_beta = cc$has_beta,
    has_se = cc$has_se,
    has_p = cc$has_p,
    has_z = cc$has_z,
    has_n = cc$has_n,
    has_eaf = cc$has_eaf,
    sumstats_column_score = cc$score,
    is_large_file = is_large,
    likely_genomewide_sumstats = likely_genomewide,
    likely_mr_only = likely_mr_only,
    notes = ifelse(likely_genomewide, "Potential LDSC-compatible genome-wide sumstats candidate.", "Not clearly LDSC-compatible genome-wide sumstats."),
    stringsAsFactors = FALSE
  )
}

inventory <- if (length(rows) > 0) do.call(rbind, rows) else data.frame()

if (nrow(inventory) > 0) {
  inventory <- inventory[order(inventory$trait_guess, -inventory$likely_genomewide_sumstats, -inventory$file_size_bytes), , drop = FALSE]
}

write.table(inventory, inventory_file, sep = "\t", quote = FALSE, row.names = FALSE)

# Decision by trait
traits_needed <- c("SBP", "GBS_IOPcomponent", "GBS_nonIOPcomponent")

decision_rows <- list()

for (tr in traits_needed) {
  sub <- inventory[inventory$trait_guess == tr, , drop = FALSE]

  n_candidates <- nrow(sub)
  n_likely <- if (n_candidates > 0) sum(sub$likely_genomewide_sumstats, na.rm = TRUE) else 0

  best_file <- if (n_likely > 0) {
    sub$file_path[which(sub$likely_genomewide_sumstats)[1]]
  } else if (n_candidates > 0) {
    sub$file_path[1]
  } else {
    NA_character_
  }

  decision_rows[[length(decision_rows) + 1]] <- data.frame(
    trait = tr,
    n_local_candidates = n_candidates,
    n_likely_genomewide_sumstats = n_likely,
    best_candidate_file = best_file,
    ldsc_ready = n_likely > 0,
    decision_note = ifelse(
      n_likely > 0,
      "Potential LDSC-compatible genome-wide summary statistics found; manual format/QC still required.",
      "No clearly LDSC-compatible genome-wide summary statistics found locally."
    ),
    stringsAsFactors = FALSE
  )
}

decision <- do.call(rbind, decision_rows)

overall_ready <- all(decision$ldsc_ready)

decision$overall_ldsc_feasibility <- ifelse(
  overall_ready,
  "POTENTIALLY_FEASIBLE_AFTER_FORMAT_QC",
  "NOT_FEASIBLE_WITH_CURRENT_LOCAL_FILES"
)

write.table(decision, decision_file, sep = "\t", quote = FALSE, row.names = FALSE)

if (overall_ready) {
  write_status(
    "PASSED_POTENTIALLY_FEASIBLE",
    "Potential local genome-wide sumstats found for SBP and both glaucoma components. Manual LDSC format/QC required next."
  )
} else {
  missing_traits <- paste(decision$trait[!decision$ldsc_ready], collapse = ";")
  write_status(
    "NOT_FEASIBLE_WITH_CURRENT_LOCAL_FILES",
    paste0("No clearly LDSC-compatible genome-wide summary statistics found for: ", missing_traits, ". Do not run LDSC unless suitable sumstats are added.")
  )
}
