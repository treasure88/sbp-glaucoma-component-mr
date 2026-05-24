options(stringsAsFactors = FALSE)

out_dir <- "../../58_component_divergence_reinforcement"

index_file <- file.path(out_dir, "phase10_5_supplementary_table_index.tsv")
summary_file <- file.path(out_dir, "phase10_5_robustness_integration_summary.txt")
claim_file <- file.path(out_dir, "phase10_5_claim_boundary.tsv")
status_file <- file.path(out_dir, "phase10_5_status.tsv")

files <- data.frame(
  supplement_item = c(
    "Contrast audit and covariance sensitivity",
    "Single-SNP/locus influence assessment",
    "Alternative SBP instrument and influence-filtered sensitivity summary",
    "Phase 10.4 manuscript-safe wording",
    "Phase 10.4 recommendation"
  ),
  phase = c("Phase10.2", "Phase10.3", "Phase10.4", "Phase10.4", "Phase10.4"),
  file_path = c(
    file.path(out_dir, "phase10_2_contrast_audit_summary.tsv"),
    file.path(out_dir, "phase10_3_status.tsv"),
    file.path(out_dir, "phase10_4_alternative_sbp_instrument_summary.tsv"),
    file.path(out_dir, "phase10_4_manuscript_safe_text.txt"),
    file.path(out_dir, "phase10_4_recommendation.txt")
  ),
  manuscript_use = c(
    "Supplementary Methods/Table: contrast definition, covariance assumptions, and FDR-family audit.",
    "Supplementary Table: influence diagnostics showing no single SNP or locus dominance.",
    "Supplementary Table: alternative SBP instrument definitions and influence-filtered sensitivity sets.",
    "Source text for Results/Methods robustness wording.",
    "Source text for final interpretation and claim boundary."
  ),
  stringsAsFactors = FALSE
)

files$file_exists <- file.exists(files$file_path)
files$file_size_bytes <- ifelse(files$file_exists, file.info(files$file_path)$size, NA)
files$status <- ifelse(files$file_exists & files$file_size_bytes > 0, "PASS", "CHECK_MISSING_OR_EMPTY")

write.table(files, index_file, sep = "\t", quote = FALSE, row.names = FALSE)

claim <- data.frame(
  claim_domain = c(
    "Primary interpretation",
    "Alternative instrument definitions",
    "Influence-filtered analyses",
    "Component-specific estimates",
    "External/causal interpretation"
  ),
  allowed_wording = c(
    "SBP showed statistically stable, hypothesis-generating component-divergence evidence.",
    "The positive IOP-dependent minus IOP-independent contrast was directionally preserved across stricter SBP instrument-quality definitions.",
    "Influence-filtered analyses did not indicate dependence on a single top SNP or locus.",
    "Component sign pattern was generally preserved, but individual component-specific estimates attenuated in some stricter subsets.",
    "Findings support robustness of the internal contrast pattern but do not establish causal confirmation."
  ),
  avoid_wording = c(
    "SBP causally affects glaucoma.",
    "Alternative instruments confirmed the mechanism.",
    "Influence-filtered analyses are independent external validation.",
    "Both component-specific estimates remained significant across every subset.",
    "External triangulation or sensitivity analyses confirmed causality."
  ),
  stringsAsFactors = FALSE
)

write.table(claim, claim_file, sep = "\t", quote = FALSE, row.names = FALSE)

txt <- c(
  "Phase 10.5 robustness integration summary",
  "=========================================",
  "",
  "Overall interpretation:",
  "The SBP component-divergence signal was supported by targeted reinforcement analyses. The formal IOP-dependent minus IOP-independent contrast was reproduced and remained directionally stable across covariance assumptions, leave-one-chromosome/locus influence checks, stricter SBP instrument-quality definitions, and influence-filtered sensitivity sets.",
  "",
  "Manuscript-safe Results wording:",
  "Targeted robustness analyses supported the stability of the SBP component-divergence pattern. The positive IOP-dependent minus IOP-independent contrast was retained across stricter instrument-quality subsets and across influence-filtered analyses excluding top contributing SNPs or loci. These analyses reduce the likelihood that the signal depends entirely on a single instrument definition or top contributor, although they do not establish causal confirmation.",
  "",
  "Manuscript-safe caveats:",
  "- Influence-filtered analyses are perturbation analyses, not independent GWAS validation.",
  "- Some stricter subsets attenuated one component-specific estimate to nominal non-significance.",
  "- Non-palindromic restrictions did not introduce additional exclusions in the locked dataset.",
  "- The evidence remains hypothesis-generating and should not be described as causal confirmation.",
  "",
  "Recommended supplement package:",
  "Use phase10_4_alternative_sbp_instrument_summary.tsv as the final supplement-ready table for alternative SBP instrument definitions and influence-filtered sensitivity sets, together with the Phase 10.2/10.3 audit outputs."
)

writeLines(txt, summary_file)

all_required_pass <- all(files$status != "CHECK_MISSING_OR_EMPTY")

status <- data.frame(
  phase = "Phase10.5",
  status = ifelse(all_required_pass, "PASSED_FINAL_ROBUSTNESS_PACKAGE_INDEX_CREATED", "CHECK_REQUIRED_MISSING_FILES"),
  n_indexed_items = nrow(files),
  n_pass = sum(files$status == "PASS"),
  n_check = sum(files$status == "CHECK_MISSING_OR_EMPTY"),
  output_directory = out_dir,
  timestamp = as.character(Sys.time()),
  stringsAsFactors = FALSE
)

write.table(status, status_file, sep = "\t", quote = FALSE, row.names = FALSE)
