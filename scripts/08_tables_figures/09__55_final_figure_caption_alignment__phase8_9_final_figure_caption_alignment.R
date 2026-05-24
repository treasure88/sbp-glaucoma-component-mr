# Phase 8.9: Final figure and caption alignment
# Purpose: align figure captions with the top-journal manuscript framing

options(stringsAsFactors = FALSE)

out_dir <- "../../55_final_figure_caption_alignment"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

caption_file <- file.path(out_dir, "phase8_9_final_figure_captions.tsv")
figure_inventory_file <- file.path(out_dir, "phase8_9_figure_file_inventory.tsv")
manuscript_alignment_file <- file.path(out_dir, "phase8_9_manuscript_figure_alignment.tsv")
status_file <- file.path(out_dir, "phase8_9_status.tsv")

current_manuscript <- "../../54_topjournal_manuscript_revision/phase8_7C_topjournal_manuscript_draft_heading_polished.md"

figure_files <- data.frame(
  figure_id = c(
    "Figure_1",
    "Figure_2A",
    "Figure_2A",
    "Figure_2B",
    "Figure_2B",
    "Figure_3",
    "Figure_3"
  ),
  file_path = c(
    "../../48_figures/Figure1_conceptual_model.pdf",
    "../../48_figures/Figure2A_internal_component_estimates.pdf",
    "../../48_figures/Figure2A_internal_component_estimates.png",
    "../../48_figures/Figure2B_internal_component_contrasts.pdf",
    "../../48_figures/Figure2B_internal_component_contrasts.png",
    "../../48_figures/Figure3_external_SBP_triangulation.pdf",
    "../../48_figures/Figure3_external_SBP_triangulation.png"
  ),
  intended_use = c(
    "conceptual schematic; optional/pending",
    "main or supplement",
    "preview",
    "main",
    "preview",
    "main",
    "preview"
  ),
  stringsAsFactors = FALSE
)

figure_files$file_exists <- file.exists(figure_files$file_path)
figure_files$file_size_bytes <- ifelse(
  figure_files$file_exists,
  file.info(figure_files$file_path)$size,
  NA
)

figure_files$status <- ifelse(
  figure_files$file_exists & figure_files$file_size_bytes > 0,
  "PASS",
  ifelse(figure_files$figure_id == "Figure_1", "PENDING_OR_OPTIONAL", "MISSING_REQUIRED")
)

write.table(
  figure_files,
  figure_inventory_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

captions <- data.frame(
  figure_id = c("Figure_1", "Figure_2", "Figure_3"),
  proposed_title = c(
    "Analytic framework for SBP component-divergence testing",
    "Internal SBP-centered component-divergence evidence across vascular traits",
    "External triangulation of the SBP component-divergence signal"
  ),
  caption = c(
    "Conceptual framework for evaluating whether systemic vascular liability maps differentially onto IOP-dependent and IOP-independent glaucoma genetic components. The primary analysis tested component-specific Mendelian randomization estimates and a formal IOP-dependent minus IOP-independent contrast. External outcomes were used for directional triangulation rather than confirmatory validation. This figure is optional unless a final conceptual schematic is generated.",
    "Internal component-specific Mendelian randomization estimates and formal IOP-dependent minus IOP-independent contrasts for the focused vascular/perfusion panel. Panel A shows component-specific estimates for IOP-independent and IOP-dependent glaucoma components. Panel B shows the formal component contrast, defined as the IOP-dependent estimate minus the IOP-independent estimate. SBP showed the strongest hypothesis-generating component-divergence signal, with a negative estimate for the IOP-independent component, a positive estimate for the IOP-dependent component, and an FDR-significant contrast retained from the original five-exposure r=0 contrast family. DBP was included as a blood-pressure comparator, and arterial stiffness was exploratory and low-powered.",
    "External triangulation of the SBP signal across measured IOP and glaucoma-related outcomes. The external evidence was directional but non-confirmatory: measured IOP showed a positive borderline estimate, whereas NTG and POAG estimates were negative and non-significant and RNFL/GCIPL analyses were non-confirmatory. HTG validation could not be completed because no verified full downloadable HTG summary-statistics file was available after structured source review. These results support an IOP-related interpretation more than a generalized positive glaucoma effect, but they do not establish causal confirmation."
  ),
  manuscript_role = c(
    "Optional main figure or graphical abstract",
    "Recommended main figure",
    "Recommended main figure"
  ),
  evidence_level = c(
    "Conceptual",
    "Primary hypothesis-generating internal evidence",
    "Directional non-confirmatory external triangulation"
  ),
  claim_safety_note = c(
    "Do not imply mechanism is proven.",
    "Do not state that SBP causally affects glaucoma or that all methods confirmed the result.",
    "Do not state that external validation confirmed the mechanism; HTG remains unavailable."
  ),
  stringsAsFactors = FALSE
)

write.table(
  captions,
  caption_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

manuscript_text <- if (file.exists(current_manuscript)) {
  paste(readLines(current_manuscript, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
} else {
  ""
}

alignment <- data.frame(
  manuscript_file = current_manuscript,
  manuscript_exists = file.exists(current_manuscript),
  mentions_Figure_1 = grepl("Figure 1|Figure_1", manuscript_text),
  mentions_Figure_2 = grepl("Figure 2|Figure_2", manuscript_text),
  mentions_Figure_3 = grepl("Figure 3|Figure_3", manuscript_text),
  has_SBP_component_divergence_language = grepl("component-divergence|directional divergence", manuscript_text),
  has_HTG_unavailable_language = grepl("HTG validation could not be completed", manuscript_text),
  has_LDSC_not_primary_language = grepl("LDSC was not used as primary evidence", manuscript_text),
  stringsAsFactors = FALSE
)

write.table(
  alignment,
  manuscript_alignment_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

required_missing <- figure_files$status == "MISSING_REQUIRED"

status <- data.frame(
  phase = "Phase8.9",
  status = ifelse(any(required_missing), "CHECK_REQUIRED_MISSING_REQUIRED_FIGURES", "PASSED_FIGURE_CAPTION_ALIGNMENT_CREATED"),
  note = ifelse(
    any(required_missing),
    paste("Missing required figure files:", paste(figure_files$file_path[required_missing], collapse = ";")),
    "Final figure captions and figure-file inventory created. Figure 1 remains optional/pending unless a conceptual schematic is generated."
  ),
  timestamp = as.character(Sys.time()),
  stringsAsFactors = FALSE
)

write.table(
  status,
  status_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
