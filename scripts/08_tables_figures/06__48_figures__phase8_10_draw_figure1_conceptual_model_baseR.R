# Phase 8.10: Draw Figure 1 conceptual schematic using base R
# Purpose: create a top-journal-style conceptual framework figure
# Language policy: English text in figures

options(stringsAsFactors = FALSE)

out_dir <- "../../48_figures"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

pdf_file <- file.path(out_dir, "Figure1_conceptual_model.pdf")
png_file <- file.path(out_dir, "Figure1_conceptual_model.png")
manifest_file <- file.path(out_dir, "phase8_10_figure1_file_manifest.tsv")
status_file <- file.path(out_dir, "phase8_10_status.tsv")

write_status <- function(status, note) {
  x <- data.frame(
    phase = "Phase8.10",
    status = status,
    note = note,
    timestamp = as.character(Sys.time()),
    stringsAsFactors = FALSE
  )
  write.table(x, status_file, sep = "\t", quote = FALSE, row.names = FALSE)
}

wrap_text <- function(x, width = 28) {
  paste(strwrap(x, width = width), collapse = "\n")
}

box <- function(xleft, ybottom, xright, ytop, label,
                cex = 0.82, font = 1, lwd = 1.4, fill = "white",
                border = "black", width = 28) {
  rect(xleft, ybottom, xright, ytop, col = fill, border = border, lwd = lwd)
  text(
    x = (xleft + xright) / 2,
    y = (ybottom + ytop) / 2,
    labels = wrap_text(label, width = width),
    cex = cex,
    font = font
  )
}

arrow <- function(x0, y0, x1, y1, lwd = 1.5) {
  arrows(x0, y0, x1, y1, length = 0.08, lwd = lwd, xpd = TRUE)
}

draw_figure <- function() {
  par(mar = c(1, 1, 2.5, 1), xpd = TRUE)
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))

  title(
    main = "Analytic framework for SBP component-divergence testing",
    cex.main = 1.15,
    font.main = 2
  )

  # Top exposure box
  box(
    0.30, 0.82, 0.70, 0.94,
    "Systemic vascular liability\nSBP primary signal\nDBP comparator\nArterial stiffness exploratory\nHypertension not analyzable",
    cex = 0.73,
    font = 2,
    fill = "gray95",
    width = 34
  )

  # Component boxes
  box(
    0.07, 0.55, 0.40, 0.72,
    "IOP-independent glaucoma component\nSBP estimate: negative\nbeta=-0.01499; p=0.00983; q=0.0632",
    cex = 0.72,
    fill = "gray98",
    width = 31
  )

  box(
    0.60, 0.55, 0.93, 0.72,
    "IOP-dependent glaucoma component\nSBP estimate: positive\nbeta=0.007721; p=0.0126; q=0.0632",
    cex = 0.72,
    fill = "gray98",
    width = 31
  )

  # Arrows from exposure to components
  arrow(0.40, 0.82, 0.25, 0.72)
  arrow(0.60, 0.82, 0.75, 0.72)

  # Contrast box
  box(
    0.28, 0.36, 0.72, 0.48,
    "Primary component contrast\nIOP-dependent minus IOP-independent\nbeta difference=0.02271\np=5.57e-04; q=0.00279",
    cex = 0.74,
    font = 2,
    fill = "gray92",
    width = 38
  )

  arrow(0.32, 0.55, 0.41, 0.48)
  arrow(0.68, 0.55, 0.59, 0.48)

  # Sensitivity audit box
  box(
    0.06, 0.18, 0.45, 0.30,
    "Targeted stability checks\nContrast audit passed\nCovariance sensitivity stable\nLeave-one-chromosome-out stable\nInstrument-quality subsets stable",
    cex = 0.70,
    fill = "white",
    width = 34
  )

  # External triangulation box
  box(
    0.55, 0.18, 0.94, 0.30,
    "External triangulation\nMeasured IOP: positive borderline\nNTG/POAG: negative non-significant\nRNFL/GCIPL: non-confirmatory\nHTG: unavailable",
    cex = 0.70,
    fill = "white",
    width = 34
  )

  arrow(0.40, 0.36, 0.27, 0.30)
  arrow(0.60, 0.36, 0.73, 0.30)

  # Bottom interpretation box
  box(
    0.20, 0.045, 0.80, 0.125,
    "Interpretation: statistically stable, hypothesis-generating component-divergence evidence; directional but non-confirmatory external support; not causal confirmation.",
    cex = 0.72,
    fill = "gray95",
    width = 65
  )

  arrow(0.27, 0.18, 0.42, 0.125)
  arrow(0.73, 0.18, 0.58, 0.125)

  # Small panel letters / labels
  text(0.02, 0.97, "Figure 1", font = 2, cex = 0.9, adj = 0)
}

tryCatch({
  pdf(pdf_file, width = 8.5, height = 6.2, onefile = TRUE)
  draw_figure()
  dev.off()

  png(png_file, width = 1800, height = 1300, res = 200)
  draw_figure()
  dev.off()

  files <- data.frame(
    figure_id = c("Figure1_conceptual_model", "Figure1_conceptual_model"),
    file_path = c(pdf_file, png_file),
    file_type = c("pdf", "png"),
    file_exists = file.exists(c(pdf_file, png_file)),
    file_size_bytes = ifelse(file.exists(c(pdf_file, png_file)), file.info(c(pdf_file, png_file))$size, NA),
    status = ifelse(file.exists(c(pdf_file, png_file)) & file.info(c(pdf_file, png_file))$size > 0, "PASS", "FAIL"),
    stringsAsFactors = FALSE
  )

  write.table(files, manifest_file, sep = "\t", quote = FALSE, row.names = FALSE)

  if (all(files$status == "PASS")) {
    write_status(
      "PASSED_FIGURE1_CREATED",
      "Figure 1 conceptual model generated successfully in PDF and PNG formats."
    )
  } else {
    write_status(
      "CHECK_REQUIRED",
      "Figure 1 script ran, but at least one output file is missing or empty."
    )
  }

}, error = function(e) {
  write_status("FAILED", paste("Error:", conditionMessage(e)))
  stop(e)
})
