# Phase 8.1R: Draw main figures using base R only
# Project: IOP-dependent vs IOP-independent glaucoma component MR
# Language policy: English labels for figures/statistical outputs

options(stringsAsFactors = FALSE)

out_dir <- "../../48_figures"
input_dir <- "../../46_figure_table_preparation"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

status_file <- file.path(out_dir, "phase8_1R_status.tsv")
manifest_file <- file.path(out_dir, "phase8_1R_figure_file_manifest.tsv")
qc_file <- file.path(out_dir, "phase8_1R_figure_input_qc.tsv")

figure2_input <- file.path(input_dir, "figure2_internal_component_contrast_plot_input.tsv")
figure3_input <- file.path(input_dir, "figure3_external_validation_plot_input.tsv")

write_status <- function(status, note) {
  x <- data.frame(
    phase = "Phase8.1R",
    status = status,
    note = note,
    timestamp = as.character(Sys.time()),
    stringsAsFactors = FALSE
  )
  write.table(x, status_file, sep = "\t", quote = FALSE, row.names = FALSE)
}

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

to_num <- function(x) {
  suppressWarnings(as.numeric(x))
}

fmt_num <- function(x, digits = 3) {
  ifelse(is.na(x), "NA", formatC(x, digits = digits, format = "f"))
}

fmt_p <- function(x) {
  if (length(x) == 0 || is.na(x)) return("NA")
  if (x < 0.001) return(formatC(x, format = "e", digits = 2))
  formatC(x, format = "f", digits = 3)
}

make_qc <- function() {
  required_f2 <- c(
    "exposure_id",
    "beta_nonIOP",
    "se_nonIOP",
    "p_nonIOP",
    "beta_IOP",
    "se_IOP",
    "p_IOP",
    "contrast_beta_diff_r0",
    "contrast_se_r0",
    "contrast_p_r0",
    "contrast_q_r0",
    "contrast_direction",
    "evidence_label"
  )

  required_f3 <- c(
    "exposure_id",
    "outcome_id",
    "validation_layer",
    "beta",
    "se",
    "pval",
    "qval",
    "direction",
    "evidence_label"
  )

  qc_rows <- list()

  for (item in list(
    list(name = "figure2_internal_component_contrast_plot_input", path = figure2_input, required = required_f2),
    list(name = "figure3_external_validation_plot_input", path = figure3_input, required = required_f3)
  )) {
    exists_flag <- file.exists(item$path)
    n_rows <- NA_integer_
    missing_cols <- "FILE_NOT_FOUND"
    status <- "FAIL"

    if (exists_flag) {
      d <- read_tsv(item$path)
      n_rows <- nrow(d)
      missing <- setdiff(item$required, names(d))
      missing_cols <- if (length(missing) == 0) "NONE" else paste(missing, collapse = ";")
      status <- if (length(missing) == 0 && n_rows > 0) "PASS" else "FAIL"
    }

    qc_rows[[length(qc_rows) + 1]] <- data.frame(
      input_name = item$name,
      path = item$path,
      file_exists = exists_flag,
      n_rows = n_rows,
      missing_required_columns = missing_cols,
      status = status,
      stringsAsFactors = FALSE
    )
  }

  qc <- do.call(rbind, qc_rows)
  write.table(qc, qc_file, sep = "\t", quote = FALSE, row.names = FALSE)
  qc
}

save_dual <- function(base_name, plot_fun, width = 8, height = 5.5, png_width = 1800, png_height = 1300, res = 220) {
  pdf_path <- file.path(out_dir, paste0(base_name, ".pdf"))
  png_path <- file.path(out_dir, paste0(base_name, ".png"))

  pdf(pdf_path, width = width, height = height, useDingbats = FALSE)
  plot_fun()
  dev.off()

  png(png_path, width = png_width, height = png_height, res = res)
  plot_fun()
  dev.off()

  c(pdf_path, png_path)
}

plot_figure2A <- function(df) {
  d <- df
  exposure_order <- c("SBP", "DBP", "ART_STIFFNESS")
  d <- d[match(intersect(exposure_order, d$exposure_id), d$exposure_id), , drop = FALSE]

  rows <- data.frame(
    exposure_id = rep(d$exposure_id, each = 2),
    component = rep(c("IOP-independent component", "IOP-dependent component"), times = nrow(d)),
    beta = as.vector(t(cbind(to_num(d$beta_nonIOP), to_num(d$beta_IOP)))),
    se = as.vector(t(cbind(to_num(d$se_nonIOP), to_num(d$se_IOP)))),
    p = as.vector(t(cbind(to_num(d$p_nonIOP), to_num(d$p_IOP)))),
    stringsAsFactors = FALSE
  )

  rows$lower <- rows$beta - 1.96 * rows$se
  rows$upper <- rows$beta + 1.96 * rows$se
  rows$label <- paste(rows$exposure_id, rows$component, sep = " | ")

  y <- rev(seq_len(nrow(rows)))
  xlim <- range(c(rows$lower, rows$upper, 0), na.rm = TRUE)
  pad <- diff(xlim) * 0.18
  xlim <- c(xlim[1] - pad, xlim[2] + pad)

  par(mar = c(5, 12, 4, 2))
  plot(
    NA,
    xlim = xlim,
    ylim = c(0.5, length(y) + 0.5),
    yaxt = "n",
    xlab = "MR estimate",
    ylab = "",
    main = "Component-specific MR estimates",
    bty = "n"
  )
  abline(v = 0, lty = 2)
  axis(2, at = y, labels = rows$label, las = 2, tick = FALSE, cex.axis = 0.75)

  for (i in seq_len(nrow(rows))) {
    segments(rows$lower[i], y[i], rows$upper[i], y[i], lwd = 1.5)
    points(rows$beta[i], y[i], pch = ifelse(grepl("IOP-dependent", rows$component[i]), 17, 16), cex = 1.1)
  }

  legend(
    "bottomright",
    legend = c("IOP-independent component", "IOP-dependent component"),
    pch = c(16, 17),
    bty = "n",
    cex = 0.85
  )
}

plot_figure2B <- function(df) {
  d <- df
  exposure_order <- c("SBP", "DBP", "ART_STIFFNESS")
  d <- d[match(intersect(exposure_order, d$exposure_id), d$exposure_id), , drop = FALSE]

  beta <- to_num(d$contrast_beta_diff_r0)
  se <- to_num(d$contrast_se_r0)
  p <- to_num(d$contrast_p_r0)
  q <- to_num(d$contrast_q_r0)

  lower <- beta - 1.96 * se
  upper <- beta + 1.96 * se

  y <- rev(seq_len(nrow(d)))
  xlim <- range(c(lower, upper, 0), na.rm = TRUE)
  pad <- diff(xlim) * 0.30
  xlim <- c(xlim[1] - pad, xlim[2] + pad)

  par(mar = c(5, 9, 4, 7))
  plot(
    NA,
    xlim = xlim,
    ylim = c(0.5, length(y) + 0.5),
    yaxt = "n",
    xlab = "Contrast estimate: IOP-dependent minus IOP-independent",
    ylab = "",
    main = "Internal component contrast",
    bty = "n"
  )
  abline(v = 0, lty = 2)
  axis(2, at = y, labels = d$exposure_id, las = 2, tick = FALSE)

  for (i in seq_len(nrow(d))) {
    segments(lower[i], y[i], upper[i], y[i], lwd = 1.5)
    points(beta[i], y[i], pch = 16, cex = ifelse(d$exposure_id[i] == "SBP", 1.4, 1.1))
    text(
      x = xlim[2],
      y = y[i],
      labels = paste0("p=", fmt_p(p[i]), "; q=", fmt_p(q[i])),
      pos = 2,
      cex = 0.78,
      xpd = TRUE
    )
  }

  mtext("Positive values indicate a more positive estimate for the IOP-dependent component", side = 1, line = 3.5, cex = 0.75)
}

plot_figure3 <- function(df) {
  d <- df[df$exposure_id == "SBP", , drop = FALSE]

  outcome_order <- c("IOP", "HTG", "NTG", "POAG", "RNFL", "GCIPL")
  d <- d[match(intersect(outcome_order, d$outcome_id), d$outcome_id), , drop = FALSE]

  beta <- to_num(d$beta)
  se <- to_num(d$se)
  p <- to_num(d$pval)
  q <- to_num(d$qval)

  lower <- beta - 1.96 * se
  upper <- beta + 1.96 * se

  finite_values <- c(lower, upper, beta, 0)
  finite_values <- finite_values[is.finite(finite_values)]
  if (length(finite_values) == 0) finite_values <- c(-1, 1)

  xlim <- range(finite_values, na.rm = TRUE)
  pad <- max(diff(xlim) * 0.35, 0.01)
  xlim <- c(xlim[1] - pad, xlim[2] + pad)

  y <- rev(seq_len(nrow(d)))

  par(mar = c(5, 8, 4, 8))
  plot(
    NA,
    xlim = xlim,
    ylim = c(0.5, length(y) + 0.5),
    yaxt = "n",
    xlab = "MR estimate",
    ylab = "",
    main = "External triangulation of the SBP signal",
    bty = "n"
  )
  abline(v = 0, lty = 2)
  axis(2, at = y, labels = d$outcome_id, las = 2, tick = FALSE)

  for (i in seq_len(nrow(d))) {
    if (is.finite(beta[i]) && is.finite(se[i])) {
      segments(lower[i], y[i], upper[i], y[i], lwd = 1.5)
      points(beta[i], y[i], pch = 16, cex = 1.1)
      label <- paste0("p=", fmt_p(p[i]), "; q=", fmt_p(q[i]))
    } else {
      text(0, y[i], labels = "Not available", cex = 0.85)
      label <- "Not available"
    }

    text(
      x = xlim[2],
      y = y[i],
      labels = label,
      pos = 2,
      cex = 0.78,
      xpd = TRUE
    )
  }

  mtext("Directional external evidence; not confirmatory", side = 1, line = 3.5, cex = 0.75)
}

run_phase <- function() {
  qc <- make_qc()
  if (any(qc$status != "PASS")) {
    write_status("FAILED_INPUT_QC", "One or more required input files or columns are missing. See phase8_1R_figure_input_qc.tsv.")
    stop("Input QC failed")
  }

  fig2 <- read_tsv(figure2_input)
  fig3 <- read_tsv(figure3_input)

  output_files <- c(
    save_dual("Figure2A_internal_component_estimates", function() plot_figure2A(fig2)),
    save_dual("Figure2B_internal_component_contrasts", function() plot_figure2B(fig2)),
    save_dual("Figure3_external_SBP_triangulation", function() plot_figure3(fig3))
  )

  manifest <- data.frame(
    figure_id = sub("\\.(pdf|png)$", "", basename(output_files)),
    file_path = output_files,
    file_type = sub("^.*\\.", "", output_files),
    file_exists = file.exists(output_files),
    file_size_bytes = ifelse(file.exists(output_files), file.info(output_files)$size, NA),
    status = ifelse(file.exists(output_files) & file.info(output_files)$size > 0, "PASS", "FAIL"),
    stringsAsFactors = FALSE
  )

  write.table(manifest, manifest_file, sep = "\t", quote = FALSE, row.names = FALSE)

  if (all(manifest$status == "PASS")) {
    write_status("PASSED", "Base R figures generated successfully: Figure2A, Figure2B, and Figure3.")
  } else {
    write_status("FAILED_OUTPUT_QC", "One or more figure files were missing or empty. See phase8_1R_figure_file_manifest.tsv.")
  }
}

tryCatch(
  run_phase(),
  error = function(e) {
    write_status("FAILED", paste("Error:", conditionMessage(e)))
    stop(e)
  }
)
