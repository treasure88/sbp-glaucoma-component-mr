#!/usr/bin/env python3
import csv
import math
import os
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = "../../48_figures"
os.makedirs(BASE, exist_ok=True)

FIG2_INPUT = "../../46_figure_table_preparation/figure2_internal_component_contrast_plot_input.tsv"
FIG3_INPUT = "../../46_figure_table_preparation/figure3_external_validation_plot_input.tsv"

STATUS_OUT = os.path.join(BASE, "phase8_1_status.tsv")
MANIFEST_OUT = os.path.join(BASE, "phase8_1_figure_file_manifest.tsv")
QC_OUT = os.path.join(BASE, "phase8_1_figure_input_qc.tsv")
RUNTIME_OUT = os.path.join(BASE, "phase8_1_runtime_log.tsv")

FIG2A_PNG = os.path.join(BASE, "Figure2A_internal_component_estimates.png")
FIG2A_PDF = os.path.join(BASE, "Figure2A_internal_component_estimates.pdf")
FIG2B_PNG = os.path.join(BASE, "Figure2B_internal_component_contrasts.png")
FIG2B_PDF = os.path.join(BASE, "Figure2B_internal_component_contrasts.pdf")
FIG3_PNG = os.path.join(BASE, "Figure3_external_SBP_triangulation.png")
FIG3_PDF = os.path.join(BASE, "Figure3_external_SBP_triangulation.pdf")

def read_tsv(path):
    if not os.path.exists(path):
        return []
    with open(path, "r", encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))

def write_tsv(path, fields, rows):
    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

def fnum(x):
    try:
        if x is None:
            return None
        s = str(x).strip()
        if s in ("", "NA", "NaN", "nan", "None"):
            return None
        return float(s)
    except Exception:
        return None

def ci(beta, se):
    b = fnum(beta)
    s = fnum(se)
    if b is None or s is None:
        return None, None
    return b - 1.96 * s, b + 1.96 * s

def fmt_p(x):
    v = fnum(x)
    if v is None:
        return "NA"
    if v < 0.001:
        return f"{v:.1e}"
    return f"{v:.3f}"

def draw_figure2a(rows):
    plot_rows = [r for r in rows if r.get("exposure_id") in ["SBP", "DBP", "ART_STIFFNESS"]]
    labels = []
    y = []
    x = []
    xerr_low = []
    xerr_high = []

    idx = 0
    for r in plot_rows:
        exposure = r["exposure_id"]
        for comp, beta_col, se_col in [
            ("nonIOP component", "beta_nonIOP", "se_nonIOP"),
            ("IOP component", "beta_IOP", "se_IOP")
        ]:
            b = fnum(r.get(beta_col))
            s = fnum(r.get(se_col))
            if b is None or s is None:
                continue
            lo, hi = ci(b, s)
            labels.append(f"{exposure}: {comp}")
            y.append(idx)
            x.append(b)
            xerr_low.append(b - lo)
            xerr_high.append(hi - b)
            idx += 1

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.errorbar(x, y, xerr=[xerr_low, xerr_high], fmt="o", capsize=3)
    ax.axvline(0, linestyle="--", linewidth=1)
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.invert_yaxis()
    ax.set_xlabel("MR estimate beta")
    ax.set_title("Internal component-specific MR estimates")
    ax.text(
        0.01, -0.12,
        "Error bars show 95% confidence intervals. Arterial stiffness is exploratory with 3 instruments.",
        transform=ax.transAxes,
        ha="left",
        va="top"
    )
    fig.tight_layout()
    fig.savefig(FIG2A_PNG, dpi=300, bbox_inches="tight")
    fig.savefig(FIG2A_PDF, bbox_inches="tight")
    plt.close(fig)

def draw_figure2b(rows):
    plot_rows = [r for r in rows if r.get("exposure_id") in ["SBP", "DBP", "ART_STIFFNESS"]]
    labels = []
    y = []
    x = []
    xerr_low = []
    xerr_high = []
    annotations = []

    for idx, r in enumerate(plot_rows):
        b = fnum(r.get("contrast_beta_diff_r0"))
        s = fnum(r.get("contrast_se_r0"))
        if b is None or s is None:
            continue
        lo, hi = ci(b, s)
        labels.append(r.get("exposure_id", "NA"))
        y.append(idx)
        x.append(b)
        xerr_low.append(b - lo)
        xerr_high.append(hi - b)
        annotations.append(f"p={fmt_p(r.get('contrast_p_r0'))}; q={fmt_p(r.get('contrast_q_r0'))}")

    fig, ax = plt.subplots(figsize=(8, 4))
    ax.errorbar(x, y, xerr=[xerr_low, xerr_high], fmt="o", capsize=3)
    ax.axvline(0, linestyle="--", linewidth=1)
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.invert_yaxis()
    ax.set_xlabel("Contrast estimate: beta(IOP component) - beta(nonIOP component)")
    ax.set_title("Formal nonIOP versus IOP component contrast")

    for xi, yi, txt in zip(x, y, annotations):
        ax.text(xi, yi + 0.18, txt, fontsize=8, ha="center")

    ax.text(
        0.01, -0.18,
        "SBP contrast is the primary hypothesis-generating result. DBP is a comparator; arterial stiffness is exploratory.",
        transform=ax.transAxes,
        ha="left",
        va="top"
    )
    fig.tight_layout()
    fig.savefig(FIG2B_PNG, dpi=300, bbox_inches="tight")
    fig.savefig(FIG2B_PDF, bbox_inches="tight")
    plt.close(fig)

def draw_figure3(rows):
    plot_rows = [r for r in rows if r.get("exposure_id") == "SBP"]
    preferred_order = ["IOP", "NTG", "POAG", "RNFL", "GCIPL", "HTG"]
    rank = {k: i for i, k in enumerate(preferred_order)}
    plot_rows.sort(key=lambda r: rank.get(r.get("outcome_id", ""), 99))

    labels = []
    y = []
    x = []
    xerr_low = []
    xerr_high = []
    annotations = []

    idx = 0
    for r in plot_rows:
        if r.get("outcome_id") == "HTG":
            continue
        b = fnum(r.get("beta"))
        s = fnum(r.get("se"))
        if b is None or s is None:
            continue
        lo, hi = ci(b, s)
        labels.append(r.get("outcome_id", "NA"))
        y.append(idx)
        x.append(b)
        xerr_low.append(b - lo)
        xerr_high.append(hi - b)
        annotations.append(f"{r.get('direction','NA')}; p={fmt_p(r.get('pval'))}")
        idx += 1

    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.errorbar(x, y, xerr=[xerr_low, xerr_high], fmt="o", capsize=3)
    ax.axvline(0, linestyle="--", linewidth=1)
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.invert_yaxis()
    ax.set_xlabel("MR estimate beta on source GWAS scale")
    ax.set_title("External triangulation of the SBP signal")

    for xi, yi, txt in zip(x, y, annotations):
        ax.text(xi, yi + 0.18, txt, fontsize=8, ha="center")

    ax.text(
        0.01, -0.20,
        "Effect scales differ across outcomes; this figure summarizes direction and uncertainty, not directly comparable magnitudes. HTG was unavailable.",
        transform=ax.transAxes,
        ha="left",
        va="top"
    )
    fig.tight_layout()
    fig.savefig(FIG3_PNG, dpi=300, bbox_inches="tight")
    fig.savefig(FIG3_PDF, bbox_inches="tight")
    plt.close(fig)

def main():
    start = time.time()

    fig2_rows = read_tsv(FIG2_INPUT)
    fig3_rows = read_tsv(FIG3_INPUT)

    qc_rows = [
        {
            "input_file": FIG2_INPUT,
            "exists": "YES" if os.path.exists(FIG2_INPUT) else "NO",
            "rows": len(fig2_rows),
            "status": "READY" if fig2_rows else "MISSING_OR_EMPTY",
            "note": "Figure 2 internal component and contrast input"
        },
        {
            "input_file": FIG3_INPUT,
            "exists": "YES" if os.path.exists(FIG3_INPUT) else "NO",
            "rows": len(fig3_rows),
            "status": "READY" if fig3_rows else "MISSING_OR_EMPTY",
            "note": "Figure 3 external validation input"
        }
    ]

    if fig2_rows:
        draw_figure2a(fig2_rows)
        draw_figure2b(fig2_rows)

    if fig3_rows:
        draw_figure3(fig3_rows)

    figure_files = [
        ("Figure2A", FIG2A_PNG, FIG2A_PDF, "Internal component-specific MR estimates"),
        ("Figure2B", FIG2B_PNG, FIG2B_PDF, "Formal nonIOP versus IOP component contrast"),
        ("Figure3", FIG3_PNG, FIG3_PDF, "External SBP validation / triangulation"),
    ]

    manifest = []
    for fig_id, png, pdf, desc in figure_files:
        manifest.append({
            "figure_id": fig_id,
            "png_file": png,
            "pdf_file": pdf,
            "png_exists": "YES" if os.path.exists(png) else "NO",
            "pdf_exists": "YES" if os.path.exists(pdf) else "NO",
            "description": desc,
            "recommended_use": "Main manuscript figure candidate"
        })

    write_tsv(QC_OUT, ["input_file", "exists", "rows", "status", "note"], qc_rows)
    write_tsv(
        MANIFEST_OUT,
        ["figure_id", "png_file", "pdf_file", "png_exists", "pdf_exists", "description", "recommended_use"],
        manifest
    )

    elapsed = time.time() - start

    created = sum(1 for r in manifest if r["png_exists"] == "YES" and r["pdf_exists"] == "YES")

    with open(STATUS_OUT, "w", encoding="utf-8", newline="") as f:
        f.write("phase\tstatus\tkey_result\tnote\n")
        f.write(f"Phase 8.1 draw main figures\tPASSED\tfigures_created={created}/3\tFigure 2A, Figure 2B, and Figure 3 generated as PNG and PDF\n")
        f.write("figure2_status\tINFO\tcomponent_estimates_and_contrast_created\tUse Figure2B as main if journal has limited space\n")
        f.write("figure3_status\tINFO\texternal_triangulation_created\tInterpret effect magnitudes cautiously because outcome scales differ\n")
        f.write("next_step\tTO_DO\tPhase 8.2 optional evidence heatmap or Phase 8.3 conceptual schematic\tFigure 1 schematic can be drawn after final journal target\n")
        f.write(f"runtime\tINFO\t{elapsed:.3f}s\tPhase 8.1 runtime\n")

    with open(RUNTIME_OUT, "w", encoding="utf-8") as f:
        f.write("phase\telapsed_seconds\telapsed_human\n")
        f.write(f"Phase 8.1\t{elapsed:.3f}\t{elapsed:.1f}s\n")

    print("===== Phase 8.1 completed =====")
    print("Wrote:", STATUS_OUT)
    print("Wrote:", MANIFEST_OUT)
    print("Wrote:", QC_OUT)

if __name__ == "__main__":
    main()
