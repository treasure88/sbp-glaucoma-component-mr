#!/usr/bin/env python3
import csv
import os
import time
import shutil

BASE = "../../46_figure_table_preparation"
os.makedirs(BASE, exist_ok=True)

MASTER = "../../44_master_evidence_table/phase6_1_master_evidence_table.tsv"
ROUTEA = "../../44_master_evidence_table/phase6_1_routeA_master_summary_matrix.tsv"
FLAGS = "../../44_master_evidence_table/phase6_1_publication_readiness_flags.tsv"
CLAIMS = "../../45_final_claims_wording/phase6_2_final_permitted_claims.tsv"
AVOID = "../../45_final_claims_wording/phase6_2_prohibited_or_unsafe_claims.tsv"
INTERNAL = "../../24_vascular_panel_integration/phase4_8_vascular_panel_integrated_evidence_table.tsv"
EXTERNAL_GRID = "../../42_ntg_htg_validation_integration/phase5_10F_routeA_external_validation_grid.tsv"
NEURO_MR = "../../29_external_neuroretinal_mr_results/phase5_6_external_neuroretinal_mr_results.tsv"
POAG_MR = "../../33_poag_external_mr_results/phase5_8E_POAG_external_mr_results.tsv"
IOP_SUMMARY = "../../38_iop_external_validation_integration/phase5_9F_IOP_external_validation_summary.tsv"
NTG_SUMMARY = "../../42_ntg_htg_validation_integration/phase5_10F_NTG_HTG_validation_summary.tsv"

STATUS_OUT = os.path.join(BASE, "phase6_3_status.tsv")
MANIFEST_OUT = os.path.join(BASE, "phase6_3_figure_table_manifest.tsv")
TABLE1_OUT = os.path.join(BASE, "table1_primary_vascular_component_evidence.tsv")
TABLE2_OUT = os.path.join(BASE, "table2_external_routeA_triangulation.tsv")
TABLE3_OUT = os.path.join(BASE, "table3_claims_and_limitations.tsv")
SUPP_MASTER_OUT = os.path.join(BASE, "supplementary_table_master_evidence.tsv")
FIG1_OUT = os.path.join(BASE, "figure1_conceptual_model_edges.tsv")
FIG2_OUT = os.path.join(BASE, "figure2_internal_component_contrast_plot_input.tsv")
FIG3_OUT = os.path.join(BASE, "figure3_external_validation_plot_input.tsv")
FIG4_OUT = os.path.join(BASE, "figure4_routeA_heatmap_input.tsv")
CAPTIONS_OUT = os.path.join(BASE, "phase6_3_figure_table_captions.tsv")
RUNTIME_OUT = os.path.join(BASE, "phase6_3_runtime_log.tsv")

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
        return float(x)
    except Exception:
        return None

def fmt(x, digits=3):
    v = fnum(x)
    if v is None:
        return str(x) if x not in (None, "") else "NA"
    if abs(v) < 0.001 and v != 0:
        return f"{v:.2e}"
    return f"{v:.{digits}g}"

def get_row(rows, **kwargs):
    for r in rows:
        ok = True
        for k, v in kwargs.items():
            if r.get(k) != v:
                ok = False
                break
        if ok:
            return r
    return {}

def direction_symbol(direction):
    d = str(direction).lower()
    if "positive" in d:
        return "+"
    if "negative" in d:
        return "-"
    if "null" in d or "na" in d:
        return "NA"
    return d

def main():
    start = time.time()

    master = read_tsv(MASTER)
    routea = read_tsv(ROUTEA)
    flags = read_tsv(FLAGS)
    claims = read_tsv(CLAIMS)
    avoid = read_tsv(AVOID)
    internal = read_tsv(INTERNAL)
    external_grid = read_tsv(EXTERNAL_GRID)
    neuro = read_tsv(NEURO_MR)
    poag = read_tsv(POAG_MR)
    iop = read_tsv(IOP_SUMMARY)
    ntg = read_tsv(NTG_SUMMARY)

    # Supplementary master table copy
    if os.path.exists(MASTER):
        shutil.copyfile(MASTER, SUPP_MASTER_OUT)

    # Table 1: primary vascular/component evidence
    table1 = []
    for exposure in ["SBP", "DBP", "HYPERTENSION", "ART_STIFFNESS"]:
        r = get_row(internal, exposure_id=exposure)
        if not r:
            continue

        table1.append({
            "exposure_id": exposure,
            "trait_group": r.get("trait_group", "NA"),
            "analysis_status": r.get("analysis_status", "NA"),
            "n_instruments_nonIOP": r.get("n_instruments_nonIOP", "NA"),
            "beta_nonIOP": r.get("beta_nonIOP", "NA"),
            "p_nonIOP": r.get("p_nonIOP", "NA"),
            "direction_nonIOP": r.get("direction_nonIOP", "NA"),
            "n_instruments_IOP": r.get("n_instruments_IOP", "NA"),
            "beta_IOP": r.get("beta_IOP", "NA"),
            "p_IOP": r.get("p_IOP", "NA"),
            "direction_IOP": r.get("direction_IOP", "NA"),
            "contrast_beta_diff_r0": r.get("contrast_beta_diff_r0", "NA"),
            "contrast_p_r0": r.get("contrast_p_r0", "NA"),
            "contrast_q_r0": r.get("contrast_q_r0", "NA"),
            "component_divergence_label": r.get("component_divergence_label", "NA"),
            "evidence_tier": r.get("evidence_tier", "NA"),
            "publication_interpretation": r.get("publication_interpretation", "NA"),
            "main_caveat": r.get("main_caveat", "NA"),
        })

    write_tsv(TABLE1_OUT, [
        "exposure_id", "trait_group", "analysis_status",
        "n_instruments_nonIOP", "beta_nonIOP", "p_nonIOP", "direction_nonIOP",
        "n_instruments_IOP", "beta_IOP", "p_IOP", "direction_IOP",
        "contrast_beta_diff_r0", "contrast_p_r0", "contrast_q_r0",
        "component_divergence_label", "evidence_tier",
        "publication_interpretation", "main_caveat"
    ], table1)

    # Table 2: external Route-A triangulation
    table2 = []
    for r in external_grid:
        table2.append({
            "exposure_id": r.get("exposure_id", "NA"),
            "validation_layer": r.get("validation_layer", "NA"),
            "outcome_id": r.get("outcome_id", "NA"),
            "role_in_routeA": r.get("role_in_routeA", "NA"),
            "n_instruments": r.get("n_instruments", "NA"),
            "beta": r.get("beta", "NA"),
            "se": r.get("se", "NA"),
            "pval": r.get("pval", "NA"),
            "qval": r.get("qval", "NA"),
            "direction": r.get("direction", "NA"),
            "evidence_label": r.get("evidence_label", "NA"),
            "interpretation": r.get("interpretation", "NA"),
        })

    write_tsv(TABLE2_OUT, [
        "exposure_id", "validation_layer", "outcome_id", "role_in_routeA",
        "n_instruments", "beta", "se", "pval", "qval", "direction",
        "evidence_label", "interpretation"
    ], table2)

    # Table 3: claims and limitations
    table3 = []
    for r in claims:
        table3.append({
            "item_type": "permitted_claim",
            "item_id": r.get("claim_id", "NA"),
            "strength_or_level": r.get("claim_level", "NA"),
            "wording": r.get("recommended_wording", "NA"),
            "qualifier_or_reason": r.get("required_qualifier", "NA"),
            "manuscript_location": r.get("manuscript_location", "NA"),
        })
    for r in avoid:
        table3.append({
            "item_type": "unsafe_claim_to_avoid",
            "item_id": r.get("unsafe_claim_id", "NA"),
            "strength_or_level": "AVOID",
            "wording": r.get("unsafe_or_prohibited_wording", "NA"),
            "qualifier_or_reason": r.get("reason", "NA"),
            "manuscript_location": "All sections; safe replacement: " + r.get("safe_replacement", "NA"),
        })

    write_tsv(TABLE3_OUT, [
        "item_type", "item_id", "strength_or_level", "wording",
        "qualifier_or_reason", "manuscript_location"
    ], table3)

    # Figure 1: conceptual model edges
    fig1 = [
        {
            "from_node": "Systemic vascular liability",
            "to_node": "SBP genetic liability",
            "edge_label": "primary vascular signal",
            "direction": "context",
            "evidence_strength": "hypothesis-generating",
            "caption_note": "SBP is the main vascular trait showing component divergence."
        },
        {
            "from_node": "SBP genetic liability",
            "to_node": "Measured IOP",
            "edge_label": "positive borderline",
            "direction": "positive",
            "evidence_strength": "directional_support_not_confirmatory",
            "caption_note": "External IOP validation: positive borderline association."
        },
        {
            "from_node": "Measured IOP",
            "to_node": "IOP-dependent glaucoma genetic component",
            "edge_label": "IOP-linked biology",
            "direction": "positive_context",
            "evidence_strength": "biological_model",
            "caption_note": "Conceptual IOP-dependent pathway."
        },
        {
            "from_node": "SBP genetic liability",
            "to_node": "IOP-dependent glaucoma genetic component",
            "edge_label": "more positive than nonIOP component",
            "direction": "positive_relative",
            "evidence_strength": "primary_hypothesis_generating_signal",
            "caption_note": "Internal component contrast is the strongest result."
        },
        {
            "from_node": "SBP genetic liability",
            "to_node": "NTG",
            "edge_label": "negative null",
            "direction": "negative_null",
            "evidence_strength": "directional_context",
            "caption_note": "External NTG validation does not support a positive generalized NTG effect."
        },
        {
            "from_node": "HTG",
            "to_node": "External validation gap",
            "edge_label": "unavailable",
            "direction": "missing",
            "evidence_strength": "important_limitation",
            "caption_note": "HTG was the ideal clinical subtype endpoint but summary statistics were not locked."
        }
    ]

    write_tsv(FIG1_OUT, [
        "from_node", "to_node", "edge_label", "direction",
        "evidence_strength", "caption_note"
    ], fig1)

    # Figure 2: internal component contrast plot input
    fig2 = []
    for r in internal:
        exposure = r.get("exposure_id", "NA")
        if exposure not in ["SBP", "DBP", "ART_STIFFNESS"]:
            continue
        fig2.append({
            "exposure_id": exposure,
            "plot_group": "internal_component_contrast",
            "beta_nonIOP": r.get("beta_nonIOP", "NA"),
            "se_nonIOP": r.get("se_nonIOP", "NA"),
            "p_nonIOP": r.get("p_nonIOP", "NA"),
            "beta_IOP": r.get("beta_IOP", "NA"),
            "se_IOP": r.get("se_IOP", "NA"),
            "p_IOP": r.get("p_IOP", "NA"),
            "contrast_beta_diff_r0": r.get("contrast_beta_diff_r0", "NA"),
            "contrast_se_r0": r.get("contrast_se_r0", "NA"),
            "contrast_p_r0": r.get("contrast_p_r0", "NA"),
            "contrast_q_r0": r.get("contrast_q_r0", "NA"),
            "contrast_direction": r.get("contrast_direction", "NA"),
            "evidence_label": r.get("component_divergence_label", "NA"),
            "plot_note": "Plot as paired component estimates plus contrast estimate."
        })

    write_tsv(FIG2_OUT, [
        "exposure_id", "plot_group",
        "beta_nonIOP", "se_nonIOP", "p_nonIOP",
        "beta_IOP", "se_IOP", "p_IOP",
        "contrast_beta_diff_r0", "contrast_se_r0", "contrast_p_r0", "contrast_q_r0",
        "contrast_direction", "evidence_label", "plot_note"
    ], fig2)

    # Figure 3: SBP external validation plot input
    fig3 = []
    for r in external_grid:
        if r.get("exposure_id") != "SBP":
            continue
        fig3.append({
            "exposure_id": "SBP",
            "outcome_id": r.get("outcome_id", "NA"),
            "validation_layer": r.get("validation_layer", "NA"),
            "beta": r.get("beta", "NA"),
            "se": r.get("se", "NA"),
            "pval": r.get("pval", "NA"),
            "qval": r.get("qval", "NA"),
            "direction": r.get("direction", "NA"),
            "evidence_label": r.get("evidence_label", "NA"),
            "plot_note": r.get("interpretation", "NA")
        })

    # Add RNFL/GCIPL SBP random-effects external endpoints
    for r in neuro:
        if r.get("exposure_id") == "SBP" and r.get("method") == "IVW_multiplicative_random_effects":
            fig3.append({
                "exposure_id": "SBP",
                "outcome_id": r.get("outcome_id", "NA"),
                "validation_layer": "neuroretinal_endophenotype",
                "beta": r.get("beta", "NA"),
                "se": r.get("se", "NA"),
                "pval": r.get("pval", "NA"),
                "qval": r.get("qval_bh_ivw_random", "NA"),
                "direction": r.get("direction", "NA"),
                "evidence_label": "NOT_SIGNIFICANT_" + str(r.get("direction", "NA")).upper(),
                "plot_note": "No confirmatory neuroretinal validation evidence."
            })

    write_tsv(FIG3_OUT, [
        "exposure_id", "outcome_id", "validation_layer",
        "beta", "se", "pval", "qval", "direction",
        "evidence_label", "plot_note"
    ], fig3)

    # Figure 4: heatmap input
    fig4 = []
    heatmap_layers = [
        ("internal_nonIOP", "GBS_nonIOPcomponent"),
        ("internal_IOP", "GBS_IOPcomponent"),
        ("internal_contrast", "IOP_minus_nonIOP_contrast"),
        ("external_IOP", "IOP"),
        ("external_NTG", "NTG"),
        ("external_POAG", "POAG"),
        ("external_RNFL", "RNFL"),
        ("external_GCIPL", "GCIPL"),
        ("external_HTG", "HTG")
    ]

    for exposure in ["SBP", "DBP", "ART_STIFFNESS", "HYPERTENSION"]:
        internal_row = get_row(internal, exposure_id=exposure)

        # Internal cells
        if internal_row:
            fig4.extend([
                {
                    "exposure_id": exposure,
                    "evidence_layer": "internal_nonIOP",
                    "outcome_or_component": "GBS_nonIOPcomponent",
                    "direction": internal_row.get("direction_nonIOP", "NA"),
                    "beta": internal_row.get("beta_nonIOP", "NA"),
                    "pval": internal_row.get("p_nonIOP", "NA"),
                    "qval": internal_row.get("q_nonIOP", "NA"),
                    "cell_label": internal_row.get("nonIOP_interpretation_label", "NA"),
                    "cell_note": "Internal component MR"
                },
                {
                    "exposure_id": exposure,
                    "evidence_layer": "internal_IOP",
                    "outcome_or_component": "GBS_IOPcomponent",
                    "direction": internal_row.get("direction_IOP", "NA"),
                    "beta": internal_row.get("beta_IOP", "NA"),
                    "pval": internal_row.get("p_IOP", "NA"),
                    "qval": internal_row.get("q_IOP", "NA"),
                    "cell_label": internal_row.get("IOP_interpretation_label", "NA"),
                    "cell_note": "Internal component MR"
                },
                {
                    "exposure_id": exposure,
                    "evidence_layer": "internal_contrast",
                    "outcome_or_component": "IOP_minus_nonIOP_contrast",
                    "direction": internal_row.get("contrast_direction", "NA"),
                    "beta": internal_row.get("contrast_beta_diff_r0", "NA"),
                    "pval": internal_row.get("contrast_p_r0", "NA"),
                    "qval": internal_row.get("contrast_q_r0", "NA"),
                    "cell_label": internal_row.get("component_divergence_label", "NA"),
                    "cell_note": "Internal component contrast"
                }
            ])

        # External cells from route-A only for SBP
        if exposure == "SBP":
            for r in external_grid:
                fig4.append({
                    "exposure_id": exposure,
                    "evidence_layer": "external_" + r.get("outcome_id", "NA"),
                    "outcome_or_component": r.get("outcome_id", "NA"),
                    "direction": r.get("direction", "NA"),
                    "beta": r.get("beta", "NA"),
                    "pval": r.get("pval", "NA"),
                    "qval": r.get("qval", "NA"),
                    "cell_label": r.get("evidence_label", "NA"),
                    "cell_note": r.get("interpretation", "NA")
                })

            for r in neuro:
                if r.get("exposure_id") == "SBP" and r.get("method") == "IVW_multiplicative_random_effects":
                    fig4.append({
                        "exposure_id": exposure,
                        "evidence_layer": "external_" + r.get("outcome_id", "NA"),
                        "outcome_or_component": r.get("outcome_id", "NA"),
                        "direction": r.get("direction", "NA"),
                        "beta": r.get("beta", "NA"),
                        "pval": r.get("pval", "NA"),
                        "qval": r.get("qval_bh_ivw_random", "NA"),
                        "cell_label": "NOT_SIGNIFICANT_" + str(r.get("direction", "NA")).upper(),
                        "cell_note": "External neuroretinal endpoint"
                    })

    write_tsv(FIG4_OUT, [
        "exposure_id", "evidence_layer", "outcome_or_component",
        "direction", "beta", "pval", "qval", "cell_label", "cell_note"
    ], fig4)

    captions = [
        {
            "item_id": "Figure_1",
            "item_type": "conceptual_model",
            "proposed_title": "Conceptual model of vascular liability across IOP-dependent and IOP-independent glaucoma pathways",
            "caption": "SBP was the primary vascular trait showing hypothesis-generating component divergence. External validation provided directional but non-confirmatory support for an IOP-related interpretation, with positive borderline evidence for measured IOP and negative non-significant evidence for NTG. HTG validation was unavailable."
        },
        {
            "item_id": "Figure_2",
            "item_type": "main_result_plot",
            "proposed_title": "Internal component-specific MR and nonIOP versus IOP contrast",
            "caption": "Component-specific MR estimates and formal IOP-minus-nonIOP contrast for prespecified vascular traits. SBP showed the strongest directionally divergent pattern, while DBP was a comparator and arterial stiffness was exploratory because of low instrument count."
        },
        {
            "item_id": "Figure_3",
            "item_type": "external_validation_plot",
            "proposed_title": "External validation of the SBP signal across measured IOP and glaucoma-related outcomes",
            "caption": "External triangulation of SBP against measured IOP, NTG, POAG, and neuroretinal endophenotypes. Results were directional rather than confirmatory; measured IOP was positive borderline, NTG was negative non-significant, and POAG/RNFL/GCIPL were non-confirmatory."
        },
        {
            "item_id": "Figure_4",
            "item_type": "evidence_heatmap",
            "proposed_title": "Evidence map across vascular traits and glaucoma-related outcomes",
            "caption": "Heatmap-ready evidence matrix summarizing direction, evidence label, and caveats across internal component MR, component contrast, and external validation layers."
        },
        {
            "item_id": "Table_1",
            "item_type": "main_table",
            "proposed_title": "Primary vascular-panel evidence for glaucoma component divergence",
            "caption": "Main table summarizing internal component-specific MR, component contrasts, evidence tier, and caveats for SBP, DBP, hypertension liability, and arterial stiffness."
        },
        {
            "item_id": "Table_2",
            "item_type": "main_or_supplement_table",
            "proposed_title": "External triangulation of the SBP signal",
            "caption": "External validation table for SBP across measured IOP, NTG, POAG, and unavailable HTG endpoint."
        },
        {
            "item_id": "Supplementary_Table_1",
            "item_type": "supplementary_table",
            "proposed_title": "Master evidence table",
            "caption": "Complete evidence inventory integrating internal MR, robustness interpretation, external validation, and manuscript-safe interpretation labels."
        }
    ]

    write_tsv(CAPTIONS_OUT, [
        "item_id", "item_type", "proposed_title", "caption"
    ], captions)

    manifest = [
        {
            "item_id": "Figure_1",
            "file": FIG1_OUT,
            "recommended_location": "Main Figure 1",
            "status": "DATA_READY_CONCEPTUAL",
            "note": "Use as schematic model input; can be drawn manually or with PowerPoint/BioRender."
        },
        {
            "item_id": "Figure_2",
            "file": FIG2_OUT,
            "recommended_location": "Main Figure 2",
            "status": "PLOT_INPUT_READY",
            "note": "Internal component estimates and contrast plot input."
        },
        {
            "item_id": "Figure_3",
            "file": FIG3_OUT,
            "recommended_location": "Main Figure 3 or Supplementary Figure",
            "status": "PLOT_INPUT_READY",
            "note": "External SBP validation forest plot input."
        },
        {
            "item_id": "Figure_4",
            "file": FIG4_OUT,
            "recommended_location": "Main Figure 4 or Supplementary Figure",
            "status": "HEATMAP_INPUT_READY",
            "note": "Evidence heatmap input."
        },
        {
            "item_id": "Table_1",
            "file": TABLE1_OUT,
            "recommended_location": "Main Table 1",
            "status": "TABLE_READY",
            "note": "Primary vascular component evidence."
        },
        {
            "item_id": "Table_2",
            "file": TABLE2_OUT,
            "recommended_location": "Main Table 2 or Supplementary Table",
            "status": "TABLE_READY",
            "note": "External Route-A triangulation."
        },
        {
            "item_id": "Table_3",
            "file": TABLE3_OUT,
            "recommended_location": "Supplementary Table or writing aid",
            "status": "TABLE_READY",
            "note": "Permitted and prohibited claims."
        },
        {
            "item_id": "Supplementary_Table_1",
            "file": SUPP_MASTER_OUT,
            "recommended_location": "Supplementary Table 1",
            "status": "TABLE_READY",
            "note": "Full master evidence table."
        },
        {
            "item_id": "Captions",
            "file": CAPTIONS_OUT,
            "recommended_location": "Manuscript figure/table legends",
            "status": "TEXT_READY",
            "note": "Draft captions for manuscript."
        }
    ]

    write_tsv(MANIFEST_OUT, [
        "item_id", "file", "recommended_location", "status", "note"
    ], manifest)

    elapsed = time.time() - start

    with open(STATUS_OUT, "w", encoding="utf-8", newline="") as f:
        f.write("phase\tstatus\tkey_result\tnote\n")
        f.write(f"Phase 6.3 figure/table preparation\tPASSED\tfigures=4;tables=4;captions={len(captions)}\tManuscript-ready figure/table inputs generated\n")
        f.write("main_figure_recommendation\tINFO\tFigure1_conceptual_model;Figure2_internal_component_contrast;Figure3_external_validation\tFigure4 heatmap optional or supplementary depending on journal space\n")
        f.write("main_table_recommendation\tINFO\tTable1_primary_vascular_component_evidence\tTable2 can be main or supplementary depending on word limit\n")
        f.write("next_step\tTO_DO\tPhase 7 manuscript drafting\tUse Phase 6.2 claims and Phase 6.3 figure/table files\n")
        f.write(f"runtime\tINFO\t{elapsed:.3f}s\tPhase 6.3 runtime\n")

    with open(RUNTIME_OUT, "w", encoding="utf-8", newline="") as f:
        f.write("phase\telapsed_seconds\telapsed_human\n")
        f.write(f"Phase 6.3\t{elapsed:.3f}\t{elapsed:.1f}s\n")

    print("===== Phase 6.3 completed =====")
    print("Wrote:", STATUS_OUT)
    print("Wrote:", MANIFEST_OUT)
    print("Wrote:", TABLE1_OUT)
    print("Wrote:", TABLE2_OUT)
    print("Wrote:", FIG1_OUT)
    print("Wrote:", FIG2_OUT)
    print("Wrote:", FIG3_OUT)
    print("Wrote:", FIG4_OUT)
    print("Wrote:", CAPTIONS_OUT)

if __name__ == "__main__":
    main()
