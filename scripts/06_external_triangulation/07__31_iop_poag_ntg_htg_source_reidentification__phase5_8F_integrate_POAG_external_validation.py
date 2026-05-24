#!/usr/bin/env python3
import csv
import os
import time

OUTDIR = "../../34_poag_external_validation_integration"
os.makedirs(OUTDIR, exist_ok=True)

PHASE4_1_SBP_CONTRAST = "../../21_sbp_component_contrast/phase4_1_SBP_component_contrast_summary.tsv"
PHASE4_2_ROBUSTNESS = "../../22_sbp_robustness/phase4_2_sbp_robustness_summary.tsv"
PHASE4_8_VASCULAR = "../../24_vascular_panel_integration/phase4_8_vascular_panel_integrated_evidence_table.tsv"
PHASE5_7_NEURO = "../../30_external_triangulation_integration/phase5_7_external_neuroretinal_triangulation_summary.tsv"
PHASE5_8E_POAG = "../../33_poag_external_mr_results/phase5_8E_POAG_external_mr_results.tsv"
PHASE5_8E_HET = "../../33_poag_external_mr_results/phase5_8E_POAG_external_heterogeneity.tsv"
PHASE5_8E_EGGER = "../../33_poag_external_mr_results/phase5_8E_POAG_external_egger_intercept.tsv"

STATUS_OUT = os.path.join(OUTDIR, "phase5_8F_status.tsv")
POAG_SUMMARY_OUT = os.path.join(OUTDIR, "phase5_8F_POAG_external_validation_summary.tsv")
MATRIX_OUT = os.path.join(OUTDIR, "phase5_8F_cross_layer_evidence_matrix.tsv")
CLAIMS_OUT = os.path.join(OUTDIR, "phase5_8F_updated_publication_claims_after_POAG.tsv")
INTERPRETATION_MD = os.path.join(OUTDIR, "phase5_8F_POAG_external_validation_interpretation.md")
RUNTIME_OUT = os.path.join(OUTDIR, "phase5_8F_runtime_log.tsv")

def read_tsv(path):
    if not os.path.exists(path):
        return []
    with open(path, "r", encoding="utf-8", errors="replace", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))

def first(rows, **criteria):
    for r in rows:
        ok = True
        for k, v in criteria.items():
            if str(r.get(k, "")) != str(v):
                ok = False
                break
        if ok:
            return r
    return {}

def val(row, key, default="NA"):
    if not row:
        return default
    x = row.get(key, default)
    if x is None or x == "":
        return default
    return x

def fnum(x):
    try:
        if x in ("", "NA", None):
            return None
        return float(x)
    except Exception:
        return None

def evidence_label(p, q, n=None):
    pnum = fnum(p)
    qnum = fnum(q)
    nnum = fnum(n) if n is not None else None

    if qnum is not None and qnum < 0.05:
        return "FDR_SIGNIFICANT"
    if pnum is not None and pnum < 0.05:
        if nnum is not None and nnum <= 3:
            return "NOMINAL_ONLY_LOW_POWER"
        return "NOMINAL_ONLY"
    return "NOT_SIGNIFICANT"

def sig_text(p, q):
    label = evidence_label(p, q)
    if label == "FDR_SIGNIFICANT":
        return "FDR significant"
    if label == "NOMINAL_ONLY":
        return "nominal only"
    if label == "NOMINAL_ONLY_LOW_POWER":
        return "nominal only; low power"
    return "not significant"

def write_tsv(path, fieldnames, rows):
    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

def main():
    start = time.time()

    sbp_contrast = read_tsv(PHASE4_1_SBP_CONTRAST)
    robustness = read_tsv(PHASE4_2_ROBUSTNESS)
    vascular = read_tsv(PHASE4_8_VASCULAR)
    neuro = read_tsv(PHASE5_7_NEURO)
    poag = read_tsv(PHASE5_8E_POAG)
    poag_het = read_tsv(PHASE5_8E_HET)
    poag_egger = read_tsv(PHASE5_8E_EGGER)

    poag_summary_rows = []

    for exposure_id in ["SBP", "ART_STIFFNESS"]:
        r = first(
            poag,
            exposure_id=exposure_id,
            outcome_id="POAG",
            method="IVW_multiplicative_random_effects",
        )
        h = first(poag_het, exposure_id=exposure_id, outcome_id="POAG")
        e = first(poag_egger, exposure_id=exposure_id, outcome_id="POAG")

        n = val(r, "n_instruments")
        p = val(r, "pval")
        q = val(r, "qval_bh_ivw_random")
        label = evidence_label(p, q, n)

        if exposure_id == "SBP":
            if label == "NOT_SIGNIFICANT":
                interpretation = "Neutral POAG external validation; does not significantly support the SBP component-divergence signal"
            else:
                interpretation = "POAG external validation provides statistical support for SBP"
        else:
            if label == "NOT_SIGNIFICANT":
                interpretation = "No clear POAG evidence for arterial stiffness; low power because only 3 instruments were available"
            else:
                interpretation = "Exploratory POAG signal for arterial stiffness; interpret cautiously due to only 3 instruments"

        poag_summary_rows.append({
            "exposure_id": exposure_id,
            "outcome_id": "POAG",
            "outcome_scale": val(r, "outcome_scale", "log_odds"),
            "method": "IVW_multiplicative_random_effects",
            "n_instruments": n,
            "beta_log_odds": val(r, "beta_log_odds"),
            "se": val(r, "se"),
            "pval": p,
            "qval_bh": q,
            "or": val(r, "or"),
            "or_ci_lower": val(r, "or_ci_lower"),
            "or_ci_upper": val(r, "or_ci_upper"),
            "direction": val(r, "direction"),
            "evidence_label": label,
            "heterogeneity_flag": val(h, "heterogeneity_flag"),
            "Q_pval": val(h, "Q_pval"),
            "egger_intercept_pval": val(e, "egger_intercept_pval"),
            "egger_pleiotropy_flag": val(e, "pleiotropy_flag"),
            "integrated_interpretation": interpretation,
        })

    write_tsv(
        POAG_SUMMARY_OUT,
        [
            "exposure_id",
            "outcome_id",
            "outcome_scale",
            "method",
            "n_instruments",
            "beta_log_odds",
            "se",
            "pval",
            "qval_bh",
            "or",
            "or_ci_lower",
            "or_ci_upper",
            "direction",
            "evidence_label",
            "heterogeneity_flag",
            "Q_pval",
            "egger_intercept_pval",
            "egger_pleiotropy_flag",
            "integrated_interpretation",
        ],
        poag_summary_rows,
    )

    sbp_contrast_r0 = first(sbp_contrast, exposure_id="SBP", assumed_component_correlation_r="0")
    sbp_vascular = first(vascular, exposure_id="SBP")
    art_vascular = first(vascular, exposure_id="ART_STIFFNESS")

    sbp_rnfl = first(neuro, exposure_id="SBP", outcome_id="RNFL")
    sbp_gcipl = first(neuro, exposure_id="SBP", outcome_id="GCIPL")
    art_rnfl = first(neuro, exposure_id="ART_STIFFNESS", outcome_id="RNFL")
    art_gcipl = first(neuro, exposure_id="ART_STIFFNESS", outcome_id="GCIPL")

    sbp_poag = first(poag_summary_rows, exposure_id="SBP")
    art_poag = first(poag_summary_rows, exposure_id="ART_STIFFNESS")

    matrix_rows = [
        {
            "evidence_layer": "Primary component MR screen",
            "main_result": "No individual component association survived FDR correction; SBP was the strongest nominal exposure",
            "direction_or_effect_pattern": "SBP negative for nonIOP component and positive for IOP component",
            "statistical_strength": "nominal individual associations only",
            "interpretation": "Primary screen is hypothesis-generating rather than confirmatory",
            "manuscript_role": "Background discovery screen",
        },
        {
            "evidence_layer": "Formal nonIOP vs IOP component contrast",
            "main_result": "SBP beta_IOPcomponent minus beta_nonIOPcomponent was significant under r=0 and across correlation sensitivity assumptions",
            "direction_or_effect_pattern": "IOP_MORE_POSITIVE_THAN_NONIOP",
            "statistical_strength": "p=" + val(sbp_contrast_r0, "p_contrast") + "; q=" + val(sbp_contrast_r0, "q_contrast_within_r"),
            "interpretation": "Strongest project-level signal",
            "manuscript_role": "Main finding",
        },
        {
            "evidence_layer": "SBP robustness package",
            "main_result": "SBP direction was broadly consistent across IVW, Egger slope, weighted median direction checks, outlier-corrected IVW, and method-specific contrasts",
            "direction_or_effect_pattern": "nonIOP mostly negative; IOP mostly positive",
            "statistical_strength": "directionally supportive, not uniformly significant",
            "interpretation": "Supports robustness of direction but does not establish confirmatory causality",
            "manuscript_role": "Sensitivity analysis",
        },
        {
            "evidence_layer": "Vascular exposure panel",
            "main_result": "ART_STIFFNESS showed the same nonIOP-negative and IOP-positive direction pattern as SBP",
            "direction_or_effect_pattern": "directionally consistent with SBP",
            "statistical_strength": val(art_vascular, "component_divergence_label", "exploratory; not significant"),
            "interpretation": "Exploratory vascular-stiffness support only because only 3 instruments were available",
            "manuscript_role": "Exploratory vascular extension",
        },
        {
            "evidence_layer": "External neuroretinal endophenotypes",
            "main_result": "RNFL and GCIPL external MR did not provide FDR-significant support for SBP",
            "direction_or_effect_pattern": "SBP positive for RNFL and GCIPL but not significant",
            "statistical_strength": "neutral; no FDR-significant evidence",
            "interpretation": "SBP component-divergence signal is not clearly mirrored by RNFL/GCIPL thickness",
            "manuscript_role": "External triangulation, neutral",
        },
        {
            "evidence_layer": "External POAG clinical endpoint",
            "main_result": "SBP -> POAG was negative but not significant under IVW random effects",
            "direction_or_effect_pattern": "SBP beta log-odds=" + val(sbp_poag, "beta_log_odds") + "; OR=" + val(sbp_poag, "or"),
            "statistical_strength": "p=" + val(sbp_poag, "pval") + "; q=" + val(sbp_poag, "qval_bh"),
            "interpretation": "Neutral POAG validation; broad POAG may mix IOP-dependent and IOP-independent biology",
            "manuscript_role": "External clinical validation, neutral",
        },
        {
            "evidence_layer": "Unresolved external pressure/subtype outcomes",
            "main_result": "IOP, NTG, and HTG outcome validation could not be completed with available downloadable sources",
            "direction_or_effect_pattern": "NA",
            "statistical_strength": "not analyzed",
            "interpretation": "This remains the most important limitation for high-impact mechanistic interpretation",
            "manuscript_role": "Limitation and future work",
        },
    ]

    write_tsv(
        MATRIX_OUT,
        [
            "evidence_layer",
            "main_result",
            "direction_or_effect_pattern",
            "statistical_strength",
            "interpretation",
            "manuscript_role",
        ],
        matrix_rows,
    )

    claim_rows = [
        {
            "claim_id": "C1",
            "claim": "SBP is the strongest project-level signal.",
            "supporting_evidence": "Formal SBP nonIOP vs IOP contrast was FDR-significant; robustness package supported direction consistency.",
            "counterbalancing_evidence": "Individual component associations did not survive FDR; RNFL/GCIPL and POAG external validation were neutral.",
            "recommended_wording": "SBP showed the strongest hypothesis-generating evidence for directional divergence between IOP-dependent and IOP-independent glaucoma genetic components.",
            "claim_strength": "HYPOTHESIS_GENERATING_STRONGEST_INTERNAL_SIGNAL",
        },
        {
            "claim_id": "C2",
            "claim": "The result should not be described as externally confirmed.",
            "supporting_evidence": "RNFL, GCIPL, and POAG external MR did not produce FDR-significant SBP associations.",
            "counterbalancing_evidence": "External POAG is a mixed clinical endpoint and may dilute opposite component-specific effects.",
            "recommended_wording": "External neuroretinal and POAG triangulation did not clearly confirm the SBP component-divergence signal.",
            "claim_strength": "CAUTION_REQUIRED",
        },
        {
            "claim_id": "C3",
            "claim": "POAG neutrality is interpretable rather than fatal.",
            "supporting_evidence": "The main finding is component divergence; overall POAG can combine pressure-dependent and pressure-independent biology.",
            "counterbalancing_evidence": "A significant POAG association would have strengthened clinical validation, but it was not observed.",
            "recommended_wording": "The neutral POAG-wide result may reflect etiologic heterogeneity within POAG rather than directly refuting the component-level divergence.",
            "claim_strength": "BALANCED_INTERPRETATION",
        },
        {
            "claim_id": "C4",
            "claim": "ART_STIFFNESS should remain exploratory.",
            "supporting_evidence": "ART_STIFFNESS showed directional consistency with SBP in components and a positive but non-significant POAG association.",
            "counterbalancing_evidence": "Only 3 instruments were available; estimates were imprecise.",
            "recommended_wording": "Arterial stiffness provided exploratory directional context but insufficient statistical evidence.",
            "claim_strength": "EXPLORATORY_ONLY",
        },
        {
            "claim_id": "C5",
            "claim": "IOP/NTG/HTG validation remains the key missing layer.",
            "supporting_evidence": "IOP, NTG, and HTG sources were not yet analyzable from public downloadable files.",
            "counterbalancing_evidence": "POAG and RNFL/GCIPL triangulation were completed but neutral.",
            "recommended_wording": "Further validation in IOP, normal-tension glaucoma, and high-tension glaucoma datasets is needed.",
            "claim_strength": "MAJOR_LIMITATION",
        },
    ]

    write_tsv(
        CLAIMS_OUT,
        [
            "claim_id",
            "claim",
            "supporting_evidence",
            "counterbalancing_evidence",
            "recommended_wording",
            "claim_strength",
        ],
        claim_rows,
    )

    with open(INTERPRETATION_MD, "w", encoding="utf-8", newline="") as f:
        f.write("# Phase 5.8F POAG External Validation Integration\n\n")
        f.write("## Current status\n\n")
        f.write("PASSED_WITH_NEUTRAL_POAG_EXTERNAL_SUPPORT\n\n")
        f.write("## Main interpretation\n\n")
        f.write("The strongest evidence remains the formal SBP component-divergence contrast between GBS_nonIOPcomponent and GBS_IOPcomponent.\n\n")
        f.write("SBP external validation using POAG did not provide statistically significant support. The IVW random-effects estimate for SBP -> POAG was negative but non-significant, with substantial heterogeneity across instruments.\n\n")
        f.write("This neutral POAG-wide result does not directly refute the component-divergence finding, because POAG is a broad clinical endpoint that may mix IOP-dependent and IOP-independent biological pathways.\n\n")
        f.write("External RNFL and GCIPL triangulation also did not provide FDR-significant support. Therefore, the current project should be framed as a hypothesis-generating vascular/hemodynamic component-divergence study rather than a fully externally confirmed causal mechanism.\n\n")
        f.write("ART_STIFFNESS remains exploratory. It showed component-level directional consistency with SBP, but only three instruments were available and POAG external MR was not significant.\n\n")
        f.write("## Recommended manuscript wording\n\n")
        f.write("Among evaluated systemic and vascular traits, SBP showed the strongest evidence for directional divergence between IOP-dependent and IOP-independent glaucoma genetic components. This divergence was supported by formal component-contrast testing and directionally consistent sensitivity analyses. However, external triangulation using RNFL, GCIPL, and POAG did not provide FDR-significant support, indicating that the SBP signal remains hypothesis-generating and requires further validation in IOP, NTG, and HTG datasets.\n\n")
        f.write("## Updated evidence hierarchy\n\n")
        f.write("1. SBP component contrast: strongest internal evidence\n")
        f.write("2. SBP robustness package: directionally supportive\n")
        f.write("3. ART_STIFFNESS: exploratory vascular-stiffness context\n")
        f.write("4. RNFL / GCIPL: neutral external neuroretinal triangulation\n")
        f.write("5. POAG: neutral external clinical endpoint validation\n")
        f.write("6. IOP / NTG / HTG: unresolved and important limitation\n\n")
        f.write("## Manuscript implication\n\n")
        f.write("The manuscript can now move toward final evidence synthesis and outline drafting, but the Results and Discussion should avoid confirmatory causal language.\n")

    elapsed = time.time() - start

    neutral_poag = 0
    for r in poag_summary_rows:
        if r["evidence_label"] == "NOT_SIGNIFICANT":
            neutral_poag += 1

    with open(STATUS_OUT, "w", encoding="utf-8", newline="") as f:
        f.write("phase\tstatus\tkey_result\tnote\n")
        f.write(f"Phase 5.8F integrate POAG external validation evidence\tPASSED_WITH_NEUTRAL_POAG_EXTERNAL_SUPPORT\tpoag_tests={len(poag_summary_rows)};poag_not_significant={neutral_poag}\tPOAG evidence integrated with component, robustness, vascular, and neuroretinal layers\n")
        f.write("main_internal_signal\tDOCUMENTED\tSBP component-divergence contrast remains strongest evidence\tInternal contrast is stronger than external validation layers\n")
        f.write("POAG_external_support\tDOCUMENTED\tNeutral; no FDR-significant POAG evidence\tPOAG-wide endpoint may mix IOP-dependent and IOP-independent biology\n")
        f.write("manuscript_framing\tUPDATED\thypothesis-generating vascular/hemodynamic component-divergence study\tAvoid externally confirmed causal language\n")
        f.write(f"runtime\tINFO\t{elapsed:.3f}s\tPhase 5.8F runtime\n")

    with open(RUNTIME_OUT, "w", encoding="utf-8", newline="") as f:
        f.write("phase\telapsed_seconds\telapsed_human\n")
        f.write(f"Phase 5.8F\t{elapsed:.3f}\t{elapsed:.1f}s\n")

    print("===== Phase 5.8F completed =====")
    print("Wrote:", STATUS_OUT)
    print("Wrote:", POAG_SUMMARY_OUT)
    print("Wrote:", MATRIX_OUT)
    print("Wrote:", CLAIMS_OUT)
    print("Wrote:", INTERPRETATION_MD)

if __name__ == "__main__":
    main()
