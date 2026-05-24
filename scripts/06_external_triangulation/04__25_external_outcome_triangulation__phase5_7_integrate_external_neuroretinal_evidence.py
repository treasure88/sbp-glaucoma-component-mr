#!/usr/bin/env python3
import csv
import os
import time

OUTDIR = "../../30_external_triangulation_integration"
os.makedirs(OUTDIR, exist_ok=True)

PHASE4_8_VASCULAR = "../../24_vascular_panel_integration/phase4_8_vascular_panel_integrated_evidence_table.tsv"
PHASE5_6_RESULTS = "../../29_external_neuroretinal_mr_results/phase5_6_external_neuroretinal_mr_results.tsv"
PHASE5_6_HET = "../../29_external_neuroretinal_mr_results/phase5_6_external_neuroretinal_heterogeneity.tsv"
PHASE5_6_EGGER = "../../29_external_neuroretinal_mr_results/phase5_6_external_neuroretinal_egger_intercept.tsv"

INTEGRATED_OUT = os.path.join(OUTDIR, "phase5_7_external_neuroretinal_triangulation_summary.tsv")
CLAIM_OUT = os.path.join(OUTDIR, "phase5_7_updated_publication_claims.tsv")
STATUS_OUT = os.path.join(OUTDIR, "phase5_7_status.tsv")
RUNTIME_OUT = os.path.join(OUTDIR, "phase5_7_runtime_log.tsv")
INTERPRETATION_MD = os.path.join(OUTDIR, "phase5_7_external_triangulation_interpretation.md")

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
    x = row.get(key, default) if row else default
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

def evidence_label(p, q, n):
    pnum = fnum(p)
    qnum = fnum(q)
    nnum = fnum(n)

    if qnum is not None and qnum < 0.05:
        return "FDR_SIGNIFICANT"
    if pnum is not None and pnum < 0.05:
        if nnum is not None and nnum <= 3:
            return "NOMINAL_ONLY_LOW_POWER"
        return "NOMINAL_ONLY"
    return "NOT_SIGNIFICANT"

def write_tsv(path, fieldnames, rows):
    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

def main():
    start = time.time()

    vascular = read_tsv(PHASE4_8_VASCULAR)
    mr = read_tsv(PHASE5_6_RESULTS)
    het = read_tsv(PHASE5_6_HET)
    egger = read_tsv(PHASE5_6_EGGER)

    rows = []

    for exposure_id in ["SBP", "ART_STIFFNESS"]:
        vascular_row = first(vascular, exposure_id=exposure_id)

        for outcome_id in ["RNFL", "GCIPL"]:
            mr_row = first(
                mr,
                exposure_id=exposure_id,
                outcome_id=outcome_id,
                method="IVW_multiplicative_random_effects",
            )
            het_row = first(het, exposure_id=exposure_id, outcome_id=outcome_id)
            egger_row = first(egger, exposure_id=exposure_id, outcome_id=outcome_id)

            n = val(mr_row, "n_instruments")
            beta = val(mr_row, "beta")
            se = val(mr_row, "se")
            p = val(mr_row, "pval")
            q = val(mr_row, "qval_bh_ivw_random")
            direction = val(mr_row, "direction")
            label = evidence_label(p, q, n)

            if exposure_id == "SBP":
                if label == "NOT_SIGNIFICANT":
                    interpretation = "No clear external neuroretinal support for the SBP component-divergence signal in this outcome"
                else:
                    interpretation = "External neuroretinal signal observed for SBP"
            else:
                if label == "NOMINAL_ONLY_LOW_POWER":
                    interpretation = "Exploratory nominal signal only; low power due to 3 instruments"
                elif label == "NOT_SIGNIFICANT":
                    interpretation = "No clear external neuroretinal support; low power due to 3 instruments"
                else:
                    interpretation = "External signal observed but should be interpreted cautiously due to limited instruments"

            rows.append({
                "exposure_id": exposure_id,
                "outcome_id": outcome_id,
                "external_outcome_layer": "neuroretinal_endophenotype",
                "component_story_role": val(vascular_row, "role_in_vascular_story"),
                "component_evidence_tier": val(vascular_row, "evidence_tier"),
                "component_contrast_direction": val(vascular_row, "contrast_direction"),
                "component_contrast_p_r0": val(vascular_row, "contrast_p_r0"),
                "component_contrast_q_r0": val(vascular_row, "contrast_q_r0"),
                "external_n_instruments": n,
                "external_beta": beta,
                "external_se": se,
                "external_pval": p,
                "external_qval_bh": q,
                "external_direction": direction,
                "external_evidence_label": label,
                "external_heterogeneity_flag": val(het_row, "heterogeneity_flag"),
                "external_Q_pval": val(het_row, "Q_pval"),
                "external_egger_intercept_pval": val(egger_row, "egger_intercept_pval"),
                "external_egger_pleiotropy_flag": val(egger_row, "pleiotropy_flag"),
                "integrated_interpretation": interpretation,
            })

    write_tsv(
        INTEGRATED_OUT,
        [
            "exposure_id",
            "outcome_id",
            "external_outcome_layer",
            "component_story_role",
            "component_evidence_tier",
            "component_contrast_direction",
            "component_contrast_p_r0",
            "component_contrast_q_r0",
            "external_n_instruments",
            "external_beta",
            "external_se",
            "external_pval",
            "external_qval_bh",
            "external_direction",
            "external_evidence_label",
            "external_heterogeneity_flag",
            "external_Q_pval",
            "external_egger_intercept_pval",
            "external_egger_pleiotropy_flag",
            "integrated_interpretation",
        ],
        rows,
    )

    claim_rows = [
        {
            "claim_id": "C1",
            "claim": "SBP remains the strongest project-level signal.",
            "supporting_result": "SBP has an FDR-significant formal component contrast between GBS_IOPcomponent and GBS_nonIOPcomponent.",
            "external_triangulation_update": "RNFL and GCIPL external MR did not show FDR-significant SBP associations.",
            "recommended_wording": "SBP showed the strongest evidence for component divergence, but this was not clearly supported by RNFL or GCIPL endophenotype triangulation.",
            "claim_strength": "HYPOTHESIS_GENERATING",
        },
        {
            "claim_id": "C2",
            "claim": "The SBP component-divergence signal should not be interpreted as proven neuroretinal structural mediation.",
            "supporting_result": "SBP associations with RNFL and GCIPL were positive but not statistically significant.",
            "external_triangulation_update": "No FDR-significant RNFL/GCIPL evidence.",
            "recommended_wording": "The available external neuroretinal endophenotype data did not indicate that SBP clearly acts through RNFL or GCIPL thickness.",
            "claim_strength": "CAUTION_REQUIRED",
        },
        {
            "claim_id": "C3",
            "claim": "ART_STIFFNESS provides exploratory but weak external support.",
            "supporting_result": "ART_STIFFNESS -> GCIPL was nominally positive, but q=0.104 and only 3 instruments were used.",
            "external_triangulation_update": "Nominal GCIPL result only; low power.",
            "recommended_wording": "Arterial stiffness showed a nominal positive association with GCIPL thickness, but this result was low-powered and did not survive FDR correction.",
            "claim_strength": "EXPLORATORY_ONLY",
        },
        {
            "claim_id": "C4",
            "claim": "The manuscript should remain focused on component-divergence rather than external endophenotype confirmation.",
            "supporting_result": "External RNFL/GCIPL triangulation did not add FDR-significant support.",
            "external_triangulation_update": "External validation layer is neutral rather than confirmatory.",
            "recommended_wording": "The main finding is directional divergence across glaucoma genetic components, with limited support from currently available neuroretinal endophenotype triangulation.",
            "claim_strength": "BALANCED_INTERPRETATION",
        },
    ]

    write_tsv(
        CLAIM_OUT,
        [
            "claim_id",
            "claim",
            "supporting_result",
            "external_triangulation_update",
            "recommended_wording",
            "claim_strength",
        ],
        claim_rows,
    )

    with open(INTERPRETATION_MD, "w", encoding="utf-8", newline="") as f:
        f.write("# Phase 5.7 External Neuroretinal Triangulation Interpretation\n\n")
        f.write("## Phase\n\n")
        f.write("Integrate external neuroretinal triangulation evidence\n\n")
        f.write("## Current status\n\n")
        f.write("PASSED_WITH_NEUTRAL_EXTERNAL_NEURORETINAL_TRIANGULATION\n\n")
        f.write("## Main interpretation\n\n")
        f.write("The strongest evidence in the project remains the SBP component-divergence signal between GBS_nonIOPcomponent and GBS_IOPcomponent.\n\n")
        f.write("External neuroretinal endophenotype MR using RNFL and GCIPL did not provide FDR-significant support for SBP effects on retinal nerve fibre layer or ganglion cell-inner plexiform layer thickness.\n\n")
        f.write("Therefore, the current data do not support a strong claim that the SBP component-divergence signal is mediated through measurable RNFL or GCIPL thickness.\n\n")
        f.write("ART_STIFFNESS showed a nominal positive association with GCIPL, but this result did not survive FDR correction and was based on only three instruments. It should be reported only as exploratory.\n\n")
        f.write("## Recommended manuscript wording\n\n")
        f.write("No individual external neuroretinal endophenotype association survived FDR correction. SBP showed positive but non-significant associations with RNFL and GCIPL thickness, suggesting that the SBP component-divergence signal was not clearly mirrored by these structural retinal endophenotypes. Arterial stiffness showed a nominal positive association with GCIPL, but this exploratory finding was based on only three instruments and did not survive FDR correction.\n\n")
        f.write("## Updated claim strength\n\n")
        f.write("- SBP component divergence: strongest hypothesis-generating signal\n")
        f.write("- RNFL / GCIPL triangulation: neutral, not confirmatory\n")
        f.write("- ART_STIFFNESS -> GCIPL: nominal exploratory only\n")
        f.write("- Overall article framing: vascular/hemodynamic component-divergence study with cautious external triangulation\n\n")
        f.write("## Next recommended direction\n\n")
        f.write("The most important remaining external validation layer is IOP / POAG / NTG / HTG. RNFL and GCIPL do not currently confirm the SBP signal, so pressure-related and glaucoma-subtype external outcomes remain higher priority for mechanistic interpretation.\n")

    elapsed = time.time() - start

    nominal = sum(1 for r in rows if r["external_evidence_label"].startswith("NOMINAL"))
    fdr = sum(1 for r in rows if r["external_evidence_label"] == "FDR_SIGNIFICANT")

    with open(STATUS_OUT, "w", encoding="utf-8", newline="") as f:
        f.write("phase\tstatus\tkey_result\tnote\n")
        f.write(f"Phase 5.7 external neuroretinal triangulation integration\tPASSED_WITH_NEUTRAL_EXTERNAL_SUPPORT\texternal_tests=4;nominal={nominal};fdr={fdr}\tIntegrated RNFL/GCIPL triangulation with component-divergence evidence\n")
        f.write("SBP_external_triangulation\tDOCUMENTED\tNo FDR-significant RNFL/GCIPL evidence\tDoes not clearly validate neuroretinal structural mediation\n")
        f.write("ART_STIFFNESS_external_triangulation\tDOCUMENTED\tGCIPL nominal only, q>0.05, n=3\tExploratory only\n")
        f.write("manuscript_interpretation\tUPDATED\tComponent-divergence signal remains main story\tExternal neuroretinal layer is neutral rather than confirmatory\n")
        f.write(f"runtime\tINFO\t{elapsed:.3f}s\tPhase 5.7 runtime\n")

    with open(RUNTIME_OUT, "w", encoding="utf-8", newline="") as f:
        f.write("phase\telapsed_seconds\telapsed_human\n")
        f.write(f"Phase 5.7\t{elapsed:.3f}\t{elapsed:.1f}s\n")

    print("===== Phase 5.7 completed =====")
    print("Wrote:", STATUS_OUT)
    print("Wrote:", INTEGRATED_OUT)
    print("Wrote:", CLAIM_OUT)
    print("Wrote:", INTERPRETATION_MD)

if __name__ == "__main__":
    main()
