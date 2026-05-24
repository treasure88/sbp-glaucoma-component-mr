#!/usr/bin/env python3
import csv
import os
import time

BASE = "../../38_iop_external_validation_integration"
os.makedirs(BASE, exist_ok=True)

MR_RESULTS = "../../37_iop_external_mr_results/phase5_9E_IOP_external_mr_results.tsv"
HET = "../../37_iop_external_mr_results/phase5_9E_IOP_external_heterogeneity.tsv"

SUMMARY_OUT = os.path.join(BASE, "phase5_9F_IOP_external_validation_summary.tsv")
CLAIMS_OUT = os.path.join(BASE, "phase5_9F_updated_routeA_claims.tsv")
STATUS_OUT = os.path.join(BASE, "phase5_9F_status.tsv")
RUNTIME_OUT = os.path.join(BASE, "phase5_9F_runtime_log.tsv")

def read_tsv(path):
    with open(path, "r", encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))

def get_row(rows, exposure, method):
    for r in rows:
        if r.get("exposure_id") == exposure and r.get("method") == method:
            return r
    return {}

def fnum(x):
    try:
        return float(x)
    except Exception:
        return None

def label_primary(row):
    beta = fnum(row.get("beta"))
    p = fnum(row.get("pval"))
    q = fnum(row.get("qval_bh_ivw_random"))
    direction = row.get("direction", "NA")

    if beta is None:
        return "NOT_EVALUABLE"

    if q is not None and q < 0.05:
        return "FDR_SIGNIFICANT_" + direction.upper()

    if p is not None and p < 0.05:
        return "NOMINAL_" + direction.upper()

    if p is not None and p < 0.10:
        return "BORDERLINE_DIRECTIONAL_" + direction.upper()

    return "NOT_SIGNIFICANT_" + direction.upper()

def main():
    start = time.time()

    rows = read_tsv(MR_RESULTS)
    het_rows = read_tsv(HET) if os.path.exists(HET) else []

    summary = []

    for exposure in ["SBP", "ART_STIFFNESS"]:
        fixed = get_row(rows, exposure, "IVW_fixed_effects")
        random = get_row(rows, exposure, "IVW_multiplicative_random_effects")
        median = get_row(rows, exposure, "weighted_median_ratio_descriptive")
        egger = get_row(rows, exposure, "MR_Egger_slope")

        het = {}
        for h in het_rows:
            if h.get("exposure_id") == exposure:
                het = h

        primary_label = label_primary(random)

        if exposure == "SBP":
            interpretation = (
                "SBP showed a positive association direction with measured IOP. "
                "The fixed-effect IVW estimate was nominally significant, while the primary random-effects estimate was borderline and not FDR-significant. "
                "This provides directional support for an IOP-mediated component but should not be interpreted as confirmatory evidence."
            )
            manuscript_use = (
                "Use as supportive external triangulation for the IOP-dependent side of the SBP component-divergence story."
            )
            caveat = (
                "Substantial heterogeneity; random-effects p was borderline; q was not significant."
            )
        else:
            interpretation = (
                "Arterial stiffness showed no reliable association with measured IOP. "
                "The estimate was imprecise and based on only three instruments."
            )
            manuscript_use = (
                "Do not use as strong validation; mention as exploratory and underpowered if needed."
            )
            caveat = (
                "Only 3 instruments; heterogeneous and unstable method-specific directions."
            )

        summary.append({
            "exposure_id": exposure,
            "outcome_id": "IOP",
            "outcome_scale": "continuous_IOP",
            "primary_method": "IVW_multiplicative_random_effects",
            "n_instruments": random.get("n_instruments", "NA"),
            "beta_random": random.get("beta", "NA"),
            "se_random": random.get("se", "NA"),
            "p_random": random.get("pval", "NA"),
            "q_random": random.get("qval_bh_ivw_random", "NA"),
            "direction_random": random.get("direction", "NA"),
            "beta_fixed": fixed.get("beta", "NA"),
            "p_fixed": fixed.get("pval", "NA"),
            "weighted_median_direction": median.get("direction", "NA"),
            "egger_direction": egger.get("direction", "NA"),
            "egger_p": egger.get("pval", "NA"),
            "Q": random.get("Q", "NA"),
            "Q_pval": random.get("Q_pval", "NA"),
            "phi": random.get("phi", "NA"),
            "heterogeneity_flag": het.get("heterogeneity_flag", "YES" if fnum(random.get("Q_pval")) is not None and fnum(random.get("Q_pval")) < 0.05 else "NA"),
            "evidence_label": primary_label,
            "interpretation": interpretation,
            "manuscript_use": manuscript_use,
            "main_caveat": caveat,
        })

    fields = list(summary[0].keys())
    with open(SUMMARY_OUT, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(summary)

    claims = [
        {
            "claim_id": "IOP_external_validation",
            "claim_strength": "DIRECTIONAL_SUPPORT_NOT_CONFIRMATORY",
            "recommended_wording": "SBP showed directionally positive, borderline evidence for association with measured IOP, supporting an IOP-mediated component of the vascular signal, although the random-effects estimate did not reach nominal significance and heterogeneity was substantial.",
            "avoid_wording": "SBP significantly increases IOP.",
            "evidence_source": "Phase 5.9E/5.9F",
            "note": "Use primary random-effects estimate for inference; fixed-effect nominal result can be shown as sensitivity."
        },
        {
            "claim_id": "routeA_story_update",
            "claim_strength": "STRENGTHENED_BUT_STILL_HYPOTHESIS_GENERATING",
            "recommended_wording": "Together with the internal component contrast, the positive SBP-to-IOP direction supports a model in which systemic blood pressure liability may map more strongly to IOP-related glaucoma biology.",
            "avoid_wording": "The mechanism is proven to be mediated by IOP.",
            "evidence_source": "Phase 4.8 + Phase 5.9E/5.9F",
            "note": "MVMR or formal mediation would be needed for stronger mediation language."
        },
        {
            "claim_id": "ART_STIFFNESS_IOP",
            "claim_strength": "NO_RELIABLE_SUPPORT",
            "recommended_wording": "Arterial stiffness did not provide reliable IOP validation evidence because the estimate was imprecise and based on only three instruments.",
            "avoid_wording": "Arterial stiffness decreases IOP.",
            "evidence_source": "Phase 5.9E/5.9F",
            "note": "Treat as exploratory only."
        }
    ]

    with open(CLAIMS_OUT, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(claims[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(claims)

    elapsed = time.time() - start

    with open(STATUS_OUT, "w", encoding="utf-8", newline="") as f:
        f.write("phase\tstatus\tkey_result\tnote\n")
        f.write("Phase 5.9F integrate IOP external validation evidence\tPASSED\tSBP_directional_positive_borderline;ART_STIFFNESS_no_reliable_support\tIOP validation integrated into route-A interpretation\n")
        f.write("SBP_IOP_interpretation\tDIRECTIONAL_SUPPORT_NOT_CONFIRMATORY\trandom_p=0.0559;q=0.1118;direction=positive\tSupports IOP-related mechanism but not definitive\n")
        f.write("ART_STIFFNESS_IOP_interpretation\tNO_RELIABLE_SUPPORT\tn=3;imprecise\tDo not emphasize as validation\n")
        f.write(f"runtime\tINFO\t{elapsed:.3f}s\tPhase 5.9F runtime\n")

    with open(RUNTIME_OUT, "w", encoding="utf-8", newline="") as f:
        f.write("phase\telapsed_seconds\telapsed_human\n")
        f.write(f"Phase 5.9F\t{elapsed:.3f}\t{elapsed:.1f}s\n")

    print("===== Phase 5.9F completed =====")
    print("Wrote:", STATUS_OUT)
    print("Wrote:", SUMMARY_OUT)
    print("Wrote:", CLAIMS_OUT)

if __name__ == "__main__":
    main()
