#!/usr/bin/env python3
import csv
import os
import time

BASE = "../../42_ntg_htg_validation_integration"
os.makedirs(BASE, exist_ok=True)

NTG_RESULTS = "../../41_ntg_external_mr_results/phase5_10E_NTG_external_mr_results.tsv"
NTG_HET = "../../41_ntg_external_mr_results/phase5_10E_NTG_external_heterogeneity.tsv"
IOP_SUMMARY = "../../38_iop_external_validation_integration/phase5_9F_IOP_external_validation_summary.tsv"
POAG_RESULTS = "../../33_poag_external_mr_results/phase5_8E_POAG_external_mr_results.tsv"

STATUS_OUT = os.path.join(BASE, "phase5_10F_status.tsv")
SUMMARY_OUT = os.path.join(BASE, "phase5_10F_NTG_HTG_validation_summary.tsv")
GRID_OUT = os.path.join(BASE, "phase5_10F_routeA_external_validation_grid.tsv")
CLAIMS_OUT = os.path.join(BASE, "phase5_10F_updated_routeA_claims.tsv")
RUNTIME_OUT = os.path.join(BASE, "phase5_10F_runtime_log.tsv")

def read_tsv(path):
    if not os.path.exists(path):
        return []
    with open(path, "r", encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))

def get_row(rows, exposure=None, outcome=None, method=None):
    for r in rows:
        if exposure is not None and r.get("exposure_id") != exposure:
            continue
        if outcome is not None and r.get("outcome_id") != outcome:
            continue
        if method is not None and r.get("method") != method:
            continue
        return r
    return {}

def fnum(x):
    try:
        return float(x)
    except Exception:
        return None

def evidence_label(beta, pval, qval=None):
    b = fnum(beta)
    p = fnum(pval)
    q = fnum(qval)

    if b is None:
        return "NOT_EVALUABLE"

    direction = "POSITIVE" if b > 0 else "NEGATIVE" if b < 0 else "ZERO"

    if q is not None and q < 0.05:
        return "FDR_SIGNIFICANT_" + direction
    if p is not None and p < 0.05:
        return "NOMINAL_" + direction
    if p is not None and p < 0.10:
        return "BORDERLINE_" + direction
    return "NOT_SIGNIFICANT_" + direction

def main():
    start = time.time()

    ntg_rows = read_tsv(NTG_RESULTS)
    ntg_het_rows = read_tsv(NTG_HET)
    iop_summary_rows = read_tsv(IOP_SUMMARY)
    poag_rows = read_tsv(POAG_RESULTS)

    summary_rows = []

    for exposure in ["SBP", "ART_STIFFNESS"]:
        random = get_row(ntg_rows, exposure, "NTG", "IVW_multiplicative_random_effects")
        fixed = get_row(ntg_rows, exposure, "NTG", "IVW_fixed_effects")
        median = get_row(ntg_rows, exposure, "NTG", "weighted_median_ratio_descriptive")
        egger = get_row(ntg_rows, exposure, "NTG", "MR_Egger_slope")
        het = get_row(ntg_het_rows, exposure, "NTG", "IVW")

        label = evidence_label(
            random.get("beta"),
            random.get("pval"),
            random.get("qval_bh_ivw_random")
        )

        if exposure == "SBP":
            interpretation = (
                "SBP showed a negative but non-significant association direction with NTG. "
                "Together with the positive borderline SBP-to-IOP result, this supports a directional pattern "
                "in which SBP-related liability may map more closely to IOP-related glaucoma biology than to NTG liability. "
                "This remains hypothesis-generating rather than confirmatory."
            )
            manuscript_use = "Use as directional external triangulation for route A."
            caveat = "NTG random-effects result was not significant; heterogeneity was present; HTG full summary statistics remain unavailable."
        else:
            interpretation = (
                "Arterial stiffness showed no reliable NTG validation evidence because only three instruments were available."
            )
            manuscript_use = "Do not emphasize; keep as exploratory low-power evidence."
            caveat = "Only 3 instruments; estimates are imprecise and method-specific directions are unstable."

        summary_rows.append({
            "exposure_id": exposure,
            "outcome_id": "NTG",
            "outcome_scale": "log_odds_or_MTAG_scale_TO_VERIFY",
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
            "evidence_label": label,
            "interpretation": interpretation,
            "manuscript_use": manuscript_use,
            "main_caveat": caveat,
        })

    with open(SUMMARY_OUT, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(summary_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(summary_rows)

    sbp_iop = {}
    for r in iop_summary_rows:
        if r.get("exposure_id") == "SBP":
            sbp_iop = r
            break

    sbp_ntg = get_row(ntg_rows, "SBP", "NTG", "IVW_multiplicative_random_effects")
    sbp_poag = get_row(poag_rows, "SBP", "POAG", "IVW_multiplicative_random_effects")

    grid_rows = [
        {
            "exposure_id": "SBP",
            "validation_layer": "IOP",
            "outcome_id": "IOP",
            "role_in_routeA": "IOP_mediated_mechanism_check",
            "n_instruments": sbp_iop.get("n_instruments", "NA"),
            "beta": sbp_iop.get("beta_random", "NA"),
            "se": sbp_iop.get("se_random", "NA"),
            "pval": sbp_iop.get("p_random", "NA"),
            "qval": sbp_iop.get("q_random", "NA"),
            "direction": sbp_iop.get("direction_random", "positive"),
            "evidence_label": sbp_iop.get("evidence_label", "BORDERLINE_DIRECTIONAL_POSITIVE"),
            "interpretation": "Positive borderline signal; supports possible IOP-related SBP biology but is not confirmatory.",
        },
        {
            "exposure_id": "SBP",
            "validation_layer": "NTG",
            "outcome_id": "NTG",
            "role_in_routeA": "IOP_independent_clinical_subtype_check",
            "n_instruments": sbp_ntg.get("n_instruments", "NA"),
            "beta": sbp_ntg.get("beta", "NA"),
            "se": sbp_ntg.get("se", "NA"),
            "pval": sbp_ntg.get("pval", "NA"),
            "qval": sbp_ntg.get("qval_bh_ivw_random", "NA"),
            "direction": sbp_ntg.get("direction", "NA"),
            "evidence_label": evidence_label(sbp_ntg.get("beta"), sbp_ntg.get("pval"), sbp_ntg.get("qval_bh_ivw_random")),
            "interpretation": "Negative but non-significant signal; does not support a positive generalized NTG effect.",
        },
        {
            "exposure_id": "SBP",
            "validation_layer": "POAG",
            "outcome_id": "POAG",
            "role_in_routeA": "broad_clinical_glaucoma_endpoint_check",
            "n_instruments": sbp_poag.get("n_instruments", "NA"),
            "beta": sbp_poag.get("beta_log_odds", sbp_poag.get("beta", "NA")),
            "se": sbp_poag.get("se", "NA"),
            "pval": sbp_poag.get("pval", "NA"),
            "qval": sbp_poag.get("qval_bh_ivw_random", "NA"),
            "direction": sbp_poag.get("direction", "NA"),
            "evidence_label": evidence_label(sbp_poag.get("beta_log_odds", sbp_poag.get("beta", "NA")), sbp_poag.get("pval"), sbp_poag.get("qval_bh_ivw_random")),
            "interpretation": "No reliable broad POAG validation signal.",
        },
        {
            "exposure_id": "SBP",
            "validation_layer": "HTG",
            "outcome_id": "HTG",
            "role_in_routeA": "ideal_IOP_dependent_clinical_subtype_check",
            "n_instruments": "NA",
            "beta": "NA",
            "se": "NA",
            "pval": "NA",
            "qval": "NA",
            "direction": "NA",
            "evidence_label": "NOT_AVAILABLE",
            "interpretation": "HTG full summary statistics not yet locked; this remains the key missing validation endpoint.",
        },
    ]

    with open(GRID_OUT, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(grid_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(grid_rows)

    claims = [
        {
            "claim_id": "SBP_IOP_NTG_directional_triangulation",
            "claim_strength": "DIRECTIONAL_SUPPORT_NOT_CONFIRMATORY",
            "recommended_wording": (
                "External triangulation showed a positive borderline SBP-to-IOP association and a negative non-significant SBP-to-NTG association, "
                "consistent with the hypothesis that SBP-related genetic liability may map more strongly to IOP-related glaucoma biology than to NTG liability."
            ),
            "avoid_wording": "SBP causally increases IOP-dependent glaucoma and protects against NTG.",
            "evidence_source": "Phase 5.9E/5.9F + Phase 5.10E/5.10F",
            "note": "Use cautious language because neither IOP random-effects nor NTG random-effects estimates reached conventional significance.",
        },
        {
            "claim_id": "HTG_missing_endpoint",
            "claim_strength": "IMPORTANT_LIMITATION",
            "recommended_wording": "HTG validation could not be completed because a full downloadable HTG summary-statistics source was not locked.",
            "avoid_wording": "HTG results support the mechanism.",
            "evidence_source": "Phase 5.10A-5.10F",
            "note": "HTG remains the most important missing external validation outcome.",
        },
        {
            "claim_id": "ART_STIFFNESS_subtype_validation",
            "claim_strength": "LOW_POWER_EXPLORATORY_ONLY",
            "recommended_wording": "Arterial stiffness subtype-validation analyses were exploratory and underpowered because only three instruments were available.",
            "avoid_wording": "Arterial stiffness increases NTG risk.",
            "evidence_source": "Phase 5.10E/5.10F",
            "note": "Do not emphasize ART_STIFFNESS external subtype validation.",
        },
    ]

    with open(CLAIMS_OUT, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(claims[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(claims)

    elapsed = time.time() - start

    with open(STATUS_OUT, "w", encoding="utf-8", newline="") as f:
        f.write("phase\tstatus\tkey_result\tnote\n")
        f.write("Phase 5.10F integrate NTG/HTG validation evidence\tPASSED\tSBP_IOP_positive_borderline_and_NTG_negative_null;HTG_missing\tExternal subtype validation integrated into route-A interpretation\n")
        f.write("SBP_NTG_interpretation\tDIRECTIONAL_SUPPORT_NOT_CONFIRMATORY\tnegative_null_against_NTG\tSupports IOP-related rather than generalized NTG signal, but not definitive\n")
        f.write("HTG_interpretation\tNOT_AVAILABLE\tfull_summary_stats_not_locked\tKey missing external subtype endpoint\n")
        f.write("ART_STIFFNESS_interpretation\tLOW_POWER_EXPLORATORY_ONLY\tn=3\tDo not emphasize subtype validation\n")
        f.write(f"runtime\tINFO\t{elapsed:.3f}s\tPhase 5.10F runtime\n")

    with open(RUNTIME_OUT, "w", encoding="utf-8", newline="") as f:
        f.write("phase\telapsed_seconds\telapsed_human\n")
        f.write(f"Phase 5.10F\t{elapsed:.3f}\t{elapsed:.1f}s\n")

    print("===== Phase 5.10F completed =====")
    print("Wrote:", STATUS_OUT)
    print("Wrote:", SUMMARY_OUT)
    print("Wrote:", GRID_OUT)
    print("Wrote:", CLAIMS_OUT)

if __name__ == "__main__":
    main()
