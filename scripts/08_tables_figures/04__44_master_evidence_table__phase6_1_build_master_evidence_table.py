#!/usr/bin/env python3
import csv
import os
import time

BASE = "../../44_master_evidence_table"
os.makedirs(BASE, exist_ok=True)

INPUTS = {
    "phase4_8_integrated": "../../24_vascular_panel_integration/phase4_8_vascular_panel_integrated_evidence_table.tsv",
    "phase4_2B_curated_package_interpretation": "../../22_sbp_robustness/phase4_2B_formal_packages_final/phase4_2B_final_curated_formal_package_interpretation.tsv",
    "phase4_2B_curated_mrpresso": "../../22_sbp_robustness/phase4_2B_formal_packages_final/phase4_2B_final_curated_mrpresso_results.tsv",
    "phase4_2B_curated_mrraps": "../../22_sbp_robustness/phase4_2B_formal_packages_final/phase4_2B_final_curated_mrraps_interpretation.tsv",
    "phase5_6_neuroretinal_mr": "../../29_external_neuroretinal_mr_results/phase5_6_external_neuroretinal_mr_results.tsv",
    "phase5_8E_poag_mr": "../../33_poag_external_mr_results/phase5_8E_POAG_external_mr_results.tsv",
    "phase5_9F_iop_summary": "../../38_iop_external_validation_integration/phase5_9F_IOP_external_validation_summary.tsv",
    "phase5_10F_ntg_summary": "../../42_ntg_htg_validation_integration/phase5_10F_NTG_HTG_validation_summary.tsv",
    "phase5_10F_routeA_grid": "../../42_ntg_htg_validation_integration/phase5_10F_routeA_external_validation_grid.tsv",
    "phase5_10F_claims": "../../42_ntg_htg_validation_integration/phase5_10F_updated_routeA_claims.tsv",
    "phase5_11B_htg_lock": "../../43_htg_source_final_attempt/phase5_11B_corrected_HTG_source_lock.tsv",
}

STATUS_OUT = os.path.join(BASE, "phase6_1_status.tsv")
MASTER_OUT = os.path.join(BASE, "phase6_1_master_evidence_table.tsv")
ROUTEA_OUT = os.path.join(BASE, "phase6_1_routeA_master_summary_matrix.tsv")
FLAGS_OUT = os.path.join(BASE, "phase6_1_publication_readiness_flags.tsv")
SOURCE_OUT = os.path.join(BASE, "phase6_1_source_file_inventory.tsv")
RUNTIME_OUT = os.path.join(BASE, "phase6_1_runtime_log.tsv")

MASTER_FIELDS = [
    "evidence_id",
    "evidence_layer",
    "exposure_id",
    "outcome_or_component",
    "role_in_story",
    "analysis_status",
    "n_instruments",
    "beta",
    "se",
    "pval",
    "qval",
    "direction",
    "heterogeneity_or_caveat_flag",
    "evidence_label",
    "publication_weight",
    "manuscript_safe_interpretation",
    "main_caveat",
    "source_file"
]

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

def pub_weight_from_label(label, exposure_id="", n_instruments=""):
    n = fnum(n_instruments)

    if "NOT_AVAILABLE" in label:
        return "IMPORTANT_LIMITATION"
    if exposure_id == "ART_STIFFNESS" or (n is not None and n <= 3):
        return "LOW_POWER_EXPLORATORY"
    if "FDR_SIGNIFICANT" in label:
        return "PRIMARY_OR_STRONG_SUPPORT"
    if "BORDERLINE" in label:
        return "DIRECTIONAL_SUPPORT_NOT_CONFIRMATORY"
    if "NOMINAL" in label:
        return "NOMINAL_SUPPORT"
    if "NOT_SIGNIFICANT" in label:
        return "NULL_OR_DIRECTIONAL_CONTEXT"
    return "CONTEXTUAL"

def add_master(rows, **kwargs):
    row = {k: kwargs.get(k, "NA") for k in MASTER_FIELDS}
    rows.append(row)

def get_method_row(rows, exposure, outcome, method):
    for r in rows:
        if r.get("exposure_id") == exposure and r.get("outcome_id") == outcome and r.get("method") == method:
            return r
    return {}

def main():
    start = time.time()

    source_inventory = []
    for name, path in INPUTS.items():
        exists = os.path.exists(path)
        source_inventory.append({
            "source_name": name,
            "path": path,
            "exists": "YES" if exists else "NO",
            "n_rows": str(len(read_tsv(path))) if exists else "NA",
            "note": "Input source for Phase 6.1 master evidence table"
        })

    master = []

    # 1. Internal vascular-panel integrated evidence
    internal_rows = read_tsv(INPUTS["phase4_8_integrated"])

    for r in internal_rows:
        exposure = r.get("exposure_id", "NA")
        analysis_status = r.get("analysis_status", "NA")
        evidence_tier = r.get("evidence_tier", r.get("component_divergence_label", "NA"))

        n_inst = "nonIOP={};IOP={}".format(
            r.get("n_instruments_nonIOP", "NA"),
            r.get("n_instruments_IOP", "NA")
        )

        beta = "nonIOP={};IOP={};diff_r0={}".format(
            r.get("beta_nonIOP", "NA"),
            r.get("beta_IOP", "NA"),
            r.get("contrast_beta_diff_r0", "NA")
        )

        se = "nonIOP={};IOP={};diff_r0={}".format(
            r.get("se_nonIOP", "NA"),
            r.get("se_IOP", "NA"),
            r.get("contrast_se_r0", "NA")
        )

        pval = "nonIOP={};IOP={};contrast_r0={}".format(
            r.get("p_nonIOP", "NA"),
            r.get("p_IOP", "NA"),
            r.get("contrast_p_r0", "NA")
        )

        qval = "nonIOP={};IOP={};contrast_r0={}".format(
            r.get("q_nonIOP", "NA"),
            r.get("q_IOP", "NA"),
            r.get("contrast_q_r0", "NA")
        )

        caveat = r.get("main_caveat", "NA")
        het = "nonIOP_het={};IOP_het={}".format(
            r.get("nonIOP_het_flag", "NA"),
            r.get("IOP_het_flag", "NA")
        )

        if exposure == "SBP":
            publication_weight = "PRIMARY_HYPOTHESIS_GENERATING_SIGNAL"
        elif exposure == "DBP":
            publication_weight = "BLOOD_PRESSURE_COMPARATOR"
        elif exposure == "HYPERTENSION":
            publication_weight = "NOT_ANALYZABLE"
        elif exposure == "ART_STIFFNESS":
            publication_weight = "EXPLORATORY_DIRECTIONAL_SUPPORT_LOW_POWER"
        else:
            publication_weight = "CONTEXTUAL"

        add_master(
            master,
            evidence_id="internal_component_{}".format(exposure),
            evidence_layer="internal_component_MR",
            exposure_id=exposure,
            outcome_or_component="GBS_nonIOPcomponent_vs_GBS_IOPcomponent",
            role_in_story=r.get("role_in_vascular_story", "NA"),
            analysis_status=analysis_status,
            n_instruments=n_inst,
            beta=beta,
            se=se,
            pval=pval,
            qval=qval,
            direction=r.get("contrast_direction", "NA"),
            heterogeneity_or_caveat_flag=het,
            evidence_label=evidence_tier,
            publication_weight=publication_weight,
            manuscript_safe_interpretation=r.get("publication_interpretation", "NA"),
            main_caveat=caveat,
            source_file=INPUTS["phase4_8_integrated"]
        )

    # 2. SBP formal package robustness summary
    curated_interpretation = read_tsv(INPUTS["phase4_2B_curated_package_interpretation"])
    if curated_interpretation:
        for i, r in enumerate(curated_interpretation, start=1):
            add_master(
                master,
                evidence_id="SBP_formal_package_interpretation_{}".format(i),
                evidence_layer="SBP_robustness_formal_packages",
                exposure_id="SBP",
                outcome_or_component=r.get("phase", r.get("item", "formal_package_interpretation")),
                role_in_story="SBP robustness and caveat documentation",
                analysis_status=r.get("status", "DOCUMENTED"),
                n_instruments="391",
                beta="NA",
                se="NA",
                pval="NA",
                qval="NA",
                direction="NA",
                heterogeneity_or_caveat_flag=r.get("key_result", "NA"),
                evidence_label="ROBUSTNESS_CAVEAT_DOCUMENTED",
                publication_weight="ROBUSTNESS_CONTEXT",
                manuscript_safe_interpretation=r.get("note", "NA"),
                main_caveat="Formal package robustness should be interpreted cautiously if method-specific results disagree.",
                source_file=INPUTS["phase4_2B_curated_package_interpretation"]
            )
    else:
        add_master(
            master,
            evidence_id="SBP_formal_package_interpretation_missing",
            evidence_layer="SBP_robustness_formal_packages",
            exposure_id="SBP",
            outcome_or_component="MR_PRESSO_MR_RAPS",
            role_in_story="SBP robustness and caveat documentation",
            analysis_status="SOURCE_MISSING",
            n_instruments="391",
            beta="NA",
            se="NA",
            pval="NA",
            qval="NA",
            direction="NA",
            heterogeneity_or_caveat_flag="NA",
            evidence_label="ROBUSTNESS_SOURCE_MISSING",
            publication_weight="ROBUSTNESS_CONTEXT",
            manuscript_safe_interpretation="Formal package robustness source was not found by Phase 6.1 script.",
            main_caveat="Do not overstate SBP robustness until curated formal-package results are available.",
            source_file=INPUTS["phase4_2B_curated_package_interpretation"]
        )

    # 3. Route-A external validation grid for SBP: IOP / NTG / POAG / HTG
    routeA_grid = read_tsv(INPUTS["phase5_10F_routeA_grid"])
    for r in routeA_grid:
        label = r.get("evidence_label", "NA")
        add_master(
            master,
            evidence_id="external_routeA_{}_{}".format(r.get("exposure_id", "NA"), r.get("outcome_id", "NA")),
            evidence_layer="external_routeA_validation",
            exposure_id=r.get("exposure_id", "NA"),
            outcome_or_component=r.get("outcome_id", "NA"),
            role_in_story=r.get("role_in_routeA", "NA"),
            analysis_status="EXTERNAL_VALIDATION_SUMMARIZED",
            n_instruments=r.get("n_instruments", "NA"),
            beta=r.get("beta", "NA"),
            se=r.get("se", "NA"),
            pval=r.get("pval", "NA"),
            qval=r.get("qval", "NA"),
            direction=r.get("direction", "NA"),
            heterogeneity_or_caveat_flag="See original external MR heterogeneity tables where available",
            evidence_label=label,
            publication_weight=pub_weight_from_label(label, r.get("exposure_id", ""), r.get("n_instruments", "")),
            manuscript_safe_interpretation=r.get("interpretation", "NA"),
            main_caveat="External validation is directional and not confirmatory; HTG unavailable where marked.",
            source_file=INPUTS["phase5_10F_routeA_grid"]
        )

    # 4. Neuroretinal external endpoints: RNFL / GCIPL
    neuro_rows = read_tsv(INPUTS["phase5_6_neuroretinal_mr"])
    for r in neuro_rows:
        if r.get("method") != "IVW_multiplicative_random_effects":
            continue

        exposure = r.get("exposure_id", "NA")
        outcome = r.get("outcome_id", "NA")
        label = evidence_label(r.get("beta"), r.get("pval"), r.get("qval_bh_ivw_random"))
        if exposure == "ART_STIFFNESS":
            label = "LOW_POWER_" + label

        add_master(
            master,
            evidence_id="external_neuroretinal_{}_{}".format(exposure, outcome),
            evidence_layer="external_neuroretinal_validation",
            exposure_id=exposure,
            outcome_or_component=outcome,
            role_in_story="External neuroretinal endophenotype triangulation",
            analysis_status="MR_ANALYZED",
            n_instruments=r.get("n_instruments", "NA"),
            beta=r.get("beta", "NA"),
            se=r.get("se", "NA"),
            pval=r.get("pval", "NA"),
            qval=r.get("qval_bh_ivw_random", "NA"),
            direction=r.get("direction", "NA"),
            heterogeneity_or_caveat_flag="Q_pval={};phi={}".format(r.get("Q_pval", "NA"), r.get("phi", "NA")),
            evidence_label=label,
            publication_weight=pub_weight_from_label(label, exposure, r.get("n_instruments", "")),
            manuscript_safe_interpretation="External neuroretinal endpoint did not provide confirmatory support unless otherwise noted.",
            main_caveat="Neuroretinal validation is external triangulation; ART_STIFFNESS has only 3 instruments.",
            source_file=INPUTS["phase5_6_neuroretinal_mr"]
        )

    # 5. ART_STIFFNESS external IOP and NTG summaries
    for source_key, layer_name in [
        ("phase5_9F_iop_summary", "external_IOP_validation"),
        ("phase5_10F_ntg_summary", "external_NTG_validation")
    ]:
        rows = read_tsv(INPUTS[source_key])
        for r in rows:
            exposure = r.get("exposure_id", "NA")
            if exposure != "ART_STIFFNESS":
                continue

            label = r.get("evidence_label", evidence_label(r.get("beta_random"), r.get("p_random"), r.get("q_random")))
            if not label.startswith("LOW_POWER"):
                label = "LOW_POWER_" + label

            add_master(
                master,
                evidence_id="{}_{}".format(layer_name, exposure),
                evidence_layer=layer_name,
                exposure_id=exposure,
                outcome_or_component=r.get("outcome_id", "NA"),
                role_in_story="Exploratory vascular-stiffness external validation",
                analysis_status="MR_ANALYZED_LOW_POWER",
                n_instruments=r.get("n_instruments", "NA"),
                beta=r.get("beta_random", "NA"),
                se=r.get("se_random", "NA"),
                pval=r.get("p_random", "NA"),
                qval=r.get("q_random", "NA"),
                direction=r.get("direction_random", "NA"),
                heterogeneity_or_caveat_flag="Q_pval={};phi={}".format(r.get("Q_pval", "NA"), r.get("phi", "NA")),
                evidence_label=label,
                publication_weight="LOW_POWER_EXPLORATORY",
                manuscript_safe_interpretation=r.get("interpretation", "Exploratory only; do not emphasize."),
                main_caveat=r.get("main_caveat", "Only 3 instruments."),
                source_file=INPUTS[source_key]
            )

    # 6. ART_STIFFNESS POAG external validation
    poag_rows = read_tsv(INPUTS["phase5_8E_poag_mr"])
    art_poag = get_method_row(poag_rows, "ART_STIFFNESS", "POAG", "IVW_multiplicative_random_effects")
    if art_poag:
        label = "LOW_POWER_" + evidence_label(
            art_poag.get("beta_log_odds", art_poag.get("beta")),
            art_poag.get("pval"),
            art_poag.get("qval_bh_ivw_random")
        )
        add_master(
            master,
            evidence_id="external_POAG_ART_STIFFNESS",
            evidence_layer="external_POAG_validation",
            exposure_id="ART_STIFFNESS",
            outcome_or_component="POAG",
            role_in_story="Exploratory vascular-stiffness clinical endpoint validation",
            analysis_status="MR_ANALYZED_LOW_POWER",
            n_instruments=art_poag.get("n_instruments", "NA"),
            beta=art_poag.get("beta_log_odds", art_poag.get("beta", "NA")),
            se=art_poag.get("se", "NA"),
            pval=art_poag.get("pval", "NA"),
            qval=art_poag.get("qval_bh_ivw_random", "NA"),
            direction=art_poag.get("direction", "NA"),
            heterogeneity_or_caveat_flag="Q_pval={};phi={}".format(art_poag.get("Q_pval", "NA"), art_poag.get("phi", "NA")),
            evidence_label=label,
            publication_weight="LOW_POWER_EXPLORATORY",
            manuscript_safe_interpretation="Arterial stiffness showed no reliable POAG validation evidence.",
            main_caveat="Only 3 instruments; do not emphasize as confirmatory.",
            source_file=INPUTS["phase5_8E_poag_mr"]
        )

    # 7. HTG source availability limitation
    htg_lock = read_tsv(INPUTS["phase5_11B_htg_lock"])
    if htg_lock:
        r = htg_lock[0]
        add_master(
            master,
            evidence_id="HTG_source_availability",
            evidence_layer="external_source_availability",
            exposure_id="SBP",
            outcome_or_component="HTG",
            role_in_story="Ideal IOP-dependent subtype validation endpoint",
            analysis_status=r.get("source_status", "NOT_AVAILABLE"),
            n_instruments="NA",
            beta="NA",
            se="NA",
            pval="NA",
            qval="NA",
            direction="NA",
            heterogeneity_or_caveat_flag="NA",
            evidence_label="NOT_AVAILABLE",
            publication_weight="IMPORTANT_LIMITATION",
            manuscript_safe_interpretation="HTG validation could not be completed because no verified full downloadable HTG summary-statistics source was locked.",
            main_caveat=r.get("note", "HTG full summary statistics unavailable."),
            source_file=INPUTS["phase5_11B_htg_lock"]
        )

    # Route-A master matrix
    routeA_rows = [
        {
            "routeA_layer": "Internal component contrast",
            "outcome_or_component": "GBS_IOPcomponent vs GBS_nonIOPcomponent",
            "SBP_direction_or_pattern": "IOP more positive than nonIOP",
            "evidence_status": "Primary hypothesis-generating signal",
            "manuscript_use": "Main result",
            "main_caveat": "Individual component estimates are nominal/heterogeneity-flagged; contrast is the strongest result."
        },
        {
            "routeA_layer": "External measured IOP",
            "outcome_or_component": "IOP",
            "SBP_direction_or_pattern": "Positive borderline",
            "evidence_status": "Directional support, not confirmatory",
            "manuscript_use": "Mechanistic triangulation",
            "main_caveat": "Random-effects p was borderline and not FDR-significant; heterogeneity present."
        },
        {
            "routeA_layer": "External NTG",
            "outcome_or_component": "NTG",
            "SBP_direction_or_pattern": "Negative null",
            "evidence_status": "Directional context",
            "manuscript_use": "Subtype triangulation",
            "main_caveat": "Not statistically significant; NTG scale/sample-size metadata require cautious reporting."
        },
        {
            "routeA_layer": "External POAG",
            "outcome_or_component": "POAG",
            "SBP_direction_or_pattern": "Negative null",
            "evidence_status": "No broad clinical endpoint support",
            "manuscript_use": "Contextual validation",
            "main_caveat": "No reliable broad POAG validation signal."
        },
        {
            "routeA_layer": "External neuroretinal endpoints",
            "outcome_or_component": "RNFL / GCIPL",
            "SBP_direction_or_pattern": "Null / non-confirmatory",
            "evidence_status": "No neuroretinal structural validation",
            "manuscript_use": "Contextual validation",
            "main_caveat": "No FDR-significant external neuroretinal evidence."
        },
        {
            "routeA_layer": "External HTG",
            "outcome_or_component": "HTG",
            "SBP_direction_or_pattern": "Unavailable",
            "evidence_status": "Key missing endpoint",
            "manuscript_use": "Important limitation",
            "main_caveat": "No verified full downloadable HTG summary statistics locked."
        }
    ]

    flags = [
        {
            "flag": "overall_result_strength",
            "status": "HYPOTHESIS_GENERATING_NOT_CONFIRMATORY",
            "interpretation": "The evidence supports a vascular/IOP-related directional model but does not establish confirmatory causality.",
            "manuscript_action": "Use cautious language throughout."
        },
        {
            "flag": "main_positive_feature",
            "status": "SBP_COMPONENT_DIVERGENCE",
            "interpretation": "The strongest evidence is SBP directional divergence between IOP-dependent and IOP-independent glaucoma components.",
            "manuscript_action": "Make component divergence, not single-outcome significance, the central result."
        },
        {
            "flag": "external_triangulation",
            "status": "PARTIALLY_DIRECTIONALLY_SUPPORTIVE",
            "interpretation": "SBP-to-IOP is positive borderline; SBP-to-NTG is negative null; POAG/RNFL/GCIPL are non-confirmatory.",
            "manuscript_action": "Frame as external triangulation supporting IOP-related interpretation, not proof."
        },
        {
            "flag": "major_limitation",
            "status": "HTG_NOT_AVAILABLE",
            "interpretation": "HTG would be the ideal IOP-dependent clinical subtype validation endpoint but could not be analyzed.",
            "manuscript_action": "State clearly as an important limitation."
        },
        {
            "flag": "ART_STIFFNESS",
            "status": "LOW_POWER_EXPLORATORY_ONLY",
            "interpretation": "Arterial stiffness has only 3 instruments and should not drive conclusions.",
            "manuscript_action": "Mention as exploratory directional context only."
        },
        {
            "flag": "HYPERTENSION",
            "status": "NOT_ANALYZABLE_UNDER_STRICT_GWS",
            "interpretation": "Hypertension liability lacked strict instruments and should not be interpreted as null causal evidence.",
            "manuscript_action": "Describe as instrument-availability limitation."
        },
        {
            "flag": "manuscript_readiness",
            "status": "READY_FOR_PHASE_6_2_CLAIMS_WORDING",
            "interpretation": "Master evidence table is sufficient to proceed to final claims and manuscript-safe wording.",
            "manuscript_action": "Proceed to Phase 6.2."
        }
    ]

    write_tsv(MASTER_OUT, MASTER_FIELDS, master)
    write_tsv(
        ROUTEA_OUT,
        ["routeA_layer", "outcome_or_component", "SBP_direction_or_pattern", "evidence_status", "manuscript_use", "main_caveat"],
        routeA_rows
    )
    write_tsv(
        FLAGS_OUT,
        ["flag", "status", "interpretation", "manuscript_action"],
        flags
    )
    write_tsv(
        SOURCE_OUT,
        ["source_name", "path", "exists", "n_rows", "note"],
        source_inventory
    )

    elapsed = time.time() - start

    primary_rows = sum(1 for r in master if r["publication_weight"] == "PRIMARY_HYPOTHESIS_GENERATING_SIGNAL")
    limitations = sum(1 for r in master if r["publication_weight"] == "IMPORTANT_LIMITATION")

    with open(STATUS_OUT, "w", encoding="utf-8", newline="") as f:
        f.write("phase\tstatus\tkey_result\tnote\n")
        f.write(f"Phase 6.1 master evidence table\tPASSED\tmaster_rows={len(master)};primary_rows={primary_rows};important_limitations={limitations}\tIntegrated internal, robustness, and external validation evidence\n")
        f.write("routeA_status\tHYPOTHESIS_GENERATING_NOT_CONFIRMATORY\tSBP component divergence plus IOP-positive/NTG-negative directional triangulation\tDo not overstate causality\n")
        f.write("HTG_status\tIMPORTANT_LIMITATION\tfull_summary_stats_not_locked\tHTG remains unavailable after final attempt\n")
        f.write("next_step\tTO_DO\tPhase 6.2 final claims and manuscript-safe wording\tUse master evidence table to define permitted claims\n")
        f.write(f"runtime\tINFO\t{elapsed:.3f}s\tPhase 6.1 runtime\n")

    with open(RUNTIME_OUT, "w", encoding="utf-8", newline="") as f:
        f.write("phase\telapsed_seconds\telapsed_human\n")
        f.write(f"Phase 6.1\t{elapsed:.3f}\t{elapsed:.1f}s\n")

    print("===== Phase 6.1 completed =====")
    print("Wrote:", STATUS_OUT)
    print("Wrote:", MASTER_OUT)
    print("Wrote:", ROUTEA_OUT)
    print("Wrote:", FLAGS_OUT)
    print("Wrote:", SOURCE_OUT)

if __name__ == "__main__":
    main()
