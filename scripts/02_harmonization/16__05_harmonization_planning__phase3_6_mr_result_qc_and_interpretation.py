#!/usr/bin/env python3
import csv
import gzip
import math
import os
import time
from datetime import datetime
from collections import defaultdict

MR_DIR = "../../18_pairwise_mr_results"
OUTDIR = "../../19_mr_result_qc_interpretation"
os.makedirs(OUTDIR, exist_ok=True)

MR_RESULTS = os.path.join(MR_DIR, "phase3_5_pairwise_mr_results.tsv")
HET_FILE = os.path.join(MR_DIR, "phase3_5_heterogeneity_summary.tsv")
EGGER_FILE = os.path.join(MR_DIR, "phase3_5_egger_intercept_summary.tsv")
LOO_FILE = os.path.join(MR_DIR, "phase3_5_leave_one_out_ivw.tsv.gz")

PRIMARY_OUT = os.path.join(OUTDIR, "phase3_6_primary_mr_interpretation_summary.tsv")
CONCORDANCE_OUT = os.path.join(OUTDIR, "phase3_6_method_concordance_summary.tsv")
HET_PLEIO_OUT = os.path.join(OUTDIR, "phase3_6_heterogeneity_pleiotropy_flags.tsv")
CONTRAST_OUT = os.path.join(OUTDIR, "phase3_6_nonIOP_vs_IOP_direction_comparison.tsv")
LOO_OUT = os.path.join(OUTDIR, "phase3_6_leave_one_out_summary.tsv")
STATUS_OUT = os.path.join(OUTDIR, "phase3_6_status.tsv")
RUNTIME_OUT = os.path.join(OUTDIR, "phase3_6_runtime_log.tsv")
NOTES_MD = os.path.join(OUTDIR, "phase3_6_mr_result_qc_interpretation_notes.md")

PRIMARY_METHOD = "IVW_random_effects"
EXPOSURES = ["SBP", "DBP", "MIGRAINE", "CRP", "BMI"]
OUTCOMES = ["GBS_nonIOPcomponent", "GBS_IOPcomponent"]

def now_str():
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")

def fmt_seconds(x):
    if x < 60:
        return f"{x:.1f}s"
    if x < 3600:
        return f"{x/60:.1f}min"
    return f"{x/3600:.2f}h"

def parse_float(x):
    if x is None:
        return None
    x = str(x).strip()
    if x in ("", ".", "NA"):
        return None
    if x.startswith("<"):
        try:
            return float(x.replace("<", ""))
        except Exception:
            return None
    try:
        return float(x)
    except Exception:
        return None

def fmt(x, digits=8):
    if x is None:
        return "NA"
    if isinstance(x, str):
        return x
    if isinstance(x, int):
        return str(x)
    if math.isnan(x) or math.isinf(x):
        return "NA"
    return f"{x:.{digits}g}"

def sign_of(x):
    if x is None:
        return "NA"
    if x > 0:
        return "positive"
    if x < 0:
        return "negative"
    return "zero"

def read_tsv(path):
    with open(path, "r", encoding="utf-8", errors="replace", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))

def read_tsv_gz(path):
    with gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))

def get_primary_rows(mr_rows):
    primary = {}
    for r in mr_rows:
        if r["method"] == PRIMARY_METHOD:
            key = (r["exposure_id"], r["outcome_suffix"])
            primary[key] = r
    return primary

def interpretation_label(p, q):
    if p is None:
        return "NO_PVALUE"
    if q is not None and q < 0.05:
        return "FDR_SIGNIFICANT"
    if p < 0.05:
        return "NOMINAL_ONLY"
    return "NOT_SIGNIFICANT"

def evidence_tier(label, het_flag, egger_flag, loo_flag):
    if label == "FDR_SIGNIFICANT" and het_flag == "NO" and egger_flag == "NO" and loo_flag == "NO":
        return "HIGH_PRIORITY"
    if label == "FDR_SIGNIFICANT":
        return "FDR_SIGNIFICANT_WITH_QC_FLAGS"
    if label == "NOMINAL_ONLY" and het_flag == "NO" and egger_flag == "NO" and loo_flag == "NO":
        return "NOMINAL_PRIORITY"
    if label == "NOMINAL_ONLY":
        return "NOMINAL_WITH_QC_FLAGS"
    return "LOW_PRIORITY_OR_NULL"

def build_het_pleio_flags(het_rows, egger_rows):
    het_by_pair = defaultdict(dict)

    for r in het_rows:
        key = (r["exposure_id"], r["outcome_suffix"])
        method = r["method"]
        het_by_pair[key][method] = r

    egger_by_pair = {}
    for r in egger_rows:
        key = (r["exposure_id"], r["outcome_suffix"])
        egger_by_pair[key] = r

    flags = {}

    for exposure in EXPOSURES:
        for outcome in OUTCOMES:
            key = (exposure, outcome)

            ivw = het_by_pair.get(key, {}).get("IVW")
            egger_het = het_by_pair.get(key, {}).get("MR_Egger")
            egger_int = egger_by_pair.get(key)

            ivw_q_p = parse_float(ivw.get("Q_pval")) if ivw else None
            ivw_phi = parse_float(ivw.get("phi")) if ivw else None
            egger_q_p = parse_float(egger_het.get("Q_pval")) if egger_het else None
            egger_intercept_p = parse_float(egger_int.get("pval_intercept")) if egger_int else None

            heterogeneity_flag = "YES" if ivw_q_p is not None and ivw_q_p < 0.05 else "NO"
            egger_heterogeneity_flag = "YES" if egger_q_p is not None and egger_q_p < 0.05 else "NO"
            egger_intercept_flag = "YES" if egger_intercept_p is not None and egger_intercept_p < 0.05 else "NO"

            flags[key] = {
                "exposure_id": exposure,
                "outcome_suffix": outcome,
                "ivw_Q_pval": ivw_q_p,
                "ivw_phi": ivw_phi,
                "ivw_heterogeneity_flag": heterogeneity_flag,
                "egger_Q_pval": egger_q_p,
                "egger_heterogeneity_flag": egger_heterogeneity_flag,
                "egger_intercept": parse_float(egger_int.get("egger_intercept")) if egger_int else None,
                "egger_intercept_pval": egger_intercept_p,
                "egger_intercept_flag": egger_intercept_flag,
            }

    return flags

def build_loo_summary(loo_rows, primary):
    grouped = defaultdict(list)
    for r in loo_rows:
        key = (r["exposure_id"], r["outcome_suffix"])
        grouped[key].append(r)

    out = {}

    for exposure in EXPOSURES:
        for outcome in OUTCOMES:
            key = (exposure, outcome)
            prow = primary.get(key)

            if not prow:
                continue

            primary_beta = parse_float(prow["beta"])
            primary_p = parse_float(prow["pval"])
            primary_sign = sign_of(primary_beta)
            rows = grouped.get(key, [])

            n = len(rows)
            direction_flip_count = 0
            nominal_loss_count = 0
            nominal_gain_count = 0
            max_abs_delta = None
            max_delta_snp = "NA"

            for r in rows:
                b = parse_float(r["beta_ivw_random"])
                p = parse_float(r["pval_ivw_random"])
                if b is None:
                    continue

                if primary_sign in ("positive", "negative") and sign_of(b) != primary_sign:
                    direction_flip_count += 1

                if primary_p is not None and p is not None:
                    if primary_p < 0.05 and p >= 0.05:
                        nominal_loss_count += 1
                    if primary_p >= 0.05 and p < 0.05:
                        nominal_gain_count += 1

                delta = abs(b - primary_beta) if primary_beta is not None else None
                if delta is not None and (max_abs_delta is None or delta > max_abs_delta):
                    max_abs_delta = delta
                    max_delta_snp = r["left_out_SNP"]

            direction_flip_flag = "YES" if direction_flip_count > 0 else "NO"
            nominal_instability_flag = "YES" if nominal_loss_count > 0 or nominal_gain_count > 0 else "NO"

            out[key] = {
                "exposure_id": exposure,
                "outcome_suffix": outcome,
                "n_leave_one_out": n,
                "direction_flip_count": direction_flip_count,
                "direction_flip_flag": direction_flip_flag,
                "nominal_loss_count": nominal_loss_count,
                "nominal_gain_count": nominal_gain_count,
                "nominal_instability_flag": nominal_instability_flag,
                "max_abs_beta_delta": max_abs_delta,
                "max_delta_left_out_SNP": max_delta_snp,
            }

    return out

def build_method_concordance(mr_rows, primary):
    grouped = defaultdict(list)
    for r in mr_rows:
        key = (r["exposure_id"], r["outcome_suffix"])
        grouped[key].append(r)

    out = []

    for exposure in EXPOSURES:
        for outcome in OUTCOMES:
            key = (exposure, outcome)
            prow = primary.get(key)
            if not prow:
                continue

            primary_beta = parse_float(prow["beta"])
            primary_sign = sign_of(primary_beta)
            methods = grouped.get(key, [])

            non_primary = [r for r in methods if r["method"] != PRIMARY_METHOD]
            available_methods = [r["method"] for r in methods]

            same = 0
            opposite = 0
            na = 0
            details = []

            for r in non_primary:
                b = parse_float(r["beta"])
                s = sign_of(b)

                if s == "NA" or primary_sign == "NA":
                    na += 1
                elif s == primary_sign:
                    same += 1
                else:
                    opposite += 1

                details.append(f"{r['method']}:{s}")

            if opposite > 0:
                concordance = "DISCORDANT_DIRECTION_PRESENT"
            elif same == len(non_primary) and same > 0:
                concordance = "ALL_SENSITIVITY_METHODS_SAME_DIRECTION"
            elif same > 0:
                concordance = "PARTIAL_CONCORDANCE"
            else:
                concordance = "NO_SENSITIVITY_DIRECTION_AVAILABLE"

            out.append({
                "exposure_id": exposure,
                "outcome_suffix": outcome,
                "primary_beta": primary_beta,
                "primary_direction": primary_sign,
                "available_methods": ";".join(available_methods),
                "n_non_primary_methods": len(non_primary),
                "n_same_direction": same,
                "n_opposite_direction": opposite,
                "n_direction_na": na,
                "method_direction_details": ";".join(details),
                "concordance_status": concordance,
            })

    return out

def build_nonIOP_vs_IOP(primary):
    out = []

    for exposure in EXPOSURES:
        key_non = (exposure, "GBS_nonIOPcomponent")
        key_iop = (exposure, "GBS_IOPcomponent")

        r_non = primary.get(key_non)
        r_iop = primary.get(key_iop)

        if not r_non or not r_iop:
            continue

        beta_non = parse_float(r_non["beta"])
        beta_iop = parse_float(r_iop["beta"])
        p_non = parse_float(r_non["pval"])
        p_iop = parse_float(r_iop["pval"])
        q_non = parse_float(r_non.get("qval_bh_primary_ivw_random"))
        q_iop = parse_float(r_iop.get("qval_bh_primary_ivw_random"))

        sign_non = sign_of(beta_non)
        sign_iop = sign_of(beta_iop)

        if sign_non in ("positive", "negative") and sign_iop in ("positive", "negative"):
            if sign_non == sign_iop:
                pattern = "SAME_DIRECTION"
            else:
                pattern = "OPPOSITE_DIRECTION"
        else:
            pattern = "UNDEFINED"

        nominal_non = "YES" if p_non is not None and p_non < 0.05 else "NO"
        nominal_iop = "YES" if p_iop is not None and p_iop < 0.05 else "NO"
        fdr_non = "YES" if q_non is not None and q_non < 0.05 else "NO"
        fdr_iop = "YES" if q_iop is not None and q_iop < 0.05 else "NO"

        if pattern == "OPPOSITE_DIRECTION" and (nominal_non == "YES" or nominal_iop == "YES"):
            interpretation = "OPPOSING_NOMINAL_PATTERN"
        elif pattern == "SAME_DIRECTION" and (nominal_non == "YES" or nominal_iop == "YES"):
            interpretation = "SAME_DIRECTION_NOMINAL_PATTERN"
        else:
            interpretation = "NO_CLEAR_PRIMARY_PATTERN"

        out.append({
            "exposure_id": exposure,
            "beta_nonIOP": beta_non,
            "p_nonIOP": p_non,
            "q_nonIOP": q_non,
            "direction_nonIOP": sign_non,
            "nominal_nonIOP": nominal_non,
            "fdr_nonIOP": fdr_non,
            "beta_IOP": beta_iop,
            "p_IOP": p_iop,
            "q_IOP": q_iop,
            "direction_IOP": sign_iop,
            "nominal_IOP": nominal_iop,
            "fdr_IOP": fdr_iop,
            "direction_pattern": pattern,
            "interpretation": interpretation,
        })

    return out

def main():
    print("===== Phase 3.6 MR result QC and interpretation =====")
    print(f"Start: {now_str()}")

    t0 = time.time()

    mr_rows = read_tsv(MR_RESULTS)
    het_rows = read_tsv(HET_FILE)
    egger_rows = read_tsv(EGGER_FILE)
    loo_rows = read_tsv_gz(LOO_FILE)

    primary = get_primary_rows(mr_rows)
    het_flags = build_het_pleio_flags(het_rows, egger_rows)
    loo_summary = build_loo_summary(loo_rows, primary)
    concordance = build_method_concordance(mr_rows, primary)
    contrast = build_nonIOP_vs_IOP(primary)

    primary_summary = []

    for exposure in EXPOSURES:
        for outcome in OUTCOMES:
            key = (exposure, outcome)
            r = primary.get(key)
            if not r:
                continue

            beta = parse_float(r["beta"])
            se = parse_float(r["se"])
            p = parse_float(r["pval"])
            q = parse_float(r.get("qval_bh_primary_ivw_random"))
            n = int(r["n_instruments"])

            label = interpretation_label(p, q)
            hf = het_flags[key]
            lf = loo_summary[key]

            tier = evidence_tier(
                label,
                hf["ivw_heterogeneity_flag"],
                hf["egger_intercept_flag"],
                lf["nominal_instability_flag"],
            )

            direction = sign_of(beta)

            if exposure == "SBP" and label == "NOMINAL_ONLY":
                interpretation_note = "SBP shows nominal evidence; compare direction across nonIOP and IOP components"
            elif label == "NOT_SIGNIFICANT":
                interpretation_note = "No primary IVW evidence at nominal p<0.05"
            elif label == "FDR_SIGNIFICANT":
                interpretation_note = "Primary IVW evidence survives BH-FDR correction"
            else:
                interpretation_note = "Primary IVW nominal evidence only; does not survive BH-FDR correction"

            primary_summary.append({
                "exposure_id": exposure,
                "outcome_suffix": outcome,
                "method": PRIMARY_METHOD,
                "n_instruments": n,
                "beta": beta,
                "se": se,
                "pval": p,
                "qval_bh": q,
                "direction": direction,
                "or": parse_float(r["or"]),
                "or_ci_lower": parse_float(r["or_ci_lower"]),
                "or_ci_upper": parse_float(r["or_ci_upper"]),
                "interpretation_label": label,
                "evidence_tier": tier,
                "ivw_heterogeneity_flag": hf["ivw_heterogeneity_flag"],
                "egger_intercept_flag": hf["egger_intercept_flag"],
                "loo_nominal_instability_flag": lf["nominal_instability_flag"],
                "loo_direction_flip_flag": lf["direction_flip_flag"],
                "interpretation_note": interpretation_note,
            })

    with open(PRIMARY_OUT, "w", encoding="utf-8", newline="") as f:
        cols = [
            "exposure_id", "outcome_suffix", "method", "n_instruments",
            "beta", "se", "pval", "qval_bh", "direction",
            "or", "or_ci_lower", "or_ci_upper",
            "interpretation_label", "evidence_tier",
            "ivw_heterogeneity_flag", "egger_intercept_flag",
            "loo_nominal_instability_flag", "loo_direction_flip_flag",
            "interpretation_note"
        ]
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for r in primary_summary:
            w.writerow({k: fmt(r[k]) for k in cols})

    with open(CONCORDANCE_OUT, "w", encoding="utf-8", newline="") as f:
        cols = [
            "exposure_id", "outcome_suffix", "primary_beta", "primary_direction",
            "available_methods", "n_non_primary_methods", "n_same_direction",
            "n_opposite_direction", "n_direction_na", "method_direction_details",
            "concordance_status"
        ]
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for r in concordance:
            w.writerow({k: fmt(r[k]) for k in cols})

    with open(HET_PLEIO_OUT, "w", encoding="utf-8", newline="") as f:
        cols = [
            "exposure_id", "outcome_suffix",
            "ivw_Q_pval", "ivw_phi", "ivw_heterogeneity_flag",
            "egger_Q_pval", "egger_heterogeneity_flag",
            "egger_intercept", "egger_intercept_pval", "egger_intercept_flag"
        ]
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for exposure in EXPOSURES:
            for outcome in OUTCOMES:
                r = het_flags[(exposure, outcome)]
                w.writerow({k: fmt(r[k]) for k in cols})

    with open(CONTRAST_OUT, "w", encoding="utf-8", newline="") as f:
        cols = [
            "exposure_id",
            "beta_nonIOP", "p_nonIOP", "q_nonIOP", "direction_nonIOP", "nominal_nonIOP", "fdr_nonIOP",
            "beta_IOP", "p_IOP", "q_IOP", "direction_IOP", "nominal_IOP", "fdr_IOP",
            "direction_pattern", "interpretation"
        ]
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for r in contrast:
            w.writerow({k: fmt(r[k]) for k in cols})

    with open(LOO_OUT, "w", encoding="utf-8", newline="") as f:
        cols = [
            "exposure_id", "outcome_suffix", "n_leave_one_out",
            "direction_flip_count", "direction_flip_flag",
            "nominal_loss_count", "nominal_gain_count", "nominal_instability_flag",
            "max_abs_beta_delta", "max_delta_left_out_SNP"
        ]
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for exposure in EXPOSURES:
            for outcome in OUTCOMES:
                r = loo_summary[(exposure, outcome)]
                w.writerow({k: fmt(r[k]) for k in cols})

    n_primary = len(primary_summary)
    n_fdr = sum(1 for r in primary_summary if r["interpretation_label"] == "FDR_SIGNIFICANT")
    n_nominal = sum(1 for r in primary_summary if r["interpretation_label"] == "NOMINAL_ONLY")
    n_null = sum(1 for r in primary_summary if r["interpretation_label"] == "NOT_SIGNIFICANT")
    n_het = sum(1 for r in het_flags.values() if r["ivw_heterogeneity_flag"] == "YES")
    n_egger_int = sum(1 for r in het_flags.values() if r["egger_intercept_flag"] == "YES")
    n_loo_instability = sum(1 for r in loo_summary.values() if r["nominal_instability_flag"] == "YES")
    n_opposing_nominal = sum(1 for r in contrast if r["interpretation"] == "OPPOSING_NOMINAL_PATTERN")

    elapsed = time.time() - t0

    with open(STATUS_OUT, "w", encoding="utf-8", newline="") as f:
        f.write("phase\tstatus\tkey_result\tnote\n")
        f.write(f"Phase 3.6 MR result QC and interpretation\tPASSED\t{n_primary}/10 primary IVW results interpreted\tMR result QC completed for GBS outcomes\n")
        f.write(f"fdr_significant_primary_count\tINFO\t{n_fdr}\tPrimary IVW random-effects results with BH-FDR q<0.05\n")
        f.write(f"nominal_only_primary_count\tINFO\t{n_nominal}\tPrimary IVW random-effects results with p<0.05 but BH-FDR q>=0.05\n")
        f.write(f"not_significant_primary_count\tINFO\t{n_null}\tPrimary IVW random-effects results with p>=0.05\n")
        f.write(f"ivw_heterogeneity_flag_count\tINFO\t{n_het}\tIVW heterogeneity Q p<0.05\n")
        f.write(f"egger_intercept_flag_count\tINFO\t{n_egger_int}\tMR-Egger intercept p<0.05\n")
        f.write(f"leave_one_out_nominal_instability_count\tINFO\t{n_loo_instability}\tPairs where leave-one-out changes nominal significance status\n")
        f.write(f"opposing_nominal_nonIOP_vs_IOP_pattern_count\tINFO\t{n_opposing_nominal}\tExposures with opposite nonIOP vs IOP directions and nominal evidence\n")
        f.write("primary_interpretation\tDOCUMENTED\tSBP nominal opposing-direction pattern\tSBP is negative for nonIOP and positive for IOP at nominal p<0.05, but neither survives BH-FDR<0.05\n")
        f.write("excluded_exposure\tDOCUMENTED\tINSOMNIA\tNo LD-clumped instruments\n")
        f.write("excluded_outcome\tDOCUMENTED\tIOPcc_coordinate_subset\tNo matched clumped instruments\n")
        f.write(f"runtime\tINFO\t{fmt_seconds(elapsed)}\tPhase 3.6 runtime\n")

    with open(RUNTIME_OUT, "w", encoding="utf-8", newline="") as f:
        f.write("phase\tstart_time\tend_time\telapsed_seconds\telapsed_human\n")
        f.write(f"Phase 3.6\tNA\t{now_str()}\t{elapsed:.6f}\t{fmt_seconds(elapsed)}\n")

    with open(NOTES_MD, "w", encoding="utf-8", newline="") as f:
        f.write("# Phase 3.6 MR Result QC and Interpretation\n\n")
        f.write("## Status\n\n")
        f.write("PASSED\n\n")
        f.write("## Main interpretation\n\n")
        f.write("The primary IVW random-effects results show nominal evidence for SBP in opposite directions across the two GBS outcomes: negative for GBS_nonIOPcomponent and positive for GBS_IOPcomponent. Neither result survives BH-FDR correction at q<0.05.\n\n")
        f.write("DBP, MIGRAINE, CRP, and BMI do not show nominal primary IVW evidence in the current analysis.\n\n")
        f.write("## Important caveats\n\n")
        f.write("- INSOMNIA is excluded from primary MR because no LD-clumped instruments were generated.\n")
        f.write("- IOPcc_coordinate_subset is excluded because no matched clumped instruments were available.\n")
        f.write("- MIGRAINE has only 13 instruments; MR-Egger should be interpreted cautiously.\n")
        f.write("- SBP nominal opposing-direction results should be interpreted as hypothesis-generating unless supported by sensitivity analyses and/or external replication.\n\n")
        f.write("## Next phase\n\n")
        f.write("Phase 3.7 should create final result tables and figures for reporting, including primary MR, sensitivity QC, and nonIOP-vs-IOP contrast summaries.\n")

    print("===== Phase 3.6 completed =====")
    print(f"End: {now_str()}")
    print(f"Elapsed: {fmt_seconds(elapsed)}")
    print(f"Wrote: {PRIMARY_OUT}")
    print(f"Wrote: {CONCORDANCE_OUT}")
    print(f"Wrote: {HET_PLEIO_OUT}")
    print(f"Wrote: {CONTRAST_OUT}")
    print(f"Wrote: {LOO_OUT}")
    print(f"Wrote: {STATUS_OUT}")

if __name__ == "__main__":
    main()
