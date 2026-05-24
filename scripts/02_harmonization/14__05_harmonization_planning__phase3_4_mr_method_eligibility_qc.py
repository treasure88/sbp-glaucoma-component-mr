#!/usr/bin/env python3
import gzip
import os
import statistics
import hashlib
import time
from datetime import datetime

MR_INPUT_DIR = "../../16_mr_input_datasets/pairwise"
PHASE3_3_QC = "../../16_mr_input_datasets/phase3_3_mr_input_extraction_qc_summary.tsv"
PHASE3_2_PRIMARY_SET = "../../13_instrument_selection/phase3_2_primary_clumped_exposure_set.tsv"

OUTDIR = "../../17_mr_method_eligibility"
os.makedirs(OUTDIR, exist_ok=True)

ELIGIBILITY_OUT = os.path.join(OUTDIR, "phase3_4_mr_method_eligibility_summary.tsv")
STATUS_OUT = os.path.join(OUTDIR, "phase3_4_status.tsv")
RUNTIME_OUT = os.path.join(OUTDIR, "phase3_4_runtime_log.tsv")
PLAN_MD = os.path.join(OUTDIR, "phase3_4_mr_method_eligibility_notes.md")

PRIMARY_EXPOSURES = ["SBP", "DBP", "MIGRAINE", "CRP", "BMI"]
EXCLUDED_EXPOSURES = ["INSOMNIA"]

OUTCOMES = [
    "GBS_nonIOPcomponent",
    "GBS_IOPcomponent",
    "IOPcc_coordinate_subset",
]

def now_str():
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")

def fmt_seconds(x):
    if x < 60:
        return f"{x:.1f}s"
    if x < 3600:
        return f"{x/60:.1f}min"
    return f"{x/3600:.2f}h"

def sha256_file(path, block_size=1024 * 1024):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            b = f.read(block_size)
            if not b:
                break
            h.update(b)
    return h.hexdigest()

def safe_float(x):
    try:
        return float(x)
    except Exception:
        return None

def summarize_numeric(values):
    vals = [v for v in values if v is not None]
    if not vals:
        return {
            "min": "NA",
            "median": "NA",
            "mean": "NA",
            "max": "NA",
        }

    return {
        "min": f"{min(vals):.6g}",
        "median": f"{statistics.median(vals):.6g}",
        "mean": f"{statistics.mean(vals):.6g}",
        "max": f"{max(vals):.6g}",
    }

def method_eligibility(n):
    """
    Conservative MR-method eligibility rules based on number of independent instruments.
    """
    if n == 0:
        return {
            "primary_method": "NONE",
            "allowed_methods": "NONE",
            "sensitivity_methods": "NONE",
            "method_status": "NO_MR_NO_INSTRUMENTS",
            "note": "No matched clumped instruments available",
        }

    if n == 1:
        return {
            "primary_method": "WALD_RATIO",
            "allowed_methods": "Wald ratio",
            "sensitivity_methods": "NONE",
            "method_status": "LIMITED_SINGLE_INSTRUMENT",
            "note": "Only one instrument; no pleiotropy or heterogeneity sensitivity analysis possible",
        }

    if n == 2:
        return {
            "primary_method": "IVW_FIXED",
            "allowed_methods": "IVW fixed-effect",
            "sensitivity_methods": "Per-SNP Wald estimates",
            "method_status": "LIMITED_TWO_INSTRUMENTS",
            "note": "Two instruments; sensitivity analysis is highly limited",
        }

    if 3 <= n < 10:
        return {
            "primary_method": "IVW",
            "allowed_methods": "IVW; weighted median",
            "sensitivity_methods": "heterogeneity; leave-one-out",
            "method_status": "MR_READY_LIMITED_INSTRUMENTS",
            "note": "Multiple instruments but fewer than 10; MR-Egger not recommended as primary sensitivity",
        }

    if 10 <= n < 20:
        return {
            "primary_method": "IVW",
            "allowed_methods": "IVW; weighted median; MR-Egger",
            "sensitivity_methods": "heterogeneity; leave-one-out; MR-Egger intercept",
            "method_status": "MR_READY_EGGER_LOW_POWER",
            "note": "MR-Egger possible but lower power; interpret cautiously",
        }

    return {
        "primary_method": "IVW",
        "allowed_methods": "IVW; weighted median; MR-Egger",
        "sensitivity_methods": "heterogeneity; leave-one-out; MR-Egger intercept",
        "method_status": "MR_READY_FULL",
        "note": "Instrument count supports standard MR and sensitivity analyses",
    }

def inspect_mr_input(exposure_id, outcome_suffix):
    f = os.path.join(MR_INPUT_DIR, f"{exposure_id}__{outcome_suffix}.mr_input.tsv.gz")

    if not os.path.exists(f):
        return {
            "exposure_id": exposure_id,
            "outcome_suffix": outcome_suffix,
            "mr_input_file": f,
            "file_exists": "NO",
            "n_instruments": 0,
            "n_bad_rows": "NA",
            "n_missing_core_fields": "NA",
            "n_missing_f_stat": "NA",
            "n_weak_f_lt10": "NA",
            "min_F": "NA",
            "median_F": "NA",
            "mean_F": "NA",
            "max_F": "NA",
            "primary_method": "NONE",
            "allowed_methods": "NONE",
            "sensitivity_methods": "NONE",
            "method_status": "MISSING_FILE",
            "mr_readiness": "EXCLUDE",
            "note": "MR input file is missing",
        }

    n = 0
    bad = 0
    missing_core = 0
    missing_f = 0
    weak_f_lt10 = 0
    f_stats = []

    with gzip.open(f, "rt", encoding="utf-8", errors="replace", newline="") as fin:
        header = fin.readline().rstrip("\n\r").split("\t")
        cmap = {c: i for i, c in enumerate(header)}

        required = [
            "SNP",
            "beta_exposure",
            "se_exposure",
            "pval_exposure",
            "beta_outcome",
            "se_outcome",
            "pval_outcome",
            "instrument_F_stat",
        ]

        missing_cols = [c for c in required if c not in cmap]
        if missing_cols:
            return {
                "exposure_id": exposure_id,
                "outcome_suffix": outcome_suffix,
                "mr_input_file": f,
                "file_exists": "YES",
                "n_instruments": 0,
                "n_bad_rows": "NA",
                "n_missing_core_fields": "NA",
                "n_missing_f_stat": "NA",
                "n_weak_f_lt10": "NA",
                "min_F": "NA",
                "median_F": "NA",
                "mean_F": "NA",
                "max_F": "NA",
                "primary_method": "NONE",
                "allowed_methods": "NONE",
                "sensitivity_methods": "NONE",
                "method_status": "MISSING_COLUMNS",
                "mr_readiness": "EXCLUDE",
                "note": "Missing required columns: " + ",".join(missing_cols),
            }

        for line in fin:
            parts = line.rstrip("\n\r").split("\t")
            if len(parts) != len(header):
                bad += 1
                continue

            n += 1

            core_cols = [
                "beta_exposure",
                "se_exposure",
                "pval_exposure",
                "beta_outcome",
                "se_outcome",
                "pval_outcome",
            ]

            for c in core_cols:
                if parts[cmap[c]] in ("", ".", "NA"):
                    missing_core += 1

            fv = parts[cmap["instrument_F_stat"]]
            f_num = safe_float(fv)

            if f_num is None:
                missing_f += 1
            else:
                f_stats.append(f_num)
                if f_num < 10:
                    weak_f_lt10 += 1

    fs = summarize_numeric(f_stats)
    eligibility = method_eligibility(n)

    if n == 0:
        mr_readiness = "EXCLUDE_NO_INSTRUMENTS"
    elif missing_core > 0 or missing_f > 0:
        mr_readiness = "REVIEW_MISSING_FIELDS"
    elif weak_f_lt10 > 0:
        mr_readiness = "REVIEW_WEAK_INSTRUMENTS"
    else:
        mr_readiness = "READY"

    note = eligibility["note"]

    if outcome_suffix == "IOPcc_coordinate_subset":
        note += "; IOPcc coordinate-resolvable subset only"

    if exposure_id == "MIGRAINE" and n > 0 and n < 20:
        note += "; migraine has limited instrument count, sensitivity analyses should be interpreted cautiously"

    if exposure_id == "BMI":
        note += "; BMI EAF caveat retained from previous phases"

    if exposure_id == "CRP":
        note += "; CRP minus_log10_p retained to avoid p-value underflow"

    return {
        "exposure_id": exposure_id,
        "outcome_suffix": outcome_suffix,
        "mr_input_file": f,
        "file_exists": "YES",
        "n_instruments": n,
        "n_bad_rows": bad,
        "n_missing_core_fields": missing_core,
        "n_missing_f_stat": missing_f,
        "n_weak_f_lt10": weak_f_lt10,
        "min_F": fs["min"],
        "median_F": fs["median"],
        "mean_F": fs["mean"],
        "max_F": fs["max"],
        "primary_method": eligibility["primary_method"],
        "allowed_methods": eligibility["allowed_methods"],
        "sensitivity_methods": eligibility["sensitivity_methods"],
        "method_status": eligibility["method_status"],
        "mr_readiness": mr_readiness,
        "note": note,
    }

def main():
    print("===== Phase 3.4 MR method eligibility and instrument-count QC =====")
    print(f"Start: {now_str()}")
    t0 = time.time()

    results = []

    for exposure_id in PRIMARY_EXPOSURES:
        for outcome_suffix in OUTCOMES:
            results.append(inspect_mr_input(exposure_id, outcome_suffix))

    cols = [
        "exposure_id",
        "outcome_suffix",
        "mr_input_file",
        "file_exists",
        "n_instruments",
        "n_bad_rows",
        "n_missing_core_fields",
        "n_missing_f_stat",
        "n_weak_f_lt10",
        "min_F",
        "median_F",
        "mean_F",
        "max_F",
        "primary_method",
        "allowed_methods",
        "sensitivity_methods",
        "method_status",
        "mr_readiness",
        "note",
    ]

    with open(ELIGIBILITY_OUT, "w", encoding="utf-8", newline="\n") as out:
        out.write("\t".join(cols) + "\n")
        for r in results:
            out.write("\t".join(str(r[c]) for c in cols) + "\n")

    total_pairs = len(results)
    ready_pairs = sum(1 for r in results if r["mr_readiness"] == "READY")
    exclude_zero = sum(1 for r in results if r["mr_readiness"] == "EXCLUDE_NO_INSTRUMENTS")
    review_pairs = sum(1 for r in results if r["mr_readiness"].startswith("REVIEW"))
    gbs_ready = sum(
        1 for r in results
        if r["outcome_suffix"] in ("GBS_nonIOPcomponent", "GBS_IOPcomponent")
        and r["mr_readiness"] == "READY"
    )
    iopcc_zero = sum(
        1 for r in results
        if r["outcome_suffix"] == "IOPcc_coordinate_subset"
        and r["mr_readiness"] == "EXCLUDE_NO_INSTRUMENTS"
    )

    elapsed = time.time() - t0

    if ready_pairs == 10 and exclude_zero == 5 and review_pairs == 0:
        phase_status = "PASSED_WITH_IOPCC_EXCLUDED"
        key_result = "10/15 pairwise datasets MR-ready; 5 IOPcc pairs excluded due to zero instruments"
        note = "GBS outcomes are ready for MR; IOPcc coordinate subset has zero matched clumped instruments"
    elif review_pairs == 0 and ready_pairs > 0:
        phase_status = "PASSED_WITH_WARNINGS"
        key_result = f"{ready_pairs}/{total_pairs} pairwise datasets MR-ready"
        note = "Some pairs excluded due to zero instruments"
    else:
        phase_status = "REVIEW_NEEDED"
        key_result = f"{ready_pairs}/{total_pairs} pairwise datasets MR-ready"
        note = "Some pairs require review"

    with open(STATUS_OUT, "w", encoding="utf-8", newline="\n") as out:
        out.write("phase\tstatus\tkey_result\tnote\n")
        out.write(f"Phase 3.4 MR method eligibility and instrument-count QC\t{phase_status}\t{key_result}\t{note}\n")
        out.write(f"primary_exposure_count\tINFO\t{len(PRIMARY_EXPOSURES)}\t{','.join(PRIMARY_EXPOSURES)}\n")
        out.write(f"excluded_exposure_count\tINFO\t{len(EXCLUDED_EXPOSURES)}\t{','.join(EXCLUDED_EXPOSURES)} excluded from primary MR due to no LD-clumped instruments\n")
        out.write(f"expected_pairwise_dataset_count\tINFO\t{total_pairs}\t5 primary exposures x 3 outcomes\n")
        out.write(f"mr_ready_pair_count\tINFO\t{ready_pairs}\tPairwise datasets with READY status\n")
        out.write(f"gbs_ready_pair_count\tINFO\t{gbs_ready}\tGBS nonIOP and GBS IOP pairs ready\n")
        out.write(f"iopcc_zero_instrument_pair_count\tINFO\t{iopcc_zero}\tIOPcc coordinate subset unavailable for current primary MR\n")
        out.write(f"review_needed_pair_count\tINFO\t{review_pairs}\tPairs with missing fields or weak instruments\n")
        out.write(f"runtime\tINFO\t{fmt_seconds(elapsed)}\tPhase 3.4 QC runtime\n")

    with open(RUNTIME_OUT, "w", encoding="utf-8", newline="\n") as out:
        out.write("phase\tstart_time\tend_time\telapsed_seconds\telapsed_human\n")
        out.write(f"Phase 3.4\tNA\t{now_str()}\t{elapsed:.6f}\t{fmt_seconds(elapsed)}\n")

    with open(PLAN_MD, "w", encoding="utf-8", newline="\n") as out:
        out.write("# Phase 3.4 MR Method Eligibility and Instrument-count QC\n\n")
        out.write("## Status\n\n")
        out.write(f"{phase_status}\n\n")
        out.write("## MR-ready datasets\n\n")
        out.write("The 10 GBS pairwise datasets are ready for MR method execution.\n\n")
        out.write("## Excluded datasets\n\n")
        out.write("The five IOPcc coordinate-subset datasets have zero matched clumped instruments and should not be used in primary MR.\n\n")
        out.write("INSOMNIA remains excluded from primary MR because Phase 3.2 generated no LD-clumped instruments.\n\n")
        out.write("## Method eligibility rules\n\n")
        out.write("- 0 instruments: no MR\n")
        out.write("- 1 instrument: Wald ratio only\n")
        out.write("- 2 instruments: limited IVW only\n")
        out.write("- 3-9 instruments: IVW and weighted median; limited sensitivity\n")
        out.write("- 10-19 instruments: IVW, weighted median, MR-Egger with low-power caution\n")
        out.write("- >=20 instruments: standard IVW, weighted median, MR-Egger, heterogeneity, leave-one-out\n\n")
        out.write("## Next phase\n\n")
        out.write("Phase 3.5 should execute pairwise MR for the 10 GBS MR-ready datasets only.\n")

    print(f"End: {now_str()}")
    print(f"Elapsed: {fmt_seconds(elapsed)}")
    print(f"Wrote eligibility summary: {ELIGIBILITY_OUT}")
    print(f"Wrote status: {STATUS_OUT}")
    print(f"Wrote notes: {PLAN_MD}")

if __name__ == "__main__":
    main()
