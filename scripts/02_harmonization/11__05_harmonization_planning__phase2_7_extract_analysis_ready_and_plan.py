#!/usr/bin/env python3
import gzip
import os
import time
import hashlib
from datetime import datetime
from collections import Counter

LOCK_DIR = "../../11_locked_harmonized_datasets"
HARM_DIR = "../../10_harmonized_datasets"
OUTDIR = "../../12_analysis_ready_datasets"
PAIRWISE_DIR = os.path.join(OUTDIR, "pairwise")
PLAN_DIR = "../../05_harmonization_planning"

os.makedirs(OUTDIR, exist_ok=True)
os.makedirs(PAIRWISE_DIR, exist_ok=True)

LOCK_MANIFEST = os.path.join(LOCK_DIR, "phase2_6_locked_harmonized_dataset_manifest.tsv")

MANIFEST_OUT = os.path.join(OUTDIR, "phase2_7_analysis_ready_manifest.tsv")
QC_OUT = os.path.join(OUTDIR, "phase2_7_analysis_ready_qc_summary.tsv")
RUNTIME_OUT = os.path.join(OUTDIR, "phase2_7_runtime_log.tsv")
STATUS_OUT = os.path.join(OUTDIR, "phase2_7_status.tsv")
PLAN_MD = os.path.join(OUTDIR, "phase2_7_mr_mvmr_planning.md")

ANALYSIS_COLUMNS = [
    "exposure_id",
    "outcome_id",
    "match_mode",
    "match_key",
    "chr",
    "pos",
    "SNP",
    "exposure_effect_allele",
    "exposure_other_allele",
    "outcome_effect_allele_original",
    "outcome_other_allele_original",
    "beta_exposure",
    "se_exposure",
    "pval_exposure",
    "minus_log10_p_exposure",
    "eaf_exposure",
    "n_exposure",
    "beta_outcome_harmonized",
    "se_outcome",
    "pval_outcome",
    "pval_for_log10_outcome",
    "harmonization_action",
    "outcome_source_dataset",
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

def normalize_path(path):
    return path.replace("\\", "/")

def read_tsv(path):
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        header = f.readline().rstrip("\n\r").split("\t")
        rows = []
        for line in f:
            parts = line.rstrip("\n\r").split("\t")
            if len(parts) == len(header):
                rows.append(dict(zip(header, parts)))
        return header, rows

def estimate_runtime_seconds(file_size_bytes):
    """
    Very rough estimate based on compressed file size.
    Assumes 20-80 MB/s effective gzip read/write throughput on Windows.
    Uses conservative 30 MB/s plus overhead.
    """
    mb = file_size_bytes / (1024 * 1024)
    base = mb / 30.0
    overhead = 3.0
    return base + overhead

def parse_outfile_name(exposure_id, outcome_suffix):
    return os.path.join(PAIRWISE_DIR, f"{exposure_id}__{outcome_suffix}.analysis_ready.tsv.gz")

def extract_one(row):
    exposure_id = row["exposure_id"]
    outcome_suffix = row["outcome_suffix"]
    in_file = normalize_path(row["locked_harmonized_file"])
    out_file = parse_outfile_name(exposure_id, outcome_suffix)
    tmp_file = out_file + ".tmp"

    if not os.path.exists(in_file):
        raise FileNotFoundError(f"Input harmonized file not found: {in_file}")

    in_size = os.path.getsize(in_file)
    estimated_sec = estimate_runtime_seconds(in_size)

    print()
    print(f"===== Extracting {exposure_id} vs {outcome_suffix} =====")
    print(f"Input: {in_file}")
    print(f"Output: {out_file}")
    print(f"Compressed input size: {in_size / (1024*1024):.1f} MB")
    print(f"Approx runtime for this file: {fmt_seconds(estimated_sec)}")
    print(f"Start: {now_str()}")

    t0 = time.time()

    n_input = 0
    n_include = 0
    n_bad = 0
    missing_core = Counter()
    action_counter = Counter()
    match_mode_counter = Counter()

    with gzip.open(in_file, "rt", encoding="utf-8", errors="replace", newline="") as fin, \
         gzip.open(tmp_file, "wt", encoding="utf-8", newline="\n") as fout:

        header = fin.readline().rstrip("\n\r").split("\t")
        cmap = {c: i for i, c in enumerate(header)}

        required = ANALYSIS_COLUMNS + ["include_in_main"]
        missing = [c for c in required if c not in cmap]
        if missing:
            raise RuntimeError(f"Missing required columns in {in_file}: {missing}")

        fout.write("\t".join(ANALYSIS_COLUMNS) + "\n")

        for line in fin:
            parts = line.rstrip("\n\r").split("\t")
            n_input += 1

            if len(parts) != len(header):
                n_bad += 1
                continue

            if parts[cmap["include_in_main"]] != "YES":
                continue

            out = [parts[cmap[c]] for c in ANALYSIS_COLUMNS]

            core_checks = {
                "beta_exposure": parts[cmap["beta_exposure"]],
                "se_exposure": parts[cmap["se_exposure"]],
                "pval_exposure": parts[cmap["pval_exposure"]],
                "beta_outcome_harmonized": parts[cmap["beta_outcome_harmonized"]],
                "se_outcome": parts[cmap["se_outcome"]],
                "pval_outcome": parts[cmap["pval_outcome"]],
            }

            for k, v in core_checks.items():
                if v in ("", ".", "NA"):
                    missing_core[k] += 1

            action_counter[parts[cmap["harmonization_action"]]] += 1
            match_mode_counter[parts[cmap["match_mode"]]] += 1

            fout.write("\t".join(out) + "\n")
            n_include += 1

    os.replace(tmp_file, out_file)

    elapsed = time.time() - t0
    out_size = os.path.getsize(out_file)
    out_sha = sha256_file(out_file)

    print(f"End: {now_str()}")
    print(f"Elapsed: {fmt_seconds(elapsed)}")
    print(f"Included rows written: {n_include}")
    print(f"Output size: {out_size / (1024*1024):.1f} MB")

    action_summary = ";".join(f"{k}={v}" for k, v in sorted(action_counter.items()))
    match_mode_summary = ";".join(f"{k}={v}" for k, v in sorted(match_mode_counter.items()))
    missing_core_summary = ";".join(f"{k}={v}" for k, v in sorted(missing_core.items())) if missing_core else "none"

    status = "ANALYSIS_READY"
    note = "include_in_main=YES rows only"
    if missing_core:
        status = "REVIEW_CORE_MISSING"
        note = "Included rows have missing core MR fields"

    if outcome_suffix == "IOPcc_coordinate_subset":
        note += "; IOPcc coordinate-resolvable subset only"
    if exposure_id == "BMI":
        note += "; BMI EAF caveat retained from Phase 2.2"
    if exposure_id == "CRP":
        note += "; retain minus_log10_p to avoid p-value underflow"

    return {
        "exposure_id": exposure_id,
        "outcome_suffix": outcome_suffix,
        "input_file": in_file,
        "analysis_ready_file": out_file,
        "input_file_size_bytes": in_size,
        "output_file_size_bytes": out_size,
        "output_sha256": out_sha,
        "n_input_harmonized_rows": n_input,
        "n_analysis_ready_rows": n_include,
        "n_bad_input_rows": n_bad,
        "action_summary": action_summary,
        "match_mode_summary": match_mode_summary,
        "missing_core_summary": missing_core_summary,
        "status": status,
        "note": note,
        "estimated_runtime_seconds": estimated_sec,
        "elapsed_seconds": elapsed,
    }

def write_outputs(results):
    manifest_cols = [
        "exposure_id",
        "outcome_suffix",
        "analysis_ready_file",
        "output_sha256",
        "n_analysis_ready_rows",
        "status",
        "note",
    ]

    qc_cols = [
        "exposure_id",
        "outcome_suffix",
        "input_file",
        "analysis_ready_file",
        "n_input_harmonized_rows",
        "n_analysis_ready_rows",
        "n_bad_input_rows",
        "action_summary",
        "match_mode_summary",
        "missing_core_summary",
        "status",
        "note",
    ]

    runtime_cols = [
        "exposure_id",
        "outcome_suffix",
        "input_file_size_bytes",
        "output_file_size_bytes",
        "estimated_runtime_seconds",
        "elapsed_seconds",
        "estimated_runtime_human",
        "elapsed_human",
    ]

    with open(MANIFEST_OUT, "w", encoding="utf-8", newline="\n") as f:
        f.write("\t".join(manifest_cols) + "\n")
        for r in results:
            f.write("\t".join(str(r[c]) for c in manifest_cols) + "\n")

    with open(QC_OUT, "w", encoding="utf-8", newline="\n") as f:
        f.write("\t".join(qc_cols) + "\n")
        for r in results:
            f.write("\t".join(str(r[c]) for c in qc_cols) + "\n")

    with open(RUNTIME_OUT, "w", encoding="utf-8", newline="\n") as f:
        f.write("\t".join(runtime_cols) + "\n")
        for r in results:
            row = dict(r)
            row["estimated_runtime_human"] = fmt_seconds(r["estimated_runtime_seconds"])
            row["elapsed_human"] = fmt_seconds(r["elapsed_seconds"])
            f.write("\t".join(str(row[c]) for c in runtime_cols) + "\n")

    n_total = len(results)
    n_ready = sum(1 for r in results if r["status"] == "ANALYSIS_READY")
    n_review = n_total - n_ready
    total_rows = sum(r["n_analysis_ready_rows"] for r in results)
    total_elapsed = sum(r["elapsed_seconds"] for r in results)

    with open(STATUS_OUT, "w", encoding="utf-8", newline="\n") as f:
        f.write("phase\tstatus\tkey_result\tnote\n")
        phase_status = "PASSED" if n_total == 18 and n_ready == 18 else "REVIEW_NEEDED"
        f.write(f"Phase 2.7 analysis-ready dataset extraction\t{phase_status}\t{n_ready}/{n_total} datasets analysis-ready\tinclude_in_main=YES rows extracted\n")
        f.write(f"analysis_ready_dataset_count\tINFO\t{n_ready}\tDatasets with status=ANALYSIS_READY\n")
        f.write(f"review_needed_dataset_count\tINFO\t{n_review}\tDatasets with missing core fields or other issues\n")
        f.write(f"total_analysis_ready_rows\tINFO\t{total_rows}\tSum across all analysis-ready pairwise datasets\n")
        f.write(f"total_elapsed_runtime\tINFO\t{fmt_seconds(total_elapsed)}\tSum of per-file extraction times\n")
        f.write("IOPcc_limitation\tDOCUMENTED\tcoordinate-resolvable subset only\tAffx-style IOPcc variants remain excluded without external mapping\n")
        f.write("BMI_EAF_caveat\tDOCUMENTED\tEAF caveat retained\tPalindromic variants remain excluded\n")
        f.write("CRP_pvalue_caveat\tDOCUMENTED\tminus_log10_p retained\tAvoid p-value underflow\n")

    with open(PLAN_MD, "w", encoding="utf-8", newline="\n") as f:
        f.write("# Phase 2.7 MR / MVMR Planning\n\n")
        f.write("## Current status\n\n")
        f.write("Analysis-ready pairwise datasets have been extracted from locked harmonized datasets.\n\n")
        f.write("## Analysis-ready extraction rule\n\n")
        f.write("Only rows with `include_in_main=YES` were retained.\n\n")
        f.write("## Output directory\n\n")
        f.write("`12_analysis_ready_datasets/pairwise`\n\n")
        f.write("## Analysis-ready columns\n\n")
        for c in ANALYSIS_COLUMNS:
            f.write(f"- {c}\n")
        f.write("\n")
        f.write("## Recommended next steps before MR\n\n")
        f.write("1. Decide instrument selection threshold per exposure, e.g. genome-wide significant variants or clumped instruments.\n")
        f.write("2. Perform LD clumping against a matched ancestry reference panel.\n")
        f.write("3. Generate per-exposure instrument sets.\n")
        f.write("4. Run pairwise MR for each exposure against:\n")
        f.write("   - GBS_nonIOPcomponent\n")
        f.write("   - GBS_IOPcomponent\n")
        f.write("   - IOPcc coordinate-resolvable subset, with documented limitation\n")
        f.write("5. Run sensitivity analyses: MR-Egger, weighted median, weighted mode, heterogeneity, leave-one-out where appropriate.\n")
        f.write("6. Plan MVMR only after finalizing instrument overlap and covariance assumptions.\n\n")
        f.write("## Do not do yet\n\n")
        f.write("- Do not interpret causal estimates before instrument selection and LD clumping.\n")
        f.write("- Do not use IOPcc Affx-style variants without external validated mapping.\n")
        f.write("- Do not include palindromic variants that were excluded during harmonization.\n\n")
        f.write("## Caveats\n\n")
        f.write("1. IOPcc is limited to coordinate-resolvable allele-aware matches.\n")
        f.write("2. BMI has EAF caveat from Phase 2.2.\n")
        f.write("3. CRP requires retention of minus_log10_p to avoid p-value underflow.\n")

def main():
    if not os.path.exists(LOCK_MANIFEST):
        raise FileNotFoundError(LOCK_MANIFEST)

    header, lock_rows = read_tsv(LOCK_MANIFEST)
    rows = [r for r in lock_rows if r.get("lock_status") == "LOCKED"]

    print("===== Phase 2.7 analysis-ready extraction =====")
    print(f"Start: {now_str()}")
    print(f"Locked manifest: {LOCK_MANIFEST}")
    print(f"Locked datasets found: {len(rows)}")
    print(f"Output directory: {OUTDIR}")

    total_estimated = sum(estimate_runtime_seconds(os.path.getsize(normalize_path(r["locked_harmonized_file"]))) for r in rows if os.path.exists(normalize_path(r["locked_harmonized_file"])))
    print(f"Approx total runtime: {fmt_seconds(total_estimated)}")
    print("Note: actual runtime depends on disk speed, gzip throughput, and antivirus scanning.")

    t_all = time.time()

    results = []
    for r in rows:
        results.append(extract_one(r))

    write_outputs(results)

    elapsed_all = time.time() - t_all

    print()
    print("===== Phase 2.7 completed =====")
    print(f"End: {now_str()}")
    print(f"Total elapsed: {fmt_seconds(elapsed_all)}")
    print(f"Wrote manifest: {MANIFEST_OUT}")
    print(f"Wrote QC summary: {QC_OUT}")
    print(f"Wrote runtime log: {RUNTIME_OUT}")
    print(f"Wrote final status: {STATUS_OUT}")
    print(f"Wrote planning file: {PLAN_MD}")

if __name__ == "__main__":
    main()
