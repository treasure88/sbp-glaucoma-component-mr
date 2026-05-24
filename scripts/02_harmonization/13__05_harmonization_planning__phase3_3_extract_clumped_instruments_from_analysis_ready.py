#!/usr/bin/env python3
import gzip
import os
import time
import hashlib
from datetime import datetime
from collections import Counter

CLUMPED_DIR = "../../13_instrument_selection/clumped_instruments"
ANALYSIS_READY_DIR = "../../12_analysis_ready_datasets/pairwise"
PRIMARY_SET_FILE = "../../13_instrument_selection/phase3_2_primary_clumped_exposure_set.tsv"

OUTDIR = "../../16_mr_input_datasets"
PAIRWISE_OUTDIR = os.path.join(OUTDIR, "pairwise")
MISSING_OUTDIR = os.path.join(OUTDIR, "missing_instruments")

os.makedirs(OUTDIR, exist_ok=True)
os.makedirs(PAIRWISE_OUTDIR, exist_ok=True)
os.makedirs(MISSING_OUTDIR, exist_ok=True)

OUTCOMES = [
    "GBS_nonIOPcomponent",
    "GBS_IOPcomponent",
    "IOPcc_coordinate_subset",
]

DATASET_IDS = {
    "SBP": "ieu-b-38",
    "DBP": "ieu-b-39",
    "MIGRAINE": "ebi-a-GCST90038646",
    "INSOMNIA": "ebi-a-GCST90018869",
    "CRP": "ebi-a-GCST90018950",
    "BMI": "ieu-a-2",
}

QC_OUT = os.path.join(OUTDIR, "phase3_3_mr_input_extraction_qc_summary.tsv")
MANIFEST_OUT = os.path.join(OUTDIR, "phase3_3_mr_input_manifest.tsv")
RUNTIME_OUT = os.path.join(OUTDIR, "phase3_3_runtime_log.tsv")
STATUS_OUT = os.path.join(OUTDIR, "phase3_3_status.tsv")
NOTES_MD = os.path.join(OUTDIR, "phase3_3_mr_input_extraction_notes.md")

MR_INPUT_COLUMNS = [
    "exposure_id",
    "outcome_id",
    "SNP",
    "chr",
    "pos",
    "effect_allele",
    "other_allele",
    "beta_exposure",
    "se_exposure",
    "pval_exposure",
    "minus_log10_p_exposure",
    "eaf_exposure",
    "n_exposure",
    "beta_outcome",
    "se_outcome",
    "pval_outcome",
    "pval_for_log10_outcome",
    "harmonization_action",
    "match_mode",
    "match_key",
    "clump_rank",
    "instrument_F_stat",
    "source_analysis_ready_file",
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

def is_missing(x):
    return x in ("", ".", "NA")

def safe_float(x):
    try:
        return float(x)
    except Exception:
        return None

def f_stat(beta, se):
    b = safe_float(beta)
    s = safe_float(se)
    if b is None or s is None or s == 0:
        return "NA"
    return f"{(b / s) ** 2:.8g}"

def read_primary_exposures():
    exposures = []

    if os.path.exists(PRIMARY_SET_FILE):
        with open(PRIMARY_SET_FILE, "r", encoding="utf-8", errors="replace") as f:
            header = f.readline().rstrip("\n\r").split("\t")
            cmap = {c: i for i, c in enumerate(header)}

            for line in f:
                parts = line.rstrip("\n\r").split("\t")
                if len(parts) != len(header):
                    continue

                exposure_id = parts[cmap["exposure_id"]]
                status = parts[cmap["primary_mr_status"]]

                if status == "INCLUDE":
                    exposures.append(exposure_id)

    if not exposures:
        exposures = ["SBP", "DBP", "MIGRAINE", "CRP", "BMI"]

    return exposures

def load_clumped_instruments(exposure_id):
    dataset_id = DATASET_IDS[exposure_id]
    f = os.path.join(CLUMPED_DIR, f"{exposure_id}__{dataset_id}.clumped_instruments.tsv.gz")

    if not os.path.exists(f):
        raise FileNotFoundError(f)

    instruments = {}
    order = []

    with gzip.open(f, "rt", encoding="utf-8", errors="replace", newline="") as fin:
        header = fin.readline().rstrip("\n\r").split("\t")
        cmap = {c: i for i, c in enumerate(header)}

        required = ["SNP", "clump_rank"]
        missing = [c for c in required if c not in cmap]
        if missing:
            raise RuntimeError(f"Missing required clumped instrument columns in {f}: {missing}")

        for line in fin:
            parts = line.rstrip("\n\r").split("\t")
            if len(parts) != len(header):
                continue

            snp = parts[cmap["SNP"]]
            if snp in instruments:
                continue

            instruments[snp] = {
                "clump_rank": parts[cmap["clump_rank"]],
                "raw_parts": parts,
            }
            order.append(snp)

    return f, instruments, order

def extract_for_pair(exposure_id, outcome_suffix, instruments, instrument_order):
    analysis_file = os.path.join(
        ANALYSIS_READY_DIR,
        f"{exposure_id}__{outcome_suffix}.analysis_ready.tsv.gz"
    )

    if not os.path.exists(analysis_file):
        raise FileNotFoundError(analysis_file)

    out_file = os.path.join(
        PAIRWISE_OUTDIR,
        f"{exposure_id}__{outcome_suffix}.mr_input.tsv.gz"
    )

    missing_file = os.path.join(
        MISSING_OUTDIR,
        f"{exposure_id}__{outcome_suffix}.missing_clumped_instruments.tsv"
    )

    tmp_out = out_file + ".tmp"
    tmp_missing = missing_file + ".tmp"

    n_analysis_rows = 0
    n_bad_rows = 0
    n_matched_rows = 0
    n_duplicate_matched_snp = 0
    matched_snps = set()
    action_counter = Counter()
    missing_core_counter = Counter()

    t0 = time.time()

    print()
    print(f"===== Extract MR input: {exposure_id} vs {outcome_suffix} =====")
    print(f"Analysis-ready input: {analysis_file}")
    print(f"Clumped instruments: {len(instruments)}")
    print(f"Approx runtime: {fmt_seconds(estimate_runtime_seconds(analysis_file))}")
    print(f"Start: {now_str()}")

    with gzip.open(analysis_file, "rt", encoding="utf-8", errors="replace", newline="") as fin, \
         gzip.open(tmp_out, "wt", encoding="utf-8", newline="\n") as fout:

        header = fin.readline().rstrip("\n\r").split("\t")
        cmap = {c: i for i, c in enumerate(header)}

        required = [
            "exposure_id",
            "outcome_id",
            "match_mode",
            "match_key",
            "chr",
            "pos",
            "SNP",
            "exposure_effect_allele",
            "exposure_other_allele",
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
        ]

        missing = [c for c in required if c not in cmap]
        if missing:
            raise RuntimeError(f"Missing required analysis-ready columns in {analysis_file}: {missing}")

        fout.write("\t".join(MR_INPUT_COLUMNS) + "\n")

        for line in fin:
            parts = line.rstrip("\n\r").split("\t")
            n_analysis_rows += 1

            if len(parts) != len(header):
                n_bad_rows += 1
                continue

            snp = parts[cmap["SNP"]]
            if snp not in instruments:
                continue

            if snp in matched_snps:
                n_duplicate_matched_snp += 1
                continue

            matched_snps.add(snp)
            n_matched_rows += 1

            beta_exposure = parts[cmap["beta_exposure"]]
            se_exposure = parts[cmap["se_exposure"]]
            beta_outcome = parts[cmap["beta_outcome_harmonized"]]
            se_outcome = parts[cmap["se_outcome"]]
            pval_exposure = parts[cmap["pval_exposure"]]
            pval_outcome = parts[cmap["pval_outcome"]]

            core = {
                "beta_exposure": beta_exposure,
                "se_exposure": se_exposure,
                "beta_outcome": beta_outcome,
                "se_outcome": se_outcome,
                "pval_exposure": pval_exposure,
                "pval_outcome": pval_outcome,
            }

            for k, v in core.items():
                if is_missing(v):
                    missing_core_counter[k] += 1

            action = parts[cmap["harmonization_action"]]
            action_counter[action] += 1

            out = [
                parts[cmap["exposure_id"]],
                parts[cmap["outcome_id"]],
                snp,
                parts[cmap["chr"]],
                parts[cmap["pos"]],
                parts[cmap["exposure_effect_allele"]],
                parts[cmap["exposure_other_allele"]],
                beta_exposure,
                se_exposure,
                pval_exposure,
                parts[cmap["minus_log10_p_exposure"]],
                parts[cmap["eaf_exposure"]],
                parts[cmap["n_exposure"]],
                beta_outcome,
                se_outcome,
                pval_outcome,
                parts[cmap["pval_for_log10_outcome"]],
                action,
                parts[cmap["match_mode"]],
                parts[cmap["match_key"]],
                instruments[snp]["clump_rank"],
                f_stat(beta_exposure, se_exposure),
                analysis_file,
            ]

            fout.write("\t".join(out) + "\n")

    missing_snps = [snp for snp in instrument_order if snp not in matched_snps]

    with open(tmp_missing, "w", encoding="utf-8", newline="\n") as f:
        f.write("exposure_id\toutcome_suffix\tmissing_SNP\tclump_rank\tnote\n")
        for snp in missing_snps:
            f.write(
                f"{exposure_id}\t{outcome_suffix}\t{snp}\t"
                f"{instruments[snp]['clump_rank']}\t"
                f"Clumped exposure instrument not found in analysis-ready harmonized dataset\n"
            )

    os.replace(tmp_out, out_file)
    os.replace(tmp_missing, missing_file)

    elapsed = time.time() - t0

    print(f"End: {now_str()}")
    print(f"Elapsed: {fmt_seconds(elapsed)}")
    print(f"MR input rows: {n_matched_rows}")
    print(f"Missing clumped instruments: {len(missing_snps)}")

    missing_core_summary = "none"
    if missing_core_counter:
        missing_core_summary = ";".join(f"{k}={v}" for k, v in sorted(missing_core_counter.items()))

    action_summary = "none"
    if action_counter:
        action_summary = ";".join(f"{k}={v}" for k, v in sorted(action_counter.items()))

    status = "MR_INPUT_READY"
    note = "All extracted rows are include_in_main=YES from Phase 2.7 analysis-ready data"

    if n_matched_rows == 0:
        status = "NO_MATCHED_INSTRUMENTS"
        note = "No clumped instruments found in this outcome analysis-ready dataset"

    if missing_core_counter:
        status = "REVIEW_MISSING_CORE_FIELDS"
        note = "Some extracted instruments have missing core MR fields"

    if outcome_suffix == "IOPcc_coordinate_subset":
        note += "; IOPcc coordinate-resolvable subset only"

    if exposure_id == "BMI":
        note += "; BMI EAF caveat retained"

    if exposure_id == "CRP":
        note += "; CRP minus_log10_p retained to avoid p-value underflow"

    return {
        "exposure_id": exposure_id,
        "outcome_suffix": outcome_suffix,
        "analysis_ready_file": analysis_file,
        "mr_input_file": out_file,
        "missing_instruments_file": missing_file,
        "n_clumped_instruments": len(instruments),
        "n_analysis_ready_rows": n_analysis_rows,
        "n_mr_input_rows": n_matched_rows,
        "n_missing_clumped_instruments": len(missing_snps),
        "n_bad_analysis_ready_rows": n_bad_rows,
        "n_duplicate_matched_snp": n_duplicate_matched_snp,
        "missing_core_summary": missing_core_summary,
        "harmonization_action_summary": action_summary,
        "status": status,
        "note": note,
        "mr_input_sha256": sha256_file(out_file),
        "elapsed_seconds": elapsed,
    }

def estimate_runtime_seconds(path):
    size = os.path.getsize(path)
    mb = size / (1024 * 1024)
    return mb / 20.0 + 2.0

def write_outputs(results, primary_exposures):
    qc_cols = [
        "exposure_id",
        "outcome_suffix",
        "analysis_ready_file",
        "mr_input_file",
        "missing_instruments_file",
        "n_clumped_instruments",
        "n_analysis_ready_rows",
        "n_mr_input_rows",
        "n_missing_clumped_instruments",
        "n_bad_analysis_ready_rows",
        "n_duplicate_matched_snp",
        "missing_core_summary",
        "harmonization_action_summary",
        "status",
        "note",
    ]

    manifest_cols = [
        "exposure_id",
        "outcome_suffix",
        "mr_input_file",
        "mr_input_sha256",
        "n_mr_input_rows",
        "status",
        "note",
    ]

    runtime_cols = [
        "exposure_id",
        "outcome_suffix",
        "elapsed_seconds",
        "elapsed_human",
    ]

    with open(QC_OUT, "w", encoding="utf-8", newline="\n") as f:
        f.write("\t".join(qc_cols) + "\n")
        for r in results:
            f.write("\t".join(str(r[c]) for c in qc_cols) + "\n")

    with open(MANIFEST_OUT, "w", encoding="utf-8", newline="\n") as f:
        f.write("\t".join(manifest_cols) + "\n")
        for r in results:
            f.write("\t".join(str(r[c]) for c in manifest_cols) + "\n")

    with open(RUNTIME_OUT, "w", encoding="utf-8", newline="\n") as f:
        f.write("\t".join(runtime_cols) + "\n")
        for r in results:
            row = dict(r)
            row["elapsed_human"] = fmt_seconds(r["elapsed_seconds"])
            f.write("\t".join(str(row[c]) for c in runtime_cols) + "\n")

    n_total = len(results)
    n_ready = sum(1 for r in results if r["status"] == "MR_INPUT_READY")
    n_zero = sum(1 for r in results if r["status"] == "NO_MATCHED_INSTRUMENTS")
    n_review = sum(1 for r in results if r["status"] not in ("MR_INPUT_READY", "NO_MATCHED_INSTRUMENTS"))
    total_mr_rows = sum(r["n_mr_input_rows"] for r in results)
    total_missing = sum(r["n_missing_clumped_instruments"] for r in results)
    total_elapsed = sum(r["elapsed_seconds"] for r in results)

    with open(STATUS_OUT, "w", encoding="utf-8", newline="\n") as f:
        f.write("phase\tstatus\tkey_result\tnote\n")

        if n_ready == n_total:
            phase_status = "PASSED"
            key_result = f"{n_ready}/{n_total} MR input datasets ready"
            note = "All primary exposure-outcome pairs have at least one matched clumped instrument"
        elif n_ready > 0 and n_zero > 0 and n_review == 0:
            phase_status = "PASSED_WITH_ZERO_MATCH_WARNINGS"
            key_result = f"{n_ready}/{n_total} MR input datasets ready; {n_zero} have zero matched instruments"
            note = "Some outcome pairs have no matched clumped instruments, likely due to IOPcc coordinate subset limitations"
        else:
            phase_status = "REVIEW_NEEDED"
            key_result = f"{n_ready}/{n_total} MR input datasets ready"
            note = "Some datasets require review"

        f.write(f"Phase 3.3 extract clumped instruments from analysis-ready datasets\t{phase_status}\t{key_result}\t{note}\n")
        f.write(f"primary_exposure_count\tINFO\t{len(primary_exposures)}\t{','.join(primary_exposures)}\n")
        f.write(f"expected_pairwise_dataset_count\tINFO\t{len(primary_exposures) * len(OUTCOMES)}\tPrimary exposures x 3 outcomes\n")
        f.write(f"mr_input_ready_count\tINFO\t{n_ready}\tDatasets with status=MR_INPUT_READY\n")
        f.write(f"zero_matched_dataset_count\tINFO\t{n_zero}\tDatasets with status=NO_MATCHED_INSTRUMENTS\n")
        f.write(f"review_needed_dataset_count\tINFO\t{n_review}\tDatasets with missing core fields or other issues\n")
        f.write(f"total_mr_input_rows\tINFO\t{total_mr_rows}\tTotal extracted clumped instruments across all pairwise MR inputs\n")
        f.write(f"total_missing_clumped_instruments\tINFO\t{total_missing}\tClumped instruments not found in pairwise analysis-ready datasets\n")
        f.write(f"total_elapsed_runtime\tINFO\t{fmt_seconds(total_elapsed)}\tSum of pairwise extraction times\n")
        f.write("INSOMNIA_status\tDOCUMENTED\tEXCLUDED_FROM_PRIMARY_MR\tNo LD-clumped instruments in Phase 3.2\n")
        f.write("IOPcc_limitation\tDOCUMENTED\tcoordinate-resolvable subset only\tSome clumped instruments may not be available for IOPcc\n")

    with open(NOTES_MD, "w", encoding="utf-8", newline="\n") as f:
        f.write("# Phase 3.3 MR Input Extraction Notes\n\n")
        f.write("## Status\n\n")
        f.write("Clumped instruments were extracted from Phase 2.7 analysis-ready pairwise harmonized datasets.\n\n")
        f.write("## Primary exposure set\n\n")
        for e in primary_exposures:
            f.write(f"- {e}\n")
        f.write("\nINSOMNIA is excluded from primary MR because no LD-clumped instruments were generated in Phase 3.2.\n\n")
        f.write("## Outcomes\n\n")
        for o in OUTCOMES:
            f.write(f"- {o}\n")
        f.write("\n")
        f.write("## MR input rule\n\n")
        f.write("Only LD-clumped exposure instruments were extracted. All rows come from Phase 2.7 analysis-ready datasets, which already retain only include_in_main=YES harmonized rows.\n\n")
        f.write("## Important caveat\n\n")
        f.write("IOPcc uses only coordinate-resolvable allele-aware variants. Therefore, some clumped exposure instruments may not be available for IOPcc.\n\n")
        f.write("## Next phase\n\n")
        f.write("Phase 3.4 should finalize MR method eligibility per exposure-outcome pair based on instrument counts.\n")

def main():
    print("===== Phase 3.3 extract clumped instruments from analysis-ready datasets =====")
    print(f"Start: {now_str()}")

    primary_exposures = read_primary_exposures()
    print("Primary exposures:", ", ".join(primary_exposures))
    print("Output directory:", OUTDIR)

    total_t0 = time.time()

    results = []

    for exposure_id in primary_exposures:
        clumped_file, instruments, order = load_clumped_instruments(exposure_id)

        print()
        print(f"### Exposure: {exposure_id}")
        print(f"Clumped file: {clumped_file}")
        print(f"Clumped instruments: {len(instruments)}")

        if len(instruments) == 0:
            print(f"Skipping {exposure_id}: no clumped instruments")
            continue

        for outcome_suffix in OUTCOMES:
            results.append(extract_for_pair(exposure_id, outcome_suffix, instruments, order))

    write_outputs(results, primary_exposures)

    total_elapsed = time.time() - total_t0

    print()
    print("===== Phase 3.3 completed =====")
    print(f"End: {now_str()}")
    print(f"Total elapsed: {fmt_seconds(total_elapsed)}")
    print(f"Wrote QC summary: {QC_OUT}")
    print(f"Wrote manifest: {MANIFEST_OUT}")
    print(f"Wrote runtime log: {RUNTIME_OUT}")
    print(f"Wrote status: {STATUS_OUT}")
    print(f"Wrote notes: {NOTES_MD}")

if __name__ == "__main__":
    main()
