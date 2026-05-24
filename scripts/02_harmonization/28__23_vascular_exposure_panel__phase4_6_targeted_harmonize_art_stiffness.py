#!/usr/bin/env python3
import csv
import gzip
import math
import os
import time

BASE = "../../23_vascular_exposure_panel"
CLUMPED = os.path.join(
    BASE,
    "phase4_5_instrument_selection",
    "clumped_instruments",
    "ART_STIFFNESS__ukb-b-11971.clumped_instruments.tsv.gz"
)

OUTCOME_INDEXES = {
    "GBS_nonIOPcomponent": "../../05_harmonization_planning/outcome_indexes/GBS_nonIOP_by_rsid.tsv.gz",
    "GBS_IOPcomponent": "../../05_harmonization_planning/outcome_indexes/GBS_IOP_by_rsid.tsv.gz",
}

OUTDIR = os.path.join(BASE, "phase4_6_targeted_harmonization")
PAIRWISE_DIR = os.path.join(OUTDIR, "pairwise_harmonized")
MR_DIR = os.path.join(OUTDIR, "mr_input")
MISS_DIR = os.path.join(OUTDIR, "missing_instruments")

for d in [OUTDIR, PAIRWISE_DIR, MR_DIR, MISS_DIR]:
    os.makedirs(d, exist_ok=True)

STATUS_OUT = os.path.join(OUTDIR, "phase4_6_status.tsv")
QC_OUT = os.path.join(OUTDIR, "phase4_6_targeted_harmonization_qc_summary.tsv")
MANIFEST_OUT = os.path.join(OUTDIR, "phase4_6_targeted_harmonization_manifest.tsv")
RUNTIME_OUT = os.path.join(OUTDIR, "phase4_6_runtime_log.tsv")

HARMONIZED_HEADER = [
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
    "beta_outcome_original",
    "beta_outcome_harmonized",
    "se_outcome",
    "pval_outcome",
    "pval_for_log10_outcome",
    "harmonization_action",
    "include_in_main",
    "outcome_source_dataset",
]

MR_HEADER = [
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
    "outcome_source_dataset",
    "F_stat",
    "instrument_note",
]

MISSING_HEADER = [
    "exposure_id",
    "outcome_id",
    "SNP",
    "chr",
    "pos",
    "variant_id",
    "reason",
]

def safe_float(x):
    try:
        if x in ("", ".", "NA", None):
            return None
        return float(x)
    except Exception:
        return None

def fmt_float(x):
    if x is None:
        return "NA"
    return f"{x:.12g}"

def is_palindromic(a, b):
    s = {a.upper(), b.upper()}
    return s == {"A", "T"} or s == {"C", "G"}

def load_clumped():
    rows = []
    with gzip.open(CLUMPED, "rt", encoding="utf-8", errors="replace", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            rows.append(row)
    return rows

def load_outcome_matches(index_file, wanted_snps):
    matches = {}
    with gzip.open(index_file, "rt", encoding="utf-8", errors="replace", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            snp = row.get("SNP", "")
            key = row.get("match_key", "")
            if snp in wanted_snps or key in wanted_snps:
                # Prefer SNP as key, because instruments are rsID based.
                matches[snp] = row
                matches[key] = row
    return matches

def harmonize_one(exp, out, outcome_id):
    exp_ea = exp["effect_allele"]
    exp_oa = exp["other_allele"]
    out_ea = out["effect_allele"]
    out_oa = out["other_allele"]

    beta_out = safe_float(out["beta"])

    pal = is_palindromic(exp_ea, exp_oa)

    if exp_ea == out_ea and exp_oa == out_oa:
        action = "aligned"
        beta_h = beta_out
    elif exp_ea == out_oa and exp_oa == out_ea:
        action = "flipped"
        beta_h = -beta_out if beta_out is not None else None
    else:
        action = "allele_mismatch_excluded"
        beta_h = None

    if pal and action == "aligned":
        action = "palindromic_aligned_excluded"
    elif pal and action == "flipped":
        action = "palindromic_flipped_excluded"

    include = "YES" if action in ("aligned", "flipped") else "NO"

    f_stat = "NA"
    bx = safe_float(exp["beta"])
    bxse = safe_float(exp["se"])
    if bx is not None and bxse is not None and bxse != 0:
        f_stat = fmt_float((bx / bxse) ** 2)

    harmonized = {
        "exposure_id": exp["exposure_id"],
        "outcome_id": outcome_id,
        "match_mode": "rsid",
        "match_key": exp["SNP"],
        "chr": exp["chr"],
        "pos": exp["pos"],
        "SNP": exp["SNP"],
        "exposure_effect_allele": exp_ea,
        "exposure_other_allele": exp_oa,
        "outcome_effect_allele_original": out_ea,
        "outcome_other_allele_original": out_oa,
        "beta_exposure": exp["beta"],
        "se_exposure": exp["se"],
        "pval_exposure": exp["pval"],
        "minus_log10_p_exposure": exp["minus_log10_p"],
        "eaf_exposure": exp["eaf"],
        "n_exposure": exp["n"],
        "beta_outcome_original": out["beta"],
        "beta_outcome_harmonized": fmt_float(beta_h),
        "se_outcome": out["se"],
        "pval_outcome": out["pval"],
        "pval_for_log10_outcome": out["pval_for_log10"],
        "harmonization_action": action,
        "include_in_main": include,
        "outcome_source_dataset": out["source_dataset"],
    }

    mr = None
    if include == "YES":
        mr = {
            "exposure_id": exp["exposure_id"],
            "outcome_id": outcome_id,
            "SNP": exp["SNP"],
            "chr": exp["chr"],
            "pos": exp["pos"],
            "effect_allele": exp_ea,
            "other_allele": exp_oa,
            "beta_exposure": exp["beta"],
            "se_exposure": exp["se"],
            "pval_exposure": exp["pval"],
            "minus_log10_p_exposure": exp["minus_log10_p"],
            "eaf_exposure": exp["eaf"],
            "n_exposure": exp["n"],
            "beta_outcome": fmt_float(beta_h),
            "se_outcome": out["se"],
            "pval_outcome": out["pval"],
            "pval_for_log10_outcome": out["pval_for_log10"],
            "harmonization_action": action,
            "match_mode": "rsid",
            "match_key": exp["SNP"],
            "outcome_source_dataset": out["source_dataset"],
            "F_stat": f_stat,
            "instrument_note": "ART_STIFFNESS strict p<5e-8 LD-clumped instrument",
        }

    return harmonized, mr

def write_gz_tsv(path, header, rows):
    with gzip.open(path, "wt", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=header, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

def write_tsv(path, header, rows):
    with open(path, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=header, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

def main():
    total_start = time.time()

    clumped = load_clumped()
    wanted_snps = set(row["SNP"] for row in clumped)

    qc_rows = []
    manifest_rows = []
    runtime_rows = []

    for outcome_id, index_file in OUTCOME_INDEXES.items():
        t0 = time.time()
        print(f"Targeted harmonization: ART_STIFFNESS -> {outcome_id}", flush=True)

        matches = load_outcome_matches(index_file, wanted_snps)

        harmonized_rows = []
        mr_rows = []
        missing_rows = []

        action_counts = {}
        n_matched = 0

        for exp in clumped:
            snp = exp["SNP"]
            out = matches.get(snp)

            if out is None:
                missing_rows.append({
                    "exposure_id": exp["exposure_id"],
                    "outcome_id": outcome_id,
                    "SNP": snp,
                    "chr": exp["chr"],
                    "pos": exp["pos"],
                    "variant_id": exp["variant_id"],
                    "reason": "SNP_NOT_FOUND_IN_OUTCOME_RSID_INDEX",
                })
                continue

            n_matched += 1
            h, mr = harmonize_one(exp, out, outcome_id)
            harmonized_rows.append(h)
            action = h["harmonization_action"]
            action_counts[action] = action_counts.get(action, 0) + 1
            if mr is not None:
                mr_rows.append(mr)

        pair_file = os.path.join(PAIRWISE_DIR, f"ART_STIFFNESS__{outcome_id}.targeted_harmonized.tsv.gz")
        mr_file = os.path.join(MR_DIR, f"ART_STIFFNESS__{outcome_id}.mr_input.tsv.gz")
        missing_file = os.path.join(MISS_DIR, f"ART_STIFFNESS__{outcome_id}.missing_clumped_instruments.tsv")

        write_gz_tsv(pair_file, HARMONIZED_HEADER, harmonized_rows)
        write_gz_tsv(mr_file, MR_HEADER, mr_rows)
        write_tsv(missing_file, MISSING_HEADER, missing_rows)

        elapsed = time.time() - t0

        action_summary = ";".join(f"{k}={v}" for k, v in sorted(action_counts.items())) if action_counts else "none"

        status = "MR_INPUT_READY" if len(mr_rows) > 0 else "NO_MR_INPUT_ROWS"

        qc_rows.append({
            "exposure_id": "ART_STIFFNESS",
            "outcome_id": outcome_id,
            "n_clumped_instruments": str(len(clumped)),
            "n_matched_in_outcome": str(n_matched),
            "n_missing_in_outcome": str(len(missing_rows)),
            "n_harmonized_rows": str(len(harmonized_rows)),
            "n_mr_input_rows": str(len(mr_rows)),
            "harmonization_action_summary": action_summary,
            "status": status,
            "note": "Palindromic instruments excluded from main analysis; only aligned/flipped non-palindromic instruments retained",
        })

        manifest_rows.append({
            "exposure_id": "ART_STIFFNESS",
            "outcome_id": outcome_id,
            "harmonized_file": pair_file,
            "mr_input_file": mr_file,
            "missing_instruments_file": missing_file,
            "status": status,
        })

        runtime_rows.append({
            "exposure_id": "ART_STIFFNESS",
            "outcome_id": outcome_id,
            "elapsed_seconds": f"{elapsed:.3f}",
            "elapsed_human": f"{elapsed:.1f}s",
        })

    write_tsv(QC_OUT, list(qc_rows[0].keys()), qc_rows)
    write_tsv(MANIFEST_OUT, list(manifest_rows[0].keys()), manifest_rows)
    write_tsv(RUNTIME_OUT, list(runtime_rows[0].keys()), runtime_rows)

    ready_count = sum(1 for r in qc_rows if r["status"] == "MR_INPUT_READY")
    total_mr_rows = sum(int(r["n_mr_input_rows"]) for r in qc_rows)
    elapsed_total = time.time() - total_start

    with open(STATUS_OUT, "w", encoding="utf-8", newline="") as f:
        f.write("phase\tstatus\tkey_result\tnote\n")
        phase_status = "PASSED" if ready_count > 0 else "REVIEW_NEEDED_NO_MR_INPUT"
        f.write(
            f"Phase 4.6 targeted ART_STIFFNESS harmonization\t{phase_status}\t"
            f"ready_pairs={ready_count}/2;total_mr_input_rows={total_mr_rows}\t"
            f"Targeted rsID lookup against GBS outcome indexes completed\n"
        )
        f.write(
            "HYPERTENSION_status\tDOCUMENTED\tEXCLUDED_NO_GWS_INSTRUMENTS\t"
            "No strict p<5e-8 instruments from Phase 4.5\n"
        )
        f.write(
            "ART_STIFFNESS_instrument_count_caveat\tDOCUMENTED\t3 clumped instruments before outcome matching\t"
            "Downstream MR should be exploratory and low-power\n"
        )
        f.write(
            f"runtime\tINFO\t{elapsed_total:.3f}s\tPhase 4.6 runtime\n"
        )

    print("===== Phase 4.6 completed =====")
    print("Wrote:", STATUS_OUT)
    print("Wrote:", QC_OUT)
    print("Wrote:", MANIFEST_OUT)
    print("Wrote:", RUNTIME_OUT)

if __name__ == "__main__":
    main()
