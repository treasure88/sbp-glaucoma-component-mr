#!/usr/bin/env python3
import csv
import gzip
import math
import os
import time

BASE = "../../25_external_outcome_triangulation"
OUTDIR = "../../28_external_outcome_triangulation_inputs"
HARM_DIR = os.path.join(OUTDIR, "pairwise_harmonized")
MR_DIR = os.path.join(OUTDIR, "mr_input")
MISS_DIR = os.path.join(OUTDIR, "missing_instruments")
os.makedirs(HARM_DIR, exist_ok=True)
os.makedirs(MR_DIR, exist_ok=True)
os.makedirs(MISS_DIR, exist_ok=True)

EXPOSURES = {
    "SBP": {
        "dataset_id": "ieu-b-38",
        "trait_name": "Systolic blood pressure",
        "file": "../../13_instrument_selection/clumped_instruments/SBP__ieu-b-38.clumped_instruments.tsv.gz",
        "note": "SBP Phase 3.2 LD-clumped instrument",
    },
    "ART_STIFFNESS": {
        "dataset_id": "ukb-b-11971",
        "trait_name": "Pulse wave Arterial Stiffness index",
        "file": "../../23_vascular_exposure_panel/phase4_5_instrument_selection/clumped_instruments/ART_STIFFNESS__ukb-b-11971.clumped_instruments.tsv.gz",
        "note": "ART_STIFFNESS Phase 4.5 LD-clumped instrument",
    },
}

OUTCOMES = {
    "RNFL": {
        "dataset_id": "GCST90014266",
        "trait_name": "Retinal nerve fibre layer thickness",
        "file": "../../27_external_outcome_standardized/RNFL__GCST90014266.standardized.tsv.gz",
    },
    "GCIPL": {
        "dataset_id": "GCST90014267",
        "trait_name": "Ganglion cell inner plexiform layer thickness",
        "file": "../../27_external_outcome_standardized/GCIPL__GCST90014267.standardized.tsv.gz",
    },
}

STATUS_OUT = os.path.join(OUTDIR, "phase5_5_status.tsv")
QC_OUT = os.path.join(OUTDIR, "phase5_5_external_harmonization_qc_summary.tsv")
MANIFEST_OUT = os.path.join(OUTDIR, "phase5_5_external_mr_input_manifest.tsv")
RUNTIME_OUT = os.path.join(OUTDIR, "phase5_5_runtime_log.tsv")

HARM_HEADER = [
    "exposure_id", "exposure_dataset_id", "outcome_id", "outcome_dataset_id",
    "SNP", "exposure_chr", "exposure_pos", "outcome_chr", "outcome_pos",
    "exposure_effect_allele", "exposure_other_allele",
    "outcome_effect_allele_raw", "outcome_other_allele_raw",
    "beta_exposure", "se_exposure", "pval_exposure", "minus_log10_p_exposure",
    "eaf_exposure", "n_exposure",
    "beta_outcome_raw", "beta_outcome_aligned", "se_outcome", "pval_outcome",
    "eaf_outcome", "n_outcome",
    "harmonization_action", "match_mode", "match_key",
    "include_in_external_mr", "exclusion_reason", "F_stat", "instrument_note",
]

MR_HEADER = [
    "exposure_id", "outcome_id", "SNP", "chr", "pos",
    "effect_allele", "other_allele",
    "beta_exposure", "se_exposure", "pval_exposure", "minus_log10_p_exposure",
    "eaf_exposure", "n_exposure",
    "beta_outcome", "se_outcome", "pval_outcome", "pval_for_log10_outcome",
    "harmonization_action", "match_mode", "match_key",
    "outcome_source_dataset", "F_stat", "instrument_note",
]

MISSING_HEADER = [
    "exposure_id", "outcome_id", "SNP", "chr", "pos",
    "effect_allele", "other_allele", "reason",
]

def clean(x):
    if x is None:
        return "NA"
    x = str(x).strip()
    return x if x else "NA"

def is_missing(x):
    return x is None or str(x).strip() in ("", ".", "NA", "na", "NaN")

def safe_float(x):
    try:
        if is_missing(x):
            return None
        return float(x)
    except Exception:
        return None

def pick(row, names, default="NA"):
    for n in names:
        if n in row and not is_missing(row[n]):
            return clean(row[n])
    return default

def upper_allele(x):
    return clean(x).upper()

def complement_allele(a):
    table = str.maketrans("ACGT", "TGCA")
    a = upper_allele(a)
    if any(ch not in "ACGT" for ch in a):
        return None
    return a.translate(table)

def f_stat(beta, se):
    b = safe_float(beta)
    s = safe_float(se)
    if b is None or s is None or s == 0:
        return "NA"
    return f"{(b / s) ** 2:.12g}"

def read_exposure_instruments(exposure_id, cfg):
    rows = []
    path = cfg["file"]
    with gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            snp = pick(row, ["SNP", "rsid", "variant_id"])
            if is_missing(snp):
                continue
            beta = pick(row, ["beta", "beta_exposure"])
            se = pick(row, ["se", "se_exposure"])
            pval = pick(row, ["pval", "pval_exposure", "pval_for_clumping"])
            mlp = pick(row, ["minus_log10_p", "minus_log10_p_exposure"])
            if is_missing(mlp):
                p = safe_float(pval)
                mlp = f"{-math.log10(p):.12g}" if p is not None and p > 0 else "NA"

            out = {
                "exposure_id": exposure_id,
                "exposure_dataset_id": cfg["dataset_id"],
                "trait_name": cfg["trait_name"],
                "chr": pick(row, ["chr", "chromosome", "#CHROM"]),
                "pos": pick(row, ["pos", "base_pair_location", "position"]),
                "SNP": snp,
                "variant_id": pick(row, ["variant_id"]),
                "effect_allele": upper_allele(pick(row, ["effect_allele", "ea", "EA"])),
                "other_allele": upper_allele(pick(row, ["other_allele", "oa", "OA"])),
                "beta": beta,
                "se": se,
                "pval": pval,
                "minus_log10_p": mlp,
                "eaf": pick(row, ["eaf", "effect_allele_frequency"]),
                "n": pick(row, ["n", "sample_size", "N"]),
                "F_stat": f_stat(beta, se),
                "instrument_note": cfg["note"],
            }
            rows.append(out)
    return rows

def collect_target_snps(all_exposure_rows):
    s = set()
    for rows in all_exposure_rows.values():
        for r in rows:
            s.add(r["SNP"])
    return s

def read_outcome_matches(outcome_id, cfg, target_snps):
    matches = {}
    path = cfg["file"]

    with gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            snp = pick(row, ["SNP", "rsid", "variant_id"])
            if snp not in target_snps:
                continue
            matches.setdefault(snp, []).append({
                "outcome_id": outcome_id,
                "outcome_dataset_id": cfg["dataset_id"],
                "trait_name": cfg["trait_name"],
                "chr": pick(row, ["chr"]),
                "pos": pick(row, ["pos"]),
                "SNP": snp,
                "variant_id": pick(row, ["variant_id"]),
                "effect_allele": upper_allele(pick(row, ["effect_allele"])),
                "other_allele": upper_allele(pick(row, ["other_allele"])),
                "beta": pick(row, ["beta"]),
                "se": pick(row, ["se"]),
                "pval": pick(row, ["pval"]),
                "eaf": pick(row, ["eaf"]),
                "n": pick(row, ["n"]),
                "source_dataset": pick(row, ["source_dataset"]),
            })

    return matches

def harmonize_one(exp, outcome_candidates):
    if not outcome_candidates:
        return None, "NO_OUTCOME_SNP_MATCH"

    ea = exp["effect_allele"]
    oa = exp["other_allele"]

    best_mismatch_reason = "ALLELE_MISMATCH"

    for out in outcome_candidates:
        oe = out["effect_allele"]
        oo = out["other_allele"]

        beta_o = safe_float(out["beta"])
        if beta_o is None:
            best_mismatch_reason = "MISSING_OUTCOME_BETA"
            continue

        if oe == ea and oo == oa:
            return (out, beta_o, "aligned"), None

        if oe == oa and oo == ea:
            return (out, -beta_o, "flipped"), None

        ce = complement_allele(oe)
        co = complement_allele(oo)

        if ce is not None and co is not None:
            if ce == ea and co == oa:
                return (out, beta_o, "aligned_complement"), None
            if ce == oa and co == ea:
                return (out, -beta_o, "flipped_complement"), None

    return None, best_mismatch_reason

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
    start_all = time.time()

    exposure_rows = {}
    for exp_id, cfg in EXPOSURES.items():
        exposure_rows[exp_id] = read_exposure_instruments(exp_id, cfg)

    target_snps = collect_target_snps(exposure_rows)

    outcome_matches = {}
    outcome_runtime = {}
    for out_id, cfg in OUTCOMES.items():
        t0 = time.time()
        print(f"Scanning external outcome {out_id}: {cfg['file']}", flush=True)
        outcome_matches[out_id] = read_outcome_matches(out_id, cfg, target_snps)
        outcome_runtime[out_id] = time.time() - t0

    qc_rows = []
    manifest_rows = []
    runtime_rows = []

    for exp_id, exp_list in exposure_rows.items():
        for out_id, out_cfg in OUTCOMES.items():
            t0 = time.time()

            harmonized_rows = []
            mr_rows = []
            missing_rows = []

            n_snp_match = 0
            n_duplicate_outcome_records = 0
            n_allele_mismatch = 0
            n_missing = 0

            matches = outcome_matches[out_id]

            for exp in exp_list:
                snp = exp["SNP"]
                candidates = matches.get(snp, [])
                if candidates:
                    n_snp_match += 1
                    if len(candidates) > 1:
                        n_duplicate_outcome_records += len(candidates) - 1

                h, reason = harmonize_one(exp, candidates)

                if h is None:
                    if reason == "NO_OUTCOME_SNP_MATCH":
                        n_missing += 1
                    else:
                        n_allele_mismatch += 1

                    missing_rows.append({
                        "exposure_id": exp_id,
                        "outcome_id": out_id,
                        "SNP": snp,
                        "chr": exp["chr"],
                        "pos": exp["pos"],
                        "effect_allele": exp["effect_allele"],
                        "other_allele": exp["other_allele"],
                        "reason": reason,
                    })

                    harmonized_rows.append({
                        "exposure_id": exp_id,
                        "exposure_dataset_id": EXPOSURES[exp_id]["dataset_id"],
                        "outcome_id": out_id,
                        "outcome_dataset_id": out_cfg["dataset_id"],
                        "SNP": snp,
                        "exposure_chr": exp["chr"],
                        "exposure_pos": exp["pos"],
                        "outcome_chr": "NA",
                        "outcome_pos": "NA",
                        "exposure_effect_allele": exp["effect_allele"],
                        "exposure_other_allele": exp["other_allele"],
                        "outcome_effect_allele_raw": "NA",
                        "outcome_other_allele_raw": "NA",
                        "beta_exposure": exp["beta"],
                        "se_exposure": exp["se"],
                        "pval_exposure": exp["pval"],
                        "minus_log10_p_exposure": exp["minus_log10_p"],
                        "eaf_exposure": exp["eaf"],
                        "n_exposure": exp["n"],
                        "beta_outcome_raw": "NA",
                        "beta_outcome_aligned": "NA",
                        "se_outcome": "NA",
                        "pval_outcome": "NA",
                        "eaf_outcome": "NA",
                        "n_outcome": "NA",
                        "harmonization_action": "not_harmonized",
                        "match_mode": "none",
                        "match_key": snp,
                        "include_in_external_mr": "NO",
                        "exclusion_reason": reason,
                        "F_stat": exp["F_stat"],
                        "instrument_note": exp["instrument_note"],
                    })
                    continue

                out, beta_aligned, action = h
                include = "YES"

                row_h = {
                    "exposure_id": exp_id,
                    "exposure_dataset_id": EXPOSURES[exp_id]["dataset_id"],
                    "outcome_id": out_id,
                    "outcome_dataset_id": out_cfg["dataset_id"],
                    "SNP": snp,
                    "exposure_chr": exp["chr"],
                    "exposure_pos": exp["pos"],
                    "outcome_chr": out["chr"],
                    "outcome_pos": out["pos"],
                    "exposure_effect_allele": exp["effect_allele"],
                    "exposure_other_allele": exp["other_allele"],
                    "outcome_effect_allele_raw": out["effect_allele"],
                    "outcome_other_allele_raw": out["other_allele"],
                    "beta_exposure": exp["beta"],
                    "se_exposure": exp["se"],
                    "pval_exposure": exp["pval"],
                    "minus_log10_p_exposure": exp["minus_log10_p"],
                    "eaf_exposure": exp["eaf"],
                    "n_exposure": exp["n"],
                    "beta_outcome_raw": out["beta"],
                    "beta_outcome_aligned": f"{beta_aligned:.12g}",
                    "se_outcome": out["se"],
                    "pval_outcome": out["pval"],
                    "eaf_outcome": out["eaf"],
                    "n_outcome": out["n"],
                    "harmonization_action": action,
                    "match_mode": "rsid",
                    "match_key": snp,
                    "include_in_external_mr": include,
                    "exclusion_reason": "none",
                    "F_stat": exp["F_stat"],
                    "instrument_note": exp["instrument_note"],
                }
                harmonized_rows.append(row_h)

                mr_rows.append({
                    "exposure_id": exp_id,
                    "outcome_id": out_id,
                    "SNP": snp,
                    "chr": exp["chr"],
                    "pos": exp["pos"],
                    "effect_allele": exp["effect_allele"],
                    "other_allele": exp["other_allele"],
                    "beta_exposure": exp["beta"],
                    "se_exposure": exp["se"],
                    "pval_exposure": exp["pval"],
                    "minus_log10_p_exposure": exp["minus_log10_p"],
                    "eaf_exposure": exp["eaf"],
                    "n_exposure": exp["n"],
                    "beta_outcome": f"{beta_aligned:.12g}",
                    "se_outcome": out["se"],
                    "pval_outcome": out["pval"],
                    "pval_for_log10_outcome": out["pval"],
                    "harmonization_action": action,
                    "match_mode": "rsid",
                    "match_key": snp,
                    "outcome_source_dataset": out_cfg["dataset_id"],
                    "F_stat": exp["F_stat"],
                    "instrument_note": exp["instrument_note"],
                })

            harm_file = os.path.join(HARM_DIR, f"{exp_id}__{out_id}.targeted_external_harmonized.tsv.gz")
            mr_file = os.path.join(MR_DIR, f"{exp_id}__{out_id}.external_mr_input.tsv.gz")
            miss_file = os.path.join(MISS_DIR, f"{exp_id}__{out_id}.missing_external_instruments.tsv")

            write_gz_tsv(harm_file, HARM_HEADER, harmonized_rows)
            write_gz_tsv(mr_file, MR_HEADER, mr_rows)
            write_tsv(miss_file, MISSING_HEADER, missing_rows)

            elapsed = time.time() - t0

            qc_rows.append({
                "exposure_id": exp_id,
                "outcome_id": out_id,
                "n_clumped_instruments": str(len(exp_list)),
                "n_snp_matched_in_outcome": str(n_snp_match),
                "n_harmonized_for_mr": str(len(mr_rows)),
                "n_missing_from_outcome": str(n_missing),
                "n_allele_mismatch_or_unusable": str(n_allele_mismatch),
                "n_duplicate_outcome_records_extra": str(n_duplicate_outcome_records),
                "harmonized_file": harm_file,
                "mr_input_file": mr_file,
                "missing_file": miss_file,
                "status": "MR_INPUT_READY" if len(mr_rows) > 0 else "NO_MATCHED_INSTRUMENTS",
                "note": "External outcome targeted harmonization against RNFL/GCIPL",
            })

            manifest_rows.append({
                "exposure_id": exp_id,
                "outcome_id": out_id,
                "mr_input_file": mr_file,
                "n_mr_input_rows": str(len(mr_rows)),
                "status": "MR_INPUT_READY" if len(mr_rows) > 0 else "NO_MATCHED_INSTRUMENTS",
            })

            runtime_rows.append({
                "exposure_id": exp_id,
                "outcome_id": out_id,
                "elapsed_seconds": f"{elapsed:.3f}",
                "elapsed_human": f"{elapsed:.1f}s",
            })

    write_tsv(QC_OUT, [
        "exposure_id", "outcome_id", "n_clumped_instruments",
        "n_snp_matched_in_outcome", "n_harmonized_for_mr",
        "n_missing_from_outcome", "n_allele_mismatch_or_unusable",
        "n_duplicate_outcome_records_extra",
        "harmonized_file", "mr_input_file", "missing_file",
        "status", "note",
    ], qc_rows)

    write_tsv(MANIFEST_OUT, [
        "exposure_id", "outcome_id", "mr_input_file", "n_mr_input_rows", "status"
    ], manifest_rows)

    write_tsv(RUNTIME_OUT, [
        "exposure_id", "outcome_id", "elapsed_seconds", "elapsed_human"
    ], runtime_rows)

    ready = sum(1 for r in qc_rows if r["status"] == "MR_INPUT_READY")
    total_rows = sum(int(r["n_harmonized_for_mr"]) for r in qc_rows)
    elapsed_all = time.time() - start_all

    with open(STATUS_OUT, "w", encoding="utf-8", newline="") as f:
        f.write("phase\tstatus\tkey_result\tnote\n")
        f.write(f"Phase 5.5 targeted external outcome harmonization\tPASSED\t{ready}/4 pairwise external MR inputs ready; total_rows={total_rows}\tSBP and ART_STIFFNESS instruments harmonized against RNFL and GCIPL\n")
        f.write("external_outcomes\tINFO\tRNFL;GCIPL\tNeuroretinal endophenotype triangulation\n")
        f.write("primary_exposure\tINFO\tSBP\tPrimary vascular/hemodynamic signal\n")
        f.write("secondary_exposure\tINFO\tART_STIFFNESS\tExploratory vascular-stiffness support\n")
        f.write(f"runtime\tINFO\t{elapsed_all:.3f}s\tPhase 5.5 runtime\n")

    print("===== Phase 5.5 completed =====")
    print("Wrote:", STATUS_OUT)
    print("Wrote:", QC_OUT)
    print("Wrote:", MANIFEST_OUT)

if __name__ == "__main__":
    main()
