#!/usr/bin/env python3
import csv
import gzip
import math
import os
import time

OUTDIR = "../../32_poag_external_validation_inputs"
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

OUTCOME = {
    "outcome_id": "POAG",
    "dataset_id": "GCST90011766",
    "trait_name": "Primary open-angle glaucoma",
    "file": "../../27_external_outcome_standardized/POAG__GCST90011766.standardized.tsv.gz",
    "scale": "log_odds",
}

STATUS_OUT = os.path.join(OUTDIR, "phase5_8D_status.tsv")
QC_OUT = os.path.join(OUTDIR, "phase5_8D_POAG_harmonization_qc_summary.tsv")
MANIFEST_OUT = os.path.join(OUTDIR, "phase5_8D_POAG_mr_input_manifest.tsv")
RUNTIME_OUT = os.path.join(OUTDIR, "phase5_8D_runtime_log.tsv")

HARM_HEADER = [
    "exposure_id", "exposure_dataset_id", "outcome_id", "outcome_dataset_id",
    "SNP", "exposure_chr", "exposure_pos", "outcome_chr", "outcome_pos",
    "exposure_effect_allele", "exposure_other_allele",
    "outcome_effect_allele_raw", "outcome_other_allele_raw",
    "beta_exposure", "se_exposure", "pval_exposure", "minus_log10_p_exposure",
    "eaf_exposure", "n_exposure",
    "beta_outcome_raw", "beta_outcome_aligned", "se_outcome", "pval_outcome",
    "eaf_outcome", "n_outcome", "outcome_scale",
    "harmonization_action", "match_mode", "match_key",
    "include_in_mr", "exclusion_reason", "F_stat", "instrument_note",
]

MR_HEADER = [
    "exposure_id", "outcome_id", "SNP", "chr", "pos",
    "effect_allele", "other_allele",
    "beta_exposure", "se_exposure", "pval_exposure", "minus_log10_p_exposure",
    "eaf_exposure", "n_exposure",
    "beta_outcome", "se_outcome", "pval_outcome", "pval_for_log10_outcome",
    "harmonization_action", "match_mode", "match_key",
    "outcome_source_dataset", "outcome_scale", "F_stat", "instrument_note",
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
    return x is None or str(x).strip() in ("", ".", "NA", "na", "NaN", "nan")

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
    a = upper_allele(a)
    if any(ch not in "ACGT" for ch in a):
        return None
    return a.translate(str.maketrans("ACGT", "TGCA"))

def f_stat(beta, se):
    b = safe_float(beta)
    s = safe_float(se)
    if b is None or s is None or s == 0:
        return "NA"
    return f"{(b / s) ** 2:.12g}"

def read_exposure_instruments(exposure_id, cfg):
    rows = []
    with gzip.open(cfg["file"], "rt", encoding="utf-8", errors="replace", newline="") as f:
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

            rows.append({
                "exposure_id": exposure_id,
                "exposure_dataset_id": cfg["dataset_id"],
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
            })
    return rows

def collect_target_snps(exposure_rows):
    s = set()
    for rows in exposure_rows.values():
        for r in rows:
            s.add(r["SNP"])
    return s

def read_poag_matches(target_snps):
    matches = {}

    with gzip.open(OUTCOME["file"], "rt", encoding="utf-8", errors="replace", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            snp = pick(row, ["SNP", "rsid", "variant_id"])
            if snp not in target_snps:
                continue

            matches.setdefault(snp, []).append({
                "outcome_id": OUTCOME["outcome_id"],
                "outcome_dataset_id": OUTCOME["dataset_id"],
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
                "outcome_scale": pick(row, ["outcome_scale"]),
            })

    return matches

def harmonize_one(exp, candidates):
    if not candidates:
        return None, "NO_OUTCOME_SNP_MATCH"

    ea = exp["effect_allele"]
    oa = exp["other_allele"]

    best_reason = "ALLELE_MISMATCH"

    for out in candidates:
        oe = out["effect_allele"]
        oo = out["other_allele"]

        beta_o = safe_float(out["beta"])
        if beta_o is None:
            best_reason = "MISSING_OUTCOME_BETA"
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

    return None, best_reason

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
    start = time.time()

    exposure_rows = {}
    for exp_id, cfg in EXPOSURES.items():
        exposure_rows[exp_id] = read_exposure_instruments(exp_id, cfg)

    target_snps = collect_target_snps(exposure_rows)

    print(f"Scanning POAG standardized file for {len(target_snps)} target SNPs", flush=True)
    scan_start = time.time()
    poag_matches = read_poag_matches(target_snps)
    scan_elapsed = time.time() - scan_start

    qc_rows = []
    manifest_rows = []
    runtime_rows = []

    for exp_id, exp_list in exposure_rows.items():
        t0 = time.time()

        harmonized_rows = []
        mr_rows = []
        missing_rows = []

        n_snp_match = 0
        n_duplicate_extra = 0
        n_missing = 0
        n_allele_mismatch = 0
        n_flipped = 0
        n_aligned = 0

        for exp in exp_list:
            snp = exp["SNP"]
            candidates = poag_matches.get(snp, [])

            if candidates:
                n_snp_match += 1
                if len(candidates) > 1:
                    n_duplicate_extra += len(candidates) - 1

            h, reason = harmonize_one(exp, candidates)

            if h is None:
                if reason == "NO_OUTCOME_SNP_MATCH":
                    n_missing += 1
                else:
                    n_allele_mismatch += 1

                missing_rows.append({
                    "exposure_id": exp_id,
                    "outcome_id": "POAG",
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
                    "outcome_id": "POAG",
                    "outcome_dataset_id": OUTCOME["dataset_id"],
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
                    "outcome_scale": OUTCOME["scale"],
                    "harmonization_action": "not_harmonized",
                    "match_mode": "none",
                    "match_key": snp,
                    "include_in_mr": "NO",
                    "exclusion_reason": reason,
                    "F_stat": exp["F_stat"],
                    "instrument_note": exp["instrument_note"],
                })
                continue

            out, beta_aligned, action = h

            if action.startswith("flipped"):
                n_flipped += 1
            else:
                n_aligned += 1

            hrow = {
                "exposure_id": exp_id,
                "exposure_dataset_id": EXPOSURES[exp_id]["dataset_id"],
                "outcome_id": "POAG",
                "outcome_dataset_id": OUTCOME["dataset_id"],
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
                "outcome_scale": OUTCOME["scale"],
                "harmonization_action": action,
                "match_mode": "rsid",
                "match_key": snp,
                "include_in_mr": "YES",
                "exclusion_reason": "none",
                "F_stat": exp["F_stat"],
                "instrument_note": exp["instrument_note"],
            }
            harmonized_rows.append(hrow)

            mr_rows.append({
                "exposure_id": exp_id,
                "outcome_id": "POAG",
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
                "outcome_source_dataset": OUTCOME["dataset_id"],
                "outcome_scale": OUTCOME["scale"],
                "F_stat": exp["F_stat"],
                "instrument_note": exp["instrument_note"],
            })

        harm_file = os.path.join(HARM_DIR, f"{exp_id}__POAG.targeted_external_harmonized.tsv.gz")
        mr_file = os.path.join(MR_DIR, f"{exp_id}__POAG.external_mr_input.tsv.gz")
        miss_file = os.path.join(MISS_DIR, f"{exp_id}__POAG.missing_external_instruments.tsv")

        write_gz_tsv(harm_file, HARM_HEADER, harmonized_rows)
        write_gz_tsv(mr_file, MR_HEADER, mr_rows)
        write_tsv(miss_file, MISSING_HEADER, missing_rows)

        elapsed = time.time() - t0

        qc_rows.append({
            "exposure_id": exp_id,
            "outcome_id": "POAG",
            "outcome_scale": OUTCOME["scale"],
            "n_clumped_instruments": str(len(exp_list)),
            "n_snp_matched_in_outcome": str(n_snp_match),
            "n_harmonized_for_mr": str(len(mr_rows)),
            "n_missing_from_outcome": str(n_missing),
            "n_allele_mismatch_or_unusable": str(n_allele_mismatch),
            "n_duplicate_outcome_records_extra": str(n_duplicate_extra),
            "n_aligned": str(n_aligned),
            "n_flipped": str(n_flipped),
            "harmonized_file": harm_file,
            "mr_input_file": mr_file,
            "missing_file": miss_file,
            "status": "MR_INPUT_READY" if len(mr_rows) > 0 else "NO_MATCHED_INSTRUMENTS",
            "note": "POAG external validation targeted harmonization",
        })

        manifest_rows.append({
            "exposure_id": exp_id,
            "outcome_id": "POAG",
            "outcome_scale": OUTCOME["scale"],
            "mr_input_file": mr_file,
            "n_mr_input_rows": str(len(mr_rows)),
            "status": "MR_INPUT_READY" if len(mr_rows) > 0 else "NO_MATCHED_INSTRUMENTS",
        })

        runtime_rows.append({
            "exposure_id": exp_id,
            "outcome_id": "POAG",
            "elapsed_seconds": f"{elapsed:.3f}",
            "elapsed_human": f"{elapsed:.1f}s",
        })

    write_tsv(QC_OUT, [
        "exposure_id", "outcome_id", "outcome_scale",
        "n_clumped_instruments", "n_snp_matched_in_outcome", "n_harmonized_for_mr",
        "n_missing_from_outcome", "n_allele_mismatch_or_unusable",
        "n_duplicate_outcome_records_extra", "n_aligned", "n_flipped",
        "harmonized_file", "mr_input_file", "missing_file", "status", "note",
    ], qc_rows)

    write_tsv(MANIFEST_OUT, [
        "exposure_id", "outcome_id", "outcome_scale", "mr_input_file",
        "n_mr_input_rows", "status",
    ], manifest_rows)

    write_tsv(RUNTIME_OUT, [
        "exposure_id", "outcome_id", "elapsed_seconds", "elapsed_human",
    ], runtime_rows)

    ready = sum(1 for r in qc_rows if r["status"] == "MR_INPUT_READY")
    total_rows = sum(int(r["n_harmonized_for_mr"]) for r in qc_rows)
    elapsed_all = time.time() - start

    with open(STATUS_OUT, "w", encoding="utf-8", newline="") as f:
        f.write("phase\tstatus\tkey_result\tnote\n")
        f.write(f"Phase 5.8D targeted POAG harmonization\tPASSED\t{ready}/2 pairwise POAG MR inputs ready; total_rows={total_rows}\tSBP and ART_STIFFNESS instruments harmonized against POAG\n")
        f.write("outcome\tINFO\tPOAG\tExternal clinical glaucoma validation outcome\n")
        f.write("outcome_scale\tDOCUMENTED\tlog_odds\tPOAG is a binary disease outcome\n")
        f.write(f"poag_scan_runtime\tINFO\t{scan_elapsed:.3f}s\tTime used to scan standardized POAG file\n")
        f.write(f"runtime\tINFO\t{elapsed_all:.3f}s\tPhase 5.8D runtime\n")

    print("===== Phase 5.8D completed =====")
    print("Wrote:", STATUS_OUT)
    print("Wrote:", QC_OUT)
    print("Wrote:", MANIFEST_OUT)

if __name__ == "__main__":
    main()
