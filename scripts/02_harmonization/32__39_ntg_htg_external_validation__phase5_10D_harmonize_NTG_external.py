#!/usr/bin/env python3
import csv
import gzip
import os
import time

BASE = "../../40_ntg_external_validation_inputs"
PAIRWISE_DIR = os.path.join(BASE, "pairwise_harmonized")
MR_INPUT_DIR = os.path.join(BASE, "mr_input")
MISSING_DIR = os.path.join(BASE, "missing_instruments")
os.makedirs(PAIRWISE_DIR, exist_ok=True)
os.makedirs(MR_INPUT_DIR, exist_ok=True)
os.makedirs(MISSING_DIR, exist_ok=True)

NTG_STD = "../../27_external_outcome_standardized/NTG__Zenodo14010557.standardized.tsv.gz"

INSTRUMENT_SOURCES = [
    {
        "exposure_id": "SBP",
        "source_file": "../../28_external_outcome_triangulation_inputs/mr_input/SBP__RNFL.external_mr_input.tsv.gz",
        "instrument_note": "SBP LD-clumped instrument from prior external harmonization source"
    },
    {
        "exposure_id": "ART_STIFFNESS",
        "source_file": "../../28_external_outcome_triangulation_inputs/mr_input/ART_STIFFNESS__RNFL.external_mr_input.tsv.gz",
        "instrument_note": "ART_STIFFNESS strict p<5e-8 LD-clumped instrument"
    },
]

STATUS_OUT = os.path.join(BASE, "phase5_10D_status.tsv")
QC_OUT = os.path.join(BASE, "phase5_10D_NTG_harmonization_qc_summary.tsv")
RUNTIME_OUT = os.path.join(BASE, "phase5_10D_runtime_log.tsv")

def is_missing(x):
    return x is None or str(x).strip() in ("", ".", "NA", "NaN", "nan")

def fnum(x):
    try:
        if is_missing(x):
            return None
        return float(x)
    except Exception:
        return None

def clean_allele(x):
    if is_missing(x):
        return "NA"
    return str(x).strip().upper()

def complement(a):
    comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
    a = clean_allele(a)
    if len(a) != 1:
        return a
    return comp.get(a, a)

def is_palindromic(a1, a2):
    a1 = clean_allele(a1)
    a2 = clean_allele(a2)
    return (a1, a2) in {("A", "T"), ("T", "A"), ("C", "G"), ("G", "C")}

def flip_beta(beta):
    b = fnum(beta)
    if b is None:
        return "NA"
    return f"{-b:.12g}"

def f_stat(beta, se):
    b = fnum(beta)
    s = fnum(se)
    if b is None or s is None or s <= 0:
        return "NA"
    return f"{(b / s) ** 2:.12g}"

def read_instruments(path, exposure_id):
    rows = []
    with gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for r in reader:
            snp = r.get("SNP", "NA")
            if is_missing(snp):
                continue
            rows.append({
                "exposure_id": exposure_id,
                "SNP": snp,
                "chr": r.get("chr", "NA"),
                "pos": r.get("pos", "NA"),
                "effect_allele": clean_allele(r.get("effect_allele", "NA")),
                "other_allele": clean_allele(r.get("other_allele", "NA")),
                "beta_exposure": r.get("beta_exposure", "NA"),
                "se_exposure": r.get("se_exposure", "NA"),
                "pval_exposure": r.get("pval_exposure", "NA"),
                "minus_log10_p_exposure": r.get("minus_log10_p_exposure", "NA"),
                "eaf_exposure": r.get("eaf_exposure", "NA"),
                "n_exposure": r.get("n_exposure", "NA"),
                "F_stat": r.get("F_stat", f_stat(r.get("beta_exposure", "NA"), r.get("se_exposure", "NA"))),
            })
    return rows

def read_ntg_for_snps(target_snps):
    matched = {}
    n_scan = 0
    with gzip.open(NTG_STD, "rt", encoding="utf-8", errors="replace", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for r in reader:
            n_scan += 1
            snp = r.get("SNP", "")
            if snp in target_snps:
                matched.setdefault(snp, []).append(r)
    return matched, n_scan

def evaluate_match(inst, out):
    exp_ea = clean_allele(inst["effect_allele"])
    exp_oa = clean_allele(inst["other_allele"])
    out_ea = clean_allele(out.get("effect_allele"))
    out_oa = clean_allele(out.get("other_allele"))
    beta_out = out.get("beta", "NA")

    if exp_ea == out_ea and exp_oa == out_oa:
        return 1, "aligned", beta_out
    if exp_ea == out_oa and exp_oa == out_ea:
        return 2, "flipped", flip_beta(beta_out)

    comp_ea = complement(out_ea)
    comp_oa = complement(out_oa)

    if exp_ea == comp_ea and exp_oa == comp_oa:
        return 3, "strand_complement_aligned", beta_out
    if exp_ea == comp_oa and exp_oa == comp_ea:
        return 4, "strand_complement_flipped", flip_beta(beta_out)

    return 99, "allele_mismatch", "NA"

def choose_best_match(inst, outcome_rows):
    candidates = []
    for out in outcome_rows:
        priority, action, beta_out = evaluate_match(inst, out)
        p = fnum(out.get("pval"))
        candidates.append((priority, p if p is not None else 1.0, action, beta_out, out))
    candidates.sort(key=lambda x: (x[0], x[1]))
    best = candidates[0]
    return {
        "priority": best[0],
        "harmonization_action": best[2],
        "beta_outcome": best[3],
        "outcome_row": best[4],
        "n_duplicate_outcome_records_extra": max(0, len(outcome_rows) - 1),
    }

def write_gz_tsv(path, fieldnames, rows):
    with gzip.open(path, "wt", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

def write_tsv(path, fieldnames, rows):
    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

def main():
    start_all = time.time()

    all_instruments = {}
    all_snps = set()

    for spec in INSTRUMENT_SOURCES:
        exposure_id = spec["exposure_id"]
        rows = read_instruments(spec["source_file"], exposure_id)
        all_instruments[exposure_id] = rows
        for r in rows:
            all_snps.add(r["SNP"])

    print(f"Scanning NTG standardized file for {len(all_snps)} target SNPs", flush=True)
    ntg_matches, n_scan = read_ntg_for_snps(all_snps)

    qc_rows = []
    runtime_rows = []

    harmonized_fields = [
        "exposure_id", "outcome_id", "outcome_scale",
        "SNP", "chr", "pos",
        "effect_allele_exposure", "other_allele_exposure",
        "beta_exposure", "se_exposure", "pval_exposure", "minus_log10_p_exposure",
        "eaf_exposure", "n_exposure",
        "effect_allele_outcome", "other_allele_outcome",
        "beta_outcome", "se_outcome", "pval_outcome", "minus_log10_p_outcome",
        "outcome_eaf", "outcome_n",
        "harmonization_action", "match_mode", "match_key",
        "palindromic_possible",
        "n_duplicate_outcome_records_extra",
        "outcome_source_dataset", "F_stat", "instrument_note"
    ]

    mr_fields = [
        "exposure_id", "outcome_id", "SNP", "chr", "pos",
        "effect_allele", "other_allele",
        "beta_exposure", "se_exposure", "pval_exposure", "minus_log10_p_exposure",
        "eaf_exposure", "n_exposure",
        "beta_outcome", "se_outcome", "pval_outcome", "pval_for_log10_outcome",
        "harmonization_action", "match_mode", "match_key",
        "outcome_source_dataset", "F_stat", "instrument_note",
        "outcome_scale", "outcome_eaf", "outcome_n", "palindromic_possible"
    ]

    missing_fields = ["exposure_id", "outcome_id", "SNP", "reason", "effect_allele", "other_allele"]

    for spec in INSTRUMENT_SOURCES:
        t0 = time.time()

        exposure_id = spec["exposure_id"]
        instruments = all_instruments[exposure_id]

        harmonized = []
        mr_input = []
        missing = []

        n_matched = 0
        n_harmonized = 0
        n_missing = 0
        n_allele_mismatch = 0
        n_duplicate_extra = 0
        n_aligned = 0
        n_flipped = 0
        n_strand_rescue = 0
        n_palindromic = 0

        for inst in instruments:
            snp = inst["SNP"]
            outcome_rows = ntg_matches.get(snp, [])

            if not outcome_rows:
                n_missing += 1
                missing.append({
                    "exposure_id": exposure_id,
                    "outcome_id": "NTG",
                    "SNP": snp,
                    "reason": "missing_from_NTG",
                    "effect_allele": inst["effect_allele"],
                    "other_allele": inst["other_allele"],
                })
                continue

            n_matched += 1
            best = choose_best_match(inst, outcome_rows)

            if best["priority"] == 99:
                n_allele_mismatch += 1
                missing.append({
                    "exposure_id": exposure_id,
                    "outcome_id": "NTG",
                    "SNP": snp,
                    "reason": "allele_mismatch",
                    "effect_allele": inst["effect_allele"],
                    "other_allele": inst["other_allele"],
                })
                continue

            out = best["outcome_row"]
            action = best["harmonization_action"]

            if action == "aligned":
                n_aligned += 1
            elif action == "flipped":
                n_flipped += 1
            elif action.startswith("strand_complement"):
                n_strand_rescue += 1

            n_duplicate_extra += best["n_duplicate_outcome_records_extra"]

            pal = "YES" if is_palindromic(inst["effect_allele"], inst["other_allele"]) else "NO"
            if pal == "YES":
                n_palindromic += 1

            row = {
                "exposure_id": exposure_id,
                "outcome_id": "NTG",
                "outcome_scale": "log_odds_or_MTAG_scale_TO_VERIFY",
                "SNP": snp,
                "chr": inst["chr"],
                "pos": inst["pos"],
                "effect_allele_exposure": inst["effect_allele"],
                "other_allele_exposure": inst["other_allele"],
                "beta_exposure": inst["beta_exposure"],
                "se_exposure": inst["se_exposure"],
                "pval_exposure": inst["pval_exposure"],
                "minus_log10_p_exposure": inst["minus_log10_p_exposure"],
                "eaf_exposure": inst["eaf_exposure"],
                "n_exposure": inst["n_exposure"],
                "effect_allele_outcome": clean_allele(out.get("effect_allele")),
                "other_allele_outcome": clean_allele(out.get("other_allele")),
                "beta_outcome": best["beta_outcome"],
                "se_outcome": out.get("se", "NA"),
                "pval_outcome": out.get("pval", "NA"),
                "minus_log10_p_outcome": out.get("minus_log10_p", "NA"),
                "outcome_eaf": out.get("eaf", "NA"),
                "outcome_n": out.get("n", "NA"),
                "harmonization_action": action,
                "match_mode": "rsid",
                "match_key": snp,
                "palindromic_possible": pal,
                "n_duplicate_outcome_records_extra": str(best["n_duplicate_outcome_records_extra"]),
                "outcome_source_dataset": "NTG_Zenodo14010557",
                "F_stat": inst["F_stat"],
                "instrument_note": spec["instrument_note"],
            }

            harmonized.append(row)

            mr_input.append({
                "exposure_id": exposure_id,
                "outcome_id": "NTG",
                "SNP": snp,
                "chr": inst["chr"],
                "pos": inst["pos"],
                "effect_allele": inst["effect_allele"],
                "other_allele": inst["other_allele"],
                "beta_exposure": inst["beta_exposure"],
                "se_exposure": inst["se_exposure"],
                "pval_exposure": inst["pval_exposure"],
                "minus_log10_p_exposure": inst["minus_log10_p_exposure"],
                "eaf_exposure": inst["eaf_exposure"],
                "n_exposure": inst["n_exposure"],
                "beta_outcome": best["beta_outcome"],
                "se_outcome": out.get("se", "NA"),
                "pval_outcome": out.get("pval", "NA"),
                "pval_for_log10_outcome": out.get("pval", "NA"),
                "harmonization_action": action,
                "match_mode": "rsid",
                "match_key": snp,
                "outcome_source_dataset": "NTG_Zenodo14010557",
                "F_stat": inst["F_stat"],
                "instrument_note": spec["instrument_note"],
                "outcome_scale": "log_odds_or_MTAG_scale_TO_VERIFY",
                "outcome_eaf": out.get("eaf", "NA"),
                "outcome_n": out.get("n", "NA"),
                "palindromic_possible": pal,
            })

            n_harmonized += 1

        hfile = os.path.join(PAIRWISE_DIR, f"{exposure_id}__NTG.targeted_external_harmonized.tsv.gz")
        mfile = os.path.join(MR_INPUT_DIR, f"{exposure_id}__NTG.external_mr_input.tsv.gz")
        missfile = os.path.join(MISSING_DIR, f"{exposure_id}__NTG.missing_external_instruments.tsv")

        write_gz_tsv(hfile, harmonized_fields, harmonized)
        write_gz_tsv(mfile, mr_fields, mr_input)
        write_tsv(missfile, missing_fields, missing)

        status = "MR_INPUT_READY" if n_harmonized > 0 else "NO_MR_INPUT"

        qc_rows.append({
            "exposure_id": exposure_id,
            "outcome_id": "NTG",
            "outcome_scale": "log_odds_or_MTAG_scale_TO_VERIFY",
            "n_clumped_instruments": str(len(instruments)),
            "n_snp_matched_in_outcome": str(n_matched),
            "n_harmonized_for_mr": str(n_harmonized),
            "n_missing_from_outcome": str(n_missing),
            "n_allele_mismatch_or_unusable": str(n_allele_mismatch),
            "n_duplicate_outcome_records_extra": str(n_duplicate_extra),
            "n_aligned": str(n_aligned),
            "n_flipped": str(n_flipped),
            "n_strand_complement_rescued": str(n_strand_rescue),
            "n_palindromic_possible": str(n_palindromic),
            "harmonized_file": hfile,
            "mr_input_file": mfile,
            "missing_file": missfile,
            "status": status,
            "note": "NTG external validation targeted harmonization using Zenodo14010557",
        })

        elapsed = time.time() - t0
        runtime_rows.append({
            "exposure_id": exposure_id,
            "outcome_id": "NTG",
            "elapsed_seconds": f"{elapsed:.3f}",
            "elapsed_human": f"{elapsed:.1f}s",
        })

    write_tsv(
        QC_OUT,
        [
            "exposure_id", "outcome_id", "outcome_scale",
            "n_clumped_instruments", "n_snp_matched_in_outcome", "n_harmonized_for_mr",
            "n_missing_from_outcome", "n_allele_mismatch_or_unusable",
            "n_duplicate_outcome_records_extra", "n_aligned", "n_flipped",
            "n_strand_complement_rescued", "n_palindromic_possible",
            "harmonized_file", "mr_input_file", "missing_file",
            "status", "note"
        ],
        qc_rows
    )

    write_tsv(RUNTIME_OUT, ["exposure_id", "outcome_id", "elapsed_seconds", "elapsed_human"], runtime_rows)

    total_rows = sum(int(r["n_harmonized_for_mr"]) for r in qc_rows)
    ready_count = sum(1 for r in qc_rows if r["status"] == "MR_INPUT_READY")
    elapsed_all = time.time() - start_all

    with open(STATUS_OUT, "w", encoding="utf-8", newline="") as f:
        f.write("phase\tstatus\tkey_result\tnote\n")
        f.write(f"Phase 5.10D targeted NTG harmonization\tPASSED\t{ready_count}/2 pairwise NTG MR inputs ready; total_rows={total_rows}\tSBP and ART_STIFFNESS instruments harmonized against NTG\n")
        f.write("external_outcome\tINFO\tNTG\tNormal tension glaucoma external validation layer\n")
        f.write("primary_exposure\tINFO\tSBP\tMain vascular/hemodynamic signal\n")
        f.write("secondary_exposure\tINFO\tART_STIFFNESS\tExploratory vascular-stiffness signal\n")
        f.write("outcome_caveat\tDOCUMENTED\tNTG_EAF_N_missing;scale_TO_VERIFY\tRaw NTG file lacks EAF/N; BETA scale from MTAG file should be described cautiously\n")
        f.write("HTG_status\tTO_DO\tfull_summary_stats_not_locked\tHTG remains unresolved\n")
        f.write(f"runtime\tINFO\t{elapsed_all:.3f}s\tPhase 5.10D runtime\n")

    print("===== Phase 5.10D completed =====")
    print("Wrote:", STATUS_OUT)
    print("Wrote:", QC_OUT)

if __name__ == "__main__":
    main()
