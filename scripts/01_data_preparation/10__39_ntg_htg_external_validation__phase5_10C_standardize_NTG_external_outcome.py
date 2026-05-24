#!/usr/bin/env python3
import csv
import gzip
import math
import os
import time

BASE = "../../39_ntg_htg_external_validation"
INFILE = "../../26_external_outcome_raw/NTG__Zenodo14010557__MTAG_NTG_IGGC_STAGE2.tab"
OUTDIR = "../../27_external_outcome_standardized"
os.makedirs(BASE, exist_ok=True)
os.makedirs(OUTDIR, exist_ok=True)

OUTFILE = os.path.join(OUTDIR, "NTG__Zenodo14010557.standardized.tsv.gz")
QC_OUT = os.path.join(BASE, "phase5_10C_standardized_NTG_qc_summary.tsv")
STATUS_OUT = os.path.join(BASE, "phase5_10C_status.tsv")
RUNTIME_OUT = os.path.join(BASE, "phase5_10C_runtime_log.tsv")

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

def minus_log10_p(p):
    x = fnum(p)
    if x is None:
        return "NA"
    if x <= 0:
        return "Inf"
    return f"{-math.log10(x):.12g}"

def main():
    start = time.time()

    n_in = 0
    n_written = 0
    n_skipped = 0
    n_missing_snp = 0
    n_missing_allele = 0
    n_missing_beta = 0
    n_missing_se = 0
    n_missing_pval = 0
    n_pval_zero = 0
    n_duplicate_variant_id = 0
    seen = set()
    max_lp = None

    out_fields = [
        "chr",
        "pos",
        "SNP",
        "variant_id",
        "effect_allele",
        "other_allele",
        "beta",
        "se",
        "pval",
        "minus_log10_p",
        "eaf",
        "n",
        "genome_build",
        "source_dataset",
        "outcome_id",
        "trait_name",
        "population",
        "outcome_scale",
        "column_source"
    ]

    with open(INFILE, "r", encoding="utf-8", errors="replace", newline="") as fin, \
         gzip.open(OUTFILE, "wt", encoding="utf-8", newline="") as fout:

        reader = csv.DictReader(fin, delimiter="\t")
        writer = csv.DictWriter(fout, fieldnames=out_fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()

        for r in reader:
            n_in += 1

            snp = r.get("SNP", "NA")
            chr_ = r.get("CHR", "NA")
            pos = r.get("BP", "NA")
            ea = clean_allele(r.get("A1"))
            oa = clean_allele(r.get("A2"))
            beta = r.get("BETA")
            se = r.get("SE")
            pval = r.get("P")

            if is_missing(snp):
                n_missing_snp += 1
            if ea == "NA" or oa == "NA":
                n_missing_allele += 1
            if is_missing(beta):
                n_missing_beta += 1
            if is_missing(se):
                n_missing_se += 1
            if is_missing(pval):
                n_missing_pval += 1
            if fnum(pval) == 0:
                n_pval_zero += 1

            if is_missing(snp) or ea == "NA" or oa == "NA" or is_missing(beta) or is_missing(se) or is_missing(pval):
                n_skipped += 1
                continue

            variant_id = str(snp).strip()

            if variant_id in seen:
                n_duplicate_variant_id += 1
            else:
                seen.add(variant_id)

            lp = minus_log10_p(pval)
            if lp == "Inf":
                lp_num = float("inf")
            else:
                lp_num = fnum(lp)

            if lp_num is not None:
                if max_lp is None or lp_num > max_lp:
                    max_lp = lp_num

            writer.writerow({
                "chr": chr_,
                "pos": pos,
                "SNP": snp,
                "variant_id": variant_id,
                "effect_allele": ea,
                "other_allele": oa,
                "beta": beta,
                "se": se,
                "pval": pval,
                "minus_log10_p": lp,
                "eaf": "NA",
                "n": "NA",
                "genome_build": "GRCh37_TO_VERIFY",
                "source_dataset": "Zenodo14010557",
                "outcome_id": "NTG",
                "trait_name": "Normal tension glaucoma",
                "population": "European_TO_VERIFY",
                "outcome_scale": "log_odds_or_MTAG_scale_TO_VERIFY",
                "column_source": "MTAG_NTG_IGGC_STAGE2_COLUMNS"
            })

            n_written += 1

    elapsed = time.time() - start

    qc = {
        "outcome_id": "NTG",
        "dataset_id": "Zenodo14010557",
        "input_file": INFILE,
        "standardized_file": OUTFILE,
        "n_input_rows": str(n_in),
        "n_written_rows": str(n_written),
        "n_skipped_missing_core": str(n_skipped),
        "n_missing_snp": str(n_missing_snp),
        "n_missing_allele": str(n_missing_allele),
        "n_missing_beta": str(n_missing_beta),
        "n_missing_se": str(n_missing_se),
        "n_missing_pval": str(n_missing_pval),
        "n_pval_zero": str(n_pval_zero),
        "n_missing_eaf": str(n_written),
        "n_missing_n": str(n_written),
        "n_duplicate_variant_id": str(n_duplicate_variant_id),
        "max_minus_log10_p": "Inf" if max_lp == float("inf") else (f"{max_lp:.12g}" if max_lp is not None else "NA"),
        "sample_size_used": "NA",
        "outcome_scale": "log_odds_or_MTAG_scale_TO_VERIFY",
        "status": "STANDARDIZED",
        "elapsed_seconds": f"{elapsed:.3f}",
    }

    with open(QC_OUT, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(qc.keys()), delimiter="\t")
        w.writeheader()
        w.writerow(qc)

    with open(STATUS_OUT, "w", encoding="utf-8", newline="") as f:
        f.write("phase\tstatus\tkey_result\tnote\n")
        f.write(f"Phase 5.10C standardize NTG external outcome\tPASSED\tn_written={n_written};n_skipped={n_skipped}\tNTG standardized from Zenodo14010557\n")
        f.write("eaf_handling\tDOCUMENTED\tEAF_set_to_NA\tRaw NTG file has no EAF column\n")
        f.write("n_handling\tDOCUMENTED\tN_set_to_NA\tRaw NTG file has no sample-size column\n")
        f.write("outcome_scale\tDOCUMENTED\tlog_odds_or_MTAG_scale_TO_VERIFY\tBETA from MTAG_NTG_IGGC_STAGE2.tab\n")
        f.write(f"runtime\tINFO\t{elapsed:.3f}s\tPhase 5.10C runtime\n")

    with open(RUNTIME_OUT, "w", encoding="utf-8", newline="") as f:
        f.write("phase\telapsed_seconds\telapsed_human\n")
        f.write(f"Phase 5.10C\t{elapsed:.3f}\t{elapsed:.1f}s\n")

    print("===== Phase 5.10C completed =====")
    print("Wrote:", STATUS_OUT)
    print("Wrote:", QC_OUT)
    print("Wrote:", OUTFILE)

if __name__ == "__main__":
    main()
