#!/usr/bin/env python3
import csv
import gzip
import math
import os
import time

BASE = "../../31_iop_poag_ntg_htg_source_reidentification"
OUTDIR = "../../27_external_outcome_standardized"
os.makedirs(OUTDIR, exist_ok=True)

INPUT_FILE = "../../26_external_outcome_raw/POAG__GCST90011766__GCST90011766.h.tsv.gz"
OUTPUT_FILE = os.path.join(OUTDIR, "POAG__GCST90011766.standardized.tsv.gz")

STATUS_OUT = os.path.join(BASE, "phase5_8C_status.tsv")
QC_OUT = os.path.join(BASE, "phase5_8C_standardized_POAG_qc_summary.tsv")
MANIFEST_OUT = os.path.join(BASE, "phase5_8C_standardized_POAG_manifest.tsv")
RUNTIME_OUT = os.path.join(BASE, "phase5_8C_runtime_log.tsv")

OUT_HEADER = [
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
    "column_source",
]

def is_missing(x):
    return x is None or str(x).strip() in ("", ".", "NA", "na", "NaN", "nan")

def clean(x):
    if x is None:
        return "NA"
    x = str(x).strip()
    return x if x else "NA"

def safe_float(x):
    try:
        if is_missing(x):
            return None
        return float(x)
    except Exception:
        return None

def minus_log10_p(p):
    p = safe_float(p)
    if p is None:
        return "NA"
    if p <= 0:
        return "NA"
    return f"{-math.log10(p):.12g}"

def valid_core(row):
    for k in ["SNP", "effect_allele", "other_allele", "beta", "se", "pval"]:
        if is_missing(row.get(k)):
            return False
    return True

def main():
    t0 = time.time()

    n_input = 0
    n_written = 0
    n_skipped = 0
    n_missing_snp = 0
    n_missing_allele = 0
    n_missing_beta = 0
    n_missing_se = 0
    n_missing_pval = 0
    n_missing_eaf = 0
    n_duplicate_variant = 0
    n_pval_zero = 0
    max_mlp = None
    seen_variant = set()

    with gzip.open(INPUT_FILE, "rt", encoding="utf-8", errors="replace", newline="") as fin, \
         gzip.open(OUTPUT_FILE, "wt", encoding="utf-8", newline="") as fout:

        reader = csv.DictReader(fin, delimiter="\t")
        writer = csv.DictWriter(fout, fieldnames=OUT_HEADER, delimiter="\t", extrasaction="ignore")
        writer.writeheader()

        for row in reader:
            n_input += 1

            chr_ = clean(row.get("chromosome"))
            pos = clean(row.get("base_pair_location"))
            snp = clean(row.get("rsid"))
            variant_id = clean(row.get("variant_id"))
            ea = clean(row.get("effect_allele")).upper()
            oa = clean(row.get("other_allele")).upper()
            beta = clean(row.get("beta"))
            se = clean(row.get("standard_error"))
            pval = clean(row.get("p_value"))
            eaf = clean(row.get("effect_allele_frequency"))

            if is_missing(snp):
                n_missing_snp += 1
            if is_missing(ea) or is_missing(oa):
                n_missing_allele += 1
            if is_missing(beta):
                n_missing_beta += 1
            if is_missing(se):
                n_missing_se += 1
            if is_missing(pval):
                n_missing_pval += 1
            if is_missing(eaf):
                n_missing_eaf += 1

            p = safe_float(pval)
            if p == 0:
                n_pval_zero += 1

            mlp = minus_log10_p(pval)
            mlp_num = safe_float(mlp)
            if mlp_num is not None and (max_mlp is None or mlp_num > max_mlp):
                max_mlp = mlp_num

            out = {
                "chr": chr_,
                "pos": pos,
                "SNP": snp,
                "variant_id": variant_id,
                "effect_allele": ea,
                "other_allele": oa,
                "beta": beta,
                "se": se,
                "pval": pval,
                "minus_log10_p": mlp,
                "eaf": eaf,
                "n": "NA",
                "genome_build": "GRCh37",
                "source_dataset": "GCST90011766",
                "outcome_id": "POAG",
                "trait_name": "Primary open-angle glaucoma",
                "population": "European",
                "outcome_scale": "log_odds",
                "column_source": "GWAS_CATALOG_COLUMNS",
            }

            if not valid_core(out):
                n_skipped += 1
                continue

            vkey = variant_id if not is_missing(variant_id) else f"{chr_}:{pos}_{ea}_{oa}"
            if vkey in seen_variant:
                n_duplicate_variant += 1
            else:
                seen_variant.add(vkey)

            writer.writerow(out)
            n_written += 1

    elapsed = time.time() - t0

    qc = {
        "outcome_id": "POAG",
        "dataset_id": "GCST90011766",
        "input_file": INPUT_FILE,
        "standardized_file": OUTPUT_FILE,
        "n_input_rows": str(n_input),
        "n_written_rows": str(n_written),
        "n_skipped_missing_core": str(n_skipped),
        "n_missing_snp": str(n_missing_snp),
        "n_missing_allele": str(n_missing_allele),
        "n_missing_beta": str(n_missing_beta),
        "n_missing_se": str(n_missing_se),
        "n_missing_pval": str(n_missing_pval),
        "n_pval_zero": str(n_pval_zero),
        "n_missing_eaf": str(n_missing_eaf),
        "n_duplicate_variant_id": str(n_duplicate_variant),
        "max_minus_log10_p": f"{max_mlp:.12g}" if max_mlp is not None else "NA",
        "sample_size_used": "NA",
        "outcome_scale": "log_odds",
        "status": "STANDARDIZED" if n_written > 0 else "NO_ROWS_WRITTEN",
        "elapsed_seconds": f"{elapsed:.3f}",
    }

    with open(QC_OUT, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(qc.keys()), delimiter="\t")
        writer.writeheader()
        writer.writerow(qc)

    manifest = {
        "outcome_id": "POAG",
        "dataset_id": "GCST90011766",
        "trait_name": "Primary open-angle glaucoma",
        "standardized_file": OUTPUT_FILE,
        "n_written_rows": str(n_written),
        "genome_build": "GRCh37",
        "population": "European",
        "sample_size": "NA",
        "outcome_scale": "log_odds",
        "status": qc["status"],
    }

    with open(MANIFEST_OUT, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(manifest.keys()), delimiter="\t")
        writer.writeheader()
        writer.writerow(manifest)

    with open(RUNTIME_OUT, "w", encoding="utf-8", newline="") as f:
        f.write("phase\telapsed_seconds\telapsed_human\n")
        f.write(f"Phase 5.8C\t{elapsed:.3f}\t{elapsed:.1f}s\n")

    with open(STATUS_OUT, "w", encoding="utf-8", newline="") as f:
        f.write("phase\tstatus\tkey_result\tnote\n")
        status = "PASSED" if qc["status"] == "STANDARDIZED" else "REVIEW_NEEDED"
        f.write(f"Phase 5.8C standardize POAG external outcome\t{status}\tn_written={n_written};n_skipped={n_skipped}\tPOAG standardized on log-odds scale\n")
        f.write("outcome\tINFO\tPOAG\tPrimary open-angle glaucoma external validation outcome\n")
        f.write("scale\tDOCUMENTED\tlog_odds\tBinary disease GWAS; beta interpreted as log odds\n")
        f.write(f"runtime\tINFO\t{elapsed:.3f}s\tPhase 5.8C runtime\n")

    print("===== Phase 5.8C completed =====")
    print("Wrote:", STATUS_OUT)
    print("Wrote:", QC_OUT)
    print("Wrote:", MANIFEST_OUT)
    print("Wrote:", OUTPUT_FILE)

if __name__ == "__main__":
    main()
