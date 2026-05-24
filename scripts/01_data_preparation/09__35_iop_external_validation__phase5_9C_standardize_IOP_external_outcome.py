#!/usr/bin/env python3
import csv
import gzip
import math
import os
import time

BASE = "../../35_iop_external_validation"
INFILE = "../../26_external_outcome_raw/IOP__GCST009413__Bonnemaijer_Intraocular_pressure_pmid.txt.gz"
OUTDIR = "../../27_external_outcome_standardized"
os.makedirs(BASE, exist_ok=True)
os.makedirs(OUTDIR, exist_ok=True)

OUTFILE = os.path.join(OUTDIR, "IOP__GCST009413.standardized.tsv.gz")
QC_OUT = os.path.join(BASE, "phase5_9C_standardized_IOP_qc_summary.tsv")
STATUS_OUT = os.path.join(BASE, "phase5_9C_status.tsv")
RUNTIME_OUT = os.path.join(BASE, "phase5_9C_runtime_log.tsv")

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

def clean_snp(snp, marker):
    if not is_missing(snp):
        return str(snp).strip()
    if not is_missing(marker):
        return str(marker).strip()
    return "NA"

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
    n_missing_n = 0
    n_missing_maf = 0
    max_lp = None
    seen_variant = set()
    n_dup_variant = 0

    out_fields = [
        "chr", "pos", "SNP", "variant_id",
        "effect_allele", "other_allele",
        "beta", "se", "pval", "minus_log10_p",
        "eaf", "maf", "n",
        "genome_build", "source_dataset", "outcome_id",
        "trait_name", "population", "outcome_scale", "column_source"
    ]

    with gzip.open(INFILE, "rt", encoding="utf-8", errors="replace", newline="") as fin, \
         gzip.open(OUTFILE, "wt", encoding="utf-8", newline="") as fout:

        reader = csv.DictReader(fin, delimiter="\t")
        writer = csv.DictWriter(fout, fieldnames=out_fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()

        for r in reader:
            n_in += 1

            marker = r.get("MarkerName")
            snp = clean_snp(r.get("SNP"), marker)
            chr_ = r.get("CHR", "NA")
            pos = r.get("POS", "NA")
            ea = clean_allele(r.get("EffectAllele"))
            oa = clean_allele(r.get("NonEffectAllele"))
            beta = r.get("Effect")
            se = r.get("StdErr")
            pval = r.get("P-value")
            n = r.get("Ntot")
            maf = r.get("MAF")

            if snp == "NA":
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
            if is_missing(n):
                n_missing_n += 1
            if is_missing(maf):
                n_missing_maf += 1

            if snp == "NA" or ea == "NA" or oa == "NA" or is_missing(beta) or is_missing(se) or is_missing(pval):
                n_skipped += 1
                continue

            variant_id = str(marker).strip() if not is_missing(marker) else f"{chr_}:{pos}_{ea}_{oa}"

            if variant_id in seen_variant:
                n_dup_variant += 1
            else:
                seen_variant.add(variant_id)

            lp = minus_log10_p(pval)
            lp_num = None
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
                "maf": maf if not is_missing(maf) else "NA",
                "n": n if not is_missing(n) else "NA",
                "genome_build": "GRCh37",
                "source_dataset": "GCST009413",
                "outcome_id": "IOP",
                "trait_name": "Intraocular pressure",
                "population": "European",
                "outcome_scale": "continuous_IOP",
                "column_source": "Bonnemaijer_GCST009413_COLUMNS"
            })
            n_written += 1

    elapsed = time.time() - start

    qc = {
        "outcome_id": "IOP",
        "dataset_id": "GCST009413",
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
        "n_missing_maf": str(n_missing_maf),
        "n_missing_n": str(n_missing_n),
        "n_duplicate_variant_id": str(n_dup_variant),
        "max_minus_log10_p": "Inf" if max_lp == float("inf") else (f"{max_lp:.12g}" if max_lp is not None else "NA"),
        "sample_size_used": "per_variant_Ntot",
        "outcome_scale": "continuous_IOP",
        "status": "STANDARDIZED",
        "elapsed_seconds": f"{elapsed:.3f}",
    }

    qc_fields = list(qc.keys())
    with open(QC_OUT, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=qc_fields, delimiter="\t")
        w.writeheader()
        w.writerow(qc)

    with open(STATUS_OUT, "w", encoding="utf-8", newline="") as f:
        f.write("phase\tstatus\tkey_result\tnote\n")
        f.write(f"Phase 5.9C standardize IOP external outcome\tPASSED\tn_written={n_written};n_skipped={n_skipped}\tIOP standardized from GCST009413\n")
        f.write("eaf_handling\tDOCUMENTED\tEAF_set_to_NA;MAF_retained\tMAF is minor allele frequency and not treated as effect allele frequency\n")
        f.write("outcome_scale\tDOCUMENTED\tcontinuous_IOP\tBeta interpreted as effect on measured intraocular pressure\n")
        f.write(f"runtime\tINFO\t{elapsed:.3f}s\tPhase 5.9C runtime\n")

    with open(RUNTIME_OUT, "w", encoding="utf-8", newline="") as f:
        f.write("phase\telapsed_seconds\telapsed_human\n")
        f.write(f"Phase 5.9C\t{elapsed:.3f}\t{elapsed:.1f}s\n")

    print("===== Phase 5.9C completed =====")
    print("Wrote:", STATUS_OUT)
    print("Wrote:", QC_OUT)
    print("Wrote:", OUTFILE)

if __name__ == "__main__":
    main()
