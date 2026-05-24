#!/usr/bin/env python3
import csv
import gzip
import math
import os
import time

BASE = "../../25_external_outcome_triangulation"
REGISTER = os.path.join(BASE, "phase5_4A_external_outcome_registration_rows.tsv")

OUTDIR = "../../27_external_outcome_standardized"
os.makedirs(OUTDIR, exist_ok=True)

STATUS_OUT = os.path.join(BASE, "phase5_4B_status.tsv")
QC_OUT = os.path.join(BASE, "phase5_4B_standardized_RNFL_GCIPL_qc_summary.tsv")
MANIFEST_OUT = os.path.join(BASE, "phase5_4B_standardized_external_outcome_manifest.tsv")
RUNTIME_OUT = os.path.join(BASE, "phase5_4B_runtime_log.tsv")

SAMPLE_SIZE = {
    "RNFL": "31434",
    "GCIPL": "31434",
}

TRAIT_NAME = {
    "RNFL": "Retinal nerve fibre layer thickness",
    "GCIPL": "Ganglion cell inner plexiform layer thickness",
}

DATASET_ID = {
    "RNFL": "GCST90014266",
    "GCIPL": "GCST90014267",
}

POPULATION = {
    "RNFL": "European",
    "GCIPL": "European",
}

GENOME_BUILD = "GRCh37"

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
    "column_source",
]

def is_missing(x):
    return x is None or x == "" or x == "." or str(x).upper() == "NA"

def clean(x):
    if x is None:
        return "NA"
    x = str(x).strip()
    return x if x != "" else "NA"

def get(row, key):
    return clean(row.get(key, "NA"))

def safe_float(x):
    try:
        if is_missing(x):
            return None
        return float(x)
    except Exception:
        return None

def minus_log10_p(p):
    p = safe_float(p)
    if p is None or p <= 0:
        return "NA"
    return f"{-math.log10(p):.12g}"

def valid_basic(row):
    required = ["SNP", "effect_allele", "other_allele", "beta", "se", "pval"]
    for k in required:
        if is_missing(row.get(k)):
            return False
    return True

def choose_columns(row):
    """
    Prefer harmonized GWAS Catalog fields when present.
    Fallback to original submitted fields when harmonized fields are missing.
    """

    hm_snp = get(row, "hm_rsid")
    hm_chr = get(row, "hm_chrom")
    hm_pos = get(row, "hm_pos")
    hm_ea = get(row, "hm_effect_allele")
    hm_oa = get(row, "hm_other_allele")
    hm_beta = get(row, "hm_beta")
    hm_eaf = get(row, "hm_effect_allele_frequency")
    hm_var = get(row, "hm_variant_id")

    use_hm = not any(is_missing(x) for x in [hm_snp, hm_chr, hm_pos, hm_ea, hm_oa, hm_beta])

    if use_hm:
        snp = hm_snp
        chr_ = hm_chr
        pos = hm_pos
        ea = hm_ea
        oa = hm_oa
        beta = hm_beta
        eaf = hm_eaf
        variant_id = hm_var if not is_missing(hm_var) else f"{chr_}:{pos}_{ea}_{oa}"
        source = "GWAS_CATALOG_HARMONIZED_COLUMNS"
    else:
        snp = get(row, "variant_id")
        chr_ = get(row, "chromosome")
        pos = get(row, "base_pair_location")
        ea = get(row, "effect_allele")
        oa = get(row, "other_allele")
        beta = get(row, "beta")
        eaf = get(row, "effect_allele_frequency")
        variant_id = f"{chr_}:{pos}_{ea}_{oa}" if not any(is_missing(x) for x in [chr_, pos, ea, oa]) else "NA"
        source = "ORIGINAL_COLUMNS_FALLBACK"

    se = get(row, "standard_error")
    pval = get(row, "p_value")

    return {
        "chr": chr_,
        "pos": pos,
        "SNP": snp,
        "variant_id": variant_id,
        "effect_allele": ea,
        "other_allele": oa,
        "beta": beta,
        "se": se,
        "pval": pval,
        "minus_log10_p": minus_log10_p(pval),
        "eaf": eaf,
        "column_source": source,
    }

def read_registration():
    rows = []
    with open(REGISTER, "r", encoding="utf-8", errors="replace", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            outcome = row.get("outcome_trait", "")
            if outcome in ("RNFL", "GCIPL"):
                rows.append(row)
    return rows

def standardize_one(reg):
    outcome_id = reg["outcome_trait"]
    dataset_id = reg["dataset_id"]
    input_file = reg["file_path"]

    output_file = os.path.join(
        OUTDIR,
        f"{outcome_id}__{dataset_id}.standardized.tsv.gz"
    )

    t0 = time.time()

    n_input = 0
    n_written = 0
    n_skipped = 0
    n_hm = 0
    n_fallback = 0
    n_missing_snp = 0
    n_missing_allele = 0
    n_missing_beta = 0
    n_missing_se = 0
    n_missing_p = 0
    n_p_zero = 0
    n_missing_eaf = 0
    max_mlp = None

    seen_variant = set()
    duplicate_variant_count = 0

    with gzip.open(input_file, "rt", encoding="utf-8", errors="replace", newline="") as fin, \
         gzip.open(output_file, "wt", encoding="utf-8", newline="") as fout:

        reader = csv.DictReader(fin, delimiter="\t")
        writer = csv.DictWriter(fout, fieldnames=OUT_HEADER, delimiter="\t", extrasaction="ignore")
        writer.writeheader()

        for raw in reader:
            n_input += 1

            std = choose_columns(raw)

            if std["column_source"] == "GWAS_CATALOG_HARMONIZED_COLUMNS":
                n_hm += 1
            else:
                n_fallback += 1

            if is_missing(std["SNP"]):
                n_missing_snp += 1
            if is_missing(std["effect_allele"]) or is_missing(std["other_allele"]):
                n_missing_allele += 1
            if is_missing(std["beta"]):
                n_missing_beta += 1
            if is_missing(std["se"]):
                n_missing_se += 1
            if is_missing(std["pval"]):
                n_missing_p += 1
            if is_missing(std["eaf"]):
                n_missing_eaf += 1

            p = safe_float(std["pval"])
            if p == 0:
                n_p_zero += 1

            mlp = safe_float(std["minus_log10_p"])
            if mlp is not None and (max_mlp is None or mlp > max_mlp):
                max_mlp = mlp

            if not valid_basic(std):
                n_skipped += 1
                continue

            vkey = std["variant_id"]
            if vkey in seen_variant:
                duplicate_variant_count += 1
            else:
                seen_variant.add(vkey)

            std.update({
                "n": SAMPLE_SIZE.get(outcome_id, "NA"),
                "genome_build": GENOME_BUILD,
                "source_dataset": dataset_id,
                "outcome_id": outcome_id,
                "trait_name": TRAIT_NAME.get(outcome_id, outcome_id),
                "population": POPULATION.get(outcome_id, "European"),
            })

            writer.writerow(std)
            n_written += 1

    elapsed = time.time() - t0

    qc = {
        "outcome_id": outcome_id,
        "dataset_id": dataset_id,
        "input_file": input_file,
        "standardized_file": output_file,
        "n_input_rows": str(n_input),
        "n_written_rows": str(n_written),
        "n_skipped_missing_core": str(n_skipped),
        "n_harmonized_column_rows": str(n_hm),
        "n_original_fallback_rows": str(n_fallback),
        "n_missing_snp": str(n_missing_snp),
        "n_missing_allele": str(n_missing_allele),
        "n_missing_beta": str(n_missing_beta),
        "n_missing_se": str(n_missing_se),
        "n_missing_pval": str(n_missing_p),
        "n_pval_zero": str(n_p_zero),
        "n_missing_eaf": str(n_missing_eaf),
        "n_duplicate_variant_id": str(duplicate_variant_count),
        "max_minus_log10_p": f"{max_mlp:.12g}" if max_mlp is not None else "NA",
        "sample_size_used": SAMPLE_SIZE.get(outcome_id, "NA"),
        "status": "STANDARDIZED" if n_written > 0 else "NO_ROWS_WRITTEN",
        "elapsed_seconds": f"{elapsed:.3f}",
    }

    manifest = {
        "outcome_id": outcome_id,
        "dataset_id": dataset_id,
        "trait_name": TRAIT_NAME.get(outcome_id, outcome_id),
        "standardized_file": output_file,
        "n_written_rows": str(n_written),
        "genome_build": GENOME_BUILD,
        "population": POPULATION.get(outcome_id, "European"),
        "sample_size": SAMPLE_SIZE.get(outcome_id, "NA"),
        "status": qc["status"],
    }

    runtime = {
        "outcome_id": outcome_id,
        "dataset_id": dataset_id,
        "elapsed_seconds": f"{elapsed:.3f}",
        "elapsed_human": f"{elapsed:.1f}s",
    }

    return qc, manifest, runtime

def write_tsv(path, fieldnames, rows):
    with open(path, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

def main():
    total_start = time.time()

    regs = read_registration()
    qc_rows = []
    manifest_rows = []
    runtime_rows = []

    for reg in regs:
        print(f"Standardizing {reg['outcome_trait']} from {reg['file_path']}", flush=True)
        qc, manifest, runtime = standardize_one(reg)
        qc_rows.append(qc)
        manifest_rows.append(manifest)
        runtime_rows.append(runtime)

    write_tsv(
        QC_OUT,
        [
            "outcome_id",
            "dataset_id",
            "input_file",
            "standardized_file",
            "n_input_rows",
            "n_written_rows",
            "n_skipped_missing_core",
            "n_harmonized_column_rows",
            "n_original_fallback_rows",
            "n_missing_snp",
            "n_missing_allele",
            "n_missing_beta",
            "n_missing_se",
            "n_missing_pval",
            "n_pval_zero",
            "n_missing_eaf",
            "n_duplicate_variant_id",
            "max_minus_log10_p",
            "sample_size_used",
            "status",
            "elapsed_seconds",
        ],
        qc_rows
    )

    write_tsv(
        MANIFEST_OUT,
        [
            "outcome_id",
            "dataset_id",
            "trait_name",
            "standardized_file",
            "n_written_rows",
            "genome_build",
            "population",
            "sample_size",
            "status",
        ],
        manifest_rows
    )

    write_tsv(
        RUNTIME_OUT,
        ["outcome_id", "dataset_id", "elapsed_seconds", "elapsed_human"],
        runtime_rows
    )

    n_standardized = sum(1 for r in qc_rows if r["status"] == "STANDARDIZED")
    total_rows = sum(int(r["n_written_rows"]) for r in qc_rows)
    elapsed_total = time.time() - total_start

    with open(STATUS_OUT, "w", encoding="utf-8", newline="") as f:
        f.write("phase\tstatus\tkey_result\tnote\n")
        status = "PASSED" if n_standardized == len(qc_rows) and n_standardized > 0 else "REVIEW_NEEDED"
        f.write(
            f"Phase 5.4B standardize RNFL and GCIPL external outcomes\t{status}\t"
            f"{n_standardized}/{len(qc_rows)} outcomes standardized; total_rows={total_rows}\t"
            f"GWAS Catalog harmonized columns preferred with original-column fallback\n"
        )
        f.write("target_outcomes\tINFO\tRNFL;GCIPL\tExternal neuroretinal endophenotype triangulation\n")
        f.write("sample_size_note\tDOCUMENTED\tN=31434 for both outcomes\tFrom Phase 5.2 OpenGWAS metadata for GCST90014266 and GCST90014267\n")
        f.write(f"runtime\tINFO\t{elapsed_total:.3f}s\tPhase 5.4B runtime\n")

    print("===== Phase 5.4B completed =====")
    print("Wrote:", STATUS_OUT)
    print("Wrote:", QC_OUT)
    print("Wrote:", MANIFEST_OUT)
    print("Wrote:", RUNTIME_OUT)

if __name__ == "__main__":
    main()
