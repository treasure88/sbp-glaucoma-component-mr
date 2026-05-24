#!/usr/bin/env python3
import csv
import gzip
import math
import os
import time
from collections import Counter

RAW = "../../07_exposure_gwas_raw"
OUTDIR = "../../08_standardized_exposures"
PANEL = "../../23_vascular_exposure_panel"
LOGDIR = "../../05_harmonization_planning/exposure_qc_logs"

os.makedirs(OUTDIR, exist_ok=True)
os.makedirs(PANEL, exist_ok=True)
os.makedirs(LOGDIR, exist_ok=True)

MANIFEST_OUT = os.path.join(PANEL, "phase4_4_standardized_new_vascular_exposure_manifest.tsv")
QC_OUT = os.path.join(PANEL, "phase4_4_standardized_new_vascular_exposure_qc_summary.tsv")
STATUS_OUT = os.path.join(PANEL, "phase4_4_status.tsv")
RUNTIME_OUT = os.path.join(PANEL, "phase4_4_runtime_log.tsv")

TARGETS = [
    {
        "exposure_id": "HYPERTENSION",
        "dataset_id": "ukb-d-I9_HYPERTENSION",
        "trait_name": "Hypertensive diseases",
        "population": "European",
        "study_level_n": "361194",
        "input_vcf": os.path.join(RAW, "HYPERTENSION__ukb-d-I9_HYPERTENSION.vcf.gz"),
        "output_file": os.path.join(OUTDIR, "HYPERTENSION__ukb-d-I9_HYPERTENSION.standardized.tsv.gz"),
        "qc_file": os.path.join(LOGDIR, "HYPERTENSION__ukb-d-I9_HYPERTENSION_phase4_4_standardization_qc.tsv"),
        "n_strategy": "FORMAT_SS",
    },
    {
        "exposure_id": "ART_STIFFNESS",
        "dataset_id": "ukb-b-11971",
        "trait_name": "Pulse wave Arterial Stiffness index",
        "population": "European",
        "study_level_n": "151053",
        "input_vcf": os.path.join(RAW, "ART_STIFFNESS__ukb-b-11971.vcf.gz"),
        "output_file": os.path.join(OUTDIR, "ART_STIFFNESS__ukb-b-11971.standardized.tsv.gz"),
        "qc_file": os.path.join(LOGDIR, "ART_STIFFNESS__ukb-b-11971_phase4_4_standardization_qc.tsv"),
        "n_strategy": "STUDY_LEVEL_N",
    },
]

STD_HEADER = [
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
    "exposure_id",
    "trait_name",
    "population",
]

def safe_float(x):
    try:
        if x in ("", ".", "NA", None):
            return None
        return float(x)
    except Exception:
        return None

def pval_from_lp(lp_value):
    lp = safe_float(lp_value)
    if lp is None:
        return "NA"
    if lp < 0:
        return "NA"
    if lp < 300:
        return f"{10 ** (-lp):.8g}"
    # Avoid underflow for extremely small p-values.
    return f"1e-{lp:.6g}"

def parse_info(info):
    out = {}
    if info in ("", ".", "NA"):
        return out
    for part in info.split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            out[k] = v
    return out

def choose_snp(vcf_id, fmap, chrom, pos, ea, oa):
    fmt_id = fmap.get("ID", "")
    if vcf_id.startswith("rs"):
        return vcf_id
    if fmt_id.startswith("rs"):
        return fmt_id
    if vcf_id not in ("", ".", "NA"):
        return vcf_id
    if fmt_id not in ("", ".", "NA"):
        return fmt_id
    return f"{chrom}:{pos}:{ea}:{oa}"

def standardize_one(target):
    t0 = time.time()

    input_vcf = target["input_vcf"]
    output_file = target["output_file"]

    n_input_data_rows = 0
    n_written = 0
    n_bad_or_skipped = 0
    n_missing_eaf = 0
    n_missing_n = 0
    n_duplicate_chr_pos_ref_alt = 0
    max_lp = None
    min_pval_string = "NA"
    format_counter = Counter()
    seen_variant_keys = set()
    first_header = ""

    if not os.path.exists(input_vcf):
        raise FileNotFoundError(input_vcf)

    with gzip.open(input_vcf, "rt", encoding="utf-8", errors="replace", newline="") as fin, \
         gzip.open(output_file, "wt", encoding="utf-8", newline="") as fout:

        writer = csv.writer(fout, delimiter="\t", lineterminator="\n")
        writer.writerow(STD_HEADER)

        for line in fin:
            line = line.rstrip("\n\r")

            if line.startswith("#CHROM"):
                first_header = line
                break

        if not first_header:
            raise RuntimeError(f"No #CHROM header found in {input_vcf}")

        for line in fin:
            line = line.rstrip("\n\r")
            if not line or line.startswith("#"):
                continue

            n_input_data_rows += 1
            parts = line.split("\t")

            if len(parts) != 10:
                n_bad_or_skipped += 1
                continue

            chrom, pos, vcf_id, ref, alt, qual, filt, info, fmt, sample = parts

            if not chrom or not pos or not ref or not alt:
                n_bad_or_skipped += 1
                continue

            # Skip multi-allelic ALT if present; OpenGWAS VCFs used here are expected to be biallelic.
            if "," in alt:
                n_bad_or_skipped += 1
                continue

            fmt_keys = fmt.split(":")
            fmt_vals = sample.split(":")
            fmap = dict(zip(fmt_keys, fmt_vals))
            imap = parse_info(info)
            format_counter[fmt] += 1

            beta = fmap.get("ES", "NA")
            se = fmap.get("SE", "NA")
            lp = fmap.get("LP", "NA")

            if beta in ("", ".", "NA") or se in ("", ".", "NA") or lp in ("", ".", "NA"):
                n_bad_or_skipped += 1
                continue

            lp_float = safe_float(lp)
            if lp_float is not None:
                if max_lp is None or lp_float > max_lp:
                    max_lp = lp_float
                    min_pval_string = pval_from_lp(lp)

            eaf = fmap.get("AF", "")
            if eaf in ("", ".", "NA"):
                eaf = imap.get("AF", "NA")
            if eaf in ("", ".", "NA"):
                n_missing_eaf += 1
                eaf = "NA"

            if target["n_strategy"] == "FORMAT_SS":
                n_value = fmap.get("SS", "")
                if n_value in ("", ".", "NA"):
                    n_missing_n += 1
                    n_value = target["study_level_n"]
            else:
                n_value = target["study_level_n"]
                # Count missing FORMAT SS for transparency, but still write study-level N.
                if fmap.get("SS", "") in ("", ".", "NA"):
                    n_missing_n += 1

            effect_allele = alt
            other_allele = ref
            variant_id = f"{chrom}:{pos}_{effect_allele}_{other_allele}"
            snp = choose_snp(vcf_id, fmap, chrom, pos, effect_allele, other_allele)

            variant_key = (chrom, pos, ref, alt)
            if variant_key in seen_variant_keys:
                n_duplicate_chr_pos_ref_alt += 1
            else:
                seen_variant_keys.add(variant_key)

            pval = pval_from_lp(lp)

            writer.writerow([
                chrom,
                pos,
                snp,
                variant_id,
                effect_allele,
                other_allele,
                beta,
                se,
                pval,
                lp,
                eaf,
                n_value,
                "GRCh37",
                target["dataset_id"],
                target["exposure_id"],
                target["trait_name"],
                target["population"],
            ])

            n_written += 1

    elapsed = time.time() - t0

    format_summary = ";".join(f"{k}={v}" for k, v in sorted(format_counter.items()))

    qc_row = {
        "exposure_id": target["exposure_id"],
        "dataset_id": target["dataset_id"],
        "n_input_data_rows": str(n_input_data_rows),
        "n_written_standardized_rows": str(n_written),
        "n_bad_or_skipped_rows": str(n_bad_or_skipped),
        "n_missing_eaf": str(n_missing_eaf),
        "n_missing_n": str(n_missing_n),
        "n_duplicate_chr_pos_ref_alt": str(n_duplicate_chr_pos_ref_alt),
        "max_minus_log10_p": "NA" if max_lp is None else f"{max_lp:.8g}",
        "min_pval_string": min_pval_string,
        "format_summary": format_summary,
        "n_strategy": target["n_strategy"],
        "status": "STANDARDIZED" if n_written > 0 else "FAILED_NO_ROWS",
    }

    with open(target["qc_file"], "w", encoding="utf-8", newline="") as f:
        fieldnames = list(qc_row.keys())
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerow(qc_row)

    runtime_row = {
        "exposure_id": target["exposure_id"],
        "dataset_id": target["dataset_id"],
        "input_file_size_bytes": str(os.path.getsize(input_vcf)),
        "output_file_size_bytes": str(os.path.getsize(output_file)),
        "elapsed_seconds": f"{elapsed:.3f}",
        "elapsed_human": f"{elapsed:.1f}s",
    }

    manifest_row = {
        "exposure_id": target["exposure_id"],
        "dataset_id": target["dataset_id"],
        "input_vcf": input_vcf,
        "output_file": output_file,
        "qc_file": target["qc_file"],
        "status": qc_row["status"],
    }

    return manifest_row, qc_row, runtime_row

def write_table(path, fieldnames, rows):
    with open(path, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

def main():
    manifest_rows = []
    qc_rows = []
    runtime_rows = []

    for target in TARGETS:
        print(f"Standardizing {target['exposure_id']} ...", flush=True)
        manifest_row, qc_row, runtime_row = standardize_one(target)
        manifest_rows.append(manifest_row)
        qc_rows.append(qc_row)
        runtime_rows.append(runtime_row)

    write_table(
        MANIFEST_OUT,
        ["exposure_id", "dataset_id", "input_vcf", "output_file", "qc_file", "status"],
        manifest_rows
    )

    write_table(
        QC_OUT,
        [
            "exposure_id",
            "dataset_id",
            "n_input_data_rows",
            "n_written_standardized_rows",
            "n_bad_or_skipped_rows",
            "n_missing_eaf",
            "n_missing_n",
            "n_duplicate_chr_pos_ref_alt",
            "max_minus_log10_p",
            "min_pval_string",
            "format_summary",
            "n_strategy",
            "status",
        ],
        qc_rows
    )

    write_table(
        RUNTIME_OUT,
        [
            "exposure_id",
            "dataset_id",
            "input_file_size_bytes",
            "output_file_size_bytes",
            "elapsed_seconds",
            "elapsed_human",
        ],
        runtime_rows
    )

    standardized_count = sum(1 for r in qc_rows if r["status"] == "STANDARDIZED")
    total_written = sum(int(r["n_written_standardized_rows"]) for r in qc_rows)

    with open(STATUS_OUT, "w", encoding="utf-8", newline="") as f:
        f.write("phase\tstatus\tkey_result\tnote\n")
        f.write(
            f"Phase 4.4 standardize newly added vascular exposures\t"
            f"{'PASSED' if standardized_count == len(TARGETS) else 'REVIEW_NEEDED'}\t"
            f"{standardized_count}/{len(TARGETS)} exposures standardized\t"
            f"Total standardized rows={total_written}\n"
        )
        f.write(
            f"HYPERTENSION_instrument_caveat\tDOCUMENTED\tmax_LP_may_be_below_5e-8\t"
            f"Instrument selection under strict p<5e-8 may produce zero SNPs\n"
        )
        f.write(
            f"ART_STIFFNESS_n_caveat\tDOCUMENTED\tstudy_level_N_used\t"
            f"FORMAT SS absent; standardized n column filled with study-level N=151053\n"
        )

    print("===== Phase 4.4 standardization completed =====")
    print("Wrote:", MANIFEST_OUT)
    print("Wrote:", QC_OUT)
    print("Wrote:", STATUS_OUT)
    print("Wrote:", RUNTIME_OUT)

if __name__ == "__main__":
    main()
