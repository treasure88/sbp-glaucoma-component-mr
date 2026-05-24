#!/usr/bin/env python3
import gzip
import math
import os
from collections import Counter

OUTDIR = "../../08_standardized_exposures"
LOGDIR = "../../05_harmonization_planning/exposure_qc_logs"

os.makedirs(OUTDIR, exist_ok=True)
os.makedirs(LOGDIR, exist_ok=True)

EXPOSURES = [
    {
        "exposure_id": "SBP",
        "dataset_id": "ieu-b-38",
        "trait_name": "systolic blood pressure",
        "population": "Mixed/European-dominant ICBP; verify from source metadata if needed",
        "input": "../../07_exposure_gwas_raw/SBP__ieu-b-38.vcf.gz",
        "study_n": "757601",
        "n_rule": "FORMAT_SS",
    },
    {
        "exposure_id": "DBP",
        "dataset_id": "ieu-b-39",
        "trait_name": "diastolic blood pressure",
        "population": "Mixed/European-dominant ICBP; verify from source metadata if needed",
        "input": "../../07_exposure_gwas_raw/DBP__ieu-b-39.vcf.gz",
        "study_n": "757601",
        "n_rule": "FORMAT_SS",
    },
    {
        "exposure_id": "MIGRAINE",
        "dataset_id": "ebi-a-GCST90038646",
        "trait_name": "migraine",
        "population": "GWAS Catalog/OpenGWAS migraine; ancestry/population to verify from source paper",
        "input": "../../07_exposure_gwas_raw/MIGRAINE__ebi-a-GCST90038646.vcf.gz",
        "study_n": "484598",
        "n_rule": "STUDY_LEVEL_N",
    },
    {
        "exposure_id": "INSOMNIA",
        "dataset_id": "ebi-a-GCST90018869",
        "trait_name": "insomnia",
        "population": "GWAS Catalog/OpenGWAS insomnia; ancestry/population to verify from source paper",
        "input": "../../07_exposure_gwas_raw/INSOMNIA__ebi-a-GCST90018869.vcf.gz",
        "study_n": "486627",
        "n_rule": "STUDY_LEVEL_N",
    },
    {
        "exposure_id": "CRP",
        "dataset_id": "ebi-a-GCST90018950",
        "trait_name": "C-reactive protein",
        "population": "GWAS Catalog/OpenGWAS C-reactive protein; ancestry/population to verify from source paper",
        "input": "../../07_exposure_gwas_raw/CRP__ebi-a-GCST90018950.vcf.gz",
        "study_n": "353466",
        "n_rule": "STUDY_LEVEL_N",
    },
    {
        "exposure_id": "BMI",
        "dataset_id": "ieu-a-2",
        "trait_name": "body mass index",
        "population": "GIANT BMI; ancestry/population to verify from source paper",
        "input": "../../07_exposure_gwas_raw/BMI__ieu-a-2.vcf.gz",
        "study_n": "339224",
        "n_rule": "FORMAT_SS",
    },
]

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
    "exposure_id",
    "trait_name",
    "population",
]


def parse_info_af(info):
    for item in info.split(";"):
        if item.startswith("AF="):
            value = item.split("=", 1)[1]
            if value not in ("", "."):
                return value
    return "NA"


def pval_from_lp(lp_str):
    """
    Convert LP = -log10(P) to a p-value string without floating-point underflow.

    Example:
    LP=1419.91 -> approximately 1.2302688e-1420
    """
    if lp_str in ("", ".", "NA"):
        return "NA"

    try:
        lp = float(lp_str)
    except Exception:
        return "NA"

    if not math.isfinite(lp) or lp < 0:
        return "NA"

    if lp == 0:
        return "1"

    k = math.floor(lp)
    frac = lp - k

    if abs(frac) < 1e-15:
        return f"1e-{k}"

    mantissa = 10 ** (-frac)
    mantissa *= 10.0
    exponent = -(k + 1)

    return f"{mantissa:.8g}e{exponent}"


def is_palindromic(a1, a2):
    a1 = a1.upper()
    a2 = a2.upper()
    if len(a1) != 1 or len(a2) != 1:
        return False
    return {a1, a2} in ({"A", "T"}, {"C", "G"})


def standardize_one(cfg):
    exposure_id = cfg["exposure_id"]
    dataset_id = cfg["dataset_id"]
    vcf_path = cfg["input"]
    out_path = os.path.join(OUTDIR, f"{exposure_id}__{dataset_id}.standardized.tsv.gz")
    tmp_path = out_path + ".tmp"
    qc_path = os.path.join(LOGDIR, f"{exposure_id}__{dataset_id}_phase2_2_standardization_qc.tsv")

    if not os.path.exists(vcf_path):
        raise FileNotFoundError(vcf_path)

    n_in_data = 0
    n_written = 0
    n_bad = 0
    n_missing_eaf = 0
    n_missing_n = 0
    n_palindromic = 0
    n_duplicate_chr_pos_ref_alt = 0
    n_pval_na = 0
    n_pval_zero_string = 0
    max_lp = None
    min_pval_string = "NA"

    format_counter = Counter()
    seen = set()

    genome_build = "GRCh37"
    sample_meta = "NA"

    with gzip.open(vcf_path, "rt", encoding="utf-8", errors="replace") as fin, \
         gzip.open(tmp_path, "wt", encoding="utf-8", newline="\n") as fout:

        fout.write("\t".join(OUT_HEADER) + "\n")

        for line in fin:
            line = line.rstrip("\n")

            if line.startswith("##SAMPLE="):
                sample_meta = line

            if line.startswith("##contig=<ID=1"):
                if "HG19" in line or "GRCh37" in line:
                    genome_build = "GRCh37"
                else:
                    genome_build = "UNKNOWN"

            if line.startswith("#"):
                continue

            n_in_data += 1
            parts = line.split("\t")

            if len(parts) != 10:
                n_bad += 1
                continue

            chrom, pos, snp, ref, alt, qual, filt, info, fmt, sample = parts
            fmt_keys = fmt.split(":")
            sample_vals = sample.split(":")
            format_counter[fmt] += 1

            if len(fmt_keys) != len(sample_vals):
                n_bad += 1
                continue

            d = dict(zip(fmt_keys, sample_vals))

            required = ["ES", "SE", "LP", "ID"]
            if any(k not in d for k in required):
                n_bad += 1
                continue

            beta = d.get("ES", "NA")
            se = d.get("SE", "NA")
            lp = d.get("LP", "NA")
            pval = pval_from_lp(lp)

            if pval == "NA":
                n_pval_na += 1
            if pval == "0":
                n_pval_zero_string += 1

            try:
                lp_float = float(lp)
                if math.isfinite(lp_float):
                    if max_lp is None or lp_float > max_lp:
                        max_lp = lp_float
                        min_pval_string = pval_from_lp(lp)
            except Exception:
                pass

            # EAF rule:
            # 1. Prefer FORMAT AF.
            # 2. If FORMAT AF missing, try INFO AF.
            # 3. If still missing, NA.
            if "AF" in d and d["AF"] not in ("", "."):
                eaf = d["AF"]
            else:
                eaf = parse_info_af(info)

            if eaf == "NA":
                n_missing_eaf += 1

            # N rule:
            # SBP/DBP/BMI use FORMAT SS when available.
            # MIGRAINE/INSOMNIA/CRP use study-level N.
            if cfg["n_rule"] == "FORMAT_SS":
                n_val = d.get("SS", "NA")
                if n_val in ("", ".", "NA"):
                    n_val = cfg["study_n"]
            else:
                n_val = cfg["study_n"]

            if n_val in ("", ".", "NA"):
                n_missing_n += 1

            key = (chrom, pos, ref, alt)
            if key in seen:
                n_duplicate_chr_pos_ref_alt += 1
            else:
                seen.add(key)

            if is_palindromic(ref, alt):
                n_palindromic += 1

            variant_id = f"{chrom}:{pos}_{alt}_{ref}"

            out_row = [
                chrom,
                pos,
                snp,
                variant_id,
                alt,
                ref,
                beta,
                se,
                pval,
                lp,
                eaf,
                n_val,
                genome_build,
                dataset_id,
                exposure_id,
                cfg["trait_name"],
                cfg["population"],
            ]

            fout.write("\t".join(out_row) + "\n")
            n_written += 1

    os.replace(tmp_path, out_path)

    format_summary = ";".join(f"{k}={v}" for k, v in sorted(format_counter.items()))

    with open(qc_path, "w", encoding="utf-8", newline="\n") as q:
        q.write("metric\tvalue\n")
        q.write(f"exposure_id\t{exposure_id}\n")
        q.write(f"dataset_id\t{dataset_id}\n")
        q.write(f"input_vcf\t{vcf_path}\n")
        q.write(f"output_file\t{out_path}\n")
        q.write(f"genome_build\t{genome_build}\n")
        q.write(f"n_input_data_rows\t{n_in_data}\n")
        q.write(f"n_written_standardized_rows\t{n_written}\n")
        q.write(f"n_bad_or_skipped_rows\t{n_bad}\n")
        q.write(f"n_missing_eaf\t{n_missing_eaf}\n")
        q.write(f"n_missing_n\t{n_missing_n}\n")
        q.write(f"n_palindromic_possible\t{n_palindromic}\n")
        q.write(f"n_duplicate_chr_pos_ref_alt\t{n_duplicate_chr_pos_ref_alt}\n")
        q.write(f"n_pval_na\t{n_pval_na}\n")
        q.write(f"n_pval_zero_string\t{n_pval_zero_string}\n")
        q.write(f"max_minus_log10_p\t{max_lp if max_lp is not None else 'NA'}\n")
        q.write(f"min_pval_string\t{min_pval_string}\n")
        q.write(f"format_summary\t{format_summary}\n")
        q.write(f"n_rule\t{cfg['n_rule']}\n")
        q.write(f"study_n\t{cfg['study_n']}\n")
        q.write(f"sample_meta\t{sample_meta}\n")

    return {
        "exposure_id": exposure_id,
        "dataset_id": dataset_id,
        "input_vcf": vcf_path,
        "output_file": out_path,
        "qc_file": qc_path,
        "n_input_data_rows": n_in_data,
        "n_written_standardized_rows": n_written,
        "n_bad_or_skipped_rows": n_bad,
        "n_missing_eaf": n_missing_eaf,
        "n_missing_n": n_missing_n,
        "n_duplicate_chr_pos_ref_alt": n_duplicate_chr_pos_ref_alt,
        "max_minus_log10_p": max_lp if max_lp is not None else "NA",
        "min_pval_string": min_pval_string,
        "format_summary": format_summary,
        "status": "STANDARDIZED",
    }


def main():
    manifest_path = os.path.join(OUTDIR, "phase2_2_standardized_exposure_manifest.tsv")
    summary_qc_path = os.path.join(OUTDIR, "phase2_2_standardized_exposure_qc_summary.tsv")

    rows = []

    for cfg in EXPOSURES:
        print(f"Standardizing {cfg['exposure_id']} from {cfg['input']}")
        rows.append(standardize_one(cfg))
        print(f"Done: {cfg['exposure_id']}")

    manifest_cols = [
        "exposure_id",
        "dataset_id",
        "input_vcf",
        "output_file",
        "qc_file",
        "status",
    ]

    qc_cols = [
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
        "status",
    ]

    with open(manifest_path, "w", encoding="utf-8", newline="\n") as f:
        f.write("\t".join(manifest_cols) + "\n")
        for r in rows:
            f.write("\t".join(str(r[c]) for c in manifest_cols) + "\n")

    with open(summary_qc_path, "w", encoding="utf-8", newline="\n") as f:
        f.write("\t".join(qc_cols) + "\n")
        for r in rows:
            f.write("\t".join(str(r[c]) for c in qc_cols) + "\n")

    print("")
    print("Wrote manifest:", manifest_path)
    print("Wrote QC summary:", summary_qc_path)
    print("")
    print("Phase 2.2 exposure standardization completed.")


if __name__ == "__main__":
    main()
