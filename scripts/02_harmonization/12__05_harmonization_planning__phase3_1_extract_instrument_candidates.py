#!/usr/bin/env python3
import gzip
import os
import time
import hashlib
from datetime import datetime
from collections import Counter

EXPOSURE_DIR = "../../08_standardized_exposures"
OUTDIR = "../../13_instrument_selection"
CANDIDATE_DIR = os.path.join(OUTDIR, "candidate_instruments")
CLUMP_INPUT_DIR = os.path.join(OUTDIR, "plink_clump_inputs")
PLAN_DIR = "../../05_harmonization_planning"

os.makedirs(OUTDIR, exist_ok=True)
os.makedirs(CANDIDATE_DIR, exist_ok=True)
os.makedirs(CLUMP_INPUT_DIR, exist_ok=True)

P_THRESHOLD = 5e-8
MINUS_LOG10_P_THRESHOLD = 7.301029995663981

EXPOSURES = [
    {
        "exposure_id": "SBP",
        "dataset_id": "ieu-b-38",
        "file": os.path.join(EXPOSURE_DIR, "SBP__ieu-b-38.standardized.tsv.gz"),
        "trait_name": "systolic blood pressure",
    },
    {
        "exposure_id": "DBP",
        "dataset_id": "ieu-b-39",
        "file": os.path.join(EXPOSURE_DIR, "DBP__ieu-b-39.standardized.tsv.gz"),
        "trait_name": "diastolic blood pressure",
    },
    {
        "exposure_id": "MIGRAINE",
        "dataset_id": "ebi-a-GCST90038646",
        "file": os.path.join(EXPOSURE_DIR, "MIGRAINE__ebi-a-GCST90038646.standardized.tsv.gz"),
        "trait_name": "migraine",
    },
    {
        "exposure_id": "INSOMNIA",
        "dataset_id": "ebi-a-GCST90018869",
        "file": os.path.join(EXPOSURE_DIR, "INSOMNIA__ebi-a-GCST90018869.standardized.tsv.gz"),
        "trait_name": "insomnia",
    },
    {
        "exposure_id": "CRP",
        "dataset_id": "ebi-a-GCST90018950",
        "file": os.path.join(EXPOSURE_DIR, "CRP__ebi-a-GCST90018950.standardized.tsv.gz"),
        "trait_name": "C-reactive protein",
    },
    {
        "exposure_id": "BMI",
        "dataset_id": "ieu-a-2",
        "file": os.path.join(EXPOSURE_DIR, "BMI__ieu-a-2.standardized.tsv.gz"),
        "trait_name": "body mass index",
    },
]

CANDIDATE_COLUMNS = [
    "exposure_id",
    "dataset_id",
    "trait_name",
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
    "pval_for_clumping",
    "eaf",
    "n",
    "genome_build",
    "source_dataset",
    "is_palindromic",
    "is_indel_or_multibase",
    "selection_rule",
]

QC_OUT = os.path.join(OUTDIR, "phase3_1_instrument_candidate_qc_summary.tsv")
RUNTIME_OUT = os.path.join(OUTDIR, "phase3_1_runtime_log.tsv")
MANIFEST_OUT = os.path.join(OUTDIR, "phase3_1_candidate_instrument_manifest.tsv")
CLUMP_CONFIG_OUT = os.path.join(OUTDIR, "phase3_1_ld_clumping_config.tsv")
STATUS_OUT = os.path.join(OUTDIR, "phase3_1_status.tsv")
PLAN_MD = os.path.join(OUTDIR, "phase3_1_instrument_selection_and_clumping_plan.md")

def now_str():
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")

def fmt_seconds(x):
    if x < 60:
        return f"{x:.1f}s"
    if x < 3600:
        return f"{x/60:.1f}min"
    return f"{x/3600:.2f}h"

def estimate_runtime_seconds(file_size_bytes):
    mb = file_size_bytes / (1024 * 1024)
    # Candidate extraction does gzip read only and write much smaller output.
    # Conservative Windows estimate: 25 MB/s + overhead.
    return mb / 25.0 + 5.0

def sha256_file(path, block_size=1024 * 1024):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            b = f.read(block_size)
            if not b:
                break
            h.update(b)
    return h.hexdigest()

def open_gz_text(path):
    return gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="")

def open_gz_write(path):
    return gzip.open(path, "wt", encoding="utf-8", newline="\n")

def is_missing(x):
    return x in ("", ".", "NA")

def is_palindromic(a1, a2):
    a1 = a1.upper()
    a2 = a2.upper()
    if len(a1) != 1 or len(a2) != 1:
        return False
    return {a1, a2} in ({"A", "T"}, {"C", "G"})

def is_indel_or_multibase(a1, a2):
    return len(a1) != 1 or len(a2) != 1

def safe_float(x):
    try:
        return float(x)
    except Exception:
        return None

def pval_for_plink_from_lp(lp):
    """
    PLINK clump requires a numeric p-value.
    For extremely small p-values, use 1e-300 floor to avoid underflow/zero.
    Original pval and minus_log10_p are retained in candidate file.
    """
    x = safe_float(lp)
    if x is None:
        return "NA"
    if x >= 300:
        return "1e-300"
    return f"{10 ** (-x):.8e}"

def extract_candidates_one(cfg):
    exposure_id = cfg["exposure_id"]
    dataset_id = cfg["dataset_id"]
    trait_name = cfg["trait_name"]
    f = cfg["file"]

    if not os.path.exists(f):
        raise FileNotFoundError(f)

    in_size = os.path.getsize(f)
    estimated = estimate_runtime_seconds(in_size)

    candidate_out = os.path.join(CANDIDATE_DIR, f"{exposure_id}__{dataset_id}.candidate_instruments.preclump.tsv.gz")
    clump_out = os.path.join(CLUMP_INPUT_DIR, f"{exposure_id}__{dataset_id}.plink_clump_input.tsv")

    tmp_candidate = candidate_out + ".tmp"
    tmp_clump = clump_out + ".tmp"

    print()
    print(f"===== Extracting candidate instruments: {exposure_id} =====")
    print(f"Input: {f}")
    print(f"Input compressed size: {in_size / (1024*1024):.1f} MB")
    print(f"Approx runtime: {fmt_seconds(estimated)}")
    print(f"Start: {now_str()}")

    t0 = time.time()

    n_rows = 0
    n_gws = 0
    n_missing_core = 0
    n_missing_snp = 0
    n_missing_lp = 0
    n_duplicate_snp = 0
    n_palindromic = 0
    n_indel = 0
    n_eaf_missing = 0

    # Deduplicate by SNP, keep row with largest minus_log10_p.
    best_by_snp = {}

    with open_gz_text(f) as fin:
        header = fin.readline().rstrip("\n\r").split("\t")
        cmap = {c: i for i, c in enumerate(header)}

        required = [
            "chr", "pos", "SNP", "variant_id",
            "effect_allele", "other_allele",
            "beta", "se", "pval", "minus_log10_p",
            "eaf", "n", "genome_build", "source_dataset",
        ]

        missing_cols = [c for c in required if c not in cmap]
        if missing_cols:
            raise RuntimeError(f"Missing columns in {f}: {missing_cols}")

        for line in fin:
            parts = line.rstrip("\n\r").split("\t")
            if len(parts) != len(header):
                continue

            n_rows += 1

            snp = parts[cmap["SNP"]]
            lp = parts[cmap["minus_log10_p"]]

            if is_missing(snp):
                n_missing_snp += 1
                continue

            lp_float = safe_float(lp)
            if lp_float is None:
                n_missing_lp += 1
                continue

            if lp_float < MINUS_LOG10_P_THRESHOLD:
                continue

            n_gws += 1

            core_names = ["effect_allele", "other_allele", "beta", "se", "pval", "minus_log10_p"]
            if any(is_missing(parts[cmap[c]]) for c in core_names):
                n_missing_core += 1
                continue

            ea = parts[cmap["effect_allele"]]
            oa = parts[cmap["other_allele"]]

            pal = is_palindromic(ea, oa)
            indel = is_indel_or_multibase(ea, oa)

            if pal:
                n_palindromic += 1
            if indel:
                n_indel += 1
            if is_missing(parts[cmap["eaf"]]):
                n_eaf_missing += 1

            rec = {
                "exposure_id": exposure_id,
                "dataset_id": dataset_id,
                "trait_name": trait_name,
                "chr": parts[cmap["chr"]],
                "pos": parts[cmap["pos"]],
                "SNP": snp,
                "variant_id": parts[cmap["variant_id"]],
                "effect_allele": ea,
                "other_allele": oa,
                "beta": parts[cmap["beta"]],
                "se": parts[cmap["se"]],
                "pval": parts[cmap["pval"]],
                "minus_log10_p": lp,
                "pval_for_clumping": pval_for_plink_from_lp(lp),
                "eaf": parts[cmap["eaf"]],
                "n": parts[cmap["n"]],
                "genome_build": parts[cmap["genome_build"]],
                "source_dataset": parts[cmap["source_dataset"]],
                "is_palindromic": "YES" if pal else "NO",
                "is_indel_or_multibase": "YES" if indel else "NO",
                "selection_rule": f"minus_log10_p>={MINUS_LOG10_P_THRESHOLD:.12g}; equivalent to p<{P_THRESHOLD}",
                "_lp_float": lp_float,
            }

            old = best_by_snp.get(snp)
            if old is not None:
                n_duplicate_snp += 1
                if lp_float > old["_lp_float"]:
                    best_by_snp[snp] = rec
            else:
                best_by_snp[snp] = rec

    n_unique = len(best_by_snp)

    # Sort by descending minus_log10_p, then chr/pos as strings.
    records = sorted(best_by_snp.values(), key=lambda r: (-r["_lp_float"], r["chr"], r["pos"], r["SNP"]))

    with open_gz_write(tmp_candidate) as fout:
        fout.write("\t".join(CANDIDATE_COLUMNS) + "\n")
        for r in records:
            fout.write("\t".join(str(r[c]) for c in CANDIDATE_COLUMNS) + "\n")

    with open(tmp_clump, "w", encoding="utf-8", newline="\n") as fout:
        fout.write("SNP\tP\n")
        for r in records:
            if r["pval_for_clumping"] != "NA":
                fout.write(f"{r['SNP']}\t{r['pval_for_clumping']}\n")

    os.replace(tmp_candidate, candidate_out)
    os.replace(tmp_clump, clump_out)

    elapsed = time.time() - t0

    print(f"End: {now_str()}")
    print(f"Elapsed: {fmt_seconds(elapsed)}")
    print(f"Genome-wide significant rows before dedup: {n_gws}")
    print(f"Unique SNP candidate instruments: {n_unique}")
    print(f"Candidate file: {candidate_out}")
    print(f"PLINK clump input: {clump_out}")

    return {
        "exposure_id": exposure_id,
        "dataset_id": dataset_id,
        "trait_name": trait_name,
        "input_file": f,
        "candidate_file": candidate_out,
        "plink_clump_input": clump_out,
        "input_file_size_bytes": in_size,
        "candidate_file_size_bytes": os.path.getsize(candidate_out),
        "candidate_sha256": sha256_file(candidate_out),
        "n_standardized_rows": n_rows,
        "n_gws_rows_before_dedup": n_gws,
        "n_unique_candidate_snps": n_unique,
        "n_duplicate_snp_records_removed": n_duplicate_snp,
        "n_missing_snp": n_missing_snp,
        "n_missing_lp": n_missing_lp,
        "n_missing_core_among_gws": n_missing_core,
        "n_palindromic_candidate_snps": sum(1 for r in records if r["is_palindromic"] == "YES"),
        "n_indel_or_multibase_candidate_snps": sum(1 for r in records if r["is_indel_or_multibase"] == "YES"),
        "n_missing_eaf_candidate_snps": sum(1 for r in records if r["eaf"] in ("", ".", "NA")),
        "p_threshold": P_THRESHOLD,
        "minus_log10_p_threshold": MINUS_LOG10_P_THRESHOLD,
        "status": "CANDIDATES_EXTRACTED" if n_unique > 0 else "NO_GWS_CANDIDATES",
        "estimated_runtime_seconds": estimated,
        "elapsed_seconds": elapsed,
    }

def write_outputs(results):
    qc_cols = [
        "exposure_id",
        "dataset_id",
        "trait_name",
        "input_file",
        "candidate_file",
        "plink_clump_input",
        "n_standardized_rows",
        "n_gws_rows_before_dedup",
        "n_unique_candidate_snps",
        "n_duplicate_snp_records_removed",
        "n_missing_snp",
        "n_missing_lp",
        "n_missing_core_among_gws",
        "n_palindromic_candidate_snps",
        "n_indel_or_multibase_candidate_snps",
        "n_missing_eaf_candidate_snps",
        "p_threshold",
        "minus_log10_p_threshold",
        "status",
    ]

    manifest_cols = [
        "exposure_id",
        "dataset_id",
        "candidate_file",
        "candidate_sha256",
        "plink_clump_input",
        "n_unique_candidate_snps",
        "status",
    ]

    runtime_cols = [
        "exposure_id",
        "dataset_id",
        "input_file_size_bytes",
        "candidate_file_size_bytes",
        "estimated_runtime_seconds",
        "elapsed_seconds",
        "estimated_runtime_human",
        "elapsed_human",
    ]

    with open(QC_OUT, "w", encoding="utf-8", newline="\n") as f:
        f.write("\t".join(qc_cols) + "\n")
        for r in results:
            f.write("\t".join(str(r[c]) for c in qc_cols) + "\n")

    with open(MANIFEST_OUT, "w", encoding="utf-8", newline="\n") as f:
        f.write("\t".join(manifest_cols) + "\n")
        for r in results:
            f.write("\t".join(str(r[c]) for c in manifest_cols) + "\n")

    with open(RUNTIME_OUT, "w", encoding="utf-8", newline="\n") as f:
        f.write("\t".join(runtime_cols) + "\n")
        for r in results:
            rr = dict(r)
            rr["estimated_runtime_human"] = fmt_seconds(r["estimated_runtime_seconds"])
            rr["elapsed_human"] = fmt_seconds(r["elapsed_seconds"])
            f.write("\t".join(str(rr[c]) for c in runtime_cols) + "\n")

    with open(CLUMP_CONFIG_OUT, "w", encoding="utf-8", newline="\n") as f:
        f.write("parameter\tvalue\tnote\n")
        f.write("genome_build\tGRCh37\tMust match outcomes and standardized exposures\n")
        f.write("ld_reference_ancestry\tEUR\tUse ancestry-matched European reference if possible\n")
        f.write("plink_binary\tplink\tReplace with actual PLINK executable if needed\n")
        f.write("reference_bfile\t/path/to/GRCh37_EUR_reference_panel\tPlaceholder; must be replaced before clumping\n")
        f.write("clump_p1\t5e-8\tIndex SNP threshold\n")
        f.write("clump_p2\t1\tSecondary SNP threshold\n")
        f.write("clump_r2\t0.001\tStrict MR clumping threshold\n")
        f.write("clump_kb\t10000\t10 Mb window\n")
        f.write("pvalue_column\tP\tPLINK clump input column\n")
        f.write("snp_column\tSNP\tPLINK clump input column\n")
        f.write("important_note\tDO_NOT_RUN_UNTIL_REFERENCE_PANEL_SET\tSet reference_bfile first\n")

    n_ok = sum(1 for r in results if r["status"] == "CANDIDATES_EXTRACTED")
    n_no_gws = sum(1 for r in results if r["status"] == "NO_GWS_CANDIDATES")
    total_candidates = sum(r["n_unique_candidate_snps"] for r in results)
    total_elapsed = sum(r["elapsed_seconds"] for r in results)

    with open(STATUS_OUT, "w", encoding="utf-8", newline="\n") as f:
        f.write("phase\tstatus\tkey_result\tnote\n")
        phase_status = "PASSED" if n_ok == len(results) else "REVIEW_NEEDED"
        f.write(f"Phase 3.1 instrument candidate extraction and LD clumping plan\t{phase_status}\t{n_ok}/{len(results)} exposures have genome-wide significant candidates\tCandidate extraction completed using p<5e-8 equivalent threshold\n")
        f.write(f"candidate_exposure_count\tINFO\t{n_ok}\tExposures with candidate instruments\n")
        f.write(f"no_gws_candidate_count\tINFO\t{n_no_gws}\tExposures with no genome-wide significant candidates\n")
        f.write(f"total_unique_candidate_snps\tINFO\t{total_candidates}\tBefore LD clumping\n")
        f.write(f"total_elapsed_runtime\tINFO\t{fmt_seconds(total_elapsed)}\tSum of exposure extraction times\n")
        f.write("ld_clumping_status\tPLANNED_NOT_RUN\tReference panel not configured\tSet reference_bfile in phase3_1_ld_clumping_config.tsv before clumping\n")
        f.write("selection_rule\tDOCUMENTED\tminus_log10_p>=7.30102999566\tEquivalent to p<5e-8\n")

    with open(PLAN_MD, "w", encoding="utf-8", newline="\n") as f:
        f.write("# Phase 3.1 Instrument Selection and LD Clumping Plan\n\n")
        f.write("## Current status\n\n")
        f.write("Genome-wide significant candidate instruments were extracted from standardized exposure GWAS files.\n\n")
        f.write("## Instrument selection rule\n\n")
        f.write("- Use exposure GWAS only.\n")
        f.write("- Do not select instruments based on outcome availability.\n")
        f.write("- Threshold: p < 5e-8.\n")
        f.write("- Operational threshold: minus_log10_p >= 7.30102999566.\n")
        f.write("- Candidate SNPs are deduplicated by rsID, keeping the smallest p-value / largest minus_log10_p.\n\n")
        f.write("## LD clumping plan\n\n")
        f.write("Default clumping parameters:\n\n")
        f.write("- clump_p1 = 5e-8\n")
        f.write("- clump_p2 = 1\n")
        f.write("- clump_r2 = 0.001\n")
        f.write("- clump_kb = 10000\n")
        f.write("- reference panel = GRCh37 European ancestry LD reference\n\n")
        f.write("## PLINK command template\n\n")
        f.write("```bash\n")
        f.write("plink \\\n")
        f.write("  --bfile /path/to/GRCh37_EUR_reference_panel \\\n")
        f.write("  --clump ../../13_instrument_selection/plink_clump_inputs/<EXPOSURE>.plink_clump_input.tsv \\\n")
        f.write("  --clump-snp-field SNP \\\n")
        f.write("  --clump-field P \\\n")
        f.write("  --clump-p1 5e-8 \\\n")
        f.write("  --clump-p2 1 \\\n")
        f.write("  --clump-r2 0.001 \\\n")
        f.write("  --clump-kb 10000 \\\n")
        f.write("  --out ../../13_instrument_selection/clumped_instruments/<EXPOSURE>\n")
        f.write("```\n\n")
        f.write("## Required before running clumping\n\n")
        f.write("1. Provide a GRCh37 LD reference panel in PLINK bed/bim/fam format.\n")
        f.write("2. Confirm reference ancestry is appropriate for the exposure/outcome data.\n")
        f.write("3. Confirm chromosome naming is compatible with candidate SNPs.\n")
        f.write("4. Confirm PLINK is installed and callable from Git Bash.\n\n")
        f.write("## Do not do yet\n\n")
        f.write("- Do not run MR before LD clumping.\n")
        f.write("- Do not use outcome matching rate to select instruments.\n")
        f.write("- Do not include palindromic variants that were excluded in harmonized analysis-ready files.\n\n")
        f.write("## Next phase\n\n")
        f.write("Phase 3.2 should run LD clumping after the reference panel is configured.\n")

def main():
    print("===== Phase 3.1 instrument candidate extraction =====")
    print(f"Start: {now_str()}")
    print(f"P threshold: {P_THRESHOLD}")
    print(f"minus_log10_p threshold: {MINUS_LOG10_P_THRESHOLD}")
    print(f"Output directory: {OUTDIR}")

    total_est = 0
    for cfg in EXPOSURES:
        if os.path.exists(cfg["file"]):
            total_est += estimate_runtime_seconds(os.path.getsize(cfg["file"]))

    print(f"Approx total runtime: {fmt_seconds(total_est)}")
    print("Note: actual runtime depends on gzip throughput, disk speed, and Windows antivirus scanning.")

    t_all = time.time()

    results = []
    for cfg in EXPOSURES:
        results.append(extract_candidates_one(cfg))

    write_outputs(results)

    elapsed = time.time() - t_all

    print()
    print("===== Phase 3.1 completed =====")
    print(f"End: {now_str()}")
    print(f"Total elapsed: {fmt_seconds(elapsed)}")
    print(f"Wrote QC summary: {QC_OUT}")
    print(f"Wrote manifest: {MANIFEST_OUT}")
    print(f"Wrote runtime log: {RUNTIME_OUT}")
    print(f"Wrote clumping config: {CLUMP_CONFIG_OUT}")
    print(f"Wrote status: {STATUS_OUT}")
    print(f"Wrote plan: {PLAN_MD}")

if __name__ == "__main__":
    main()
