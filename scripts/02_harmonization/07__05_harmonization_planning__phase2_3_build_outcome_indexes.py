#!/usr/bin/env python3
import gzip
import os

OUTDIR = "../../05_harmonization_planning/outcome_indexes"
PLAN = "../../05_harmonization_planning"

os.makedirs(OUTDIR, exist_ok=True)

OUTCOMES = [
    {
        "outcome_id": "GBS_nonIOPcomponent",
        "source_file": "../../04_standardized_gwas/core_outcomes/GBS_nonIOPcomponent_main_Figure2.standardized.tsv.gz",
        "index_file": OUTDIR + "/GBS_nonIOP_by_rsid.tsv.gz",
        "match_type": "SNP/rsID",
        "kind": "GBS",
        "note": "rsID index for IOP-independent POAG component",
    },
    {
        "outcome_id": "GBS_IOPcomponent",
        "source_file": "../../04_standardized_gwas/core_outcomes/GBS_IOPcomponent_Supplementary_Figure_4b.standardized.tsv.gz",
        "index_file": OUTDIR + "/GBS_IOP_by_rsid.tsv.gz",
        "match_type": "SNP/rsID",
        "kind": "GBS",
        "note": "rsID index for IOP-dependent POAG component",
    },
    {
        "outcome_id": "IOPcc",
        "source_file": "../../04_standardized_gwas/core_outcomes/IOPcc_Supplementary_Figure_4a.standardized.tsv.gz",
        "index_file": OUTDIR + "/IOPcc_by_chr_pos_alleles.tsv.gz",
        "match_type": "chr:pos:effect_allele:other_allele",
        "kind": "IOPCC",
        "note": "Allele-aware coordinate index for IOPcc; do not use position-only final matching",
    },
]

GBS_OUT_HEADER = [
    "outcome_id",
    "match_key",
    "SNP",
    "chr",
    "pos",
    "effect_allele",
    "other_allele",
    "eaf_or_maf",
    "beta",
    "se",
    "pval",
    "pval_for_log10",
    "genome_build",
    "source_dataset",
]

IOPCC_OUT_HEADER = [
    "outcome_id",
    "match_key",
    "reverse_match_key",
    "chr",
    "pos",
    "variant_id",
    "SNP",
    "effect_allele",
    "other_allele",
    "info",
    "beta",
    "se",
    "pval",
    "pval_for_log10",
    "genome_build",
    "source_dataset",
]


def read_header_map(header_line):
    cols = header_line.rstrip("\n").split("\t")
    return {c: i for i, c in enumerate(cols)}, cols


def require_columns(colmap, required, source_file):
    missing = [c for c in required if c not in colmap]
    if missing:
        raise RuntimeError(f"Missing required columns in {source_file}: {missing}")


def build_gbs_index(cfg):
    source_file = cfg["source_file"]
    index_file = cfg["index_file"]
    tmp_file = index_file + ".tmp"

    required = [
        "chr",
        "pos",
        "SNP",
        "effect_allele",
        "other_allele",
        "eaf_or_maf",
        "beta",
        "se",
        "pval",
        "pval_for_log10",
        "genome_build",
        "source_dataset",
    ]

    n_source = 0
    n_index = 0
    n_bad = 0
    n_missing_key = 0
    n_duplicate_key = 0
    seen = set()

    with gzip.open(source_file, "rt", encoding="utf-8", errors="replace") as fin, \
         gzip.open(tmp_file, "wt", encoding="utf-8", newline="\n") as fout:

        header = fin.readline()
        colmap, cols = read_header_map(header)
        require_columns(colmap, required, source_file)

        fout.write("\t".join(GBS_OUT_HEADER) + "\n")

        for line in fin:
            line = line.rstrip("\n")
            if not line:
                continue

            n_source += 1
            parts = line.split("\t")

            if len(parts) != len(cols):
                n_bad += 1
                continue

            snp = parts[colmap["SNP"]]
            if snp in ("", ".", "NA"):
                n_missing_key += 1
                continue

            match_key = snp

            if match_key in seen:
                n_duplicate_key += 1
            else:
                seen.add(match_key)

            out = [
                cfg["outcome_id"],
                match_key,
                snp,
                parts[colmap["chr"]],
                parts[colmap["pos"]],
                parts[colmap["effect_allele"]],
                parts[colmap["other_allele"]],
                parts[colmap["eaf_or_maf"]],
                parts[colmap["beta"]],
                parts[colmap["se"]],
                parts[colmap["pval"]],
                parts[colmap["pval_for_log10"]],
                parts[colmap["genome_build"]],
                parts[colmap["source_dataset"]],
            ]

            fout.write("\t".join(out) + "\n")
            n_index += 1

    os.replace(tmp_file, index_file)

    return {
        "outcome_id": cfg["outcome_id"],
        "source_file": source_file,
        "index_file": index_file,
        "match_type": cfg["match_type"],
        "n_source_rows": n_source,
        "n_index_rows": n_index,
        "n_bad_rows": n_bad,
        "n_missing_key": n_missing_key,
        "n_duplicate_match_keys": n_duplicate_key,
        "status": "INDEXED",
        "note": cfg["note"],
    }


def build_iopcc_index(cfg):
    source_file = cfg["source_file"]
    index_file = cfg["index_file"]
    tmp_file = index_file + ".tmp"

    required = [
        "chr",
        "pos",
        "variant_id",
        "SNP",
        "effect_allele",
        "other_allele",
        "info",
        "beta",
        "se",
        "pval",
        "pval_for_log10",
        "genome_build",
        "source_dataset",
    ]

    n_source = 0
    n_index = 0
    n_bad = 0
    n_missing_key = 0
    n_duplicate_key = 0
    seen = set()

    with gzip.open(source_file, "rt", encoding="utf-8", errors="replace") as fin, \
         gzip.open(tmp_file, "wt", encoding="utf-8", newline="\n") as fout:

        header = fin.readline()
        colmap, cols = read_header_map(header)
        require_columns(colmap, required, source_file)

        fout.write("\t".join(IOPCC_OUT_HEADER) + "\n")

        for line in fin:
            line = line.rstrip("\n")
            if not line:
                continue

            n_source += 1
            parts = line.split("\t")

            if len(parts) != len(cols):
                n_bad += 1
                continue

            chrom = parts[colmap["chr"]]
            pos = parts[colmap["pos"]]
            ea = parts[colmap["effect_allele"]]
            oa = parts[colmap["other_allele"]]

            if chrom in ("", ".", "NA") or pos in ("", ".", "NA") or ea in ("", ".", "NA") or oa in ("", ".", "NA"):
                n_missing_key += 1
                continue

            match_key = f"{chrom}:{pos}:{ea}:{oa}"
            reverse_match_key = f"{chrom}:{pos}:{oa}:{ea}"

            if match_key in seen:
                n_duplicate_key += 1
            else:
                seen.add(match_key)

            out = [
                cfg["outcome_id"],
                match_key,
                reverse_match_key,
                chrom,
                pos,
                parts[colmap["variant_id"]],
                parts[colmap["SNP"]],
                ea,
                oa,
                parts[colmap["info"]],
                parts[colmap["beta"]],
                parts[colmap["se"]],
                parts[colmap["pval"]],
                parts[colmap["pval_for_log10"]],
                parts[colmap["genome_build"]],
                parts[colmap["source_dataset"]],
            ]

            fout.write("\t".join(out) + "\n")
            n_index += 1

    os.replace(tmp_file, index_file)

    return {
        "outcome_id": cfg["outcome_id"],
        "source_file": source_file,
        "index_file": index_file,
        "match_type": cfg["match_type"],
        "n_source_rows": n_source,
        "n_index_rows": n_index,
        "n_bad_rows": n_bad,
        "n_missing_key": n_missing_key,
        "n_duplicate_match_keys": n_duplicate_key,
        "status": "INDEXED",
        "note": cfg["note"],
    }


def main():
    rows = []

    for cfg in OUTCOMES:
        print(f"Building index for {cfg['outcome_id']}")
        if cfg["kind"] == "GBS":
            rows.append(build_gbs_index(cfg))
        elif cfg["kind"] == "IOPCC":
            rows.append(build_iopcc_index(cfg))
        else:
            raise RuntimeError(f"Unknown outcome kind: {cfg['kind']}")
        print(f"Done: {cfg['outcome_id']}")

    manifest = PLAN + "/outcome_index_manifest.tsv"
    qc_summary = OUTDIR + "/phase2_3_outcome_index_qc_summary.tsv"

    cols = [
        "outcome_id",
        "source_file",
        "index_file",
        "match_type",
        "n_source_rows",
        "n_index_rows",
        "n_bad_rows",
        "n_missing_key",
        "n_duplicate_match_keys",
        "status",
        "note",
    ]

    with open(manifest, "w", encoding="utf-8", newline="\n") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(str(r[c]) for c in cols) + "\n")

    with open(qc_summary, "w", encoding="utf-8", newline="\n") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(str(r[c]) for c in cols) + "\n")

    print("")
    print("Wrote manifest:", manifest)
    print("Wrote QC summary:", qc_summary)
    print("Phase 2.3 outcome index construction completed.")


if __name__ == "__main__":
    main()
