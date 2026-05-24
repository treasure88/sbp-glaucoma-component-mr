#!/usr/bin/env python3
import gzip
import os
from collections import Counter

OUTDIR = "../../10_harmonized_datasets"
PLAN = "../../05_harmonization_planning"

os.makedirs(OUTDIR, exist_ok=True)

EXPOSURES = [
    {
        "exposure_id": "SBP",
        "file": "../../08_standardized_exposures/SBP__ieu-b-38.standardized.tsv.gz",
    },
    {
        "exposure_id": "DBP",
        "file": "../../08_standardized_exposures/DBP__ieu-b-39.standardized.tsv.gz",
    },
    {
        "exposure_id": "MIGRAINE",
        "file": "../../08_standardized_exposures/MIGRAINE__ebi-a-GCST90038646.standardized.tsv.gz",
    },
    {
        "exposure_id": "INSOMNIA",
        "file": "../../08_standardized_exposures/INSOMNIA__ebi-a-GCST90018869.standardized.tsv.gz",
    },
    {
        "exposure_id": "CRP",
        "file": "../../08_standardized_exposures/CRP__ebi-a-GCST90018950.standardized.tsv.gz",
    },
    {
        "exposure_id": "BMI",
        "file": "../../08_standardized_exposures/BMI__ieu-a-2.standardized.tsv.gz",
    },
]

OUTCOMES = [
    {
        "outcome_id": "GBS_nonIOPcomponent",
        "index_file": "../../05_harmonization_planning/outcome_indexes/GBS_nonIOP_by_rsid.tsv.gz",
        "match_mode": "rsid",
        "suffix": "GBS_nonIOPcomponent",
    },
    {
        "outcome_id": "GBS_IOPcomponent",
        "index_file": "../../05_harmonization_planning/outcome_indexes/GBS_IOP_by_rsid.tsv.gz",
        "match_mode": "rsid",
        "suffix": "GBS_IOPcomponent",
    },
    {
        "outcome_id": "IOPcc",
        "index_file": "../../05_harmonization_planning/outcome_indexes/IOPcc_by_chr_pos_alleles.tsv.gz",
        "match_mode": "chr_pos_allele",
        "suffix": "IOPcc_coordinate_subset",
    },
]

OUT_HEADER = [
    "exposure_id",
    "outcome_id",
    "match_mode",
    "match_key",
    "chr",
    "pos",
    "SNP",
    "exposure_effect_allele",
    "exposure_other_allele",
    "outcome_effect_allele_original",
    "outcome_other_allele_original",
    "beta_exposure",
    "se_exposure",
    "pval_exposure",
    "minus_log10_p_exposure",
    "eaf_exposure",
    "n_exposure",
    "beta_outcome_original",
    "beta_outcome_harmonized",
    "se_outcome",
    "pval_outcome",
    "pval_for_log10_outcome",
    "harmonization_action",
    "include_in_main",
    "outcome_source_dataset",
]


def open_gz_text(path):
    return gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="")


def open_gz_write(path):
    return gzip.open(path, "wt", encoding="utf-8", newline="\n")


def header_map(header_line):
    cols = header_line.rstrip("\n\r").split("\t")
    return cols, {c: i for i, c in enumerate(cols)}


def is_palindromic(a1, a2):
    a1 = a1.upper()
    a2 = a2.upper()
    if len(a1) != 1 or len(a2) != 1:
        return False
    return {a1, a2} in ({"A", "T"}, {"C", "G"})


def is_indel_or_multibase(a1, a2):
    return len(a1) != 1 or len(a2) != 1


def flip_beta(beta):
    try:
        return f"{-float(beta):.12g}"
    except Exception:
        return "NA"


def decide_harmonization(exp_ea, exp_oa, out_ea, out_oa, out_beta):
    exp_ea = exp_ea.upper()
    exp_oa = exp_oa.upper()
    out_ea = out_ea.upper()
    out_oa = out_oa.upper()

    pal = is_palindromic(exp_ea, exp_oa)
    indel = is_indel_or_multibase(exp_ea, exp_oa) or is_indel_or_multibase(out_ea, out_oa)

    if exp_ea == out_ea and exp_oa == out_oa:
        if pal:
            return "palindromic_aligned_excluded", "NO", "NA"
        if indel:
            return "indel_aligned", "YES", out_beta
        return "aligned", "YES", out_beta

    if exp_ea == out_oa and exp_oa == out_ea:
        if pal:
            return "palindromic_flipped_excluded", "NO", "NA"
        if indel:
            return "indel_flipped", "YES", flip_beta(out_beta)
        return "flipped", "YES", flip_beta(out_beta)

    return "allele_mismatch_excluded", "NO", "NA"


def load_gbs_index(index_file):
    lookup = {}
    dup = 0

    with open_gz_text(index_file) as f:
        header = f.readline()
        cols, cmap = header_map(header)

        for line in f:
            parts = line.rstrip("\n\r").split("\t")
            if len(parts) != len(cols):
                continue

            key = parts[cmap["match_key"]]
            if key in lookup:
                dup += 1
                continue

            lookup[key] = {
                "chr": parts[cmap["chr"]],
                "pos": parts[cmap["pos"]],
                "SNP": parts[cmap["SNP"]],
                "ea": parts[cmap["effect_allele"]],
                "oa": parts[cmap["other_allele"]],
                "beta": parts[cmap["beta"]],
                "se": parts[cmap["se"]],
                "pval": parts[cmap["pval"]],
                "pval_for_log10": parts[cmap["pval_for_log10"]],
                "source_dataset": parts[cmap["source_dataset"]],
            }

    return lookup, dup


def load_iopcc_index(index_file):
    direct_lookup = {}
    reverse_lookup = {}
    dup_direct = 0
    dup_reverse = 0

    with open_gz_text(index_file) as f:
        header = f.readline()
        cols, cmap = header_map(header)

        for line in f:
            parts = line.rstrip("\n\r").split("\t")
            if len(parts) != len(cols):
                continue

            direct_key = parts[cmap["match_key"]]
            reverse_key = parts[cmap["reverse_match_key"]]

            rec = {
                "chr": parts[cmap["chr"]],
                "pos": parts[cmap["pos"]],
                "SNP": parts[cmap["variant_id"]],
                "ea": parts[cmap["effect_allele"]],
                "oa": parts[cmap["other_allele"]],
                "beta": parts[cmap["beta"]],
                "se": parts[cmap["se"]],
                "pval": parts[cmap["pval"]],
                "pval_for_log10": parts[cmap["pval_for_log10"]],
                "source_dataset": parts[cmap["source_dataset"]],
            }

            if direct_key in direct_lookup:
                dup_direct += 1
            else:
                direct_lookup[direct_key] = rec

            if reverse_key in reverse_lookup:
                dup_reverse += 1
            else:
                reverse_lookup[reverse_key] = rec

    return direct_lookup, reverse_lookup, dup_direct, dup_reverse


def harmonize_one_exposure_to_outcome(exposure_cfg, outcome_cfg, direct_lookup, reverse_lookup=None):
    exposure_id = exposure_cfg["exposure_id"]
    exposure_file = exposure_cfg["file"]
    outcome_id = outcome_cfg["outcome_id"]
    match_mode = outcome_cfg["match_mode"]

    exposure_outdir = os.path.join(OUTDIR, exposure_id)
    os.makedirs(exposure_outdir, exist_ok=True)

    out_file = os.path.join(
        exposure_outdir,
        f"{exposure_id}__{outcome_cfg['suffix']}.harmonized.tsv.gz"
    )
    tmp_file = out_file + ".tmp"

    n_exposure = 0
    n_matched = 0
    n_written = 0
    action_counter = Counter()

    with open_gz_text(exposure_file) as fin, open_gz_write(tmp_file) as fout:
        header = fin.readline()
        exp_cols, exp_cmap = header_map(header)

        fout.write("\t".join(OUT_HEADER) + "\n")

        for line in fin:
            parts = line.rstrip("\n\r").split("\t")
            if len(parts) != len(exp_cols):
                continue

            n_exposure += 1

            exp_chr = parts[exp_cmap["chr"]]
            exp_pos = parts[exp_cmap["pos"]]
            exp_snp = parts[exp_cmap["SNP"]]
            exp_ea = parts[exp_cmap["effect_allele"]]
            exp_oa = parts[exp_cmap["other_allele"]]
            exp_beta = parts[exp_cmap["beta"]]
            exp_se = parts[exp_cmap["se"]]
            exp_pval = parts[exp_cmap["pval"]]
            exp_lp = parts[exp_cmap["minus_log10_p"]]
            exp_eaf = parts[exp_cmap["eaf"]]
            exp_n = parts[exp_cmap["n"]]

            if match_mode == "rsid":
                if exp_snp in ("", ".", "NA"):
                    continue
                match_key = exp_snp
                rec = direct_lookup.get(match_key)
                if rec is None:
                    continue

            else:
                if exp_chr in ("", ".", "NA") or exp_pos in ("", ".", "NA") or exp_ea in ("", ".", "NA") or exp_oa in ("", ".", "NA"):
                    continue

                match_key = f"{exp_chr}:{exp_pos}:{exp_ea}:{exp_oa}"

                rec = direct_lookup.get(match_key)
                if rec is None:
                    rec = reverse_lookup.get(match_key)
                    if rec is None:
                        continue

            action, include, beta_h = decide_harmonization(
                exp_ea,
                exp_oa,
                rec["ea"],
                rec["oa"],
                rec["beta"],
            )

            n_matched += 1
            action_counter[action] += 1

            out = [
                exposure_id,
                outcome_id,
                match_mode,
                match_key,
                exp_chr,
                exp_pos,
                exp_snp,
                exp_ea,
                exp_oa,
                rec["ea"],
                rec["oa"],
                exp_beta,
                exp_se,
                exp_pval,
                exp_lp,
                exp_eaf,
                exp_n,
                rec["beta"],
                beta_h,
                rec["se"],
                rec["pval"],
                rec["pval_for_log10"],
                action,
                include,
                rec["source_dataset"],
            ]

            fout.write("\t".join(out) + "\n")
            n_written += 1

    os.replace(tmp_file, out_file)

    n_final = sum(
        v for k, v in action_counter.items()
        if not k.endswith("_excluded") and k != "allele_mismatch_excluded"
    )

    return {
        "exposure_id": exposure_id,
        "outcome_id": outcome_id,
        "match_mode": match_mode,
        "exposure_file": exposure_file,
        "index_file": outcome_cfg["index_file"],
        "output_file": out_file,
        "n_exposure_variants": n_exposure,
        "n_matched_before_allele_filter": n_matched,
        "n_written_rows": n_written,
        "n_final_include_in_main": n_final,
        "match_rate": n_matched / n_exposure if n_exposure else 0,
        "final_rate": n_final / n_exposure if n_exposure else 0,
        "actions": dict(action_counter),
    }


def main():
    all_qc_rows = []

    for outcome_cfg in OUTCOMES:
        print()
        print(f"===== Loading outcome index: {outcome_cfg['outcome_id']} =====")

        if outcome_cfg["match_mode"] == "rsid":
            direct_lookup, dup = load_gbs_index(outcome_cfg["index_file"])
            reverse_lookup = None
            print(f"Loaded {len(direct_lookup)} rsID records; duplicated keys skipped={dup}")
        else:
            direct_lookup, reverse_lookup, dup_d, dup_r = load_iopcc_index(outcome_cfg["index_file"])
            print(f"Loaded {len(direct_lookup)} direct IOPcc keys; {len(reverse_lookup)} reverse keys; duplicated direct={dup_d}; duplicated reverse={dup_r}")

        for exposure_cfg in EXPOSURES:
            print(f"Harmonizing {exposure_cfg['exposure_id']} vs {outcome_cfg['outcome_id']}")
            result = harmonize_one_exposure_to_outcome(
                exposure_cfg,
                outcome_cfg,
                direct_lookup,
                reverse_lookup,
            )
            all_qc_rows.append(result)

    global_qc = os.path.join(OUTDIR, "phase2_5_all_pilot_harmonization_qc.tsv")
    manifest = os.path.join(OUTDIR, "phase2_5_harmonized_dataset_manifest.tsv")

    with open(global_qc, "w", encoding="utf-8", newline="\n") as f:
        f.write(
            "exposure_id\toutcome_id\tmatch_mode\texposure_file\tindex_file\toutput_file\t"
            "n_exposure_variants\tn_matched_before_allele_filter\tn_written_rows\t"
            "n_final_include_in_main\tmatch_rate\tfinal_rate\tharmonization_action_counts\n"
        )

        for r in all_qc_rows:
            action_str = ";".join(f"{k}={v}" for k, v in sorted(r["actions"].items()))
            f.write(
                f"{r['exposure_id']}\t{r['outcome_id']}\t{r['match_mode']}\t"
                f"{r['exposure_file']}\t{r['index_file']}\t{r['output_file']}\t"
                f"{r['n_exposure_variants']}\t{r['n_matched_before_allele_filter']}\t"
                f"{r['n_written_rows']}\t{r['n_final_include_in_main']}\t"
                f"{r['match_rate']:.8f}\t{r['final_rate']:.8f}\t{action_str}\n"
            )

    with open(manifest, "w", encoding="utf-8", newline="\n") as f:
        f.write("exposure_id\toutcome_id\tmatch_mode\tharmonized_file\tstatus\n")
        for r in all_qc_rows:
            f.write(
                f"{r['exposure_id']}\t{r['outcome_id']}\t{r['match_mode']}\t"
                f"{r['output_file']}\tHARMONIZED\n"
            )

    print()
    print("Wrote global QC:", global_qc)
    print("Wrote manifest:", manifest)
    print("Phase 2.5 all pilot exposure harmonization completed.")


if __name__ == "__main__":
    main()
