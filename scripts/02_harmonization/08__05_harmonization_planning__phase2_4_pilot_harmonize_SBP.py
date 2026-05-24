#!/usr/bin/env python3
import gzip
import os
from collections import Counter

EXPOSURE_ID = "SBP"
EXPOSURE_FILE = "../../08_standardized_exposures/SBP__ieu-b-38.standardized.tsv.gz"

OUTDIR = "../../09_harmonized_pilot/SBP"
os.makedirs(OUTDIR, exist_ok=True)

OUTCOMES = [
    {
        "outcome_id": "GBS_nonIOPcomponent",
        "index_file": "../../05_harmonization_planning/outcome_indexes/GBS_nonIOP_by_rsid.tsv.gz",
        "match_mode": "rsid",
        "output_file": OUTDIR + "/SBP__GBS_nonIOPcomponent.pilot_harmonized.tsv.gz",
    },
    {
        "outcome_id": "GBS_IOPcomponent",
        "index_file": "../../05_harmonization_planning/outcome_indexes/GBS_IOP_by_rsid.tsv.gz",
        "match_mode": "rsid",
        "output_file": OUTDIR + "/SBP__GBS_IOPcomponent.pilot_harmonized.tsv.gz",
    },
    {
        "outcome_id": "IOPcc",
        "index_file": "../../05_harmonization_planning/outcome_indexes/IOPcc_by_chr_pos_alleles.tsv.gz",
        "match_mode": "chr_pos_allele",
        "output_file": OUTDIR + "/SBP__IOPcc_coordinate_subset.pilot_harmonized.tsv.gz",
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
        x = float(beta)
        y = -x
        return f"{y:.12g}"
    except Exception:
        return "NA"


def load_gbs_index(index_file):
    """
    Load GBS outcome rsID index:
    columns:
    outcome_id, match_key, SNP, chr, pos, effect_allele, other_allele,
    eaf_or_maf, beta, se, pval, pval_for_log10, genome_build, source_dataset
    """
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
    """
    Load IOPcc coordinate-resolvable allele-aware index.

    Direct key:
      outcome match_key = chr:pos:outcome_EA:outcome_OA

    Reverse key:
      outcome reverse_match_key = chr:pos:outcome_OA:outcome_EA

    During exposure scan:
      exposure_key = chr:pos:exposure_EA:exposure_OA

    If exposure_key in direct_lookup:
      aligned
    If exposure_key in reverse_lookup:
      flipped
    """
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


def decide_harmonization(exp_ea, exp_oa, out_ea, out_oa, out_beta):
    """
    Returns:
      action, include_in_main, beta_outcome_harmonized
    """
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


def harmonize_against_outcome(cfg):
    outcome_id = cfg["outcome_id"]
    match_mode = cfg["match_mode"]
    index_file = cfg["index_file"]
    out_file = cfg["output_file"]
    tmp_file = out_file + ".tmp"

    print(f"Loading outcome index: {outcome_id}")

    if match_mode == "rsid":
        gbs_lookup, dup = load_gbs_index(index_file)
        print(f"Loaded {len(gbs_lookup)} rsID outcome records; duplicated outcome keys skipped={dup}")
        direct_lookup = gbs_lookup
        reverse_lookup = None
    elif match_mode == "chr_pos_allele":
        direct_lookup, reverse_lookup, dup_d, dup_r = load_iopcc_index(index_file)
        print(f"Loaded {len(direct_lookup)} direct IOPcc keys and {len(reverse_lookup)} reverse IOPcc keys; duplicated direct={dup_d}; duplicated reverse={dup_r}")
    else:
        raise RuntimeError(f"Unknown match_mode: {match_mode}")

    n_exposure = 0
    n_matched = 0
    n_written = 0
    action_counter = Counter()

    with open_gz_text(EXPOSURE_FILE) as fin, open_gz_write(tmp_file) as fout:
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

                action, include, beta_h = decide_harmonization(
                    exp_ea, exp_oa, rec["ea"], rec["oa"], rec["beta"]
                )

            else:
                if exp_chr in ("", ".", "NA") or exp_pos in ("", ".", "NA") or exp_ea in ("", ".", "NA") or exp_oa in ("", ".", "NA"):
                    continue

                match_key = f"{exp_chr}:{exp_pos}:{exp_ea}:{exp_oa}"

                rec = direct_lookup.get(match_key)
                if rec is not None:
                    action, include, beta_h = decide_harmonization(
                        exp_ea, exp_oa, rec["ea"], rec["oa"], rec["beta"]
                    )
                else:
                    rec = reverse_lookup.get(match_key)
                    if rec is None:
                        continue
                    action, include, beta_h = decide_harmonization(
                        exp_ea, exp_oa, rec["ea"], rec["oa"], rec["beta"]
                    )

            n_matched += 1
            action_counter[action] += 1

            out = [
                EXPOSURE_ID,
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

    n_final = sum(v for k, v in action_counter.items() if not k.endswith("_excluded") and k != "allele_mismatch_excluded")

    return {
        "exposure_id": EXPOSURE_ID,
        "outcome_id": outcome_id,
        "match_mode": match_mode,
        "index_file": index_file,
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
    qc_out = OUTDIR + "/SBP__pilot_harmonization_qc.tsv"

    qc_rows = []

    for cfg in OUTCOMES:
        print()
        print(f"===== Harmonizing SBP against {cfg['outcome_id']} =====")
        result = harmonize_against_outcome(cfg)
        qc_rows.append(result)

    with open(qc_out, "w", encoding="utf-8", newline="\n") as f:
        f.write(
            "exposure_id\toutcome_id\tmatch_mode\tindex_file\toutput_file\t"
            "n_exposure_variants\tn_matched_before_allele_filter\tn_written_rows\t"
            "n_final_include_in_main\tmatch_rate\tfinal_rate\tharmonization_action_counts\n"
        )

        for r in qc_rows:
            action_str = ";".join(f"{k}={v}" for k, v in sorted(r["actions"].items()))
            f.write(
                f"{r['exposure_id']}\t{r['outcome_id']}\t{r['match_mode']}\t"
                f"{r['index_file']}\t{r['output_file']}\t"
                f"{r['n_exposure_variants']}\t{r['n_matched_before_allele_filter']}\t"
                f"{r['n_written_rows']}\t{r['n_final_include_in_main']}\t"
                f"{r['match_rate']:.8f}\t{r['final_rate']:.8f}\t{action_str}\n"
            )

    print()
    print("Wrote pilot QC:", qc_out)
    print("Phase 2.4 SBP pilot harmonization completed.")


if __name__ == "__main__":
    main()
