#!/usr/bin/env python3
import csv
import math

NONIOP_FILE = "GBS_nonIOP_no_chr9p21_no_chr14SIX6.txt"
BED_FILE = "Tier1_clean_nonIOP_loci_500kb_windows.bed"
NTG_FILE = "../../../02_raw_downloads/MTAG_NTG_IGGC_STAGE2.tab"

ALL_OUT = "Tier1_NTG_locus_shared_variants.tsv"
BEST_NONIOP_OUT = "Tier1_NTG_best_shared_proxy_by_nonIOP_P.tsv"
BEST_NTG_OUT = "Tier1_NTG_best_shared_proxy_by_NTG_P.tsv"

def to_float(x):
    try:
        return float(x)
    except Exception:
        return None

def sign(x):
    x = to_float(x)
    if x is None:
        return "NA"
    if x > 0:
        return "+"
    if x < 0:
        return "-"
    return "0"

def align_beta(ref_a1, ref_a2, target_a1, target_a2, target_beta):
    target_beta = to_float(target_beta)
    if target_beta is None:
        return "missing_beta", None

    ref_a1 = ref_a1.upper()
    ref_a2 = ref_a2.upper()
    target_a1 = target_a1.upper()
    target_a2 = target_a2.upper()

    if ref_a1 == target_a1 and ref_a2 == target_a2:
        return "same_orientation", target_beta
    if ref_a1 == target_a2 and ref_a2 == target_a1:
        return "flipped_orientation", -target_beta
    return "allele_mismatch", None

# Read windows
windows = []
with open(BED_FILE, "r", encoding="utf-8", errors="replace") as f:
    for line in f:
        if not line.strip():
            continue
        chrom, start, end, lead_snp = line.strip().split()[:4]
        windows.append({
            "lead_snp": lead_snp,
            "chr": chrom,
            "start": int(start),
            "end": int(end)
        })

# Read nonIOP variants in Tier1 windows
non_variants = {}
snp_to_lead = {}

with open(NONIOP_FILE, "r", encoding="utf-8", errors="replace") as f:
    header = f.readline().strip().split()
    idx = {h: i for i, h in enumerate(header)}

    for line in f:
        if not line.strip():
            continue
        p = line.strip().split()
        chrom = p[idx["CHR"]]
        pos = int(float(p[idx["POS"]]))
        snp = p[idx["SNP"]]

        for w in windows:
            if chrom == w["chr"] and w["start"] <= pos <= w["end"]:
                key = (w["lead_snp"], snp)
                non_variants[key] = {
                    "lead_snp": w["lead_snp"],
                    "CHR": chrom,
                    "POS": pos,
                    "SNP": snp,
                    "A1": p[idx["A1"]],
                    "A2": p[idx["A2"]],
                    "BETA": p[idx["BETA"]],
                    "SE": p[idx["SE"]],
                    "P": p[idx["P"]],
                    "MAF": p[idx["MAF"]]
                }
                snp_to_lead.setdefault(snp, set()).add(w["lead_snp"])

print("nonIOP variants in Tier1 windows:", len(non_variants))
print("unique SNPs to search in NTG:", len(snp_to_lead))

# Scan NTG and keep shared variants
shared_rows = []

with open(NTG_FILE, "r", encoding="utf-8", errors="replace") as f:
    header = f.readline().strip().split()
    idx = {h: i for i, h in enumerate(header)}

    for line in f:
        if not line.strip():
            continue
        p = line.strip().split()
        snp = p[idx["SNP"]]
        if snp not in snp_to_lead:
            continue

        for lead_snp in snp_to_lead[snp]:
            n = non_variants[(lead_snp, snp)]

            orientation, ntg_beta_aligned = align_beta(
                n["A1"], n["A2"],
                p[idx["A1"]], p[idx["A2"]],
                p[idx["BETA"]]
            )

            non_beta = to_float(n["BETA"])
            if non_beta is not None and ntg_beta_aligned is not None:
                if non_beta * ntg_beta_aligned > 0:
                    concordance = "same_direction"
                elif non_beta * ntg_beta_aligned < 0:
                    concordance = "opposite_direction"
                else:
                    concordance = "zero_or_unclear"
            else:
                concordance = "NA"

            shared_rows.append({
                "lead_snp": lead_snp,
                "shared_snp": snp,
                "CHR": n["CHR"],
                "POS": n["POS"],
                "nonIOP_A1": n["A1"],
                "nonIOP_A2": n["A2"],
                "nonIOP_BETA": n["BETA"],
                "nonIOP_SE": n["SE"],
                "nonIOP_P": n["P"],
                "NTG_CHR": p[idx["CHR"]],
                "NTG_BP": p[idx["BP"]],
                "NTG_A1": p[idx["A1"]],
                "NTG_A2": p[idx["A2"]],
                "NTG_BETA_raw": p[idx["BETA"]],
                "NTG_BETA_aligned_to_nonIOP_A1": ntg_beta_aligned if ntg_beta_aligned is not None else "NA",
                "NTG_SE": p[idx["SE"]],
                "NTG_P": p[idx["P"]],
                "allele_orientation": orientation,
                "nonIOP_sign": sign(n["BETA"]),
                "NTG_aligned_sign": sign(ntg_beta_aligned),
                "direction_concordance": concordance
            })

fields = [
    "lead_snp", "shared_snp", "CHR", "POS",
    "nonIOP_A1", "nonIOP_A2", "nonIOP_BETA", "nonIOP_SE", "nonIOP_P",
    "NTG_CHR", "NTG_BP", "NTG_A1", "NTG_A2", "NTG_BETA_raw",
    "NTG_BETA_aligned_to_nonIOP_A1", "NTG_SE", "NTG_P",
    "allele_orientation", "nonIOP_sign", "NTG_aligned_sign", "direction_concordance"
]

with open(ALL_OUT, "w", newline="", encoding="utf-8") as out:
    w = csv.DictWriter(out, fieldnames=fields, delimiter="\t")
    w.writeheader()
    for row in shared_rows:
        w.writerow(row)

def best_by_nonIOP_p(rows):
    best = {}
    for r in rows:
        p = to_float(r["nonIOP_P"])
        if p is None:
            continue
        lead = r["lead_snp"]
        if lead not in best or p < to_float(best[lead]["nonIOP_P"]):
            best[lead] = r
    return best

def best_by_ntg_p(rows):
    best = {}
    for r in rows:
        p = to_float(r["NTG_P"])
        if p is None:
            continue
        lead = r["lead_snp"]
        if lead not in best or p < to_float(best[lead]["NTG_P"]):
            best[lead] = r
    return best

for outfile, bestdict in [
    (BEST_NONIOP_OUT, best_by_nonIOP_p(shared_rows)),
    (BEST_NTG_OUT, best_by_ntg_p(shared_rows))
]:
    with open(outfile, "w", newline="", encoding="utf-8") as out:
        w = csv.DictWriter(out, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for lead in sorted(bestdict):
            w.writerow(bestdict[lead])

print("shared variants:", len(shared_rows))
print("Wrote:", ALL_OUT)
print("Wrote:", BEST_NONIOP_OUT)
print("Wrote:", BEST_NTG_OUT)
