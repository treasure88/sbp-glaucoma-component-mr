#!/usr/bin/env python3
import csv

FILES = {
    "SRSF3": "SRSF3_formal_coloc_input.tsv",
    "BICC1": "BICC1_formal_coloc_input.tsv"
}

def fnum(x):
    try:
        if x in ["NA", "", None]:
            return None
        return float(x)
    except Exception:
        return None

def effect_allele_freq(effect_allele, allele_a, allele_b, allele_b_freq):
    fb = fnum(allele_b_freq)
    if fb is None:
        return "NA"

    effect_allele = str(effect_allele).upper()
    allele_a = str(allele_a).upper()
    allele_b = str(allele_b).upper()

    if effect_allele == allele_b:
        return fb
    if effect_allele == allele_a:
        return 1.0 - fb
    return "NA"

for gene, infile in FILES.items():
    out_non = f"{gene}_nonIOP_for_SMR_locus.ma"
    out_ntg = f"{gene}_NTG_for_SMR_locus.ma"

    with open(infile, "r", encoding="utf-8", errors="replace") as f, \
         open(out_non, "w", newline="", encoding="utf-8") as on, \
         open(out_ntg, "w", newline="", encoding="utf-8") as ot:

        reader = csv.DictReader(f, delimiter="\t")

        fields = ["SNP", "A1", "A2", "freq", "b", "se", "p", "N"]
        wn = csv.DictWriter(on, fieldnames=fields, delimiter="\t")
        wt = csv.DictWriter(ot, fieldnames=fields, delimiter="\t")
        wn.writeheader()
        wt.writeheader()

        for r in reader:
            # Use nonIOP A1/A2 as the aligned allele frame.
            a1 = r["nonIOP_A1"]
            a2 = r["nonIOP_A2"]

            freq = effect_allele_freq(
                a1,
                r["eQTL_AF_AlleleA"],
                r["eQTL_AF_AlleleB"],
                r["eQTL_AlleleB_freq"]
            )

            wn.writerow({
                "SNP": r["snp"],
                "A1": a1,
                "A2": a2,
                "freq": freq,
                "b": r["nonIOP_beta"],
                "se": r["nonIOP_se"],
                "p": r["nonIOP_p"],
                "N": "NA"
            })

            wt.writerow({
                "SNP": r["snp"],
                "A1": a1,
                "A2": a2,
                "freq": freq,
                "b": r["NTG_beta_aligned"],
                "se": r["NTG_se"],
                "p": r["NTG_p"],
                "N": "NA"
            })

    print("Wrote", out_non)
    print("Wrote", out_ntg)
