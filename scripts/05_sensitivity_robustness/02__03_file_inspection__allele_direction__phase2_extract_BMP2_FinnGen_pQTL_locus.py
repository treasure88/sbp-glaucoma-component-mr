#!/usr/bin/env python3
import csv
import gzip
import math
from pathlib import Path

IN = Path("Phase2_free_pQTL_FinnGen/Somascan/pQTL/BMP2__SomaScan_Batch2_seq.15666.21.txt.gz")
OUT = Path("Phase2_free_pQTL_FinnGen/BMP2_FinnGen_SomaScan_locus_pQTL.tsv")
LOG = Path("Phase2_free_pQTL_FinnGen/BMP2_FinnGen_SomaScan_locus_pQTL_extract_log.txt")

CHR = "20"
START = 6007717
END = 7007717

if not IN.exists():
    raise SystemExit(f"Missing input: {IN}")

def get(row, names):
    for n in names:
        if n in row:
            return row[n]
    lower = {k.lower(): k for k in row.keys()}
    for n in names:
        if n.lower() in lower:
            return row[lower[n.lower()]]
    return ""

def to_float(x):
    try:
        return float(x)
    except Exception:
        return None

n = 0
best = None

with gzip.open(IN, "rt", encoding="utf-8", errors="replace") as f, \
     open(OUT, "w", newline="", encoding="utf-8") as out:

    reader = csv.DictReader(f, delimiter="\t")
    fields = ["platform", "gene_symbol", "probe"] + reader.fieldnames
    writer = csv.DictWriter(out, fieldnames=fields, delimiter="\t")
    writer.writeheader()

    for r in reader:
        chrom = str(get(r, ["CHR", "CHROM", "chr", "chrom"])).replace("chr", "").replace("CHR", "")
        pos_s = get(r, ["POS", "BP", "pos", "position"])

        try:
            pos = int(float(pos_s))
        except Exception:
            continue

        if chrom == CHR and START <= pos <= END:
            row = {
                "platform": "FinnGen_SomaScan",
                "gene_symbol": "BMP2",
                "probe": "seq.15666.21",
            }
            row.update(r)
            writer.writerow(row)
            n += 1

            log10p = to_float(get(r, ["log10_P", "LOG10P", "LOGP", "logp"]))
            pval = to_float(get(r, ["P", "p", "PVAL", "pval"]))

            if log10p is None and pval is not None and pval > 0:
                log10p = -math.log10(pval)

            if log10p is not None:
                if best is None or log10p > best["log10p"]:
                    best = {
                        "log10p": log10p,
                        "p": pval,
                        "id": get(r, ["ID", "SNP", "rsid", "variant"]),
                        "pos": pos,
                        "ref": get(r, ["REF", "A2", "OtherAllele"]),
                        "alt": get(r, ["ALT", "A1", "EffectAllele"]),
                        "beta": get(r, ["BETA", "beta", "b"]),
                        "se": get(r, ["SE", "se"]),
                    }

with open(LOG, "w", encoding="utf-8") as log:
    log.write(f"input={IN}\n")
    log.write(f"window=chr{CHR}:{START}-{END}\n")
    log.write(f"n_locus_rows={n}\n")
    if best:
        log.write(f"best_variant={best['id']}\n")
        log.write(f"best_pos={best['pos']}\n")
        log.write(f"best_REF={best['ref']}\n")
        log.write(f"best_ALT_effect_allele={best['alt']}\n")
        log.write(f"best_beta={best['beta']}\n")
        log.write(f"best_se={best['se']}\n")
        log.write(f"best_P={best['p']}\n")
        log.write(f"best_log10_P={best['log10p']}\n")
    else:
        log.write("best_variant=NA\n")
        log.write("best_P=NA\n")
        log.write("best_log10_P=NA\n")

print("Wrote", OUT)
print("Wrote", LOG)
