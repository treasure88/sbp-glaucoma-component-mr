#!/usr/bin/env python3
import csv

VALIDATION_FILE = "Tier1_nonIOP_NTG_direction_validation.tsv"
BEST_NTG_PROXY_FILE = "Tier1_NTG_best_shared_proxy_by_NTG_P.tsv"
GENE_FILE = "Tier1_nearest_protein_coding_genes.tsv"
OUT_FILE = "Tier1_final_priority_table.tsv"

# read exact NTG validation
exact = {}
with open(VALIDATION_FILE, "r", encoding="utf-8", errors="replace") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        exact[row["SNP"]] = row

# read best NTG proxy by NTG P
proxy = {}
with open(BEST_NTG_PROXY_FILE, "r", encoding="utf-8", errors="replace") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        proxy[row["lead_snp"]] = row

# read nearest protein-coding gene
genes = {}
with open(GENE_FILE, "r", encoding="utf-8", errors="replace") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        genes[row["lead_snp"]] = row

def fnum(x):
    try:
        return float(x)
    except Exception:
        return None

def assign_priority(snp, exact_row, proxy_row, gene_row):
    gene = gene_row.get("nearest_protein_coding_gene", "NA")
    relation = gene_row.get("relation_to_lead", "NA")
    dist = fnum(gene_row.get("distance_bp", "NA"))

    exact_status = exact_row.get("validation_status", "NA")
    proxy_ntg_p = fnum(proxy_row.get("NTG_P", "NA"))

    if gene == "BICC1" and relation == "inside_gene" and exact_status == "strong_NTG_support":
        return "Tier1A_highest_priority", "lead_inside_BICC1_exact_same_direction_strong_NTG"

    if snp == "rs12208086" and proxy_ntg_p is not None and proxy_ntg_p < 5e-8:
        return "Tier1A_highest_priority", "lead_missing_but_locus_proxy_same_direction_strong_NTG"

    if proxy_ntg_p is not None and proxy_ntg_p < 5e-8:
        return "Tier1B_strong_locus_support", "proxy_same_direction_strong_NTG_but_gene_assignment_less_direct"

    if exact_status in ["suggestive_NTG_support", "nominal_NTG_support"]:
        return "Tier1C_supportive", "exact_same_direction_but_NTG_support_weaker"

    return "Tier1C_supportive", "requires_additional_eQTL_pQTL_or_coloc"

fields = [
    "lead_snp",
    "candidate_nearest_protein_coding_gene",
    "ensembl_gene_id",
    "gene_distance_bp",
    "relation_to_lead",
    "nonIOP_A1",
    "nonIOP_A2",
    "nonIOP_BETA",
    "nonIOP_P",
    "exact_NTG_available",
    "exact_NTG_orientation",
    "exact_NTG_beta_aligned",
    "exact_NTG_P",
    "exact_NTG_validation_status",
    "best_proxy_snp_by_NTG_P",
    "best_proxy_nonIOP_BETA",
    "best_proxy_nonIOP_P",
    "best_proxy_NTG_beta_aligned",
    "best_proxy_NTG_P",
    "best_proxy_direction",
    "priority_class",
    "priority_reason"
]

with open(OUT_FILE, "w", newline="", encoding="utf-8") as out:
    writer = csv.DictWriter(out, fieldnames=fields, delimiter="\t")
    writer.writeheader()

    for snp in sorted(genes.keys()):
        e = exact.get(snp, {})
        p = proxy.get(snp, {})
        g = genes.get(snp, {})

        priority, reason = assign_priority(snp, e, p, g)

        writer.writerow({
            "lead_snp": snp,
            "candidate_nearest_protein_coding_gene": g.get("nearest_protein_coding_gene", "NA"),
            "ensembl_gene_id": g.get("ensembl_gene_id", "NA"),
            "gene_distance_bp": g.get("distance_bp", "NA"),
            "relation_to_lead": g.get("relation_to_lead", "NA"),
            "nonIOP_A1": e.get("nonIOP_A1", "NA"),
            "nonIOP_A2": e.get("nonIOP_A2", "NA"),
            "nonIOP_BETA": e.get("nonIOP_BETA", "NA"),
            "nonIOP_P": e.get("nonIOP_P", "NA"),
            "exact_NTG_available": "No" if e.get("validation_status") == "missing_exact_rsID" else "Yes",
            "exact_NTG_orientation": e.get("allele_orientation", "NA"),
            "exact_NTG_beta_aligned": e.get("NTG_BETA_aligned_to_nonIOP_A1", "NA"),
            "exact_NTG_P": e.get("NTG_P", "NA"),
            "exact_NTG_validation_status": e.get("validation_status", "NA"),
            "best_proxy_snp_by_NTG_P": p.get("shared_snp", "NA"),
            "best_proxy_nonIOP_BETA": p.get("nonIOP_BETA", "NA"),
            "best_proxy_nonIOP_P": p.get("nonIOP_P", "NA"),
            "best_proxy_NTG_beta_aligned": p.get("NTG_BETA_aligned_to_nonIOP_A1", "NA"),
            "best_proxy_NTG_P": p.get("NTG_P", "NA"),
            "best_proxy_direction": p.get("direction_concordance", "NA"),
            "priority_class": priority,
            "priority_reason": reason
        })

print("Wrote", OUT_FILE)
