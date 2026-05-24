#!/usr/bin/env python3
import csv
import json
from pathlib import Path

GENES = [
    {
        "gene_symbol": "SRSF3",
        "ensembl_gene_id": "ENSG00000112081",
        "phase1_priority": "Tier1A_coloc_SMR_supported_expression_candidate_with_HEIDI_heterogeneity",
        "phase1_summary": "Strong coloc and SMR support; expression direction concordant; HEIDI heterogeneity"
    },
    {
        "gene_symbol": "BICC1",
        "ensembl_gene_id": "ENSG00000122870",
        "phase1_priority": "Tier1A_coloc_SMR_supported_inverse_expression_candidate",
        "phase1_summary": "Strong coloc and SMR support; inverse-expression pattern"
    },
    {
        "gene_symbol": "BMP2",
        "ensembl_gene_id": "ENSG00000125845",
        "phase1_priority": "Tier1B_locus_supported",
        "phase1_summary": "NTG-supported locus; gene assignment less direct"
    },
    {
        "gene_symbol": "TMTC2",
        "ensembl_gene_id": "ENSG00000179104",
        "phase1_priority": "Tier1B_locus_supported",
        "phase1_summary": "NTG-supported locus; gene assignment less direct"
    }
]

API_DIR = Path("Phase2_API_results")

def load_json(path):
    if not path.exists():
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8", errors="replace"))
    except Exception as e:
        return {"_error": str(e)}

def summarize_chembl(gene):
    path = API_DIR / f"{gene}_ChEMBL_target_search.json"
    obj = load_json(path)

    exact = []
    human_hits = []

    for t in obj.get("targets", []) if isinstance(obj, dict) else []:
        tid = t.get("target_chembl_id", "")
        pref = t.get("pref_name", "")
        org = t.get("organism", "")
        typ = t.get("target_type", "")

        if org == "Homo sapiens":
            item = f"{tid}:{pref}:{typ}"
            human_hits.append(item)

            pref_upper = pref.upper()
            gene_upper = gene.upper()
            if gene_upper == pref_upper or gene_upper in pref_upper:
                exact.append(item)

    return {
        "chembl_exact_hits": ";".join(exact) if exact else "none",
        "chembl_human_hit_count": len(human_hits),
        "chembl_top_human_hits": ";".join(human_hits[:8]) if human_hits else "none"
    }

def summarize_dgidb(gene):
    path = API_DIR / f"{gene}_DGIdb_interactions.json"
    obj = load_json(path)

    if "_error" in obj:
        return {
            "dgidb_interaction_count": "API_error",
            "dgidb_interactions": obj["_error"]
        }

    interactions = []
    try:
        nodes = obj.get("data", {}).get("genes", {}).get("nodes", [])
        for node in nodes:
            for inter in node.get("interactions", []) or []:
                drug = inter.get("drug") or {}
                drug_name = drug.get("name", "NA")
                types = inter.get("interactionTypes") or []
                type_str = ",".join([x.get("type", "") for x in types if isinstance(x, dict)])
                sources = inter.get("sources") or []
                source_str = ",".join([x.get("sourceDbName", "") for x in sources if isinstance(x, dict)])
                interactions.append(f"{drug_name}|{type_str}|{source_str}")
    except Exception as e:
        return {
            "dgidb_interaction_count": "parse_error",
            "dgidb_interactions": str(e)
        }

    return {
        "dgidb_interaction_count": len(interactions),
        "dgidb_interactions": ";".join(interactions[:12]) if interactions else "none"
    }

def read_opentargets_fixed():
    path = Path("Phase2_OpenTargets_fixed_summary.tsv")
    out = {}
    if not path.exists():
        return out
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for r in reader:
            out[r["gene_symbol"]] = r
    return out

def classify_gene(gene, phase1_priority, chembl_exact_hits, chembl_human_count, dgidb_count, ot_row):
    ot_drugs = str(ot_row.get("drug_candidate_count", "0"))
    ot_tract = str(ot_row.get("tractability_true", "none"))

    has_dgidb = False
    try:
        has_dgidb = int(dgidb_count) > 0
    except Exception:
        has_dgidb = False

    has_chembl = chembl_exact_hits != "none"
    has_ot_drugs = ot_drugs not in ["0", "NA", "none", ""]

    if gene == "SRSF3":
        return (
            "B_genetically_strong_but_direct_druggability_limited",
            "Strong genetic-expression evidence, but no current known drug candidate or DGIdb interaction; likely intracellular splicing factor."
        )

    if gene == "BICC1":
        return (
            "B_genetically_strong_druggability_uncertain",
            "Strong inverse-expression genetic evidence, but current direct druggability evidence is weak."
        )

    if gene == "BMP2":
        return (
            "C_translationally_druggable_but_gene_assignment_uncertain",
            "Best druggability/translational signal among four genes, but Phase 1 genetic gene assignment is less direct."
        )

    if gene == "TMTC2":
        return (
            "D_locus_supported_low_current_druggability",
            "Currently weak druggability evidence and Tier1B genetic support only."
        )

    return ("unclassified", "manual review needed")

ot = read_opentargets_fixed()
rows = []

for g in GENES:
    gene = g["gene_symbol"]
    chembl = summarize_chembl(gene)
    dgidb = summarize_dgidb(gene)
    ot_row = ot.get(gene, {})

    drug_class, interpretation = classify_gene(
        gene,
        g["phase1_priority"],
        chembl["chembl_exact_hits"],
        chembl["chembl_human_hit_count"],
        dgidb["dgidb_interaction_count"],
        ot_row
    )

    rows.append({
        "gene_symbol": gene,
        "ensembl_gene_id": g["ensembl_gene_id"],
        "phase1_priority": g["phase1_priority"],
        "phase1_summary": g["phase1_summary"],
        "local_pQTL_status": "No local pQTL/proteomics summary statistics found",
        "OpenTargets_status": ot_row.get("OpenTargets_status", "NA"),
        "OpenTargets_tractability_true": ot_row.get("tractability_true", "NA"),
        "OpenTargets_targetClass": ot_row.get("targetClass", "NA"),
        "OpenTargets_drug_candidate_count": ot_row.get("drug_candidate_count", "NA"),
        "ChEMBL_exact_hits": chembl["chembl_exact_hits"],
        "ChEMBL_human_hit_count": chembl["chembl_human_hit_count"],
        "ChEMBL_top_human_hits": chembl["chembl_top_human_hits"],
        "DGIdb_interaction_count": dgidb["dgidb_interaction_count"],
        "DGIdb_interactions": dgidb["dgidb_interactions"],
        "phase2_druggability_class": drug_class,
        "phase2_interpretation": interpretation,
        "phase2_next_step": (
            "Download/query UKB-PPP or other pQTL resources; then run pQTL direction/coloc if protein is available"
        )
    })

fields = [
    "gene_symbol",
    "ensembl_gene_id",
    "phase1_priority",
    "phase1_summary",
    "local_pQTL_status",
    "OpenTargets_status",
    "OpenTargets_tractability_true",
    "OpenTargets_targetClass",
    "OpenTargets_drug_candidate_count",
    "ChEMBL_exact_hits",
    "ChEMBL_human_hit_count",
    "ChEMBL_top_human_hits",
    "DGIdb_interaction_count",
    "DGIdb_interactions",
    "phase2_druggability_class",
    "phase2_interpretation",
    "phase2_next_step"
]

with open("Phase2_druggability_integrated_table.tsv", "w", newline="", encoding="utf-8") as out:
    writer = csv.DictWriter(out, fieldnames=fields, delimiter="\t")
    writer.writeheader()
    for r in rows:
        writer.writerow(r)

print("Wrote Phase2_druggability_integrated_table.tsv")
