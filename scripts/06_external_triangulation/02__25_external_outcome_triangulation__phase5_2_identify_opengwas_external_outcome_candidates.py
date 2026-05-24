#!/usr/bin/env python3
import csv
import json
import os
import re
import time
import urllib.error
import urllib.request

BASE = "../../25_external_outcome_triangulation"
SEARCH_TERMS = os.path.join(BASE, "phase5_2_external_outcome_search_terms.tsv")
RAW_JSON = os.path.join(BASE, "phase5_2_opengwas_gwasinfo_all_raw.json")
CANDIDATES_OUT = os.path.join(BASE, "phase5_2_opengwas_external_outcome_candidates.tsv")
STATUS_OUT = os.path.join(BASE, "phase5_2_status.tsv")
RUNTIME_OUT = os.path.join(BASE, "phase5_2_runtime_log.tsv")

API_URL = "https://api.opengwas.io/api/gwasinfo"
JWT = os.environ.get("OPENGWAS_JWT", "").strip()

def request_json(url):
    req = urllib.request.Request(url)
    req.add_header("Accept", "application/json")
    if JWT:
        req.add_header("Authorization", "Bearer " + JWT)
    with urllib.request.urlopen(req, timeout=900) as resp:
        return json.loads(resp.read().decode("utf-8", errors="replace")), resp.status

def norm(x):
    if x is None:
        return ""
    return str(x).lower()

def pick(d, keys):
    for k in keys:
        if isinstance(d, dict) and k in d and d[k] not in (None, ""):
            return d[k]
    return ""

def flatten_records(obj):
    if isinstance(obj, list):
        return obj
    if isinstance(obj, dict):
        for key in ["data", "results", "studies"]:
            if key in obj and isinstance(obj[key], list):
                return obj[key]
        # OpenGWAS often returns a dict keyed by id.
        if all(isinstance(v, dict) for v in obj.values()):
            rows = []
            for k, v in obj.items():
                vv = dict(v)
                if "id" not in vv:
                    vv["id"] = k
                rows.append(vv)
            return rows
    return []

def read_search_terms():
    rows = []
    with open(SEARCH_TERMS, "r", encoding="utf-8", errors="replace", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            terms = [t.strip().lower() for t in row["search_terms"].split("|") if t.strip()]
            row["_terms"] = terms
            rows.append(row)
    return rows

def score_record(record, terms):
    text_fields = [
        pick(record, ["id", "study_id"]),
        pick(record, ["trait", "phenotype", "name"]),
        pick(record, ["category"]),
        pick(record, ["subcategory"]),
        pick(record, ["author"]),
        pick(record, ["note", "description"]),
        pick(record, ["population"]),
    ]
    text = " | ".join(norm(x) for x in text_fields)

    hits = []
    for term in terms:
        if term in text:
            hits.append(term)

    return hits, text

def main():
    os.makedirs(BASE, exist_ok=True)
    start = time.time()

    status = "NOT_STARTED"
    note = ""
    api_status = "NA"
    records = []

    try:
        print("Fetching OpenGWAS metadata from /gwasinfo ...", flush=True)
        obj, api_status = request_json(API_URL)
        with open(RAW_JSON, "w", encoding="utf-8") as f:
            json.dump(obj, f, ensure_ascii=False)
        records = flatten_records(obj)
        status = "API_METADATA_FETCHED"
        note = f"records={len(records)}"
    except Exception as e:
        status = "API_FAILED"
        note = repr(e)

    search_rows = read_search_terms()
    out_rows = []

    if records:
        for search in search_rows:
            for rec in records:
                hits, text = score_record(rec, search["_terms"])
                if not hits:
                    continue

                dataset_id = pick(rec, ["id", "study_id"])
                trait = pick(rec, ["trait", "phenotype", "name"])
                population = pick(rec, ["population", "sample_ancestry"])
                sample_size = pick(rec, ["sample_size", "n", "nsample"])
                nsnp = pick(rec, ["nsnp", "nvariants"])
                build = pick(rec, ["build", "genome_build", "assembly"])
                year = pick(rec, ["year"])
                author = pick(rec, ["author"])
                category = pick(rec, ["category"])
                subcategory = pick(rec, ["subcategory"])

                # Lightweight prioritization.
                priority_flag = "REVIEW"
                pop_l = norm(population)
                build_l = norm(build)
                trait_l = norm(trait)

                if ("european" in pop_l or "eur" in pop_l or population == ""):
                    priority_flag = "HIGH_REVIEW"
                if search["priority"] == "1" and priority_flag == "HIGH_REVIEW":
                    priority_flag = "CORE_REVIEW"

                out_rows.append({
                    "outcome_group": search["outcome_group"],
                    "outcome_trait": search["outcome_trait"],
                    "priority": search["priority"],
                    "candidate_dataset_id": dataset_id,
                    "trait_from_api": trait,
                    "matched_terms": ";".join(hits),
                    "sample_size_from_api": str(sample_size),
                    "population_from_api": str(population),
                    "nsnp_from_api": str(nsnp),
                    "build_from_api": str(build),
                    "year_from_api": str(year),
                    "author_from_api": str(author),
                    "category_from_api": str(category),
                    "subcategory_from_api": str(subcategory),
                    "review_priority": priority_flag,
                    "decision": "VERIFY_FILES_BEFORE_DOWNLOAD",
                    "note": "Candidate identified by metadata keyword search; manually review phenotype, ancestry, build, and file availability",
                })

    # Deduplicate by outcome_trait + dataset_id.
    seen = set()
    dedup = []
    for r in out_rows:
        key = (r["outcome_trait"], r["candidate_dataset_id"])
        if key in seen:
            continue
        seen.add(key)
        dedup.append(r)

    fieldnames = [
        "outcome_group",
        "outcome_trait",
        "priority",
        "candidate_dataset_id",
        "trait_from_api",
        "matched_terms",
        "sample_size_from_api",
        "population_from_api",
        "nsnp_from_api",
        "build_from_api",
        "year_from_api",
        "author_from_api",
        "category_from_api",
        "subcategory_from_api",
        "review_priority",
        "decision",
        "note",
    ]

    with open(CANDIDATES_OUT, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(dedup)

    elapsed = time.time() - start

    with open(STATUS_OUT, "w", encoding="utf-8", newline="") as f:
        f.write("phase\tstatus\tkey_result\tnote\n")
        f.write(f"Phase 5.2 identify OpenGWAS external outcome candidates\t{status}\t{len(dedup)} candidate rows\t{note}\n")
        f.write(f"api_status\tINFO\t{api_status}\tOpenGWAS /gwasinfo response status if available\n")
        f.write("download_status\tNOT_RUN\t0 files downloaded\tPhase 5.2 only identifies candidate dataset IDs\n")
        f.write("next_step\tTO_DO\tManual review and file availability verification\tProceed to Phase 5.3 only for locked datasets\n")
        f.write(f"runtime\tINFO\t{elapsed:.3f}s\tPhase 5.2 runtime\n")

    with open(RUNTIME_OUT, "w", encoding="utf-8", newline="") as f:
        f.write("phase\telapsed_seconds\telapsed_human\n")
        f.write(f"Phase 5.2\t{elapsed:.3f}\t{elapsed:.1f}s\n")

    print("===== Phase 5.2 completed =====")
    print("Wrote:", STATUS_OUT)
    print("Wrote:", CANDIDATES_OUT)
    print("Wrote:", RUNTIME_OUT)

if __name__ == "__main__":
    main()
