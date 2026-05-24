#!/usr/bin/env python3
import csv
import html.parser
import json
import os
import re
import time
import urllib.request
import urllib.error

BASE = "../../31_iop_poag_ntg_htg_source_reidentification"
os.makedirs(BASE, exist_ok=True)

GWAS_CATALOG_OUT = os.path.join(BASE, "phase5_8A_gwas_catalog_ftp_verification.tsv")
ZENODO_OUT = os.path.join(BASE, "phase5_8A_zenodo_verification.tsv")
SOURCE_LOCK_OUT = os.path.join(BASE, "phase5_8A_iop_poag_ntg_htg_source_lock.tsv")
STATUS_OUT = os.path.join(BASE, "phase5_8A_status.tsv")
RUNTIME_OUT = os.path.join(BASE, "phase5_8A_runtime_log.tsv")

GWAS_TARGETS = [
    {
        "outcome_trait": "IOP",
        "gcst": "GCST004074",
        "trait_name": "Intraocular pressure",
        "priority": "1",
    },
    {
        "outcome_trait": "POAG",
        "gcst": "GCST90011766",
        "trait_name": "Primary open-angle glaucoma",
        "priority": "1",
    },
]

ZENODO_TARGETS = [
    {
        "outcome_trait": "NTG",
        "zenodo_record": "14010557",
        "trait_name": "Normal tension glaucoma MTAG / IGGC stage 2",
        "priority": "1",
    }
]

MANUAL_TARGETS = [
    {
        "outcome_trait": "NTG",
        "trait_name": "Normotensive glaucoma",
        "source_route": "FinnGen",
        "source_id": "H7_GLAUCOMA_NTG",
        "priority": "1",
        "source_status": "ACCESS_REQUIRED",
        "recommended_action": "USE_AS_BACKUP_OR_SECONDARY_AFTER_ACCESS",
        "note": "FinnGen summary statistics require access request and emailed download instructions",
    },
    {
        "outcome_trait": "POAG",
        "trait_name": "Primary open-angle glaucoma, strict",
        "source_route": "FinnGen",
        "source_id": "H7_GLAUCOMA_POAG",
        "priority": "2",
        "source_status": "ACCESS_REQUIRED",
        "recommended_action": "USE_AS_BACKUP_OR_SECONDARY_AFTER_ACCESS",
        "note": "FinnGen summary statistics require access request and emailed download instructions",
    },
    {
        "outcome_trait": "HTG",
        "trait_name": "High-tension glaucoma",
        "source_route": "publication_or_consortium",
        "source_id": "Gharahkhani_HTG_GWAS_or_IGGC",
        "priority": "1",
        "source_status": "TO_REIDENTIFY_OR_REQUEST",
        "recommended_action": "SEARCH_PUBLIC_SUMSTATS_OR_REQUEST_FROM_AUTHORS",
        "note": "Nature Communications 2024 used a previously published HTG GWAS; direct download route not locked yet",
    },
]

class LinkParser(html.parser.HTMLParser):
    def __init__(self):
        super().__init__()
        self.links = []
    def handle_starttag(self, tag, attrs):
        if tag.lower() != "a":
            return
        d = dict(attrs)
        href = d.get("href")
        if href:
            self.links.append(href)

def fetch_text(url, timeout=120):
    req = urllib.request.Request(url)
    req.add_header("User-Agent", "Mozilla/5.0")
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return resp.read().decode("utf-8", errors="replace"), resp.status

def fetch_json(url, timeout=120):
    req = urllib.request.Request(url)
    req.add_header("Accept", "application/json")
    req.add_header("User-Agent", "Mozilla/5.0")
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return json.loads(resp.read().decode("utf-8", errors="replace")), resp.status

def gcst_group_url(gcst):
    m = re.search(r"GCST(\d+)", gcst)
    if not m:
        raise ValueError("Bad GCST accession: " + gcst)
    n = int(m.group(1))
    start = ((n - 1) // 1000) * 1000 + 1
    end = start + 999
    return f"https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST{start}-GCST{end}/{gcst}/"

def list_dir(url):
    txt, status = fetch_text(url)
    p = LinkParser()
    p.feed(txt)
    links = []
    for href in p.links:
        if href in ("../", "./"):
            continue
        links.append(urllib.request.urljoin(url, href))
    return links, status

def discover_gwas_catalog_files(root_url, max_depth=2):
    found = []
    to_visit = [(root_url, 0)]
    visited = set()
    status_first = "NA"

    while to_visit:
        url, depth = to_visit.pop(0)
        if url in visited or depth > max_depth:
            continue
        visited.add(url)

        try:
            links, status = list_dir(url)
            if status_first == "NA":
                status_first = str(status)
        except Exception as e:
            if status_first == "NA":
                status_first = "ERROR: " + repr(e)
            continue

        for link in links:
            low = link.lower()
            if low.endswith("/"):
                to_visit.append((link, depth + 1))
            else:
                found.append(link)

    return sorted(set(found)), status_first

def classify_file(url):
    low = url.lower()
    if low.endswith(".vcf.gz.tbi") or low.endswith(".tbi"):
        return "tbi"
    if low.endswith(".vcf.gz"):
        return "vcf.gz"
    if low.endswith(".h.tsv.gz"):
        return "harmonised_tsv.gz"
    if low.endswith(".tsv.gz"):
        return "tsv.gz"
    if low.endswith(".txt.gz"):
        return "txt.gz"
    if low.endswith(".gz"):
        return "gz"
    if low.endswith(".yaml") or low.endswith(".yml"):
        return "metadata_yaml"
    if low.endswith(".json"):
        return "metadata_json"
    if low.endswith(".tsv"):
        return "tsv"
    if low.endswith(".txt"):
        return "txt"
    return "other"

def is_sumstats_candidate(url):
    low = url.lower()
    if not low.endswith((".h.tsv.gz", ".tsv.gz", ".txt.gz", ".gz", ".tab", ".txt")):
        return False
    bad = ["readme", "metadata", "checksum", ".md5", ".sha", ".tbi"]
    return not any(x in low for x in bad)

def score_sumstats_url(url):
    low = url.lower()
    score = 0
    if "harmonised" in low or "harmonized" in low:
        score += 100
    if ".h.tsv.gz" in low:
        score += 90
    if "grch37" in low or "hg19" in low:
        score += 50
    if "grch38" in low or "hg38" in low:
        score += 10
    if low.endswith(".tsv.gz"):
        score += 20
    if low.endswith(".txt.gz"):
        score += 15
    if low.endswith(".tab") or low.endswith(".txt"):
        score += 5
    return score

def write_tsv(path, fieldnames, rows):
    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

def main():
    start = time.time()

    gwas_rows = []
    zenodo_rows = []
    lock_rows = []

    for target in GWAS_TARGETS:
        outcome = target["outcome_trait"]
        gcst = target["gcst"]
        root = gcst_group_url(gcst)
        print(f"Checking GWAS Catalog FTP for {outcome}: {gcst}", flush=True)

        files, http_status = discover_gwas_catalog_files(root)
        sumstats = [u for u in files if is_sumstats_candidate(u)]
        selected = sorted(sumstats, key=score_sumstats_url, reverse=True)[0] if sumstats else ""

        for u in files:
            selected_flag = "YES" if u == selected else "NO"
            gwas_rows.append({
                "outcome_trait": outcome,
                "candidate_id": gcst,
                "trait_name": target["trait_name"],
                "source_route": "GWAS_Catalog_FTP",
                "root_url": root,
                "http_status": http_status,
                "file_type": classify_file(u),
                "is_sumstats_candidate": "YES" if is_sumstats_candidate(u) else "NO",
                "selected_for_download": selected_flag,
                "url": u,
                "note": f"n_files={len(files)};n_sumstats_candidates={len(sumstats)}",
            })

        if selected:
            lock_rows.append({
                "outcome_trait": outcome,
                "trait_name": target["trait_name"],
                "priority": target["priority"],
                "source_route": "GWAS_Catalog_FTP",
                "source_id": gcst,
                "source_status": "DOWNLOAD_CANDIDATE",
                "selected_file_or_resource": selected,
                "recommended_action": "DOWNLOAD_IN_PHASE_5_8B",
                "note": "GWAS Catalog FTP sumstats candidate found",
            })
        else:
            lock_rows.append({
                "outcome_trait": outcome,
                "trait_name": target["trait_name"],
                "priority": target["priority"],
                "source_route": "GWAS_Catalog_FTP",
                "source_id": gcst,
                "source_status": "NO_SUMSTATS_FILE_FOUND",
                "selected_file_or_resource": "NA",
                "recommended_action": "REIDENTIFY_OR_ALTERNATIVE_SOURCE",
                "note": f"GWAS Catalog FTP checked; n_files={len(files)}; n_sumstats_candidates={len(sumstats)}",
            })

    for target in ZENODO_TARGETS:
        outcome = target["outcome_trait"]
        rec = target["zenodo_record"]
        api = f"https://zenodo.org/api/records/{rec}"
        print(f"Checking Zenodo for {outcome}: {rec}", flush=True)

        try:
            obj, status = fetch_json(api)
            files = obj.get("files", [])
            if isinstance(files, dict):
                files = list(files.values())
            http_status = str(status)
        except Exception as e:
            files = []
            http_status = "ERROR: " + repr(e)

        selected = ""
        for f in files:
            key = f.get("key") or f.get("filename") or ""
            links = f.get("links", {})
            url = links.get("self") or links.get("download") or ""
            size = f.get("size", "NA")
            checksum = f.get("checksum", "NA")
            is_candidate = key.lower().endswith((".tab", ".txt", ".tsv", ".gz")) or "ntg" in key.lower()
            if is_candidate and not selected:
                selected = url if url else key

            zenodo_rows.append({
                "outcome_trait": outcome,
                "candidate_id": rec,
                "trait_name": target["trait_name"],
                "source_route": "Zenodo",
                "http_status": http_status,
                "file_name": key,
                "size_bytes": str(size),
                "checksum": str(checksum),
                "is_sumstats_candidate": "YES" if is_candidate else "NO",
                "selected_for_download": "YES" if (url if url else key) == selected else "NO",
                "url": url,
                "note": f"n_files={len(files)}",
            })

        if selected:
            lock_rows.append({
                "outcome_trait": outcome,
                "trait_name": target["trait_name"],
                "priority": target["priority"],
                "source_route": "Zenodo",
                "source_id": rec,
                "source_status": "DOWNLOAD_CANDIDATE",
                "selected_file_or_resource": selected,
                "recommended_action": "DOWNLOAD_IN_PHASE_5_8B",
                "note": "Zenodo summary statistics candidate found",
            })
        else:
            lock_rows.append({
                "outcome_trait": outcome,
                "trait_name": target["trait_name"],
                "priority": target["priority"],
                "source_route": "Zenodo",
                "source_id": rec,
                "source_status": "NO_FILE_FOUND",
                "selected_file_or_resource": "NA",
                "recommended_action": "REIDENTIFY_OR_MANUAL_REVIEW",
                "note": "Zenodo API did not return a usable file candidate",
            })

    for m in MANUAL_TARGETS:
        lock_rows.append({
            "outcome_trait": m["outcome_trait"],
            "trait_name": m["trait_name"],
            "priority": m["priority"],
            "source_route": m["source_route"],
            "source_id": m["source_id"],
            "source_status": m["source_status"],
            "selected_file_or_resource": "NA",
            "recommended_action": m["recommended_action"],
            "note": m["note"],
        })

    write_tsv(
        GWAS_CATALOG_OUT,
        [
            "outcome_trait",
            "candidate_id",
            "trait_name",
            "source_route",
            "root_url",
            "http_status",
            "file_type",
            "is_sumstats_candidate",
            "selected_for_download",
            "url",
            "note",
        ],
        gwas_rows,
    )

    write_tsv(
        ZENODO_OUT,
        [
            "outcome_trait",
            "candidate_id",
            "trait_name",
            "source_route",
            "http_status",
            "file_name",
            "size_bytes",
            "checksum",
            "is_sumstats_candidate",
            "selected_for_download",
            "url",
            "note",
        ],
        zenodo_rows,
    )

    write_tsv(
        SOURCE_LOCK_OUT,
        [
            "outcome_trait",
            "trait_name",
            "priority",
            "source_route",
            "source_id",
            "source_status",
            "selected_file_or_resource",
            "recommended_action",
            "note",
        ],
        lock_rows,
    )

    download_candidates = sum(1 for r in lock_rows if r["source_status"] == "DOWNLOAD_CANDIDATE")
    need_manual = sum(1 for r in lock_rows if r["source_status"] != "DOWNLOAD_CANDIDATE")

    elapsed = time.time() - start

    with open(STATUS_OUT, "w", encoding="utf-8", newline="") as f:
        f.write("phase\tstatus\tkey_result\tnote\n")
        f.write(f"Phase 5.8A re-identify downloadable IOP/POAG/NTG/HTG outcome sources\tPASSED\tdownload_candidates={download_candidates};manual_or_unresolved={need_manual}\tNo large files downloaded\n")
        f.write("source_routes_checked\tINFO\tGWAS_Catalog_FTP;Zenodo;manual_FinnGen_or_consortium_routes\tIOP/POAG/NTG/HTG external validation planning\n")
        f.write("download_status\tNOT_RUN\t0 files downloaded\tPhase 5.8A only verifies source availability\n")
        f.write("next_step\tTO_DO\tPhase 5.8B download/register confirmed public candidates\tOnly candidates marked DOWNLOAD_CANDIDATE should be downloaded\n")
        f.write(f"runtime\tINFO\t{elapsed:.3f}s\tPhase 5.8A runtime\n")

    with open(RUNTIME_OUT, "w", encoding="utf-8", newline="") as f:
        f.write("phase\telapsed_seconds\telapsed_human\n")
        f.write(f"Phase 5.8A\t{elapsed:.3f}\t{elapsed:.1f}s\n")

    print("===== Phase 5.8A completed =====")
    print("Wrote:", STATUS_OUT)
    print("Wrote:", GWAS_CATALOG_OUT)
    print("Wrote:", ZENODO_OUT)
    print("Wrote:", SOURCE_LOCK_OUT)

if __name__ == "__main__":
    main()
