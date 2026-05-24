#!/usr/bin/env python3
import csv
import json
import os
import re
import time
import urllib.error
import urllib.parse
import urllib.request

BASE = "../../39_ntg_htg_external_validation"
os.makedirs(BASE, exist_ok=True)

ZENODO_RECORD = "14010557"
ZENODO_API = f"https://zenodo.org/api/records/{ZENODO_RECORD}"
ZENODO_HTML = f"https://zenodo.org/records/{ZENODO_RECORD}"

ZENODO_OUT = os.path.join(BASE, "phase5_10A_zenodo_NTG_verification.tsv")
LOCK_OUT = os.path.join(BASE, "phase5_10A_NTG_HTG_source_lock.tsv")
STATUS_OUT = os.path.join(BASE, "phase5_10A_status.tsv")
RUNTIME_OUT = os.path.join(BASE, "phase5_10A_runtime_log.tsv")
RAW_JSON = os.path.join(BASE, "phase5_10A_zenodo_14010557_raw.json")

def fetch_json(url):
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0", "Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=90) as r:
        return r.status, json.loads(r.read().decode("utf-8", errors="replace"))

def score_zenodo_file(filename):
    low = filename.lower()
    score = 0

    if low == "mtag_ntg_iggc_stage2.tab":
        score += 200
    if "ntg" in low:
        score += 50
    if "iggc" in low:
        score += 30
    if "stage2" in low:
        score += 20

    if low.endswith(".tab"):
        score += 30
    elif low.endswith(".txt") or low.endswith(".tsv"):
        score += 20
    elif low.endswith(".gz"):
        score += 10

    return score

def human_size(n):
    try:
        n = float(n)
    except Exception:
        return "NA"
    units = ["B", "KB", "MB", "GB"]
    i = 0
    while n >= 1024 and i < len(units)-1:
        n /= 1024
        i += 1
    return f"{n:.1f}{units[i]}"

def verify_zenodo():
    rows = []
    status = "API_FAILED"
    note = "NA"

    try:
        api_status, obj = fetch_json(ZENODO_API)
        status = f"API_{api_status}"

        with open(RAW_JSON, "w", encoding="utf-8") as f:
            json.dump(obj, f, indent=2)

        files = obj.get("files", [])
        for fobj in files:
            filename = fobj.get("key") or fobj.get("filename") or ""
            size_bytes = fobj.get("size", "")
            checksum = fobj.get("checksum", "")
            links = fobj.get("links", {})
            download_url = links.get("self") or links.get("download") or ""

            rows.append({
                "outcome_trait": "NTG",
                "trait_name": "Normal tension glaucoma",
                "source_route": "Zenodo",
                "source_id": ZENODO_RECORD,
                "api_status": str(api_status),
                "filename": filename,
                "size_bytes": str(size_bytes),
                "size_human": human_size(size_bytes),
                "checksum": checksum,
                "download_url": download_url,
                "file_score": str(score_zenodo_file(filename)),
                "is_download_candidate": "YES" if score_zenodo_file(filename) >= 200 and download_url else "NO",
                "note": "Zenodo API file record",
            })

        note = f"Zenodo API returned {len(files)} file records"

    except Exception as e:
        note = repr(e)

    return rows, status, note

def write_tsv(path, fields, rows):
    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

def main():
    start = time.time()

    zenodo_rows, zenodo_status, zenodo_note = verify_zenodo()

    zenodo_fields = [
        "outcome_trait", "trait_name", "source_route", "source_id",
        "api_status", "filename", "size_bytes", "size_human",
        "checksum", "download_url", "file_score", "is_download_candidate", "note"
    ]
    write_tsv(ZENODO_OUT, zenodo_fields, zenodo_rows)

    ntg_candidates = [r for r in zenodo_rows if r["is_download_candidate"] == "YES"]
    ntg_candidates.sort(key=lambda r: int(r["file_score"]), reverse=True)

    lock_rows = []

    if ntg_candidates:
        best = ntg_candidates[0]
        lock_rows.append({
            "outcome_trait": "NTG",
            "trait_name": "Normal tension glaucoma",
            "priority": "1",
            "source_route": "Zenodo",
            "source_id": ZENODO_RECORD,
            "source_status": "DOWNLOAD_CANDIDATE",
            "selected_file_or_resource": best["download_url"],
            "selected_filename": best["filename"],
            "size_bytes": best["size_bytes"],
            "checksum": best["checksum"],
            "recommended_action": "DOWNLOAD_IN_PHASE_5_10B",
            "note": "Exact NTG Zenodo file candidate verified",
        })
    else:
        lock_rows.append({
            "outcome_trait": "NTG",
            "trait_name": "Normal tension glaucoma",
            "priority": "1",
            "source_route": "Zenodo",
            "source_id": ZENODO_RECORD,
            "source_status": "NO_FILE_CANDIDATE",
            "selected_file_or_resource": "NA",
            "selected_filename": "NA",
            "size_bytes": "NA",
            "checksum": "NA",
            "recommended_action": "REVIEW_ZENODO_PAGE_OR_MANUAL_DOWNLOAD",
            "note": zenodo_note,
        })

    lock_rows.append({
        "outcome_trait": "NTG",
        "trait_name": "Normotensive glaucoma",
        "priority": "2",
        "source_route": "FinnGen",
        "source_id": "H7_GLAUCOMA_NTG",
        "source_status": "ACCESS_REQUIRED",
        "selected_file_or_resource": "NA",
        "selected_filename": "NA",
        "size_bytes": "NA",
        "checksum": "NA",
        "recommended_action": "USE_AS_BACKUP_OR_SECONDARY_AFTER_ACCESS",
        "note": "FinnGen route requires access request/download instructions",
    })

    lock_rows.append({
        "outcome_trait": "HTG",
        "trait_name": "High-tension glaucoma",
        "priority": "1",
        "source_route": "publication_or_consortium",
        "source_id": "Gharahkhani_HTG_GWAS_or_IGGC",
        "source_status": "TO_REIDENTIFY_OR_REQUEST",
        "selected_file_or_resource": "NA",
        "selected_filename": "NA",
        "size_bytes": "NA",
        "checksum": "NA",
        "recommended_action": "SEARCH_PUBLIC_SUMSTATS_OR_REQUEST_FROM_AUTHORS",
        "note": "Direct downloadable HTG full summary statistics not locked in this phase",
    })

    lock_fields = [
        "outcome_trait", "trait_name", "priority", "source_route", "source_id",
        "source_status", "selected_file_or_resource", "selected_filename",
        "size_bytes", "checksum", "recommended_action", "note"
    ]
    write_tsv(LOCK_OUT, lock_fields, lock_rows)

    elapsed = time.time() - start

    download_candidates = sum(1 for r in lock_rows if r["source_status"] == "DOWNLOAD_CANDIDATE")
    htgs = [r for r in lock_rows if r["outcome_trait"] == "HTG"]
    htg_status = htgs[0]["source_status"] if htgs else "NA"

    with open(STATUS_OUT, "w", encoding="utf-8", newline="") as f:
        f.write("phase\tstatus\tkey_result\tnote\n")
        phase_status = "PASSED_NTG_DOWNLOAD_CANDIDATE_HTG_UNRESOLVED" if download_candidates > 0 else "PASSED_NO_DOWNLOAD_CANDIDATE"
        f.write(f"Phase 5.10A verify NTG/HTG source routes\t{phase_status}\tntg_download_candidates={download_candidates};htg_status={htg_status}\tNo large files downloaded\n")
        f.write(f"NTG_Zenodo_API\tINFO\t{zenodo_status}\t{zenodo_note}\n")
        f.write("HTG_route\tTO_DO\tpublic_full_sumstats_not_locked\tSearch/request required before HTG MR\n")
        f.write("download_status\tNOT_RUN\t0 files downloaded\tPhase 5.10A only verifies source availability\n")
        f.write(f"runtime\tINFO\t{elapsed:.3f}s\tPhase 5.10A runtime\n")

    with open(RUNTIME_OUT, "w", encoding="utf-8", newline="") as f:
        f.write("phase\telapsed_seconds\telapsed_human\n")
        f.write(f"Phase 5.10A\t{elapsed:.3f}\t{elapsed:.1f}s\n")

    print("===== Phase 5.10A completed =====")
    print("Wrote:", STATUS_OUT)
    print("Wrote:", ZENODO_OUT)
    print("Wrote:", LOCK_OUT)

if __name__ == "__main__":
    main()
