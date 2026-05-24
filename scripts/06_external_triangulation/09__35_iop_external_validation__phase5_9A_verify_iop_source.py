#!/usr/bin/env python3
import csv
import html.parser
import os
import re
import time
import urllib.error
import urllib.parse
import urllib.request

BASE = "../../35_iop_external_validation"
os.makedirs(BASE, exist_ok=True)

FTP_VERIFICATION_OUT = os.path.join(BASE, "phase5_9A_gwas_catalog_ftp_verification.tsv")
SOURCE_LOCK_OUT = os.path.join(BASE, "phase5_9A_iop_source_lock.tsv")
STATUS_OUT = os.path.join(BASE, "phase5_9A_status.tsv")
RUNTIME_OUT = os.path.join(BASE, "phase5_9A_runtime_log.tsv")

TARGET = {
    "outcome_trait": "IOP",
    "trait_name": "Intraocular pressure",
    "gcst": "GCST009413",
    "priority": "1",
}

class LinkParser(html.parser.HTMLParser):
    def __init__(self):
        super().__init__()
        self.links = []

    def handle_starttag(self, tag, attrs):
        if tag.lower() != "a":
            return
        attrs = dict(attrs)
        href = attrs.get("href")
        if href:
            self.links.append(href)

def gcst_group_dir(gcst):
    m = re.search(r"GCST0*([0-9]+)$", gcst)
    if not m:
        raise ValueError("Cannot parse GCST number: " + gcst)
    num = int(m.group(1))
    start = ((num - 1) // 1000) * 1000 + 1
    end = start + 999
    return f"GCST{start:06d}-GCST{end:06d}"

def fetch_text(url):
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return r.status, r.read().decode("utf-8", errors="replace")

def score_file(url):
    lower = url.lower()
    base = os.path.basename(urllib.parse.urlparse(url).path).lower()

    bad_terms = ["readme", "metadata", "manifest", "md5", "sha", "log", "citation"]
    if any(x in base for x in bad_terms):
        return -100

    score = 0

    if lower.endswith(".h.tsv.gz"):
        score += 100
    elif lower.endswith(".tsv.gz"):
        score += 80
    elif lower.endswith(".txt.gz"):
        score += 60
    elif lower.endswith(".csv.gz"):
        score += 50
    elif lower.endswith(".zip"):
        score += 40
    else:
        return -100

    if "/harmonised/" in lower or "/harmonized/" in lower:
        score += 30
    if "gcst009413" in lower:
        score += 20
    if "sumstats" in lower or "summary" in lower:
        score += 10

    return score

def is_candidate_file(url):
    return score_file(url) > 0

def crawl_gcst(gcst, max_depth=3):
    group = gcst_group_dir(gcst)
    root = f"https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/{group}/{gcst}/"

    seen_dirs = set()
    queue = [(root, 0)]
    records = []

    while queue:
        url, depth = queue.pop(0)
        if url in seen_dirs or depth > max_depth:
            continue
        seen_dirs.add(url)

        try:
            status, html = fetch_text(url)
            parser = LinkParser()
            parser.feed(html)
            links = parser.links
            records.append({
                "outcome_trait": TARGET["outcome_trait"],
                "trait_name": TARGET["trait_name"],
                "gcst": gcst,
                "checked_url": url,
                "depth": str(depth),
                "entry_type": "DIRECTORY",
                "filename": "",
                "http_status": str(status),
                "is_sumstats_candidate": "NO",
                "file_score": "NA",
                "note": "directory listed",
            })
        except Exception as e:
            records.append({
                "outcome_trait": TARGET["outcome_trait"],
                "trait_name": TARGET["trait_name"],
                "gcst": gcst,
                "checked_url": url,
                "depth": str(depth),
                "entry_type": "DIRECTORY_ERROR",
                "filename": "",
                "http_status": "NA",
                "is_sumstats_candidate": "NO",
                "file_score": "NA",
                "note": repr(e),
            })
            continue

        for href in links:
            if href in ("../", "./"):
                continue

            full = urllib.parse.urljoin(url, href)
            parsed = urllib.parse.urlparse(full)
            filename = os.path.basename(parsed.path.rstrip("/"))

            if href.endswith("/"):
                if depth + 1 <= max_depth:
                    queue.append((full, depth + 1))
                records.append({
                    "outcome_trait": TARGET["outcome_trait"],
                    "trait_name": TARGET["trait_name"],
                    "gcst": gcst,
                    "checked_url": full,
                    "depth": str(depth + 1),
                    "entry_type": "SUBDIRECTORY",
                    "filename": filename,
                    "http_status": "NA",
                    "is_sumstats_candidate": "NO",
                    "file_score": "NA",
                    "note": "subdirectory discovered",
                })
            else:
                sc = score_file(full)
                records.append({
                    "outcome_trait": TARGET["outcome_trait"],
                    "trait_name": TARGET["trait_name"],
                    "gcst": gcst,
                    "checked_url": full,
                    "depth": str(depth),
                    "entry_type": "FILE",
                    "filename": filename,
                    "http_status": "NA",
                    "is_sumstats_candidate": "YES" if sc > 0 else "NO",
                    "file_score": str(sc) if sc > 0 else "NA",
                    "note": "file discovered",
                })

    return records

def write_tsv(path, fieldnames, rows):
    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

def main():
    start = time.time()

    print("Checking GWAS Catalog FTP for IOP:", TARGET["gcst"], flush=True)
    records = crawl_gcst(TARGET["gcst"])

    fields = [
        "outcome_trait",
        "trait_name",
        "gcst",
        "checked_url",
        "depth",
        "entry_type",
        "filename",
        "http_status",
        "is_sumstats_candidate",
        "file_score",
        "note",
    ]

    write_tsv(FTP_VERIFICATION_OUT, fields, records)

    candidates = [
        r for r in records
        if r["entry_type"] == "FILE" and r["is_sumstats_candidate"] == "YES"
    ]
    candidates.sort(key=lambda r: int(r["file_score"]), reverse=True)

    if candidates:
        best = candidates[0]
        lock_status = "DOWNLOAD_CANDIDATE"
        selected = best["checked_url"]
        action = "DOWNLOAD_IN_PHASE_5_9B"
        note = f"GWAS Catalog FTP sumstats candidate found; n_candidates={len(candidates)}"
    else:
        lock_status = "NO_SUMSTATS_FILE_FOUND"
        selected = "NA"
        action = "REIDENTIFY_OR_ALTERNATIVE_SOURCE"
        note = "GWAS Catalog FTP checked; no usable sumstats candidate found"

    lock_rows = [{
        "outcome_trait": TARGET["outcome_trait"],
        "trait_name": TARGET["trait_name"],
        "priority": TARGET["priority"],
        "source_route": "GWAS_Catalog_FTP",
        "source_id": TARGET["gcst"],
        "source_status": lock_status,
        "selected_file_or_resource": selected,
        "recommended_action": action,
        "note": note,
    }]

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

    elapsed = time.time() - start

    with open(STATUS_OUT, "w", encoding="utf-8", newline="") as f:
        f.write("phase\tstatus\tkey_result\tnote\n")
        phase_status = "PASSED_DOWNLOAD_CANDIDATE_FOUND" if candidates else "PASSED_NO_DOWNLOAD_CANDIDATE"
        f.write(f"Phase 5.9A verify IOP downloadable source\t{phase_status}\tdownload_candidates={len(candidates)}\tNo large files downloaded\n")
        f.write("target_outcome\tINFO\tIOP\tExternal intraocular pressure validation outcome\n")
        f.write(f"candidate_dataset\tINFO\t{TARGET['gcst']}\tGWAS Catalog FTP route checked\n")
        f.write("download_status\tNOT_RUN\t0 files downloaded\tPhase 5.9A only verifies source availability\n")
        f.write(f"runtime\tINFO\t{elapsed:.3f}s\tPhase 5.9A runtime\n")

    with open(RUNTIME_OUT, "w", encoding="utf-8", newline="") as f:
        f.write("phase\telapsed_seconds\telapsed_human\n")
        f.write(f"Phase 5.9A\t{elapsed:.3f}\t{elapsed:.1f}s\n")

    print("===== Phase 5.9A completed =====")
    print("Wrote:", STATUS_OUT)
    print("Wrote:", FTP_VERIFICATION_OUT)
    print("Wrote:", SOURCE_LOCK_OUT)

if __name__ == "__main__":
    main()
