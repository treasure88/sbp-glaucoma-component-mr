#!/usr/bin/env python3
import csv
import html.parser
import os
import re
import time
import urllib.parse
import urllib.request

BASE = "../../35_iop_external_validation"
os.makedirs(BASE, exist_ok=True)

GCST = "GCST009413"
OUT = os.path.join(BASE, "phase5_9A2_strict_gwas_catalog_ftp_verification.tsv")
LOCK = os.path.join(BASE, "phase5_9A2_iop_source_lock.tsv")
STATUS = os.path.join(BASE, "phase5_9A2_status.tsv")
RUNTIME = os.path.join(BASE, "phase5_9A2_runtime_log.tsv")

class LinkParser(html.parser.HTMLParser):
    def __init__(self):
        super().__init__()
        self.links = []
    def handle_starttag(self, tag, attrs):
        if tag.lower() != "a":
            return
        d = dict(attrs)
        if d.get("href"):
            self.links.append(d["href"])

def gcst_group_dir(gcst):
    m = re.search(r"GCST0*([0-9]+)$", gcst)
    num = int(m.group(1))
    start = ((num - 1) // 1000) * 1000 + 1
    end = start + 999
    return f"GCST{start:06d}-GCST{end:06d}"

def fetch_links(url):
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=60) as r:
        html = r.read().decode("utf-8", errors="replace")
    p = LinkParser()
    p.feed(html)
    return p.links

def score(url):
    low = url.lower()
    base = os.path.basename(urllib.parse.urlparse(url).path).lower()

    if f"/{GCST.lower()}/" not in low:
        return -1000
    if any(x in base for x in ["readme", "metadata", "manifest", "md5", "sha", "log"]):
        return -100

    s = 0
    if low.endswith(".h.tsv.gz"):
        s += 100
    elif low.endswith(".tsv.gz"):
        s += 80
    elif low.endswith(".txt.gz"):
        s += 60
    elif low.endswith(".csv.gz"):
        s += 50
    elif low.endswith(".zip"):
        s += 40
    else:
        return -100

    if "/harmonised/" in low or "/harmonized/" in low:
        s += 30
    if GCST.lower() in low:
        s += 20
    return s

def main():
    t0 = time.time()
    group = gcst_group_dir(GCST)
    root = f"https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/{group}/{GCST}/"

    queue = [(root, 0)]
    seen = set()
    rows = []

    while queue:
        url, depth = queue.pop(0)
        if url in seen or depth > 3:
            continue
        seen.add(url)

        try:
            links = fetch_links(url)
            rows.append({
                "gcst": GCST,
                "checked_url": url,
                "depth": str(depth),
                "entry_type": "DIRECTORY",
                "filename": "",
                "is_exact_gcst_path": "YES" if f"/{GCST}/" in url else "NO",
                "is_sumstats_candidate": "NO",
                "file_score": "NA",
                "note": "directory listed",
            })
        except Exception as e:
            rows.append({
                "gcst": GCST,
                "checked_url": url,
                "depth": str(depth),
                "entry_type": "DIRECTORY_ERROR",
                "filename": "",
                "is_exact_gcst_path": "YES" if f"/{GCST}/" in url else "NO",
                "is_sumstats_candidate": "NO",
                "file_score": "NA",
                "note": repr(e),
            })
            continue

        for href in links:
            if href in ("../", "./"):
                continue
            full = urllib.parse.urljoin(url, href)
            filename = os.path.basename(urllib.parse.urlparse(full).path.rstrip("/"))
            exact = "YES" if f"/{GCST}/" in full else "NO"

            if href.endswith("/"):
                if exact == "YES":
                    queue.append((full, depth + 1))
                rows.append({
                    "gcst": GCST,
                    "checked_url": full,
                    "depth": str(depth + 1),
                    "entry_type": "SUBDIRECTORY",
                    "filename": filename,
                    "is_exact_gcst_path": exact,
                    "is_sumstats_candidate": "NO",
                    "file_score": "NA",
                    "note": "subdirectory discovered",
                })
            else:
                sc = score(full)
                rows.append({
                    "gcst": GCST,
                    "checked_url": full,
                    "depth": str(depth),
                    "entry_type": "FILE",
                    "filename": filename,
                    "is_exact_gcst_path": exact,
                    "is_sumstats_candidate": "YES" if sc > 0 else "NO",
                    "file_score": str(sc) if sc > 0 else "NA",
                    "note": "file discovered",
                })

    fields = [
        "gcst", "checked_url", "depth", "entry_type", "filename",
        "is_exact_gcst_path", "is_sumstats_candidate", "file_score", "note"
    ]

    with open(OUT, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    candidates = [
        r for r in rows
        if r["entry_type"] == "FILE"
        and r["is_exact_gcst_path"] == "YES"
        and r["is_sumstats_candidate"] == "YES"
    ]
    candidates.sort(key=lambda r: int(r["file_score"]), reverse=True)

    if candidates:
        status = "DOWNLOAD_CANDIDATE"
        selected = candidates[0]["checked_url"]
        action = "DOWNLOAD_IN_PHASE_5_9B"
        note = f"Strict exact-GCST candidate found; n_candidates={len(candidates)}"
    else:
        status = "NO_EXACT_GCST_SUMSTATS_FOUND"
        selected = "NA"
        action = "REIDENTIFY_OR_ALTERNATIVE_SOURCE"
        note = "No downloadable sumstats found under exact GCST009413 directory"

    lock_rows = [{
        "outcome_trait": "IOP",
        "trait_name": "Intraocular pressure",
        "source_route": "GWAS_Catalog_FTP_STRICT",
        "source_id": GCST,
        "source_status": status,
        "selected_file_or_resource": selected,
        "recommended_action": action,
        "note": note,
    }]

    with open(LOCK, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(lock_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(lock_rows)

    elapsed = time.time() - t0

    with open(STATUS, "w", encoding="utf-8") as f:
        f.write("phase\tstatus\tkey_result\tnote\n")
        phase_status = "PASSED_EXACT_DOWNLOAD_CANDIDATE_FOUND" if candidates else "PASSED_NO_EXACT_DOWNLOAD_CANDIDATE"
        f.write(f"Phase 5.9A2 strict IOP source verification\t{phase_status}\texact_download_candidates={len(candidates)}\tOnly files under /GCST009413/ are eligible\n")
        f.write(f"previous_phase5_9A\tINVALIDATED_FALSE_POSITIVE\tselected_non_target_GCST009236\tDo not download previous selected file\n")
        f.write(f"runtime\tINFO\t{elapsed:.3f}s\tPhase 5.9A2 runtime\n")

    with open(RUNTIME, "w", encoding="utf-8") as f:
        f.write("phase\telapsed_seconds\telapsed_human\n")
        f.write(f"Phase 5.9A2\t{elapsed:.3f}\t{elapsed:.1f}s\n")

    print("===== Phase 5.9A2 completed =====")
    print("Wrote:", STATUS)
    print("Wrote:", LOCK)
    print("Wrote:", OUT)

if __name__ == "__main__":
    main()
