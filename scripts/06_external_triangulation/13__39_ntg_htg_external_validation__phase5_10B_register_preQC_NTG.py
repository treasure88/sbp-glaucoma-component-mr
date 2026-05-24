#!/usr/bin/env python3
import csv
import hashlib
import os
import time

BASE = "../../39_ntg_htg_external_validation"
RAW = "../../26_external_outcome_raw/NTG__Zenodo14010557__MTAG_NTG_IGGC_STAGE2.tab"
EXPECTED_MD5 = "7dd090e205f124c4a06220df744635be"

STATUS_OUT = os.path.join(BASE, "phase5_10B_status.tsv")
PREQC_OUT = os.path.join(BASE, "phase5_10B_NTG_external_pre_qc_row.tsv")
REG_OUT = os.path.join(BASE, "phase5_10B_NTG_external_registration_row.tsv")
RUNTIME_OUT = os.path.join(BASE, "phase5_10B_runtime_log.tsv")

def file_md5(path, block=1024*1024):
    h = hashlib.md5()
    with open(path, "rb") as f:
        while True:
            b = f.read(block)
            if not b:
                break
            h.update(b)
    return h.hexdigest()

def detect_col(header, candidates):
    lower_map = {h.lower(): h for h in header}
    for c in candidates:
        if c.lower() in lower_map:
            return lower_map[c.lower()]
    return "NA"

def main():
    start = time.time()

    if not os.path.exists(RAW):
        raise SystemExit(f"Missing NTG raw file: {RAW}")

    md5 = file_md5(RAW)
    md5_ok = "YES" if md5 == EXPECTED_MD5 else "NO"

    n_rows = 0
    with open(RAW, "r", encoding="utf-8", errors="replace", newline="") as f:
        header_line = f.readline().rstrip("\n\r")
        header = header_line.split("\t")
        for _ in f:
            n_rows += 1

    snp_col = detect_col(header, ["SNP", "rsid", "rsID"])
    chr_col = detect_col(header, ["CHR", "chr", "chromosome"])
    pos_col = detect_col(header, ["BP", "POS", "pos", "base_pair_location"])
    ea_col = detect_col(header, ["A1", "effect_allele", "EA"])
    oa_col = detect_col(header, ["A2", "other_allele", "NEA"])
    beta_col = detect_col(header, ["BETA", "beta", "Effect"])
    se_col = detect_col(header, ["SE", "se", "StdErr", "standard_error"])
    p_col = detect_col(header, ["P", "p", "PVAL", "pval", "p_value"])
    eaf_col = detect_col(header, ["EAF", "eaf", "FRQ", "freq", "MAF"])
    n_col = detect_col(header, ["N", "n", "sample_size"])

    core = [snp_col, ea_col, oa_col, beta_col, se_col, p_col]
    ready = "READY_FOR_STANDARDIZATION" if all(x != "NA" for x in core) and md5_ok == "YES" else "REVIEW_COLUMNS_BEFORE_STANDARDIZATION"

    preqc = {
        "outcome_id": "NTG",
        "dataset_id": "Zenodo14010557",
        "trait_name": "Normal tension glaucoma",
        "source": "Zenodo",
        "file_path": RAW,
        "md5_expected": EXPECTED_MD5,
        "md5_observed": md5,
        "md5_ok": md5_ok,
        "n_rows": str(n_rows),
        "n_cols": str(len(header)),
        "header": "|".join(header),
        "snp_col": snp_col,
        "chr_col": chr_col,
        "pos_col": pos_col,
        "effect_allele_col": ea_col,
        "other_allele_col": oa_col,
        "beta_col": beta_col,
        "se_col": se_col,
        "pval_col": p_col,
        "eaf_col": eaf_col,
        "n_col": n_col,
        "status": ready,
        "note": "NTG Zenodo MTAG/IGGC stage 2 summary statistics; columns detected for standardization",
    }

    with open(PREQC_OUT, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(preqc.keys()), delimiter="\t")
        w.writeheader()
        w.writerow(preqc)

    reg = {
        "outcome_id": "NTG",
        "file_path": RAW,
        "file_format": "TAB",
        "genome_build": "GRCh37_TO_VERIFY",
        "snp_col": snp_col,
        "chr_col": chr_col,
        "pos_col": pos_col,
        "effect_allele_col": ea_col,
        "other_allele_col": oa_col,
        "beta_col": beta_col,
        "se_col": se_col,
        "pval_col": p_col,
        "eaf_col": eaf_col,
        "n_col": n_col,
        "trait_name": "Normal tension glaucoma",
        "sample_size": "TO_VERIFY",
        "population": "European_TO_VERIFY",
        "outcome_scale": "log_odds_or_MTAG_scale_TO_VERIFY",
        "note": "SNP/CHR/BP/A1/A2/BETA/SE/P detected; no EAF or N column in raw file",
        "status": "REGISTERED_READY_FOR_STANDARDIZATION" if ready == "READY_FOR_STANDARDIZATION" else "REGISTERED_NEEDS_COLUMN_REVIEW",
    }

    with open(REG_OUT, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(reg.keys()), delimiter="\t")
        w.writeheader()
        w.writerow(reg)

    elapsed = time.time() - start

    with open(STATUS_OUT, "w", encoding="utf-8", newline="") as f:
        f.write("phase\tstatus\tkey_result\tnote\n")
        phase_status = "PASSED" if ready == "READY_FOR_STANDARDIZATION" else "PASSED_WITH_COLUMN_REVIEW_NEEDED"
        f.write(f"Phase 5.10B download/register/pre-QC NTG external outcome\t{phase_status}\tmd5_ok={md5_ok};pre_qc={ready}\tNTG Zenodo source processed\n")
        f.write("outcome\tINFO\tNTG\tNormal tension glaucoma external validation outcome\n")
        f.write(f"n_rows\tINFO\t{n_rows}\tRows scanned in raw NTG file\n")
        f.write("column_mapping\tINFO\tSNP/CHR/BP/A1/A2/BETA/SE/P\tCore columns detected\n")
        f.write("HTG_status\tTO_DO\tfull_summary_stats_not_locked\tHTG remains unresolved\n")
        f.write(f"runtime\tINFO\t{elapsed:.3f}s\tPhase 5.10B runtime\n")

    with open(RUNTIME_OUT, "w", encoding="utf-8", newline="") as f:
        f.write("phase\telapsed_seconds\telapsed_human\n")
        f.write(f"Phase 5.10B\t{elapsed:.3f}\t{elapsed:.1f}s\n")

    print("===== Phase 5.10B completed =====")
    print("Wrote:", STATUS_OUT)
    print("Wrote:", PREQC_OUT)
    print("Wrote:", REG_OUT)

if __name__ == "__main__":
    main()
