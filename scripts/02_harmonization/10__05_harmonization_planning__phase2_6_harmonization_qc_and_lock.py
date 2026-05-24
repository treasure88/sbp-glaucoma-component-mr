#!/usr/bin/env python3
import gzip
import hashlib
import os
from collections import Counter

HARM_DIR = "../../10_harmonized_datasets"
LOCK_DIR = "../../11_locked_harmonized_datasets"
PLAN = "../../05_harmonization_planning"

os.makedirs(LOCK_DIR, exist_ok=True)

EXPECTED_FIELDS = 25

QC_OUT = os.path.join(LOCK_DIR, "phase2_6_harmonization_qc_summary.tsv")
LOCK_MANIFEST = os.path.join(LOCK_DIR, "phase2_6_locked_harmonized_dataset_manifest.tsv")
FINAL_STATUS = os.path.join(LOCK_DIR, "phase2_6_dataset_lock_final_status.tsv")

EXPECTED_EXPOSURES = ["SBP", "DBP", "MIGRAINE", "INSOMNIA", "CRP", "BMI"]
EXPECTED_OUTCOME_SUFFIXES = [
    "GBS_nonIOPcomponent",
    "GBS_IOPcomponent",
    "IOPcc_coordinate_subset",
]

def sha256_file(path, block_size=1024 * 1024):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            block = f.read(block_size)
            if not block:
                break
            h.update(block)
    return h.hexdigest()

def parse_expected_from_filename(path):
    base = os.path.basename(path)
    if not base.endswith(".harmonized.tsv.gz"):
        return ("NA", "NA")

    stem = base.replace(".harmonized.tsv.gz", "")
    if "__" not in stem:
        return ("NA", "NA")

    exposure_id, outcome_suffix = stem.split("__", 1)
    return exposure_id, outcome_suffix

def qc_one_file(path):
    exposure_from_name, outcome_suffix = parse_expected_from_filename(path)

    n_rows = 0
    bad_rows = 0
    include_yes = 0
    include_no = 0
    missing_beta_exposure_included = 0
    missing_beta_outcome_harmonized_included = 0
    missing_se_exposure_included = 0
    missing_se_outcome_included = 0
    missing_pval_exposure_included = 0
    missing_pval_outcome_included = 0
    missing_eaf_exposure_included = 0
    palindromic_included = 0
    allele_mismatch_included = 0
    action_counter = Counter()
    outcome_id_seen = Counter()
    match_mode_seen = Counter()

    with gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") as f:
        header = f.readline().rstrip("\n\r")
        cols = header.split("\t")
        header_fields = len(cols)
        cmap = {c: i for i, c in enumerate(cols)}

        required = [
            "exposure_id",
            "outcome_id",
            "match_mode",
            "match_key",
            "chr",
            "pos",
            "SNP",
            "exposure_effect_allele",
            "exposure_other_allele",
            "outcome_effect_allele_original",
            "outcome_other_allele_original",
            "beta_exposure",
            "se_exposure",
            "pval_exposure",
            "minus_log10_p_exposure",
            "eaf_exposure",
            "n_exposure",
            "beta_outcome_original",
            "beta_outcome_harmonized",
            "se_outcome",
            "pval_outcome",
            "pval_for_log10_outcome",
            "harmonization_action",
            "include_in_main",
            "outcome_source_dataset",
        ]

        missing_required = [c for c in required if c not in cmap]

        for line in f:
            parts = line.rstrip("\n\r").split("\t")
            n_rows += 1

            if len(parts) != EXPECTED_FIELDS:
                bad_rows += 1
                continue

            action = parts[cmap["harmonization_action"]]
            include = parts[cmap["include_in_main"]]
            outcome_id = parts[cmap["outcome_id"]]
            match_mode = parts[cmap["match_mode"]]

            action_counter[action] += 1
            outcome_id_seen[outcome_id] += 1
            match_mode_seen[match_mode] += 1

            if include == "YES":
                include_yes += 1

                if action.startswith("palindromic"):
                    palindromic_included += 1
                if action == "allele_mismatch_excluded":
                    allele_mismatch_included += 1

                def is_missing(x):
                    return x in ("", ".", "NA")

                if is_missing(parts[cmap["beta_exposure"]]):
                    missing_beta_exposure_included += 1
                if is_missing(parts[cmap["beta_outcome_harmonized"]]):
                    missing_beta_outcome_harmonized_included += 1
                if is_missing(parts[cmap["se_exposure"]]):
                    missing_se_exposure_included += 1
                if is_missing(parts[cmap["se_outcome"]]):
                    missing_se_outcome_included += 1
                if is_missing(parts[cmap["pval_exposure"]]):
                    missing_pval_exposure_included += 1
                if is_missing(parts[cmap["pval_outcome"]]):
                    missing_pval_outcome_included += 1
                if is_missing(parts[cmap["eaf_exposure"]]):
                    missing_eaf_exposure_included += 1

            elif include == "NO":
                include_no += 1

    action_summary = ";".join(f"{k}={v}" for k, v in sorted(action_counter.items()))
    outcome_id_summary = ";".join(f"{k}={v}" for k, v in sorted(outcome_id_seen.items()))
    match_mode_summary = ";".join(f"{k}={v}" for k, v in sorted(match_mode_seen.items()))

    structure_status = "PASS" if header_fields == EXPECTED_FIELDS and bad_rows == 0 and not missing_required else "FAIL"

    inclusion_status = "PASS"
    notes = []

    if include_yes == 0:
        inclusion_status = "FAIL"
        notes.append("No variants included in main analysis")

    if palindromic_included != 0:
        inclusion_status = "FAIL"
        notes.append("Palindromic variants were included unexpectedly")

    if allele_mismatch_included != 0:
        inclusion_status = "FAIL"
        notes.append("Allele mismatch variants were included unexpectedly")

    if missing_beta_exposure_included != 0 or missing_beta_outcome_harmonized_included != 0:
        inclusion_status = "FAIL"
        notes.append("Included rows have missing beta values")

    if missing_se_exposure_included != 0 or missing_se_outcome_included != 0:
        inclusion_status = "WARN"
        notes.append("Included rows have missing SE values")

    if missing_eaf_exposure_included != 0:
        notes.append("Included rows have missing EAF; this may affect palindromic resolution but palindromic variants are excluded")

    if outcome_suffix == "IOPcc_coordinate_subset":
        notes.append("IOPcc uses coordinate-resolvable subset only")

    if exposure_from_name == "BMI":
        notes.append("BMI has known EAF caveat from Phase 2.2")

    lock_status = "LOCKED" if structure_status == "PASS" and inclusion_status in ("PASS", "WARN") else "NOT_LOCKED"

    return {
        "exposure_id": exposure_from_name,
        "outcome_suffix": outcome_suffix,
        "file_path": path,
        "file_size_bytes": os.path.getsize(path),
        "sha256": sha256_file(path),
        "header_fields": header_fields,
        "data_rows": n_rows,
        "bad_rows_expected_25_fields": bad_rows,
        "include_in_main_yes": include_yes,
        "include_in_main_no": include_no,
        "missing_beta_exposure_included": missing_beta_exposure_included,
        "missing_beta_outcome_harmonized_included": missing_beta_outcome_harmonized_included,
        "missing_se_exposure_included": missing_se_exposure_included,
        "missing_se_outcome_included": missing_se_outcome_included,
        "missing_pval_exposure_included": missing_pval_exposure_included,
        "missing_pval_outcome_included": missing_pval_outcome_included,
        "missing_eaf_exposure_included": missing_eaf_exposure_included,
        "palindromic_included": palindromic_included,
        "allele_mismatch_included": allele_mismatch_included,
        "outcome_id_summary": outcome_id_summary,
        "match_mode_summary": match_mode_summary,
        "harmonization_action_counts": action_summary,
        "structure_status": structure_status,
        "inclusion_status": inclusion_status,
        "lock_status": lock_status,
        "qc_note": "; ".join(notes) if notes else "OK",
    }

def main():
    files = []

    for root, dirs, names in os.walk(HARM_DIR):
        for name in names:
            if name.endswith(".harmonized.tsv.gz"):
                files.append(os.path.join(root, name))

    files = sorted(files)

    expected_count = len(EXPECTED_EXPOSURES) * len(EXPECTED_OUTCOME_SUFFIXES)

    print(f"Found harmonized files: {len(files)}")
    print(f"Expected harmonized files: {expected_count}")

    rows = []

    for path in files:
        print(f"QC and locking: {path}")
        rows.append(qc_one_file(path))

    qc_cols = [
        "exposure_id",
        "outcome_suffix",
        "file_path",
        "file_size_bytes",
        "sha256",
        "header_fields",
        "data_rows",
        "bad_rows_expected_25_fields",
        "include_in_main_yes",
        "include_in_main_no",
        "missing_beta_exposure_included",
        "missing_beta_outcome_harmonized_included",
        "missing_se_exposure_included",
        "missing_se_outcome_included",
        "missing_pval_exposure_included",
        "missing_pval_outcome_included",
        "missing_eaf_exposure_included",
        "palindromic_included",
        "allele_mismatch_included",
        "outcome_id_summary",
        "match_mode_summary",
        "harmonization_action_counts",
        "structure_status",
        "inclusion_status",
        "lock_status",
        "qc_note",
    ]

    with open(QC_OUT, "w", encoding="utf-8", newline="\n") as f:
        f.write("\t".join(qc_cols) + "\n")
        for r in rows:
            f.write("\t".join(str(r[c]) for c in qc_cols) + "\n")

    manifest_cols = [
        "exposure_id",
        "outcome_suffix",
        "locked_harmonized_file",
        "sha256",
        "data_rows",
        "include_in_main_yes",
        "lock_status",
        "qc_note",
    ]

    with open(LOCK_MANIFEST, "w", encoding="utf-8", newline="\n") as f:
        f.write("\t".join(manifest_cols) + "\n")
        for r in rows:
            f.write("\t".join(str(r[c]) for c in [
                "exposure_id",
                "outcome_suffix",
                "file_path",
                "sha256",
                "data_rows",
                "include_in_main_yes",
                "lock_status",
                "qc_note",
            ]) + "\n")

    n_locked = sum(1 for r in rows if r["lock_status"] == "LOCKED")
    n_not_locked = sum(1 for r in rows if r["lock_status"] != "LOCKED")
    n_structure_fail = sum(1 for r in rows if r["structure_status"] != "PASS")
    n_inclusion_fail = sum(1 for r in rows if r["inclusion_status"] == "FAIL")
    n_pal_included = sum(r["palindromic_included"] for r in rows)
    n_mismatch_included = sum(r["allele_mismatch_included"] for r in rows)

    with open(FINAL_STATUS, "w", encoding="utf-8", newline="\n") as f:
        f.write("phase\tstatus\tkey_result\tnote\n")

        if len(files) == expected_count and n_locked == expected_count and n_not_locked == 0:
            phase_status = "PASSED"
            key_result = f"{n_locked}/{expected_count} harmonized datasets locked"
            note = "All expected harmonized datasets passed structural and inclusion QC"
        else:
            phase_status = "REVIEW_NEEDED"
            key_result = f"{n_locked}/{expected_count} harmonized datasets locked"
            note = "Some expected datasets are missing or failed QC"

        f.write(f"Phase 2.6 harmonization QC summary and dataset locking\t{phase_status}\t{key_result}\t{note}\n")
        f.write(f"expected_harmonized_file_count\tINFO\t{expected_count}\t6 exposures x 3 outcomes\n")
        f.write(f"observed_harmonized_file_count\tINFO\t{len(files)}\tDetected .harmonized.tsv.gz files\n")
        f.write(f"locked_dataset_count\tINFO\t{n_locked}\tFiles with lock_status=LOCKED\n")
        f.write(f"not_locked_dataset_count\tINFO\t{n_not_locked}\tFiles with lock_status other than LOCKED\n")
        f.write(f"structure_fail_count\tINFO\t{n_structure_fail}\tFiles failing 25-column structural QC\n")
        f.write(f"inclusion_fail_count\tINFO\t{n_inclusion_fail}\tFiles failing inclusion QC\n")
        f.write(f"palindromic_included_total\tINFO\t{n_pal_included}\tShould be 0\n")
        f.write(f"allele_mismatch_included_total\tINFO\t{n_mismatch_included}\tShould be 0\n")
        f.write("IOPcc_limitation\tDOCUMENTED\tcoordinate-resolvable subset only\tAffx-style IOPcc variants are not included without external mapping\n")
        f.write("BMI_EAF_caveat\tDOCUMENTED\tBMI has known EAF caveat\tPalindromic variants remain excluded from main analysis\n")
        f.write("CRP_pvalue_caveat\tDOCUMENTED\tCRP retains minus_log10_p\tAvoid p-value underflow in downstream analysis\n")

    print("")
    print("Wrote QC summary:", QC_OUT)
    print("Wrote locked manifest:", LOCK_MANIFEST)
    print("Wrote final status:", FINAL_STATUS)
    print("")
    print(f"Locked datasets: {n_locked}/{expected_count}")
    print(f"Not locked datasets: {n_not_locked}")
    print(f"Structure failures: {n_structure_fail}")
    print(f"Inclusion failures: {n_inclusion_fail}")
    print(f"Palindromic included total: {n_pal_included}")
    print(f"Allele mismatch included total: {n_mismatch_included}")

if __name__ == "__main__":
    main()
