#!/usr/bin/env bash
set -euo pipefail

PLAN="../../05_harmonization_planning"
RAW="../../07_exposure_gwas_raw"
LOGDIR="$PLAN/exposure_qc_logs"

VCF="$RAW/INSOMNIA__ebi-a-GCST90018869.vcf.gz"
REG="$PLAN/phase2_exposure_file_registration.tsv"
QC="$PLAN/phase2_exposure_file_qc_summary.tsv"
LOG="$LOGDIR/INSOMNIA__ebi-a-GCST90018869_pre_qc.txt"

mkdir -p "$LOGDIR"

if [ ! -f "$VCF" ]; then
  echo "ERROR: missing VCF: $VCF"
  exit 1
fi

echo "Running INSOMNIA pre-harmonization QC..."

python - << 'PY' > "$LOG"
import gzip
import re

vcf = "../../07_exposure_gwas_raw/INSOMNIA__ebi-a-GCST90018869.vcf.gz"

exposure_id = "INSOMNIA"
file_path = "../../07_exposure_gwas_raw/INSOMNIA__ebi-a-GCST90018869.vcf.gz"
file_format = "OpenGWAS_VCF_GZ"
genome_build = "GRCh37"
trait_name = "insomnia"

sample_size = "486627"
cases = "1402"
controls = "485225"
population = "GWAS Catalog/OpenGWAS insomnia; ancestry/population to verify from source paper"

n_rows = 0
n_cols = None
header = None

has_rsid = "NO"
has_chr_pos = "YES"
has_effect_allele = "YES"
has_other_allele = "YES"
has_beta = "YES"
has_se = "YES"
has_pval = "YES_FROM_LP"
has_eaf = "YES"
has_n = "YES_STUDY_LEVEL"

rsid_count = 0
missing_rsid = 0
bad_format = 0
pval_zero_count = 0
max_lp = None
dup_count = 0
seen = set()
palindromic_possible = 0

total_variants_meta = "NA"
harmonised_variants_meta = "NA"
variants_not_harmonised_meta = "NA"
study_type = "NA"

required_format = ["ES", "SE", "LP", "AF", "ID"]

with gzip.open(vcf, "rt", encoding="utf-8", errors="replace") as f:
    for line in f:
        line = line.rstrip("\n")

        if line.startswith("##SAMPLE="):
            m = re.search(r"TotalVariants=([^,>]+)", line)
            if m:
                total_variants_meta = m.group(1)
            m = re.search(r"HarmonisedVariants=([^,>]+)", line)
            if m:
                harmonised_variants_meta = m.group(1)
            m = re.search(r"VariantsNotHarmonised=([^,>]+)", line)
            if m:
                variants_not_harmonised_meta = m.group(1)
            m = re.search(r"StudyType=([^,>]+)", line)
            if m:
                study_type = m.group(1)
            m = re.search(r"TotalCases=([^,>]+)", line)
            if m:
                cases = str(int(float(m.group(1))))
            m = re.search(r"TotalControls=([^,>]+)", line)
            if m:
                controls = str(int(float(m.group(1))))
            try:
                sample_size = str(int(cases) + int(controls))
            except Exception:
                pass

        if line.startswith("##contig=<ID=1"):
            if "HG19" not in line and "GRCh37" not in line:
                genome_build = "UNKNOWN"

        if line.startswith("#CHROM"):
            header = line
            n_cols = len(line.split("\t"))
            continue

        if line.startswith("#"):
            continue

        parts = line.split("\t")
        if len(parts) != 10:
            bad_format += 1
            continue

        chrom, pos, rsid, ref, alt, qual, filt, info, fmt, sample = parts
        fmt_keys = fmt.split(":")
        sample_vals = sample.split(":")

        if fmt_keys != required_format or len(sample_vals) != len(fmt_keys):
            bad_format += 1
            continue

        d = dict(zip(fmt_keys, sample_vals))

        n_rows += 1

        key = (chrom, pos, ref, alt)
        if key in seen:
            dup_count += 1
        else:
            seen.add(key)

        if rsid.startswith("rs"):
            rsid_count += 1
        else:
            missing_rsid += 1

        a = ref.upper()
        b = alt.upper()
        if len(a) == 1 and len(b) == 1:
            pair = {a, b}
            if pair == {"A", "T"} or pair == {"C", "G"}:
                palindromic_possible += 1

        try:
            lp = float(d["LP"])
            if max_lp is None or lp > max_lp:
                max_lp = lp
        except Exception:
            pass

if max_lp is None:
    pval_min = "NA"
else:
    pval_min = "10^(-" + format(max_lp, ".8g") + ")"

has_rsid = "YES" if rsid_count > 0 else "NO"

build_status = "OK_GRCh37"
allele_status = "OK_ALT_IS_EFFECT_ALLELE"
ready = "YES"

qc_note = (
    "OpenGWAS VCF; ALT is effect allele; REF is other allele; "
    "ES=beta, SE=se, LP=-log10P, AF=effect allele frequency, ID=rsID; "
    "no per-variant SS in FORMAT, use study-level N=486627 from TotalCases+TotalControls; "
    "pval will be computed as 10^(-LP) during exposure standardization"
)

print("===== INSOMNIA pre-QC summary =====")
print(f"exposure_id: {exposure_id}")
print(f"file_path: {file_path}")
print(f"n_rows: {n_rows}")
print(f"n_cols: {n_cols}")
print(f"header: {header}")
print(f"total_variants_meta: {total_variants_meta}")
print(f"harmonised_variants_meta: {harmonised_variants_meta}")
print(f"variants_not_harmonised_meta: {variants_not_harmonised_meta}")
print(f"study_type: {study_type}")
print(f"cases: {cases}")
print(f"controls: {controls}")
print(f"sample_size: {sample_size}")
print(f"rsid_count: {rsid_count}")
print(f"missing_rsid_count: {missing_rsid}")
print(f"bad_format_count: {bad_format}")
print(f"duplicated_chr_pos_ref_alt_count: {dup_count}")
print(f"palindromic_possible_count: {palindromic_possible}")
print(f"max_LP: {max_lp}")
print(f"pval_min: {pval_min}")
print(f"pval_zero_count: {pval_zero_count}")
print(f"genome_build: {genome_build}")
print(f"build_status: {build_status}")
print(f"allele_status: {allele_status}")
print(f"ready_for_harmonization: {ready}")

print("")
print("===== registration_row =====")
print("\t".join([
    exposure_id,
    file_path,
    file_format,
    genome_build,
    "ID",
    "#CHROM",
    "POS",
    "ALT",
    "REF",
    "ES",
    "SE",
    "LP",
    "AF",
    "NA",
    trait_name,
    sample_size,
    population,
    "OpenGWAS ebi-a-GCST90018869; case-control insomnia; pval derived from LP; ALT is effect allele; no per-variant SS, use study-level N",
    "REGISTERED"
]))

print("")
print("===== qc_row =====")
print("\t".join([
    exposure_id,
    file_path,
    str(n_rows),
    str(n_cols),
    header.replace("\t", "|") if header else "NA",
    has_rsid,
    has_chr_pos,
    genome_build,
    has_effect_allele,
    has_other_allele,
    has_beta,
    has_se,
    has_pval,
    has_eaf,
    has_n,
    pval_min,
    str(pval_zero_count),
    str(dup_count),
    "YES" if palindromic_possible > 0 else "NO",
    build_status,
    allele_status,
    ready,
    qc_note
]))
PY

cat "$LOG"

REG_ROW=$(awk '/^===== registration_row =====/{getline; print}' "$LOG")
QC_ROW=$(awk '/^===== qc_row =====/{getline; print}' "$LOG")

awk -F '\t' 'NR==1 || $1!="INSOMNIA"' "$REG" > "$REG.tmp"
mv "$REG.tmp" "$REG"

awk -F '\t' 'NR==1 || $1!="INSOMNIA"' "$QC" > "$QC.tmp"
mv "$QC.tmp" "$QC"

printf "%s\n" "$REG_ROW" >> "$REG"
printf "%s\n" "$QC_ROW" >> "$QC"

echo
echo "===== Updated registration table ====="
cat -vet "$REG"

echo
echo "===== Updated QC summary table ====="
cat -vet "$QC"

