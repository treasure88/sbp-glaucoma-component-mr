#!/usr/bin/env bash
set -euo pipefail

RAWDIR="../../02_raw_downloads"
OUTDIR="../../04_standardized_gwas/core_outcomes"

mkdir -p "$OUTDIR"

F_NONIOP="$RAWDIR/GBS_nonIOPcomponent_main_Figure2.zip"
F_IOP="$RAWDIR/GBS_IOPcomponent_Supplementary_Figure_4b.zip"
F_IOPCC="$RAWDIR/IOPcc_Supplementary_Figure_4a.zip"

OUT_NONIOP="$OUTDIR/GBS_nonIOPcomponent_main_Figure2.standardized.tsv.gz"
OUT_IOP="$OUTDIR/GBS_IOPcomponent_Supplementary_Figure_4b.standardized.tsv.gz"
OUT_IOPCC="$OUTDIR/IOPcc_Supplementary_Figure_4a.standardized.tsv.gz"

echo "Checking input files..."
for f in "$F_NONIOP" "$F_IOP" "$F_IOPCC"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: file not found: $f"
    exit 1
  fi
done

echo "Standardizing GBS non-IOP component..."
unzip -p "$F_NONIOP" | awk '
BEGIN {
  OFS="\t";
  print "chr","pos","SNP","effect_allele","other_allele","eaf_or_maf","beta","se","pval","pval_for_log10","genome_build","source_dataset"
}
NR==1 {next}
{
  p_for_log = ($9 == 0 ? "NA" : $9);
  print $1,$2,$3,$5,$6,$4,$7,$8,$9,p_for_log,"GRCh37","GBS_nonIOPcomponent_main_Figure2"
}
' | gzip -c > "$OUT_NONIOP"

echo "Standardizing GBS IOP component..."
unzip -p "$F_IOP" | awk '
BEGIN {
  OFS="\t";
  print "chr","pos","SNP","effect_allele","other_allele","eaf_or_maf","beta","se","pval","pval_for_log10","genome_build","source_dataset"
}
NR==1 {next}
{
  p_for_log = ($9 == 0 ? "NA" : $9);
  print $1,$2,$3,$5,$6,$4,$7,$8,$9,p_for_log,"GRCh37","GBS_IOPcomponent_Supplementary_Figure_4b"
}
' | gzip -c > "$OUT_IOP"

echo "Standardizing IOPcc..."
unzip -p "$F_IOPCC" | awk '
BEGIN {
  OFS="\t";
  zero_p_replacement="2.9000001e-39";
  print "chr","pos","variant_id","SNP","effect_allele","other_allele","info","beta","se","pval","pval_for_log10","genome_build","source_dataset"
}
NR==1 {next}
{
  variant=$1;
  split(variant, a, ":");
  chr=a[1];
  rest=a[2];
  split(rest, b, "_");
  pos=b[1];

  p_for_log = ($7 == 0 ? zero_p_replacement : $7);

  print chr,pos,variant,"NA",$2,$3,$4,$5,$6,$7,p_for_log,"GRCh37","IOPcc_Supplementary_Figure_4a"
}
' | gzip -c > "$OUT_IOPCC"

cat > "$OUTDIR/phase1_standardization_metadata.tsv" << 'META'
dataset	raw_file	standardized_file	genome_build	allele_direction	sample_metadata	notes
GBS_nonIOPcomponent_main_Figure2	GBS_nonIOPcomponent_main_Figure2.zip	GBS_nonIOPcomponent_main_Figure2.standardized.tsv.gz	GRCh37 / hg19	A1 is effect allele; BETA is effect per copy of A1	Derived GWAS-by-subtraction trait; source POAG 14853 cases / 106544 controls; source IOPcc n=97644	MAF retained as eaf_or_maf; no per-SNP N column
GBS_IOPcomponent_Supplementary_Figure_4b	GBS_IOPcomponent_Supplementary_Figure_4b.zip	GBS_IOPcomponent_Supplementary_Figure_4b.standardized.tsv.gz	GRCh37 / hg19	A1 treated as effect allele using same pipeline convention as non-IOP component	Derived GWAS-by-subtraction trait; source POAG 14853 cases / 106544 controls; source IOPcc n=97644	MAF retained as eaf_or_maf; no per-SNP N column
IOPcc_Supplementary_Figure_4a	IOPcc_Supplementary_Figure_4a.zip	IOPcc_Supplementary_Figure_4a.standardized.tsv.gz	GRCh37 / hg19	a1 treated as tested/effect allele; beta is linear regression effect per copy of a1	UK Biobank IOPcc GWAS; n=97644	No rsID, no EAF/MAF; chr/pos parsed from variant_id; 52 P=0 values retained in pval and replaced only in pval_for_log10
META

echo "Standardization complete."
echo "Output directory: $OUTDIR"
ls -lh "$OUTDIR"
