#!/usr/bin/env bash
set -eu

OUTDIR="../../04_standardized_gwas/core_outcomes"
QCOUT="$OUTDIR/phase1_standardized_file_qc.txt"

F_NONIOP="$OUTDIR/GBS_nonIOPcomponent_main_Figure2.standardized.tsv.gz"
F_IOP="$OUTDIR/GBS_IOPcomponent_Supplementary_Figure_4b.standardized.tsv.gz"
F_IOPCC="$OUTDIR/IOPcc_Supplementary_Figure_4a.standardized.tsv.gz"

for f in "$F_NONIOP" "$F_IOP" "$F_IOPCC"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: standardized file not found: $f"
    exit 1
  fi
done

{
echo "===== Standardized files ====="
ls -lh "$F_NONIOP" "$F_IOP" "$F_IOPCC"

echo
echo "===== Standardized file headers ====="

echo
echo "GBS non-IOP header:"
gzip -cd "$F_NONIOP" | awk 'NR==1 {print; exit}' || true

echo
echo "GBS IOP header:"
gzip -cd "$F_IOP" | awk 'NR==1 {print; exit}' || true

echo
echo "IOPcc header:"
gzip -cd "$F_IOPCC" | awk 'NR==1 {print; exit}' || true

echo
echo "===== Standardized SNP counts ====="

echo -n "GBS non-IOP standardized SNP count: "
gzip -cd "$F_NONIOP" | wc -l | awk '{print $1-1}'

echo -n "GBS IOP standardized SNP count: "
gzip -cd "$F_IOP" | wc -l | awk '{print $1-1}'

echo -n "IOPcc standardized SNP count: "
gzip -cd "$F_IOPCC" | wc -l | awk '{print $1-1}'

echo
echo "===== Bad row counts by expected field number ====="

echo -n "GBS non-IOP bad rows, expecting 12 fields: "
gzip -cd "$F_NONIOP" | awk 'NR>1 && NF!=12 {bad++} END {print bad+0}'

echo -n "GBS IOP bad rows, expecting 12 fields: "
gzip -cd "$F_IOP" | awk 'NR>1 && NF!=12 {bad++} END {print bad+0}'

echo -n "IOPcc bad rows, expecting 13 fields: "
gzip -cd "$F_IOPCC" | awk 'NR>1 && NF!=13 {bad++} END {print bad+0}'

echo
echo "===== IOPcc P=0 retention check ====="

echo -n "IOPcc pval == 0 count: "
gzip -cd "$F_IOPCC" | awk 'NR>1 && $10==0 {n++} END {print n+0}'

echo -n "IOPcc pval_for_log10 == 2.9000001e-39 count: "
gzip -cd "$F_IOPCC" | awk 'NR>1 && $11=="2.9000001e-39" {n++} END {print n+0}'

echo
echo "===== First 3 data rows ====="

echo
echo "GBS non-IOP:"
gzip -cd "$F_NONIOP" | awk 'NR<=4 {print}' || true

echo
echo "GBS IOP:"
gzip -cd "$F_IOP" | awk 'NR<=4 {print}' || true

echo
echo "IOPcc:"
gzip -cd "$F_IOPCC" | awk 'NR<=4 {print}' || true

} > "$QCOUT"

echo "Wrote: $QCOUT"
cat "$QCOUT"
