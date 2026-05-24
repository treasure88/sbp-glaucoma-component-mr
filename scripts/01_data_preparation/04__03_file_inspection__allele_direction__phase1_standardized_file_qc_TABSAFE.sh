#!/usr/bin/env bash
set -euo pipefail

OUTDIR="../../04_standardized_gwas/core_outcomes"

F_NONIOP="$OUTDIR/GBS_nonIOPcomponent_main_Figure2.standardized.tsv.gz"
F_IOP="$OUTDIR/GBS_IOPcomponent_Supplementary_Figure_4b.standardized.tsv.gz"
F_IOPCC="$OUTDIR/IOPcc_Supplementary_Figure_4a.standardized.tsv.gz"

QC_OUT="$OUTDIR/phase1_standardized_file_qc_TABSAFE.txt"

{
echo "===== Standardized files ====="
ls -lh "$F_NONIOP" "$F_IOP" "$F_IOPCC"

echo
echo "===== Standardized file headers ====="

echo
echo "GBS non-IOP header:"
gzip -cd "$F_NONIOP" | sed -n '1p'

echo
echo "GBS IOP header:"
gzip -cd "$F_IOP" | sed -n '1p'

echo
echo "IOPcc header:"
gzip -cd "$F_IOPCC" | sed -n '1p'

echo
echo "===== Standardized SNP counts ====="

echo -n "GBS non-IOP standardized SNP count: "
gzip -cd "$F_NONIOP" | perl -ne '$n++ if $. > 1; END {print $n+0, "\n"}'

echo -n "GBS IOP standardized SNP count: "
gzip -cd "$F_IOP" | perl -ne '$n++ if $. > 1; END {print $n+0, "\n"}'

echo -n "IOPcc standardized SNP count: "
gzip -cd "$F_IOPCC" | perl -ne '$n++ if $. > 1; END {print $n+0, "\n"}'

echo
echo "===== Bad row counts by expected TAB-delimited field number ====="

echo -n "GBS non-IOP bad rows, expecting 12 fields: "
gzip -cd "$F_NONIOP" | perl -ne 'chomp; @f=split(/\t/,$_, -1); $bad++ if $.>1 && scalar(@f)!=12; END {print $bad+0, "\n"}'

echo -n "GBS IOP bad rows, expecting 12 fields: "
gzip -cd "$F_IOP" | perl -ne 'chomp; @f=split(/\t/,$_, -1); $bad++ if $.>1 && scalar(@f)!=12; END {print $bad+0, "\n"}'

echo -n "IOPcc bad rows, expecting 13 fields: "
gzip -cd "$F_IOPCC" | perl -ne 'chomp; @f=split(/\t/,$_, -1); $bad++ if $.>1 && scalar(@f)!=13; END {print $bad+0, "\n"}'

echo
echo "===== IOPcc P=0 retention check ====="

echo -n "IOPcc pval == 0 count: "
gzip -cd "$F_IOPCC" | perl -ne 'chomp; @f=split(/\t/,$_, -1); if ($.>1 && scalar(@f)==13 && $f[9] == 0) {$n++}; END {print $n+0, "\n"}'

echo -n "IOPcc pval_for_log10 == 2.9000001e-39 count: "
gzip -cd "$F_IOPCC" | perl -ne 'chomp; @f=split(/\t/,$_, -1); if ($.>1 && scalar(@f)==13 && $f[10] eq "2.9000001e-39") {$n++}; END {print $n+0, "\n"}'

echo
echo "===== IOPcc genome_build/source_dataset checks ====="

echo -n "IOPcc genome_build bad count: "
gzip -cd "$F_IOPCC" | perl -ne 'chomp; @f=split(/\t/,$_, -1); if ($.>1 && scalar(@f)==13 && $f[11] ne "GRCh37") {$bad++}; END {print $bad+0, "\n"}'

echo -n "IOPcc source_dataset bad count: "
gzip -cd "$F_IOPCC" | perl -ne 'chomp; @f=split(/\t/,$_, -1); if ($.>1 && scalar(@f)==13 && $f[12] ne "IOPcc_Supplementary_Figure_4a") {$bad++}; END {print $bad+0, "\n"}'

echo
echo "===== First 3 data rows ====="

echo
echo "GBS non-IOP:"
gzip -cd "$F_NONIOP" | sed -n '1,4p'

echo
echo "GBS IOP:"
gzip -cd "$F_IOP" | sed -n '1,4p'

echo
echo "IOPcc:"
gzip -cd "$F_IOPCC" | sed -n '1,4p'

} > "$QC_OUT"

echo "Wrote: $QC_OUT"
cat "$QC_OUT"
