#!/usr/bin/env bash
set -euo pipefail

RAWDIR="../../02_raw_downloads"
NONIOP="$RAWDIR/GBS_nonIOPcomponent_main_Figure2.zip"

OUT="phase0_nonIOP_table1_direction_build_check.txt"

{
echo "===== Local rows for selected Table 1 lead SNPs ====="
echo "Expected interpretation: A1 = effect allele; BETA = effect per copy of A1; POS = GRCh37."

echo
echo "Header and selected SNPs from local non-IOP file:"
unzip -p "$NONIOP" | awk '
NR==1 {print; next}
$3=="rs2790049" || \
$3=="rs6475604" || \
$3=="rs2472494" || \
$3=="rs10740731" || \
$3=="rs34935520" {print}
'

echo
echo "===== Expected rounded values from article Table 1 ====="
cat << 'TABLE'
SNP	CHR	POS	A1	A2	Expected_BETA_rounded	Expected_P_rounded
rs2790049	1	165743523	A	G	0.34	1.42E-21
rs6475604	9	22052734	T	C	-0.47	9.78E-75
rs2472494	9	107695539	T	C	0.22	7.12E-18
rs10740731	10	60348886	G	A	0.16	5.07E-10
rs34935520	14	61091401	G	A	0.28	1.12E-27
TABLE

echo
echo "Interpretation:"
echo "If CHR/POS/A1/A2 match the expected Table 1 rows and BETA/P are close to the rounded published values, this confirms local file direction and GRCh37 build for the non-IOP component."
} > "$OUT"

echo "Wrote: $OUT"
cat "$OUT"
