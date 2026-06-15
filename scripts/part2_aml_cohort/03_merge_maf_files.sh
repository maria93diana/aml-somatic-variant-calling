#!/usr/bin/env bash
# =============================================================================
# 03_merge_maf_files.sh
# Merge all per-sample MAF files into a single MAF for maftools analysis.
# Keeps header from first file, appends data rows from all others.
#
# Output: data/processed/TCGA-LAML-merged.maf
#
# Usage: bash scripts/part2_aml_cohort/03_merge_maf_files.sh
# =============================================================================

set -euo pipefail

INDIR="data/tcga_maf"
OUTDIR="data/processed"
OUTFILE="${OUTDIR}/TCGA-LAML-merged.maf"

mkdir -p "${OUTDIR}"

echo "[INFO] Merging MAF files from ${INDIR}..."

# Extract header from first MAF file
FIRST=$(find "${INDIR}" -name "*.maf.gz" | sort | head -1)

if [[ -z "${FIRST}" ]]; then
    echo "[ERROR] No MAF files found in ${INDIR}"
    exit 1
fi

echo "[INFO] Using header from: ${FIRST}"
gzcat "${FIRST}" | grep -E "^#|^Hugo_Symbol" > "${OUTFILE}"

# Append data rows from all files
N=0
find "${INDIR}" -name "*.maf.gz" | sort | while read -r f; do
    gzcat "${f}" | grep -v "^#" | grep -v "^Hugo_Symbol" >> "${OUTFILE}"
    N=$((N + 1))
done

TOTAL=$(grep -v "^#" "${OUTFILE}" | grep -v "^Hugo_Symbol" | wc -l)
echo "[INFO] Merged MAF written to ${OUTFILE}"
echo "[INFO] Total mutation records: ${TOTAL}"