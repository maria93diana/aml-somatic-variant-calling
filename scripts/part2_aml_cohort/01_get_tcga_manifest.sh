#!/usr/bin/env bash
# =============================================================================
# 01_get_tcga_manifest.sh
# Download list of TCGA-LAML open access somatic mutation file UUIDs from GDC.
#
# Output: data/tcga_maf/manifest_uuids.txt (one UUID per line, 153 files)
#
# Usage: bash scripts/part2_aml_cohort/01_get_tcga_manifest.sh
# =============================================================================

set -euo pipefail

OUTDIR="data/tcga_maf"
OUTFILE="${OUTDIR}/manifest_uuids.txt"

mkdir -p "${OUTDIR}"

echo "[INFO] Querying GDC API for TCGA-LAML masked somatic mutation files..."

curl -s "https://api.gdc.cancer.gov/files?filters=%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B%22op%22%3A%22%3D%22%2C%22content%22%3A%7B%22field%22%3A%22cases.project.project_id%22%2C%22value%22%3A%22TCGA-LAML%22%7D%7D%2C%7B%22op%22%3A%22%3D%22%2C%22content%22%3A%7B%22field%22%3A%22data_type%22%2C%22value%22%3A%22Masked+Somatic+Mutation%22%7D%7D%2C%7B%22op%22%3A%22%3D%22%2C%22content%22%3A%7B%22field%22%3A%22access%22%2C%22value%22%3A%22open%22%7D%7D%5D%7D&fields=file_id&size=200" \
| python3 -c "
import json, sys
data = json.load(sys.stdin)
for hit in data['data']['hits']:
    print(hit['file_id'])
" > "${OUTFILE}"

N=$(wc -l < "${OUTFILE}")
echo "[INFO] Found ${N} files. UUIDs saved to ${OUTFILE}"