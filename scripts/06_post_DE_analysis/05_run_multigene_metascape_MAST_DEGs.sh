#!/bin/bash
set -euo pipefail
module load singularity

# --- Paths ---
MSBIO_BASE_DIR="/tgen_labs/jfryer/kolney/tools/msbio_v3.5.20250701"
XLSX_DIR="/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/results/metascape_xlsx_staging"
OUTPUT_BASE="/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/results/metascape_multilist_MAST_DEGs"

CONTAINER_SCRIPT="/msbio/mylib/ms/smsbio2.sh"
CONTAINER_IMAGE="docker://metadocker8/msbio2"

# --- sanity check ---
shopt -s nullglob
XLSX_FILES=("${XLSX_DIR}"/*.xlsx)

if [[ ${#XLSX_FILES[@]} -eq 0 ]]; then
  echo "ERROR: No .xlsx files found in ${XLSX_DIR}"
  exit 1
fi

echo "Found ${#XLSX_FILES[@]} XLSX files"
echo "------------------------------------------------------------"

for XLSX in "${XLSX_FILES[@]}"; do

  BASE_NAME="$(basename "${XLSX}" .xlsx)"
  JOB_OUTPUT_DIR="${OUTPUT_BASE}/${BASE_NAME}"
  JOB_DATA_DIR="${JOB_OUTPUT_DIR}/msbio_data"
  TMP_DIR="${JOB_DATA_DIR}/tmp"

  mkdir -p "${JOB_OUTPUT_DIR}"
  mkdir -p "${JOB_DATA_DIR}/logs"
  mkdir -p "${TMP_DIR}"
  : > "${JOB_DATA_DIR}/log.txt"

  echo "------------------------------------------------------------"
  echo "Metascape Analysis"
  echo "Input:  ${XLSX}"
  echo "Output: ${JOB_OUTPUT_DIR}"
  echo "Tmp:    ${TMP_DIR}"
  echo "------------------------------------------------------------"

  # Optional XLSX drawing check
  if command -v unzip >/dev/null 2>&1; then
    unzip -l "${XLSX}" | grep -E "xl/drawings|drawing" || true
  fi

  singularity exec \
    --cleanenv \
    -B /home,/scratch,/tgen_labs \
    -B "${TMP_DIR}":/tmp \
    -B "${MSBIO_BASE_DIR}/license/license.txt":/license/license.txt \
    -B "${MSBIO_BASE_DIR}/license/pubkey.pem":/license/pubkey.pem \
    -B "${JOB_OUTPUT_DIR}":/data/out \
    -B "${JOB_DATA_DIR}/logs":/data/logs \
    -B "${JOB_DATA_DIR}/log.txt":/data/log.txt \
    -B "${MSBIO_BASE_DIR}":"${MSBIO_BASE_DIR}" \
    "${CONTAINER_IMAGE}" \
    "${CONTAINER_SCRIPT}" \
    -S 9606 -T 9606 \
    -o /data/out \
    "${XLSX}"

  echo "DONE: ${BASE_NAME}"
done

echo "============================================================"
echo "ALL METASCAPE JOBS COMPLETED"
echo "============================================================"
