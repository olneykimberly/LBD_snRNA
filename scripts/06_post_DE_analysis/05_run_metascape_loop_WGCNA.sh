#!/bin/bash
module load singularity
# --- 1. Define Base Paths ---
# Base directory for the msbio tool and configuration files (HOST PATH)
MSBIO_BASE_DIR="/tgen_labs/jfryer/kolney/tools/msbio_v3.5.20250701"

# Base directory where input gene list files are located (HOST PATH)
INPUT_DIR="/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/results/metascape_hdWGCNA_modules/"

# Base directory where all Metascape results will be saved (HOST PATH)
OUTPUT_BASE_DIR="/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/results/metascape_hdWGCNA_modules/"

# Name of the Metascape script inside the container (CONTAINER PATH)
CONTAINER_SCRIPT="/msbio/mylib/ms/smsbio2.sh"

# Name of the Docker image to run
CONTAINER_IMAGE="docker://metadocker8/msbio2"

# Ensure the main output directory exists on the host
mkdir -p "${OUTPUT_BASE_DIR}"

# --- 2. Define File Pattern ---
# Find all input files matching the pattern 'pigs_*.txt' in the INPUT_DIR
INPUT_FILES=$(find "${INPUT_DIR}" -maxdepth 1 -type f -name 'hdWGCNA_*.txt' -print)

# --- 3. Loop and Execute ---
echo "Starting Metascape analysis loop..."

for INPUT_FILE_PATH in ${INPUT_FILES}; do
    # Extract the filename (e.g., pigs_black.txt)
    FILENAME=$(basename "${INPUT_FILE_PATH}")
    
    # Extract the base name for the output directory (e.g., pigs_black)
    BASE_NAME="${FILENAME%.txt}"
    
    # Define the unique output directory for this job (HOST PATH)
    JOB_OUTPUT_DIR="${OUTPUT_BASE_DIR}/${BASE_NAME}"

    # Define the Host path for the log files/files that map to /data
    JOB_DATA_DIR="${JOB_OUTPUT_DIR}/msbio_data"
    
    # Create the unique output and log directories on the host
    mkdir -p "${JOB_OUTPUT_DIR}"
    mkdir -p "${JOB_DATA_DIR}/logs"

    # Create an empty log.txt file as required
    touch "${JOB_DATA_DIR}/log.txt"
    
    echo "--- Processing ${FILENAME} ---"
    echo "Output: ${JOB_OUTPUT_DIR}"
    
    # --- 4. Singularity Execution (Fixing Line Continuation) ---
    # NOTE: The backslash (\) MUST be the absolute last character on each line.
    
    singularity exec \
        --cleanenv \
        -B /home,/scratch,/tgen_labs \
        \
        -B "${MSBIO_BASE_DIR}/license/license.txt":/license/license.txt \
        -B "${MSBIO_BASE_DIR}/license/pubkey.pem":/license/pubkey.pem \
        \
        -B "${JOB_OUTPUT_DIR}":/data/out \
        -B "${JOB_DATA_DIR}/logs":/data/logs \
        -B "${JOB_DATA_DIR}/log.txt":/data/log.txt \
        \
        -B "${MSBIO_BASE_DIR}":"${MSBIO_BASE_DIR}" \
        -B "${INPUT_DIR}":"${INPUT_DIR}" \
        \
        "${CONTAINER_IMAGE}" \
        "${CONTAINER_SCRIPT}" \
        -u \
        -o /data/out \
        "${INPUT_FILE_PATH}"

    # Print a separator for clarity
    echo "--- ${FILENAME} COMPLETE ---"
done

echo "Metascape loop finished."

