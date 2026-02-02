#!/bin/sh
#SBATCH --signal=USR2
#SBATCH --ntasks=1
#SBATCH --mem=25G
#SBATCH --cpus-per-task=2
#SBATCH --time=24:00:00
#SBATCH --output=/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/06_post_DE_analysis/metascape_hdWGCNA_celltype.job.%j
#SBATCH --job-name=metascape_hdWGCNA_celltype.job # Job name
#SBATCH --mail-type=NONE     # Mail events (NONE, BEGIN, END, FAIL, ALL)


sh 05_run_metascape_loop_WGCNA.sh

