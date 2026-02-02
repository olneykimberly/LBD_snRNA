#!/bin/sh
#SBATCH --signal=USR2
#SBATCH --ntasks=1
#SBATCH --mem=25G
#SBATCH --cpus-per-task=2
#SBATCH --time=10:00:00
#SBATCH --output=/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/06_post_DE_analysis/metascape_multilist_DEGs_celltype.job.%j
#SBATCH --job-name=metascape_multilist_DEGs_celltype.job # Job name
#SBATCH --mail-type=NONE     # Mail events (NONE, BEGIN, END, FAIL, ALL)


sh 05_run_multigene_metascape_MAST_DEGs.sh

