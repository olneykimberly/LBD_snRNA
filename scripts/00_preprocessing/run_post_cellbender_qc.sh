#!/bin/bash
#SBATCH --job-name=cellranger                         
##SBATCH --nodes=1                                     
##SBATCH --tasks=32                                      
#SBATCH --time=48:00:00 # 8 hours   
#SBATCH -n 4 # threaded 
#SBATCH --mem=4G 
#SBATCH -o slurm.cellbender.out
#SBATCH -e slurm.cellbender.err

# source your bach profile to get your specific settings  
source $HOME/.bash_profile

module load python3
conda activate cellbender #LBD_sn

cd /tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts

R
setwd("/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts")
source("/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/post_cellbender_QC.R")
Rscript /tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/post_cellbender_QC.R

