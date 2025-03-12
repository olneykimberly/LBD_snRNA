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

# 1) get read information
#sh 01_sample_read_info.sh

# 2) create config
#python 02_create_10X_config.py


#cd /tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/cellbender
#cd BR_Nuclei_0385
#cellbender remove-background --input ../../cellranger/BR_Nuclei_0385/outs/raw_feature_bc_matrix.h5 --output BR_Nuclei_0385 --checkpoint ckpt.tar.gz

# 3) run snakemake - metaphlan alignment 
snakemake -s Snakefile -j 40 --nolock --latency-wait 15 --cluster "sbatch --ntasks 1 --cpus-per-task=8 --mem=64G --time=60:00:00"
#snakemake -s Snakefile -j 32 --nolock --latency-wait 15 --rerun-incomplete --cluster "sbatch --ntasks 1 --cpus-per-task=2 --mem=4G --time=00:20:00"
