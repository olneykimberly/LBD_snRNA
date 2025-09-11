#!/bin/bash
#SBATCH --job-name=cwow_sn_alignment                        
#SBATCH --time=48:00:00 # 8 hours   
#SBATCH -n 4 # threaded 
#SBATCH --mem=4G 
#SBATCH -o slurm.cwow_sn_alignment.out
#SBATCH -e slurm.cwow_sn_alignment.err

# source your bach profile to get your specific settings  
source $HOME/.bash_profile

module load python3
conda activate LBD_sn
#conda activate cellbender 

#------ Run steps 1 & 2 if you have not created the config file. 
# 1) get read information
#sh 01_sample_read_info.sh

# 2) create config
#python 02_create_10X_config.py

#------
# 3) run snakemake for cellranger and cellbender
snakemake -s Snakefile -j 40 --nolock --latency-wait 15 --cluster "sbatch --ntasks 1 --cpus-per-task=8 --mem=200000 --gres=gpu:1 -p gpu-a100 --time=8:00:00"
