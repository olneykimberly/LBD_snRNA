#!/bin/bash
#SBATCH --job-name=rsync                         
#SBATCH --time=08:00:00                               
#SBATCH --mem=24G
#SBATCH -n 6 # threaded 
##SBATCH -o slurm.MAST.out
##SBATCH -e slurm.MAST.err
##SBATCH --mail-user=olney.kimberly@mayo.edu

cd /tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts
rsync -avh --progress /tgen_labs/jfryer/projects/LBD_CWOW/snRNA/tgen_10x/fastq/ /tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/fastq/
