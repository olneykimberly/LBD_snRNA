#!/bin/bash
#SBATCH --job-name=SRA                         
#SBATCH --time=24:00:00                               
#SBATCH --mem=24G
#SBATCH -n 6 # threaded 
##SBATCH -o slurm.MAST.out
##SBATCH -e slurm.MAST.err
##SBATCH --mail-user=olney.kimberly@mayo.edu

source $HOME/.bash_profile
module load python3

conda activate SRA


ascp -i /tgen_labs/jfryer/kolney/tools/aspera.openssh -QT -l 3000m -k1  -d  /tgen_labs/jfryer/projects/LBD_CWOW/snRNA/tgen_10x/fastq/ subasp@upload.ncbi.nlm.nih.gov:uploads/olneykimberly_gmail.com_Mivffuea/LBD_sn/
