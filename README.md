# LBD_snRNA
Characterize cell-type transcriptional alterations across neuropathologies.\
Dementia pathologies elicit pronounced transcriptional responses in multiple cell types, including damaged and dying neurons, disease-associated microglia (DAM), and reactive disease-associated astrocytes (DAA). In Lewy body disease (LBD), increased SNCA (α-synuclein) expression leads to Lewy bodies' formation, particularly affecting dopaminergic neurons. In Alzheimer's disease (AD), elevated expression of amyloid precursor protein (APP) produces amyloid-β, a key component of amyloid plaques.\
The transcription patterns of cell types in cases with a combination of Aβ plaques, NFTs, and α-synuclein inclusions are not well understood. Thus this project explores:\
  Aim 1.1 generates and characterizes snRNAseq data (5 per sex per neuropathology; N=40) to test the hypothesis that each neuropathological defined disease exhibits distinct cell-type transcriptional dysregulation.\
Moreover, studies have shown female-specific upregulation of DAM genes attributed to sex differences in microglial responses in AD cases.\
Aim 1.2 will examine transcriptional sex differences among cell types, testing the hypothesis that XX females display more pronounced upregulation of DAM and DAA genes than XY males with similar neuropathological burdens.

Explore genes among the human brain samples in our published [shiny app](https://fryerlab.shinyapps.io/LBD_snRNA/)


## Set up conda environment
This workflow uses conda. For information on how to install conda [here](https://docs.conda.io/projects/conda/en/latest/user-guide/index.html)

To create the environment:
```
conda env create -n LBD_sn --file LBD_sn.yml

# To activate this environment, use
#
#     $ conda activate LBD_sn
#
# To deactivate an active environment, use
#
#     $ conda deactivate LBD_sn
```

# Set up
conda activate LBD_sn
1. install cellranger 9.0.0
Requires going to the 10X website and filling in information to get to the download page
wget -O cellranger-9.0.0.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-9.0.0.tar.gz?Expires=1733980355&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=j3k~sgM2XOA2Bv0yWJe-QzOxuWrp07rAqTvWmCCCfV8o-SHEtlFp4W3a~Gp1XYaKU1P4c6BCX9YDNb18dJOVLP5FVUxkMV~6Vr9VWhyWAe~bu7TQGNV3N4k8cjnXowWaY8qtSXubiq2SYw~95TQjL4DZx58leAFfcvPEfbUbKjOyUk8tGZhjKF7M~qGkIYrZcbz8uzB0nH~nP5p0XiLwUQgspyOhrN4xdufLVhguoUrMG25u8tWpUUfy~uxmo2z-CA4~YIcrvj7dGQzjMJ6NNQkdZDyOvfRCdcV2cVfIIvc2wNt50zAk0taqy21afWj0MWzEFOQ0eyrsUlONsiuSIA__"

2. Human reference 
wget "https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz"

3. get sample information
sh get_readgroup_info.sh


#----- rename files 
1. sequence order S1-S40. There are 40 samples
sequence_order.txt

2. create rename map file 
python create_rename_map.py
output is rename_map.txt

3. rename the fastq files using the rename_map.txt
sh rename_fastqs.sh

#----- Snakemake
1. get read group info
sh get_readgroup_info.sh

2. create config.json file 
python create_config.py

# Snakemake 
3. Alignment
salloc --ntasks=1 --cpus-per-task=8 --mem=64G --time=4:00:00

