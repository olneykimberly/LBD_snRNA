# LBD_snRNA
Lewy Body Dementia Center With Out Walls (LBD CWOW) single nucleus RNAseq processing and analysis.
Anterior cingulate cortex tissue samples from the Mayo Clinic brain bank were collected for 619 individuals. This single nuclues RNAseq data set is for a subset of those samples; N=40. 
Raw fastq files are available on SRA, PRJNA1023207.

5 per sex per disease type. 

| Disease                   | Count   |
| ------------------------- |:-------:|
| Control                     | 10    |
| Alzheimer’s disease (AD)    | 10    |
| Lewy body disease (LBD(S))  | 10    |
| Lewy body disease (LBD(AS)) | 10    |
| Lewy body disease (LBD(ATS))| 10    |

This git repo contains scripts for the following:
-   Metadata analysis
-   Per processing of single nucleus RNA-sequencing data
-   Analysis of single nucleus RNA-sequencing data 
-   Generation of manuscript figures from Olney et al. 2025 publication 
-   Generation of shiny app for exploration of the results presented in Olney et al. 2025 publication, view app [here](https://fryerlab.shinyapps.io/LBD_snRNA/)


*Characterize cell-type transcriptional alterations across neuropathologies:* Dementia pathologies elicit pronounced transcriptional responses in multiple cell types, including damaged and dying neurons, disease-associated microglia (DAM), and reactive disease-associated astrocytes (DAA). In Lewy body disease (LBD), increased SNCA (α-synuclein) expression leads to Lewy bodies' formation, particularly affecting dopaminergic neurons. In Alzheimer's disease (AD), elevated expression of amyloid precursor protein (APP) produces amyloid-β, a key component of amyloid plaques.The transcription patterns of cell types in cases with a combination of Aβ plaques, NFTs, and α-synuclein inclusions are not well understood. Moreover, studies have shown female-specific upregulation of DAM genes attributed to sex differences in microglial responses in AD cases. This project addresses two aims using human single nuclues RNAseq data: 

Aim 1.1 generates and characterizes snRNAseq data (5 per sex per neuropathology; N=40) to test the hypothesis that each neuropathological defined disease exhibits distinct cell-type transcriptional dysregulation.\

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

Additionally, the workflow uses cellranger and cellbender, which need to be installed separately.\
1. Install cellranger 9.0.0\
Requires going to the 10X website and filling in information to get to the download page.\
```
wget -O cellranger-9.0.0.tar.gz "[cellranger_link](https://cf.10xgenomics.com/releases/cell-exp/cellranger-9.0.0.tar.gz?Expires=1733980355&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=j3k~sgM2XOA2Bv0yWJe-QzOxuWrp07rAqTvWmCCCfV8o-SHEtlFp4W3a~Gp1XYaKU1P4c6BCX9YDNb18dJOVLP5FVUxkMV~6Vr9VWhyWAe~bu7TQGNV3N4k8cjnXowWaY8qtSXubiq2SYw~95TQjL4DZx58leAFfcvPEfbUbKjOyUk8tGZhjKF7M~qGkIYrZcbz8uzB0nH~nP5p0XiLwUQgspyOhrN4xdufLVhguoUrMG25u8tWpUUfy~uxmo2z-CA4~YIcrvj7dGQzjMJ6NNQkdZDyOvfRCdcV2cVfIIvc2wNt50zAk0taqy21afWj0MWzEFOQ0eyrsUlONsiuSIA__)"
```
2. Install cellbender as a seperete conda environment, as recommeneded by the developers. [cellbender_docs](https://cellbender.readthedocs.io/en/latest/installation/index.html)
```
conda env create -n cellbender --file cellbender.yml
```

3. Obtain Human reference from the 10X downloads page 
```
wget "https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz"
```

## File set up
Samples were sequenced at The Translational Genomics Research Institute (TGen) in Phoenix, Arizona. Files were named like so:
Example: MAYO_0407_1_BR_Nuclei_C1_X3SC3_A16845_22KJLTLT4_AGCAAGAAGC_L007_R1_001.fastq.gz\
[0] study; [1] patient; [2] visit (1 for frozen banked samples); [3] source (2 letter code; BR = brain); [4] fraction (Nuclie or whole); [5] subgroup - C1 is control as they don't know the disease status; [6] assay; [7] library (1 letter and 5 numbers); [8] sequence order - was [8:9] barcode; [9] lane; [10] strand or index

Because of the file naming provided by TGen, we will need to rename the files to be compatiable with cellranger which requires the sequence order to be in the name of the file. First we will obtain the sequencing order, which we have since this was defined by us. This information is the sequence_order.txt file. 
```
# inspect the file 
cat sequence_order.txt
```

We will create a rename map file to match the TGen filenames to the sequence order. 
```
# create rename map file 
python create_rename_map.py
# output is rename_map.txt
```

Finally we run the rename script to rename the fastq files using the rename_map file created in the step above. 
```
# rename the fastq files using the rename_map.txt
sh rename_fastqs.sh
```

## Snakemake for data alignment and removal of ambeint RNA
1. get read group info
To obtain sample information, run the get readgroup script.
```
sh get_readgroup_info.sh
```

2. create config.json file 
```
python create_config.py
```

3. run the snakemake file. 
