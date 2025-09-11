# LBD_snRNA
Lewy Body Dementia Center With Out Walls (LBD CWOW) single nucleus RNAseq processing and analysis.
Anterior cingulate cortex tissue samples from the Mayo Clinic brain bank were collected for 619 individuals. This single nucleus RNAseq data set is for a subset of those samples; n=35. 
Raw fastq files are available on SRA, PRJNA1023207.

5 per sex per disease type. 

| Disease                   | Count   | Sex (F, M)|
| ------------------------- |:-------:|:-------:|
| Control                     | 7    | (3,4)   |
| Alzheimer’s disease (AD)    | 8    | (4,4)   |
| Lewy body disease (LBD(S))  | 7    | (3,4)   |
| Lewy body disease (LBD(AS)) | 7    | (3,4)   |
| Lewy body disease (LBD(ATS))| 6    | (2,4)   |

This git repo contains scripts for the following:
-   Metadata analysis
-   Per processing of single nucleus RNA-sequencing data
-   Analysis of single nucleus RNA-sequencing data 
-   Generation of manuscript figures from Olney et al. 2025 publication 
-   Generation of shiny app for exploration of the results, view app [here](https://fryerlab.shinyapps.io/lbd_cwow_snrna/)


*Characterize cell-type transcriptional alterations across neuropathologies:* Dementia pathologies elicit pronounced transcriptional responses in multiple cell types, including damaged and dying neurons, disease-associated microglia (DAM), and reactive disease-associated astrocytes (DAA). In Lewy body disease (LBD), increased SNCA (α-synuclein) expression leads to Lewy bodies' formation, particularly affecting dopaminergic neurons. In Alzheimer's disease (AD), elevated expression of amyloid precursor protein (APP) produces amyloid-β, a key component of amyloid plaques.The transcription patterns of cell types in cases with a combination of Aβ plaques, NFTs, and α-synuclein inclusions are not well understood. Moreover, studies have shown female-specific upregulation of DAM genes attributed to sex differences in microglial responses in AD cases. This project addresses two aims using human single nuclues RNAseq data: 

Aim 1 generates and characterizes snRNAseq data (~4 per sex per neuropathology; n=35) to test the hypothesis that each neuropathological defined disease exhibits distinct cell-type transcriptional dysregulation.

Aim 2 examines transcriptional sex differences among cell types, testing the hypothesis that XX females display more pronounced upregulation of DAM and DAA genes than XY males with similar neuropathological burdens.

Explore bulk-level gene differential expression among the human brain samples in our published [shiny app](https://fryerlab.shinyapps.io/lbd_cwow_snrna/)


## Set up conda environment
This workflow uses conda. For information on how to install conda [here](https://docs.conda.io/projects/conda/en/latest/user-guide/index.html)

To create the conda environment:
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
The LBD_sn conda environment contains cellbender.  [cellbender_docs](https://cellbender.readthedocs.io/en/latest/installation/index.html)

Additionally, the workflow uses cellranger, which needs to be installed separately.
1. Install cellranger 9.0.0\
Requires going to the 10X website and filling in information to get to the download page.
```
wget -O cellranger-9.0.0.tar.gz "[cellranger_link](https://cf.10xgenomics.com/releases/cell-exp/cellranger-9.0.0.tar.gz?Expires=1733980355&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=j3k~sgM2XOA2Bv0yWJe-QzOxuWrp07rAqTvWmCCCfV8o-SHEtlFp4W3a~Gp1XYaKU1P4c6BCX9YDNb18dJOVLP5FVUxkMV~6Vr9VWhyWAe~bu7TQGNV3N4k8cjnXowWaY8qtSXubiq2SYw~95TQjL4DZx58leAFfcvPEfbUbKjOyUk8tGZhjKF7M~qGkIYrZcbz8uzB0nH~nP5p0XiLwUQgspyOhrN4xdufLVhguoUrMG25u8tWpUUfy~uxmo2z-CA4~YIcrvj7dGQzjMJ6NNQkdZDyOvfRCdcV2cVfIIvc2wNt50zAk0taqy21afWj0MWzEFOQ0eyrsUlONsiuSIA__)"
```

2. Obtain Human reference from the 10X downloads page 
```
wget "https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz"
```

## File name set up
Samples were sequenced at The Translational Genomics Research Institute (TGen) in Phoenix, Arizona. Files were named like so:
MAYO_0407_1_BR_Nuclei_C1_X3SC3_A16845_22KJLTLT4_AGCAAGAAGC_L007_R1_001.fastq.gz\
[0] study; [1] patient; [2] visit (1 for frozen banked samples); [3] source (2 letter code; BR = brain); [4] fraction (Nuclie or whole); [5] subgroup - C1 is control as they don't know the disease status; [6] assay; [7] library (1 letter and 5 numbers); [8] sequence order - was [8:9] barcode; [9] lane; [10] strand or index

Because of the file naming provided by TGen, we will need to rename the files to be compatible with cellranger which requires the sequence order to be in the name of the file. First we will obtain the sequencing order, which we have since this was defined by us. This information is the sequence_order.txt file. 
First change directory to be in the scripts preprocessing folder. 
```
cd scripts/00_preprocessing/
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

## Snakemake for data alignment and removal/correction of ambeint RNA
1. Obtain read group info from the fastq files.\
To obtain sample information, run the get_readgroup_info script. This script is located in the scripts preprocessing folder. 
```
sh get_readgroup_info.sh
```

2. Create config.json file 
```
python create_config.py
```

3. Run the snakemake file. 
TGen uses sbatch. This script will need to be updated based on your institution's High-Performance Computing (HPC) system and job scheduler.
```
cd scripts/01_alignment/
sbatch run_snakemake.sh
```
Post alignment will include both cellranger and cellbender outputs. 

## Set up R environment
The next steps in the workflow uses R. The R environment has been preserved using renv:  [renv_docs](https://rstudio.github.io/renv/articles/renv.html)
When cloning this git repo, renv should automatically activate. If not, you can manually call renv::activate().
You can then run renv::restore(). This command reads the renv.lock file and installs the exact versions of the R packages specified in the lockfile into your own private project library, ensuring a consistent environment.

## Quality control 
The scripts for QC are in the scripts/02_qc/directory. 

1. Create Seurat object. **Note that this workflow uses Seurat v5.  
```
# cd scripts/02_qc/
01_create_seurat_object_cellranger.R # output is robject CWOW_cellbender.rds
```

2. Quality filtering such as min and max nCount, nFeature, and percent MT
```
02a_quality_control_post_cellranger.Rmd 
# output is robject CWOW_cellbender_filtered.rds & various QC plots

02b_sex_check_cellbender.Rmd 
# output violin plot of XIST and UTY expression to confirm the sex of the samples
```

  3 - 5. Doublet finder, doublet exploration, and identifying potential doublets to keep. 
```
03_doublet_removal_post_cellbender.R 
# output is two robjects:
# 1) singlets only CWOW_cellbender_singlets.rds  
# 2) doublets only CWOW_cellbender_doublets.rds

# 04 explore just the doublets to determine if they are truly doublets
# doublets will be intergrated via harmony, and annotated via AZIMUTH
04_explore_doublets_cellbender.R 
# output is various plots to assess the doublets and CWOW_cellbender_doublets_harmony_azimuth.rds object

# 05 merges the doublets to keep with the singlets 
# ** MUST update line 23 of the subset for doublet clusters to be removed
05a_doublets_to_keep_and_merge_with_singlets_cellbender.R 
# output is CWOW_cellbender_mereged_singlets_with_kept_doublets.rds


# recaluculate qc metrics for percentage of mt, ribo, hb, choroid, ect. 
05b_explore_doublets_and_singlets_cellbender.R 
# output is various qc plots 
```

## Find markers and annotation
The scripts for find markers and annotation are in the 03_markers_and_annotation folder. First the data will be explored to determine if integration is needed, and then integrated via rpca. Post integration, a pass1 of findmarkers and cluster annotation is manually determined. Then each identified cell type will be re-clustered. This is done to determine that the nuclei within that cell type are indeed of that cell type or potentially contaminated nuclei that should be removed. Once the "clean" nuclei are identified for each cell type, all the cell types are merged back together and the Seurat object is re-processed to determine cluster markers and then pass2 manual annotations are completed. The "dirty" nuclei will also be merged and re-processed to ensure that these are indeed uninformative nuclei that should be removed so as not to contaminate the biological signal. 

1 - 2. Determine if data needs to be integrated and prepare data for integration
```
# cd scripts/03_markers_and_annotation/

# Assess the data without integration 
01_no_integration_assessment_cellbender.R 
# output plots named "CWOW_cellbender_joinlayers_azimuth"
# output rObject CWOW_cellbender_joinlayers_azimuth.rds

# Input for the 02 script is the rObject CWOW_cellbender_joinlayers_azimuth.rds obtained for script 01
# split the object by sample and run SCTransform
# SCTransform replaces NormalizeData(), ScaleData(), and FindVariableFeatures()
# after SCTransform, the script will run PCA, UMAP, FindNeighbors, FindClusters
02_preprocess_for_integration_split_SCTransform_cellbender.R 
# Output is two rObjects:
# 1) SCTransform_only before UMAP and Find nieghbors and clusters
# 2) CWOW_cellbender_SCTransform.rds with Find nieghbors and clusters, and UMAP
```

3. Data will be integrated using the rpca approach. 
Robust Principal Component Analysis (RPCA). Instead of standard PCA within the CCA framework, RPCA is designed to be more robust to outliers and noise in the data. This can be particularly beneficial when dealing with highly heterogeneous datasets or those with significant technical artifacts.
```
03_integration_rpca.R 
# output is CWOW_cellbender_RPCAIntegration.rds object. 
# **Note the umap is named integrated.rpca
```

4. Find markers
```
04_find_markers_pass1_rpca.R 
# Output is two rds objects:
# 1) all makers CWOW_cellbender_RPCAIntegration_markers.rds and 
# 2) CWOW_cellbender_RPCAIntegration_markers_log2FC1_q0.01.rds which is markers with log2FC > 1 & q-value < 0.01
```

5. Pass 1 annotations 
```
05_annotations_pass1_rpca.Rmd 
# output CWOW_cellbender_RPCAIntegration_annotated_before_recluster.rds
```

6. Recluster each cell type
```
06_recluster_pass1_rpca.R 
```

7. merge the reclusters pass1. Split by clean and dirty 
```
07_merge_recluster_clean_pass1_rpca.R
```

8. Recluster the clean nuclei and recluster the dirty nuclei
```
08_recluster_clean.R
08_recluster_dirty.R
```

9. Find markers and annotation - pass 2 clean nuclei only 
```
09a_find_markers_pass2_rpca.R
09b_post_reclustering_annotations_pass2_rpca.Rmd
```

## Differential expression 

1. Assess the percentage of variance explained by covariates and factors
```
01_variance_assessment.Rmd
```

  2 - 4. Differential expression and output tables and plots 
```
02_differential_expression.Rmd
03_make_DEG_excel_and_volcano_plots.Rmd
04_make_heatmaps_and_upset_plots.Rmd
```


## Shiny app 
Generation of shiny app for exploration of the results\
[cwow_snrna_shiny](https://fryerlab.shinyapps.io/lbd_cwow_snrna/)
```
01_make_shiny.R
```

## References
All packages used in this workflow are publicity available. If you use this workflow please cite the packages used. 
If you use the data in this workflow cite the following:
[Olney et al. 2025](https://academic.oup.com/brain/article/148/1/69/7698452)


## Contacts

| Contact | Email |
| --- | --- |
| Kimberly Olney, PhD | kolney@tgen.org |
| John Fryer, PhD | jfryer@tgen.org |