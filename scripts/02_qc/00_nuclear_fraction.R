## ----working_directory-------------------------------------------------------------------------------------------------------------------
knitr::opts_knit$set(root.dir = "/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/")
setwd("/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/")

library(DropletQC)
library(Rsamtools)
library(GenomicRanges)

## ----echo=FALSE, message=FALSE-----------------------------------------------------------------------------------------------------------
source(here::here("/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/", "file_paths_and_colours.R"))

## ----gene_info---------------------------------------------------------------------------------------------------------------------------
gtf.file <- paste0(path_ref, "/genes/genes.gtf")

## ----nuclear_fraction---------------------------------------------------------------------------------------------------------------------------
prefix <- "../cellranger/"
suffix <- "/outs/"
# Create an empty list to store Seurat objects
for (i in order_samples) {
    print(i)
    sample_name <- basename(i)
    # get nuclear fraction
    nf <- nuclear_fraction_annotation(
      annotation_path = gtf.file,
      tiles = 100, cores = 8, 
      bam = paste0(prefix,i,suffix, "possorted_genome_bam.bam"),
      barcodes = paste0(prefix,i,suffix, "filtered_feature_bc_matrix/barcodes.tsv.gz"), # or use cellbender's _cell_barcodes.csv
      verbose = FALSE
    )
    # reformat
    nf$barcode <- rownames(nf)
    # write output
    write.table(x = nf,
                file = paste0("../nuclear_fraction/", i, "_nuclear_fraction.tsv"),
                quote = FALSE,
                row.names = FALSE)
    rm(nf)
  }
