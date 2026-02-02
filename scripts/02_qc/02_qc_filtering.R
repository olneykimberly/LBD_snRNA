knitr::opts_knit$set(root.dir = "/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/")
knitr::opts_chunk$set(message = FALSE, warning = FALSE)
setwd("/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/")

source(here::here("/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/", "file_paths_and_colours.R"))
project_ID <- "CWOW_cellbender"

library(future)
plan("multicore", workers = 8)
options(future.globals.maxSize = 250 * 1024^3) # 250GB 
options(future.rng.onMisuse = "ignore") # trigger Garbage Collection

# genes
ann_path <- "../rObjects/annotation.rds"

if (file.exists(ann_path)) {
  genes <- readRDS(ann_path)
} else {
  gtf.file <- file.path(path_ref, "genes/genes.gtf")
  genes <- rtracklayer::import(gtf.file)
  genes <- as.data.frame(genes)
  genes <- genes[genes$type == "gene", ]
  saveRDS(genes, ann_path)
}

gene_type_table <- table(genes$gene_type)
write.table(gene_type_table, "gene_type_table.tsv", row.names = FALSE, quote = FALSE, sep = "\t")
write.table(genes, "genes.tsv", row.names = FALSE, quote = FALSE, sep = "\t")

mt.genes.df <- subset(genes, seqnames == "chrM")
genes_mt <- mt.genes.df$gene_name

rm(mt.genes.df, gene_type_table)

nCount.min <- 1000 
nCount.max <- 100000 
nFeature.min <- 800
nFeature.max<- 12000
cutoff_complexity <- 0.80
cutoff_mt <- 5
cutoff_hb <- 1
cutoff_ribo <- 1
cutoff_choroid <- 1

dataObject <- readRDS(paste0("../rObjects/", project_ID, ".rds"))
DefaultAssay(dataObject) <- "RNA"

# Seurat barcodes
barcodes <- colnames(dataObject)

# Extract BR_Nuclei_#### from cell names like "BR_Nuclei_0373_AAAC..."
m <- stringr::str_match(barcodes, "^(BR_Nuclei_)([0-9]+)_")
if (any(is.na(m[,1]))) {
  stop("Some cell names did not match 'BR_Nuclei_####_'. Example: ",
       barcodes[which(is.na(m[,1]))[1]])
}

# Standardize Seurat key to BR_Nuclei_#### (4-digit)
sample_br <- paste0(m[,2], sprintf("%04d", as.integer(m[,3])))
dataObject$Sample_BR <- sample_br  # keep original namespace 
metadata$sample_br <- unique(dataObject$Sample_BR) # add to metadata
# make meta key 
meta_key_raw <- as.character(metadata$sample_br)

# extract prefix and number; then pad to 4 digits
mm <- stringr::str_match(meta_key_raw, "^(BR_Nuclei_)([0-9]+)$")
if (any(is.na(mm[,1]))) {
  bad <- unique(meta_key_raw[is.na(mm[,1])])
  stop("metadata$sample_br contains values not matching 'BR_Nuclei_####':\n",
       paste(head(bad, 20), collapse = ", "),
       if (length(bad) > 20) "\n... (truncated)" else "")
}
meta_key <- paste0(mm[,2], sprintf("%04d", as.integer(mm[,3])))

# build mapping: names = BR_Nuclei_####, values = short Sample_ID ----
br_to_short <- as.character(metadata$Sample_ID)
names(br_to_short) <- meta_key
if (anyDuplicated(names(br_to_short))) {
  dup <- unique(names(br_to_short)[duplicated(names(br_to_short))])
  stop("Duplicate BR keys in metadata after standardization:\n",
       paste(dup, collapse = ", "))
}

sample_short <- unname(br_to_short[sample_br])
if (any(is.na(sample_short))) {
  missing_br <- sort(unique(sample_br[is.na(sample_short)]))
  # helpful debug: show what metadata has vs what Seurat has
  stop("These BR_Nuclei samples are in Seurat but missing from metadata mapping:\n",
       paste(missing_br, collapse = ", "),
       "\n\nQuick debug:\n",
       "Example Seurat keys: ", paste(head(sort(unique(sample_br)), 5), collapse = ", "), "\n",
       "Example metadata keys:", paste(head(sort(unique(names(br_to_short))), 5), collapse = ", "))
}
dataObject$Sample_ID <- sample_short
dataObject$Sample_ID <- factor(dataObject$Sample_ID, levels = order_samples_short)
Idents(dataObject) <- dataObject$Sample_ID

# Sanity checks
print(table(dataObject$Sample_BR))
print(table(dataObject$Sample_ID))

rm(barcodes, m, mm, meta_key_raw, meta_key, sample_br, sample_short, br_to_short)

metadata2 <- metadata
rownames(metadata2) <- metadata2$Sample_ID  # key by short IDs

# confirm  overlap
missing_in_meta <- setdiff(levels(dataObject$Sample_ID), rownames(metadata2))
if (length(missing_in_meta) > 0) {
stop("These Sample_IDs are in Seurat but missing from metadata: ", paste(missing_in_meta, collapse = ", "))
}

meta_mapped <- metadata2[as.character(dataObject$Sample_ID), , drop = FALSE]
dataObject <- AddMetaData(dataObject, metadata = meta_mapped)
dataObject$group <- factor(dataObject$TYPE, levels = c("CONTROL", "AD_AT", "LBD_S", "LBD_AS", "LBD_ATS"))

#  checks
table(dataObject$Sample_ID)
table(dataObject$Sex)
table(dataObject$group)

count_layers <- grep("^counts", Layers(dataObject), value = TRUE)

sample_ids <- unique(dataObject$Sample_ID)          # or unique(dataObject@meta.data$Sample_ID)
stopifnot(length(sample_ids) == length(count_layers))

layer_map <- data.frame(
  Sample_ID = sample_ids,
  layer     = count_layers,
  stringsAsFactors = FALSE
)

layer_map

summary(dataObject$nCount_RNA)
summary(dataObject$nFeature_RNA)

# safer complexity (avoid Inf / -Inf)
dataObject$cell.complexity <- log10(dataObject$nFeature_RNA + 1) / log10(dataObject$nCount_RNA + 1)

gene.names <- rownames(dataObject)

# ribosomal proteins
ribo.genes <- gene.names[grep("^RPS|^RPL", gene.names)]
mt.ribo    <- gene.names[grep("^MRP[SL]", gene.names)]
ribo.combined <- unique(c(mt.ribo, ribo.genes))

# hemoglobin proteins
hb.genes <- gene.names[grep("^HB[BA]", gene.names)]

# Ensure your features are present in the object
genes_mt_valid <- intersect(genes_mt, rownames(dataObject))
ribo_valid <- intersect(ribo.combined, rownames(dataObject))
hb_valid <- intersect(hb.genes, rownames(dataObject))
choroid_genes <- intersect(c("TTR","FOLR1","PRLR"), rownames(dataObject))

# Run vectorized QC (Seurat v5 will handle the layers internally)
dataObject$percent.mt <- PercentageFeatureSet(dataObject, features = genes_mt_valid)
dataObject$percent.ribo <- PercentageFeatureSet(dataObject, features = ribo_valid)
dataObject$percent.hb <- PercentageFeatureSet(dataObject, features = hb_valid)
dataObject$percent.choroid <- PercentageFeatureSet(dataObject, features = choroid_genes)


summary(dataObject$percent.mt)
summary(dataObject$percent.ribo)
summary(dataObject$percent.hb)
summary(dataObject$percent.choroid)

# Visualize the number of cell counts per sample
data <- as.data.frame(table(dataObject$Sample_ID))
colnames(data) <- c("Sample_ID","frequency")

ncells1 <- ggplot(data, aes(x = Sample_ID, y = frequency, fill = Sample_ID)) + 
  geom_col() +
  theme_classic() +
  geom_text(aes(label = frequency), 
            position=position_dodge(width=0.9), 
            #vjust=-0.25, 
            hjust = -.025,
            angle = 90) +
  #scale_fill_manual(values = sample_colors) + 
  scale_y_continuous(breaks = seq(0,30000, by = 2000), limits = c(0,30000)) +
  ggtitle("Raw: nuclei per sample") +
  theme(legend.position =  "none") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))
ncells1

rm(data)

# set graphical parameter
par(mfrow = c(7,1))
# Visualize nCount_RNA
den1 <- ggplot(dataObject@meta.data,
       aes(color = Sample_ID,
           x = nCount_RNA,
           fill = Sample_ID)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  scale_x_log10() +
  #scale_color_manual(values = sample_colors) +
  #scale_fill_manual(values = sample_colors) +
  xlab("nCount_RNA") +
  ylab("Density") +
  geom_vline(xintercept = nCount.min) +
  geom_vline(xintercept = nCount.max) + theme(legend.position="none") +
  annotate("text", x = nCount.min, y = 0.1, label = paste("Min:", nCount.min), 
           angle = 0, vjust = -0.5, hjust = -0.1, size = 3) +
  annotate("text", x = nCount.max, y = 0.1, label = paste("Max:", nCount.max), 
           angle = 0, vjust = -0.5, hjust = -0.1, size = 3)

# Visualize nFeature
den2 <- ggplot(dataObject@meta.data,
       aes(color = Sample_ID,
           x = nFeature_RNA,
           fill = Sample_ID)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  scale_x_log10() +
  #scale_color_manual(values = sample_colors) +
  #scale_fill_manual(values = sample_colors) +
  xlab("nFeature_RNA") +
  ylab("Density") +
  geom_vline(xintercept = nFeature.min) +
  geom_vline(xintercept = nFeature.max) + theme(legend.position="none") +
  annotate("text", x = nFeature.min, y = 0.1, label = paste("Min:", nFeature.min), 
           angle = 0, vjust = -0.5, hjust = -0.1, size = 3) +
  annotate("text", x = nFeature.max, y = 0.1, label = paste("Max:", nFeature.max), 
           angle = 0, vjust = -0.5, hjust = -0.1, size = 3)


# Visualize cell complexity
# Quality cells are usually above 0.85
den3 <- ggplot(dataObject@meta.data,
       aes(color = Sample_ID,
           x = cell.complexity,
           fill = Sample_ID)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  #scale_color_manual(values = sample_colors) +
  #scale_fill_manual(values = sample_colors) +
  xlab("Cell Complexity (log10(nFeature/nCount))") +
  ylab("Density") +
  geom_vline(xintercept = cutoff_complexity) + theme(legend.position="none")  +
  annotate("text", x = cutoff_complexity, y = 0.1, label = paste("Min:", cutoff_complexity), 
           angle = 0, vjust = -0.5, hjust = -0.1, size = 3)

# Visualize percent.mt
den4 <- ggplot(dataObject@meta.data,
       aes(color = Sample_ID,
           x = percent.mt,
           fill = Sample_ID)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  scale_x_continuous(n.breaks = 4) +
  geom_vline(xintercept = cutoff_mt) +
  #scale_color_manual(values = sample_colors) +
  #scale_fill_manual(values = sample_colors) +
  xlab("% Mitochondrial Genes") +
  ylab("Density")+ theme(legend.position="none")  +
  annotate("text", x = cutoff_mt, y = 0.1, label = paste("Max:", cutoff_mt), 
           angle = 0, vjust = -0.5, hjust = -0.1, size = 3)


# Visualize percent.ribo
den5 <- ggplot(dataObject@meta.data,
       aes(color = Sample_ID,
           x = percent.ribo,
           fill = Sample_ID)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  scale_x_continuous(n.breaks = 4) +
  geom_vline(xintercept = cutoff_ribo) +
  #scale_color_manual(values = sample_colors) +
  #scale_fill_manual(values = sample_colors) +
  xlab("% Ribosomal Genes") +
  ylab("Density")+ theme(legend.position="none")  +
  annotate("text", x = cutoff_ribo, y = 0.1, label = paste("Max:", cutoff_ribo), 
           angle = 0, vjust = -0.5, hjust = -0.1, size = 3)


# Visualize percent.hb
den6 <- ggplot(dataObject@meta.data,
       aes(color = Sample_ID,
           x = percent.hb,
           fill = Sample_ID)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  scale_x_continuous(n.breaks = 4) +
  geom_vline(xintercept = cutoff_hb) +
  #scale_color_manual(values = sample_colors) +
  #scale_fill_manual(values = sample_colors) +
  xlab("% Hemoglobin Genes") +
  ylab("Density")+ theme(legend.position="none") +
  annotate("text", x = cutoff_hb, y = 0.1, label = paste("Max:", cutoff_hb), 
           angle = 0, vjust = -0.5, hjust = -0.1, size = 3)


den7 <- ggplot(dataObject@meta.data,
       aes(color = Sample_ID,
           x = percent.choroid,
           fill = Sample_ID)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  scale_x_continuous(n.breaks = 4) +
  geom_vline(xintercept = cutoff_choroid) +
  #scale_color_manual(values = sample_colors) +
  #scale_fill_manual(values = sample_colors) +
  xlab("% Choroid Marker Genes") +
  ylab("Density")+ theme(legend.position="none") +
  annotate("text", x = cutoff_choroid, y = 0.1, label = paste("Max:", cutoff_choroid), 
           angle = 0, vjust = -0.5, hjust = -0.1, size = 3)

# Arrange graphs in grid
plots1 <- list(den1,den2,den3,den4,den5,den6,den7)
layout1 <- rbind(c(1),c(2),c(3),c(4),c(5),c(6),c(7))
grid1 <- grid.arrange(grobs = plots1, layout_matrix = layout1)

grid1 <- grid.arrange(grobs = plots1, layout_matrix = layout1)
path <- paste0("../results/density/",project_ID,"_density_raw")
saveToPDF(paste0(path, ".pdf"), width = 6, height = 10)
dev.off()

rm(grid1, plots1)

# nFeature, nCount, and cell.complexity violins
v1 <- VlnPlot(dataObject,
              features = c( "nCount_RNA","nFeature_RNA", "cell.complexity"),
              ncol = 1,
              group.by = 'Sample_ID',
              ##cols = sample_colors,
              pt.size = 0)
v1

#  percent violins
v2 <- VlnPlot(dataObject,
              features = c("percent.mt","percent.ribo","percent.hb", "percent.choroid"),
              ncol = 2,
              group.by = 'Sample_ID',
              #cols = sample_colors,
              pt.size = 0)
v2

# save v1
v1
path <- paste0("../results/violin/",project_ID,"_nFeature_nCount_complexity_raw")
saveToPDF(paste0(path, ".pdf"), width = 12, height = 9)
dev.off()

# save v2
v2
path <- paste0("../results/violin/",project_ID,"_percent_raw")
saveToPDF(paste0(path, ".pdf"), width = 16, height = 7)
dev.off()

# cleanup
#remove(v1,v2)

s1 <- ggplot(
  dataObject@meta.data,
  aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt)) + 
  geom_point() + 
  stat_smooth(method=lm) +
	scale_x_log10() +   	
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = nCount.min) + 
  geom_vline(xintercept = nCount.max) + 
  geom_hline(yintercept = nFeature.min) + 
  geom_hline(yintercept = nFeature.max) + 
  facet_wrap(~Sample_ID, ncol = 7) +
  scale_colour_gradient(low = "gray90", high = "black", limits =c(0,100))
s1

# save
s1
path <- paste0("../results/scatter/",project_ID,"_nFeature_vs_nCount_percentMt_raw")
saveToPDF(paste0(path, ".pdf"), width = 20, height = 15)
dev.off()

# cleanup
remove(s1)

s2 <- FeatureScatter(dataObject,
               feature1 = "nCount_RNA",
               feature2 = "percent.mt",
               group.by = 'Sample_ID',
               #cols = sample_colors,
               shuffle = TRUE)
s2

# save
s2
path <- paste0("../results/scatter/",project_ID,"_nCount_vs_percentMT_raw")
saveToPDF(paste0(path, ".pdf"), width = 6, height = 4)

# cleanup
remove(s2)

dataObject.filtered <- subset(dataObject, 
                              subset = nCount_RNA > nCount.min & 
                                nCount_RNA < nCount.max &
                                nFeature_RNA > nFeature.min & 
                                nFeature_RNA < nFeature.max &
                                cell.complexity > cutoff_complexity &
                                percent.mt <= cutoff_mt &
                                percent.hb <= cutoff_hb &
                                percent.ribo <= cutoff_ribo &
                                percent.choroid <= cutoff_choroid)
# print nuclei removed
print(paste0(dim(dataObject)[2] - dim(dataObject.filtered)[2]," nuclei removed"))
# 
print("before filter: ")
table(dataObject$Sample_ID)
print("after filter: ")
table(dataObject.filtered$Sample_ID)

saveRDS(dataObject.filtered, paste0("../../rObjects/",project_ID,"_filtered.rds"), compress = FALSE)

# Visualize the number of cell counts per sample
data <- as.data.frame(table(dataObject.filtered$Sample_ID))
colnames(data) <- c("Sample_ID","frequency")

ncells2 <- ggplot(data, aes(x = Sample_ID, y = frequency, fill = Sample_ID)) + 
  geom_col() +
  theme_classic() +
  geom_text(aes(label = frequency), 
            position=position_dodge(width=0.9), 
            #vjust=-0.25, 
            hjust = -.025,
            angle = 90) +
  #scale_fill_manual(values = sample_colors) + 
  scale_y_continuous(breaks = seq(0,30000, by = 2000), limits = c(0,30000)) +
  ggtitle("Filtered: nuclei per sample") +
  theme(legend.position =  "none") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))
ncells2
# Arrange graphs in grid
plots2 <- list(ncells1,ncells2)
layout2 <- rbind(c(1),c(2))
grid2 <- grid.arrange(grobs = plots2, layout_matrix = layout2)

mean(data$frequency)
median(data$frequency)

# save
grid2 <- grid.arrange(grobs = plots2, layout_matrix = layout2)
path <- paste0("../results/nuclei_count/",project_ID, 
               "_cells_per_sample")
saveToPDF(paste0(path, ".pdf"), width = 12, height = 8)


# cleanup
remove(ncells1,ncells2,plots2,layout2,grid2)

# set graphical parameter
par(mfrow = c(7,2))

den8 <- ggplot(dataObject.filtered@meta.data,
       aes(color = Sample_ID,
           x = nCount_RNA,
           fill = Sample_ID)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  scale_x_log10() +
  #scale_color_manual(values = sample_colors) +
  #scale_fill_manual(values = sample_colors) +
  xlab("nCount_RNA") +
  ylab("Density") +
  geom_vline(xintercept = nCount.min) +
  geom_vline(xintercept = nCount.max) + theme(legend.position="none") +
  annotate("text", x = nCount.min, y = 0.1, label = paste("Min:", nCount.min), 
           angle = 0, vjust = -0.5, hjust = -0.1, size = 3) +
  annotate("text", x = nCount.max, y = 0.1, label = paste("Max:", nCount.max), 
           angle = 0, vjust = -0.5, hjust = -0.1, size = 3)

# Visualize nFeature
den9 <- ggplot(dataObject.filtered@meta.data,
       aes(color = Sample_ID,
           x = nFeature_RNA,
           fill = Sample_ID)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  scale_x_log10() +
  #scale_color_manual(values = sample_colors) +
  #scale_fill_manual(values = sample_colors) +
  xlab("nFeature_RNA") +
  ylab("Density") +
  geom_vline(xintercept = nFeature.min) +
  geom_vline(xintercept = nFeature.max) + theme(legend.position="none") +
  annotate("text", x = nFeature.min, y = 0.1, label = paste("Min:", nFeature.min), 
           angle = 0, vjust = -0.5, hjust = -0.1, size = 3) +
  annotate("text", x = nFeature.max, y = 0.1, label = paste("Max:", nFeature.max), 
           angle = 0, vjust = -0.5, hjust = -0.1, size = 3)


# Visualize cell complexitydataObject.filtered.split.doublets
# Quality cells are usually above 0.85
den10 <- ggplot(dataObject.filtered@meta.data,
       aes(color = Sample_ID,
           x = cell.complexity,
           fill = Sample_ID)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  #scale_color_manual(values = sample_colors) +
  #scale_fill_manual(values = sample_colors) +
  xlab("Cell Complexity (log10(nFeature/nCount))") +
  ylab("Density") +
  geom_vline(xintercept = cutoff_complexity) + theme(legend.position="none")  +
  annotate("text", x = cutoff_complexity, y = 0.1, label = paste("Min:", cutoff_complexity), 
           angle = 0, vjust = -0.5, hjust = -0.1, size = 3)

# Visualize percent.mt
den11 <- ggplot(dataObject.filtered@meta.data,
       aes(color = Sample_ID,
           x = percent.mt,
           fill = Sample_ID)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  scale_x_continuous(n.breaks = 4) +
  geom_vline(xintercept = cutoff_mt) +
  #scale_color_manual(values = sample_colors) +
  #scale_fill_manual(values = sample_colors) +
  xlab("% Mitochondrial Genes") +
  ylab("Density")+ theme(legend.position="none")  +
  annotate("text", x = cutoff_mt, y = 0.1, label = paste("Max:", cutoff_mt), 
           angle = 0, vjust = -0.5, hjust = -0.1, size = 3)


# Visualize percent.ribo
den12 <- ggplot(dataObject.filtered@meta.data,
       aes(color = Sample_ID,
           x = percent.ribo,
           fill = Sample_ID)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  scale_x_continuous(n.breaks = 4) +
  geom_vline(xintercept = cutoff_ribo) +
  #scale_color_manual(values = sample_colors) +
  #scale_fill_manual(values = sample_colors) +
  xlab("% Ribosomal Genes") +
  ylab("Density")+ theme(legend.position="none")  +
  annotate("text", x = cutoff_ribo, y = 0.1, label = paste("Max:", cutoff_ribo), 
           angle = 0, vjust = -0.5, hjust = -0.1, size = 3)


# Visualize percent.hb
den13 <- ggplot(dataObject.filtered@meta.data,
       aes(color = Sample_ID,
           x = percent.hb,
           fill = Sample_ID)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  scale_x_continuous(n.breaks = 4) +
  geom_vline(xintercept = cutoff_hb) +
  #scale_color_manual(values = sample_colors) +
  #scale_fill_manual(values = sample_colors) +
  xlab("% Hemoglobin Genes") +
  ylab("Density")+ theme(legend.position="none") +
  annotate("text", x = cutoff_hb, y = 0.1, label = paste("Max:", cutoff_hb), 
           angle = 0, vjust = -0.5, hjust = -0.1, size = 3)


den14 <- ggplot(dataObject.filtered@meta.data,
       aes(color = Sample_ID,
           x = percent.choroid,
           fill = Sample_ID)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  scale_x_continuous(n.breaks = 4) +
  geom_vline(xintercept = cutoff_choroid) +
  #scale_color_manual(values = sample_colors) +
  #scale_fill_manual(values = sample_colors) +
  xlab("% Choroid Marker Genes") +
  ylab("Density")+ theme(legend.position="none") +
  annotate("text", x = cutoff_choroid, y = 0.1, label = paste("Max:", cutoff_choroid), 
           angle = 0, vjust = -0.5, hjust = -0.1, size = 3)

# Arrange graphs in grid
plots3 <- list(den1,den2,den3,den4,den5,den6,den7,den8,den9,den10,den11,den12,den13,den14)
layout3 <- rbind(c(1,8),c(2,9),c(3,10), c(4,11),c(5,12),c(6,13),c(7,14))
grid3 <- grid.arrange(grobs = plots3, layout_matrix = layout3)

# save
grid3 <- grid.arrange(grobs = plots3, layout_matrix = layout3)
path <- paste0("../results/density/",project_ID, 
               "_density_filtered")
saveToPDF(paste0(path, ".pdf"), width = 10, height = 10)

# cleanup
remove(den1,den2,den3,den4,den5,den6,den7,den8,den9,den10,den11,den12,den13,den14,plots3,layout3,grid3)

# nFeature, nCount, and cell.complexity violins
v3 <- VlnPlot(dataObject.filtered,
              features = c("nFeature_RNA", "nCount_RNA","cell.complexity"),
              ncol = 1,
              group.by = 'Sample_ID',
              #cols = sample_colors,
              pt.size = 0)
v3

#  percent violins
v4 <- VlnPlot(dataObject.filtered,
              features = c("percent.mt","percent.ribo","percent.hb", "percent.choroid"),
              ncol = 2,
              group.by = 'Sample_ID',
              #cols = sample_colors,
              pt.size = 0)
v4

# save
v3
path <- paste0("../results/violin/",project_ID, 
               "_nFeature_nCount_complexity_filtered")
saveToPDF(paste0(path, ".pdf"), width = 12, height = 11)

v4
path <- paste0("../results/violin/",project_ID, 
               "_percent_filtered")
saveToPDF(paste0(path, ".pdf"), width = 14, height = 9)

# cleanup
remove(v3,v4)

s3 <- ggplot(
  dataObject.filtered@meta.data,
  aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt)) + 
  geom_point() + 
  stat_smooth(method=lm) +
	scale_x_log10() +   	
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = nCount.min) + 
  geom_vline(xintercept = nCount.max) + 
  geom_hline(yintercept = nFeature.min) + 
  geom_hline(yintercept = nFeature.max) + 
  facet_wrap(~Sample_ID, ncol = 8) +
  scale_colour_gradient(low = "gray90", high = "black", limits =c(0,100))
s3

# save
s3
path <- paste0("../results/scatter/",project_ID,
               "_nFeature_vs_nCount_perecentMt_filtered")
saveToPDF(paste0(path, ".pdf"), width = 20, height = 15)

# cleanup
remove(s3)

s4 <- FeatureScatter(dataObject.filtered,
               feature1 = "nCount_RNA",
               feature2 = "percent.mt",
               group.by = 'Sample_ID',
               #cols = sample_colors,
               shuffle = TRUE)
s4

# save
s4
path <- paste0("../results/scatter/",project_ID,
               "_nCount_vs_percentMT_filtered")
saveToPDF(paste0(path, ".pdf"), width = 6, height = 4)

# cleanup
remove(s4)

# Visualize the distribution of genes detected per cell via boxplot
b1 <- ggplot(dataObject.filtered@meta.data,
       aes(x = Sample_ID, 
           y = log10(nFeature_RNA), 
           fill=Sample_ID)) + 
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        legend.position = "none") +
  theme(plot.title = element_text(hjust = 0.5, face="bold")) +
  ggtitle("Unique Genes / Nuclei / Sample") +
  #scale_color_manual(values = sample_colors) +
  #scale_fill_manual(values = sample_colors) +
  xlab("Sample")
b1

# save
b1
path <- paste0("../results/nuclei_count/",project_ID,
               "_nFeature_per_sample")
saveToPDF(paste0(path, ".pdf"), width = 12, height = 4)

# cleanup
remove(b1)

