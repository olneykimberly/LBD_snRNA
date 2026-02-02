setwd("/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/")
source("file_paths_and_colours.R")
project_ID <- "CWOW_cellbender"

dataObject <- readRDS(paste0("../rObjects/", project_ID, "_filtered_subclusters_pass2.rds"))

# thresholds already applied 
nCount.min <- 1000 
nCount.max <- 100000 
nFeature.min <- 800
nFeature.max<- 12000
cutoff_complexity <- 0.80
cutoff_mt <- 5
cutoff_hb <- 1
cutoff_ribo <- 1
cutoff_choroid <- 1


#  checks
table(dataObject$Sample_ID)
table(dataObject$Sex)
table(dataObject$group)
dataObject

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
  scale_y_continuous(breaks = seq(0,20000, by = 2000), limits = c(0,20000)) +
  ggtitle("Nuclei per sample") +
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
path <- paste0("../results/density/",project_ID,"_density_final_clean")
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
path <- paste0("../results/violin/",project_ID,"_nFeature_nCount_complexity_final_clean")
saveToPDF(paste0(path, ".pdf"), width = 12, height = 9)
dev.off()

# save v2
v2
path <- paste0("../results/violin/",project_ID,"_percent_final_clean")
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
path <- paste0("../results/scatter/",project_ID,"_nFeature_vs_nCount_percentMt_final_clean")
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
path <- paste0("../results/scatter/",project_ID,"_nCount_vs_percentMT_final_clean")
saveToPDF(paste0(path, ".pdf"), width = 6, height = 4)

# Visualize the number of cell counts per sample
data <- as.data.frame(table(dataObject$Sample_ID))
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
  scale_y_continuous(breaks = seq(0,20000, by = 2000), limits = c(0,20000)) +
  ggtitle("Filtered: nuclei per sample") +
  theme(legend.position =  "none") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))
ncells2

mean(data$frequency)
median(data$frequency)
path <- paste0("../results/nuclei_count/",project_ID, 
               "_cells_per_sample_final_clean")
ncells2
saveToPDF(paste0(path, ".pdf"), width = 6, height = 4)


# 1. Calculate global abundance (one row per cell type)
relative_abundance <- dataObject@meta.data %>%
  group_by(cell_type) %>% # Use your cluster/cell type column
  dplyr::count() %>%
  ungroup() %>%
  dplyr::mutate(percent = 100 * n / sum(n)) %>%
  arrange(desc(percent))

# 2. Plot with a static X-axis label
rel_abun <- ggplot(relative_abundance, aes(x = "Total Dataset", y = percent, fill = cell_type)) +
  geom_col(width = 0.5) + # Adjusted width for a single bar
  geom_text(aes(label = paste0(round(percent, 1), "%")), 
            position = position_stack(vjust = 0.5), 
            size = 3.5, 
            color = "white",
            fontface = "bold") +
  scale_fill_manual(values = color_panel) +
  labs(x = NULL, y = "Percentage (%)", title = "Global Cell Type Composition") +
  theme_bw() +
  theme(panel.grid.major.x = element_blank()) # Clean up the X-axis

rel_abun
pdf(
  paste0(
    "../results/UMAP/",
    project_ID,
    "_annotated_relative_abundance_pass2_final_clean.pdf"
  ),
  width = 4,
  height = 5
)
rel_abun
dev.off()



relative_abundance <- dataObject@meta.data %>%
  group_by(cell_type, Sample_ID) %>%
  dplyr::count() %>%
  group_by(Sample_ID) %>%
  dplyr::mutate(percent = 100 * n / sum(n)) %>%
  ungroup()


rel_abun <- ggplot(relative_abundance, aes(x = Sample_ID, y = percent, fill = cell_type)) +
  geom_col() +
  geom_text(aes(label = paste0(round(percent), "%")), 
            position = position_stack(vjust = 0.5), size = 3, color = "white") +
  scale_fill_manual(values = color_panel) +
  ggtitle("Percentage of cell type") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
rel_abun

pdf(
  paste0(
    "../results/UMAP/",
    project_ID,
    "_annotated_relative_abundance_Sample_ID_pass2_clean_final.pdf"
  ),
  width = 13,
  height = 6
)
rel_abun
dev.off()

# Group
relative_abundance <- dataObject@meta.data %>%
  group_by(cell_type, group) %>%
  dplyr::count() %>%
  group_by(group) %>%
  dplyr::mutate(percent = 100 * n / sum(n)) %>%
  ungroup()


rel_abun <- ggplot(relative_abundance, aes(x = group, y = percent, fill = cell_type)) +
  geom_col() +
  geom_text(aes(label = paste0(round(percent), "%")), 
            position = position_stack(vjust = 0.5), size = 3, color = "white") +
  scale_fill_manual(values = color_panel) +
  ggtitle("Percentage of cell type") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
rel_abun

pdf(
  paste0(
    "../results/UMAP/",
    project_ID,
    "_annotated_relative_abundance_group_pass2_clean_final.pdf"
  ),
  width = 8,
  height = 6
)
rel_abun
dev.off()

dot_ind <- DotPlot(dataObject,
                   features = genes_markers, 
                   cluster.idents = TRUE,
                   dot.scale = 8) + RotatedAxis()
dot_ind


pdf(
  paste0(
    "../results/DotPlot/",
    project_ID,
    "_DotPlot_annotation_pass2_final_clean.pdf"
  ),
  width = 14,
  height = 5
)
dot_ind
dev.off()
