## ----setup-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
knitr::opts_knit$set(root.dir = "/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/")
setwd("/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/")

## ----source------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
source(here::here("/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/", "file_paths_and_colours.R"))
projectID <- "CWOW_cellbender_mereged_singlets_with_kept_doublets"
color.panel <- dittoColors()

## ----read_object--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
dataObject <- readRDS(paste0("../rObjects/",projectID,".rds"))
dataObject # inspect

## ----metadata-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
metadata <- subset(metadata, Sample_ID != "LBD_AS_F4")
metadata$sampleID <- factor(metadata$Sample_ID, levels = c(metadata$Sample_ID))
samples <- metadata$sampleID 
order_sex <- factor(metadata$sex_inferred, levels = unique(metadata$sex_inferred))
order_disease <- factor(metadata$TYPE, levels = c("CONTROL", "AD_AT", "LBD_S", "LBD_AS", "LBD_ATS"))

metadata <- metadata %>%
  mutate(sampleID = gsub(".*_(\\d+)_.*_(BR_Nuclei).*", "\\2_\\1", Lane.Name))
samples <- metadata$sampleID 

# sampleID with order_disease
order <- metadata %>%
  arrange(order_disease) %>%
  dplyr::select(TYPE, sampleID, Sample_ID)
samples <- order$sampleID
order_disease <- order$TYPE
sample_order <- factor(order$Sample_ID, levels = order$Sample_ID)

seurat_sample_order <- as.character(dataObject$sample)
matched_metadata <- metadata[match(seurat_sample_order, as.character(metadata$sampleID)), ]
dataObject$Sample_ID <- factor(matched_metadata$Sample_ID, levels = unique(matched_metadata$Sample_ID))

genes <- readRDS("../rObjects/annotation.rds")
mt.genes.df <- subset(genes, seqnames == "chrM")
mt.genes <- mt.genes.df$gene_name
## ----QC_metrics-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
summary(dataObject$nCount_RNA)
summary(dataObject$nFeature_RNA)
# cell.complexity
dataObject$cell.complexity <- log10(dataObject$nFeature_RNA) / log10(dataObject$nCount_RNA)

# Chromosome M
gene.names <- rownames(dataObject)
dataObject$percent.mt <- PercentageFeatureSet(dataObject, features = mt.genes)
summary(dataObject$percent.mt)

"^RPS|^RPL"
# ribosomal proteins 
ribo.genes <- gene.names[grep("^RPS|^RPL", gene.names)] 
mt.ribo <- gene.names[grep("^MRP[SL]", gene.names)]
ribo.combined <- c(mt.ribo,ribo.genes)
dataObject$percent.ribo <- PercentageFeatureSet(dataObject, features = ribo.combined)
summary(dataObject$percent.ribo)

# hemoglobin proteins
hb.genes <- gene.names[grep("^HB[BA]", gene.names)]
dataObject$percent.hb <- PercentageFeatureSet(dataObject, features = hb.genes)
summary(dataObject$percent.hb)

# percent choroid plexus
dataObject$percent.choroid <- PercentageFeatureSet(dataObject, features = c("TTR","FOLR1", "PRLR"))
summary(dataObject$percent.choroid)

# percent MALAT1 
dataObject$percent.MALAT1 <- PercentageFeatureSet(dataObject, features = c("MALAT1"))
summary(dataObject$percent.MALAT1)


## ----cells_per_sample-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
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
  ggtitle("Nuclei per sample") +
  theme(legend.position =  "none") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))
ncells1
path <- paste0("../results/nuclei_count/",projectID, 
               "_cells_per_sample")
ncells1
saveToPDF(paste0(path, ".pdf"), width = 12, height = 4.5)
mean_counts <- mean(data$frequency)
median(data$frequency)
sd_counts <- sd(data$frequency)

upper_threshold <- mean_counts + 2 * sd_counts
lower_threshold <- mean_counts - 2 * sd_counts

data$is_outlier <- ifelse(data$frequency > upper_threshold | data$frequency < lower_threshold, TRUE, FALSE)

# To view the outlier samples:
outlier_samples <- data[data$is_outlier == TRUE, ]
print(outlier_samples)

## ----Density-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
df <- dataObject@meta.data
density_plot <- ggplot(df, aes(x = nCount_RNA, color = Sample_ID, fill = Sample_ID)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  scale_x_log10() +
  xlab("nCount_RNA") +
  ylab("Density") +
  ggtitle("Density nCount") +
  theme(axis.title.x = element_blank(), # Remove x-axis title for cleaner alignment
        axis.text.x = element_blank(),   # Remove x-axis text as boxplot will provide labels
        axis.ticks.x = element_blank())  # Remove x-axis ticks

# Calculate summary statistics per Sample_ID
summary_stats <- df %>%
  group_by(Sample_ID) %>%
  summarise(
    ymin = min(nCount_RNA),
    ymax = max(nCount_RNA),
    Q1 = quantile(nCount_RNA, 0.25),
    Median = quantile(nCount_RNA, 0.5),
    Q3 = quantile(nCount_RNA, 0.75),
    Mean = mean(nCount_RNA)
  )

# Convert to log10 if your x-axis in the plot is log10
summary_stats_log10 <- summary_stats %>%
  mutate(across(c(ymin, ymax, Q1, Median, Q3, Mean), ~log10(.x)))
# For multiple boxplots:
box_plot_multiple <- ggplot(df, aes(x = nCount_RNA, y = Sample_ID, color = Sample_ID, fill = Sample_ID)) +
  geom_boxplot(width = 0.5, alpha = 0.2, outlier.shape = NA) +
  theme_classic() +
  scale_x_log10() +
  ggtitle("Boxplot nCount") +
  xlab("") + # No x-label here
  ylab("") + # No y-label here, as it will be provided by the density plot
  theme(legend.position = "none",
        axis.text.y = element_blank(), # Hide y-axis text to clean up
        axis.ticks.y = element_blank(), # Hide y-axis ticks
        axis.title.y = element_blank(), # Hide y-axis title
        plot.margin = margin(b = 0)) # Reduce bottom margin
# For individual boxplots for each Sample_ID, arranged on top:
combined_plot <- box_plot_multiple / density_plot + plot_layout(heights = c(2, 3)) # Adjust heights as needed
path <- paste0("../results/density/",projectID,"_nCount")
print(combined_plot)
saveToPDF(paste0(path, ".pdf"), width = 8, height = 8)

density_plot <- ggplot(df, aes(x = nFeature_RNA, color = Sample_ID, fill = Sample_ID)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  scale_x_log10() +
  xlab("nFeature_RNA") +
  ylab("Density") +
  ggtitle("Density nFeature") +
  theme(axis.title.x = element_blank(), # Remove x-axis title for cleaner alignment
        axis.text.x = element_blank(),   # Remove x-axis text as boxplot will provide labels
        axis.ticks.x = element_blank())  # Remove x-axis ticks

# Calculate summary statistics per Sample_ID
summary_stats <- df %>%
  group_by(Sample_ID) %>%
  summarise(
    ymin = min(nFeature_RNA),
    ymax = max(nFeature_RNA),
    Q1 = quantile(nFeature_RNA, 0.25),
    Median = quantile(nFeature_RNA, 0.5),
    Q3 = quantile(nFeature_RNA, 0.75),
    Mean = mean(nFeature_RNA)
  )

# Convert to log10 if your x-axis in the plot is log10
summary_stats_log10 <- summary_stats %>%
  mutate(across(c(ymin, ymax, Q1, Median, Q3, Mean), ~log10(.x)))
# For multiple boxplots:
box_plot_multiple <- ggplot(df, aes(x = nFeature_RNA, y = Sample_ID, color = Sample_ID, fill = Sample_ID)) +
  geom_boxplot(width = 0.5, alpha = 0.2, outlier.shape = NA) +
  theme_classic() +
  scale_x_log10() +
  ggtitle("Boxplot nFeature") +
  xlab("") + # No x-label here
  ylab("") + # No y-label here, as it will be provided by the density plot
  theme(legend.position = "none",
        axis.text.y = element_blank(), # Hide y-axis text to clean up
        axis.ticks.y = element_blank(), # Hide y-axis ticks
        axis.title.y = element_blank(), # Hide y-axis title
        plot.margin = margin(b = 0)) # Reduce bottom margin
# For individual boxplots for each Sample_ID, arranged on top:
combined_plot <- box_plot_multiple / density_plot + plot_layout(heights = c(2, 3)) # Adjust heights as needed
path <- paste0("../results/density/",projectID,"_nFeature")
print(combined_plot)
saveToPDF(paste0(path, ".pdf"), width = 8, height = 8)

density_plot <- ggplot(df, aes(x = percent.mt, color = Sample_ID, fill = Sample_ID)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  xlab("percent.mt") +
  ylab("Density") +
  ggtitle("Density percent.mt") +
  theme(axis.title.x = element_blank(), # Remove x-axis title for cleaner alignment
        axis.text.x = element_blank(),   # Remove x-axis text as boxplot will provide labels
        axis.ticks.x = element_blank())  # Remove x-axis ticks

# Calculate summary statistics per Sample_ID
summary_stats <- df %>%
  group_by(Sample_ID) %>%
  summarise(
    ymin = min(percent.mt),
    ymax = max(percent.mt),
    Q1 = quantile(percent.mt, 0.25),
    Median = quantile(percent.mt, 0.5),
    Q3 = quantile(percent.mt, 0.75),
    Mean = mean(percent.mt)
  )

# Convert to log10 if your x-axis in the plot is log10
summary_stats_log10 <- summary_stats %>%
  mutate(across(c(ymin, ymax, Q1, Median, Q3, Mean),))
# For multiple boxplots:
box_plot_multiple <- ggplot(df, aes(x = percent.mt, y = Sample_ID, color = Sample_ID, fill = Sample_ID)) +
  geom_boxplot(width = 0.5, alpha = 0.2, outlier.shape = NA) +
  theme_classic() +
  ggtitle("Boxplot percent.mt") +
  xlab("") + # No x-label here
  ylab("") + # No y-label here, as it will be provided by the density plot
  theme(legend.position = "none",
        axis.text.y = element_blank(), # Hide y-axis text to clean up
        axis.ticks.y = element_blank(), # Hide y-axis ticks
        axis.title.y = element_blank(), # Hide y-axis title
        plot.margin = margin(b = 0)) # Reduce bottom margin
# For individual boxplots for each Sample_ID, arranged on top:
combined_plot <- box_plot_multiple / density_plot + plot_layout(heights = c(2, 3)) # Adjust heights as needed
path <- paste0("../results/density/",projectID,"_percent.mt")
print(combined_plot)
saveToPDF(paste0(path, ".pdf"), width = 8, height = 8)

density_plot <- ggplot(df, aes(x = percent.MALAT1, color = Sample_ID, fill = Sample_ID)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  scale_x_log10() +
  xlab("percent.MALAT1") +
  ylab("Density") +
  ggtitle("Density percent.MALAT1") +
  theme(axis.title.x = element_blank(), # Remove x-axis title for cleaner alignment
        axis.text.x = element_blank(),   # Remove x-axis text as boxplot will provide labels
        axis.ticks.x = element_blank())  # Remove x-axis ticks

# Calculate summary statistics per Sample_ID
summary_stats <- df %>%
  group_by(Sample_ID) %>%
  summarise(
    ymin = min(percent.MALAT1),
    ymax = max(percent.MALAT1),
    Q1 = quantile(percent.MALAT1, 0.25),
    Median = quantile(percent.MALAT1, 0.5),
    Q3 = quantile(percent.MALAT1, 0.75),
    Mean = mean(percent.MALAT1)
  )

# Convert to log10 if your x-axis in the plot is log10
summary_stats_log10 <- summary_stats %>%
  mutate(across(c(ymin, ymax, Q1, Median, Q3, Mean), ~log10(.x)))
# For multiple boxplots:
box_plot_multiple <- ggplot(df, aes(x = percent.MALAT1, y = Sample_ID, color = Sample_ID, fill = Sample_ID)) +
  geom_boxplot(width = 0.5, alpha = 0.2, outlier.shape = NA) +
  theme_classic() +
  scale_x_log10() +
  ggtitle("Boxplot percent.MALAT1") +
  xlab("") + # No x-label here
  ylab("") + # No y-label here, as it will be provided by the density plot
  theme(legend.position = "none",
        axis.text.y = element_blank(), # Hide y-axis text to clean up
        axis.ticks.y = element_blank(), # Hide y-axis ticks
        axis.title.y = element_blank(), # Hide y-axis title
        plot.margin = margin(b = 0)) # Reduce bottom margin
# For individual boxplots for each Sample_ID, arranged on top:
combined_plot <- box_plot_multiple / density_plot + plot_layout(heights = c(2, 3)) # Adjust heights as needed
path <- paste0("../results/density/",projectID,"_percent.MALAT1")
print(combined_plot)
saveToPDF(paste0(path, ".pdf"), width = 8, height = 8)