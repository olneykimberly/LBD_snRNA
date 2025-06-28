## ----setup-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
knitr::opts_knit$set(root.dir = "/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/")
setwd("/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/")

## ----source------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
source(here::here("/tgen_labs/jfryer/kolney/LBD_CWOW/LBD_snRNA/scripts/", "file_paths_and_colours.R"))
projectID <- "CWOW_cellbender"
color.panel <- dittoColors()


## ----read_object--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# read object
dataObject <- readRDS(file = paste0("../rObjects/", projectID, "_filtered.rds"))


## ----split_object-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# split object by sample
dataObject.split <- SplitObject(dataObject, split.by = "Sample_ID") 
dataObject.split.doublets <- SplitObject(dataObject, split.by = "Sample_ID") 

## ----doubletFinder------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
for (i in 1:length(dataObject.split)) {
  # normalize and find PCs
  print(i)
  Sample_ID <- dataObject.split[[i]]$Sample_ID[[i]]
  print(Sample_ID)
  obj_sample <- NormalizeData(dataObject.split[[i]])
  sampleID <- levels(droplevels(obj_sample@meta.data$Sample_ID))
  obj_sample <- FindVariableFeatures(obj_sample, selection.method = "vst", nfeatures = 2000)
  obj_sample <- ScaleData(obj_sample)
  obj_sample <- RunPCA(obj_sample)
  
  # get significant PCs
  stdv <- obj_sample[["pca"]]@stdev
  sum.stdv <- sum(obj_sample[["pca"]]@stdev)
  percent.stdv <- (stdv / sum.stdv) * 100
  cumulative <- cumsum(percent.stdv)
  co1 <- which(cumulative > 90 & percent.stdv < 5)[1]
  co2 <- sort(which((percent.stdv[1:length(percent.stdv) - 1] - 
                       percent.stdv[2:length(percent.stdv)]) > 0.1), 
              decreasing = T)[1] + 1
  min.pc <- min(co1, co2)
  min.pc
  
  # run umap
  obj_sample <- RunUMAP(obj_sample, dims = 1:min.pc, reduction = "pca")
  
  # cluster
  obj_sample <- FindNeighbors(object = obj_sample, dims = 1:min.pc)                           
  obj_sample <- FindClusters(object = obj_sample, resolution = 0.2)
  
  # Assign identity of clusters
  Idents(object = obj_sample) <- "seurat_clusters"

  # number of cells in each cluster
  n_cells <- FetchData(obj_sample, vars = c("ident")) %>% dplyr::count(ident) %>%tidyr::spread(ident, n)
  
  ## pK Identification (no ground-truth) 
  sweep.res.list <- paramSweep(obj_sample, PCs = 1:min.pc, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  
  # Optimal pK for any scRNA-seq data can be manually discerned as maxima in BCmvn distributions
  bcmvn_max <- bcmvn[which.max(bcmvn$BCmetric),]
  pK_value <- bcmvn_max$pK
  pK_value <- as.numeric(levels(pK_value))[pK_value]
  
  # Homotypic Doublet Proportion Estimate 
  annotations <- obj_sample@meta.data$individual_clusters
  homotypic.prop <- modelHomotypic(annotations) 
  nExp_poi <- round(pK_value*nrow(obj_sample@meta.data))
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  # Run DoubletFinder with varying classification  
  obj_sample <- doubletFinder(obj_sample, PCs = 1:min.pc, 
                pN = 0.25, pK = pK_value, nExp = nExp_poi.adj, 
                reuse.pANN = FALSE, sct = FALSE)
  
  # set DF class for calling doublets
  DF_class <- obj_sample@meta.data[, grep("DF.classifications",colnames(obj_sample@meta.data)),]
  DF_class[which(DF_class == "Doublet")] <- "Doublet"
  table(DF_class)
  

  # table showing the number of doublets and singlets
  write.table(table(DF_class), paste0("../results/DoubletFinder/",projectID, "_doubletFinder_table_",Sample_ID), sep = "\t", 
              row.names = FALSE, quote = FALSE)
  obj_sample@meta.data[,"CellTypes_DF"] <- DF_class
  
  # plot
  d2 <- DimPlot(obj_sample, group.by="CellTypes_DF", reduction="umap",
          order=c("Coll.Duct.TC","Doublet"), 
          cols=c("#66C2A5","black"))
  path <- paste0("../results/DoubletFinder/",projectID,
               "_doubletFinder_UMAP_",Sample_ID)
  pdf(paste0(path, ".pdf"), width = 5,height = 4)
  print(d2)
  dev.off()
  
  # plot
  f1 <- FeaturePlot(obj_sample, 
            reduction = "umap", 
            features = c("nFeature_RNA", "nCount_RNA", 
                         "cell.complexity", "percent.mt"),
            pt.size = 0.4, 
            order = TRUE,
            label = TRUE)
  path <- paste0("../results/DoubletFinder/",projectID,
               "_doubletFinder_FeaturePlot_",Sample_ID)
  pdf(paste0(path, ".pdf"), width = 10, height = 7)
  print(f1)
  dev.off()
  
  #only keep singlets
  obj_sample_singlets <- subset(obj_sample, subset = CellTypes_DF == "Singlet")
  
  #only keep doublets
  obj_sample_doublets <- subset(obj_sample, subset = CellTypes_DF == "Doublet")
  
  # number of cells in each cluster per and post removing doublets
  n_cells_singlets <- FetchData(obj_sample_singlets, vars = c("ident")) %>% dplyr::count(ident) %>% tidyr::spread(ident, n)
  n_cells_singlets
  # Find the column names that are in n_cells_singlets but not in n_cells
  missing_cols <- setdiff(colnames(n_cells), colnames(n_cells_singlets))

  # Add the missing column(s) to n_cells, setting their values to 0
  for (col in missing_cols) {
  n_cells_singlets[[col]] <- 0
  }

  # Now combine the two data frames
  ncells_per_cluster <- rbind(n_cells, n_cells_singlets)

  # View the combined result
  ncells_per_cluster

  row.names(ncells_per_cluster) <- c("Doublets and singlets", "Singlets only")
  ncells_per_cluster
  difference <- diff(as.matrix(ncells_per_cluster))
  difference <- as.data.frame(difference)
  row.names(difference) <- c("difference")
  cbind(difference, ncells_per_cluster)
  write.table(ncells_per_cluster, paste0(
    "../results/DoubletFinder/",projectID,
    "_doubletFinder_ncells_per_cluster_",Sample_ID, ".txt"), sep = "\t", 
    row.names = FALSE, quote = FALSE)
  # plot the number of cells in each cluster per and post doubletFinder
  ncell_matrix <- as.matrix(ncells_per_cluster)
  ncells_melt <- melt(ncell_matrix)
  colnames(ncells_melt) <- c("doublet type","cluster","number of cells")
  ncell_max <- ncells_melt[which.max(ncells_melt$`number of cells`),]
  ncell_max_value <- ncell_max$`number of cells`
  cellmax <- ncell_max_value + 800 # so that the figure doesn't cut off the text
  b1 <- ggplot(ncells_melt, aes(x = factor(cluster), y = `number of cells`,
                          fill = `doublet type`)) + 
    geom_bar(stat="identity", colour="black", width=1, position = position_dodge(width=0.8)) +
    geom_text(aes(label = `number of cells`), 
              position=position_dodge(width=0.9), vjust=-0.25, angle = 45, hjust=-.01) + 
    theme_classic() + scale_fill_manual(values = c("gray", "#66C2A5")) +
    ggtitle("Number of cells per cluster") +  xlab("cluster") +
    theme(axis.text.x = element_text(angle = 45, hjust=1)) +
    scale_y_continuous(limits = c(0,cellmax))
  path <- paste0("../results/DoubletFinder/",projectID,
               "_doubletFinder_barplot_ncells_per_cluster_",Sample_ID)
  pdf(paste0(path, ".pdf"), width = 7,height = 5)
  print(b1)
  dev.off()
  f2 <- FeaturePlot(obj_sample_singlets, 
            reduction = "umap", 
            features = c("nFeature_RNA", "nCount_RNA", 
                         "cell.complexity", "percent.mt"),
            pt.size = 0.4, 
            order = TRUE,
            label = TRUE)
  path <- paste0("../results/DoubletFinder/",projectID,
               "_doubletFinder_FeaturePlot_singlets_",Sample_ID)
  pdf(paste0(path, ".pdf"), width = 10,height = 7)
  print(f2)
  dev.off()
  
  # put the PMI together again
  dataObject.split[[i]] <- obj_sample_singlets
  dataObject.split.doublets[[i]] <- obj_sample_doublets
  
}


## -----------Merge------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# converge obj.split
dataObject.singlets <- merge(x = dataObject.split[[1]],
                             y = c(dataObject.split[[2]], dataObject.split[[3]], dataObject.split[[4]], dataObject.split[[5]], dataObject.split[[6]],
                                   dataObject.split[[7]], dataObject.split[[8]], dataObject.split[[9]], dataObject.split[[10]], dataObject.split[[11]],
                                   dataObject.split[[12]], dataObject.split[[13]], dataObject.split[[14]], dataObject.split[[15]], dataObject.split[[16]],
                                   dataObject.split[[17]], dataObject.split[[18]], dataObject.split[[19]], dataObject.split[[20]], dataObject.split[[21]], 
                                   dataObject.split[[22]], dataObject.split[[23]], dataObject.split[[24]], dataObject.split[[25]], dataObject.split[[26]],
                                   dataObject.split[[27]], dataObject.split[[28]], dataObject.split[[29]], dataObject.split[[30]], dataObject.split[[31]],
                                   dataObject.split[[32]], dataObject.split[[33]], dataObject.split[[34]], dataObject.split[[35]], dataObject.split[[36]],
                                   dataObject.split[[37]], dataObject.split[[38]], dataObject.split[[39]]),
                             project = paste0("cwow"))

# print how many cells removed
print(paste0(dim(dataObject)[2] - dim(dataObject.singlets)[2]," nuclie removed"))

dataObject.doublets <- merge(x = dataObject.split.doublets[[1]],
                             y = c(dataObject.split.doublets[[2]], dataObject.split.doublets[[3]], dataObject.split.doublets[[4]], dataObject.split.doublets[[5]], dataObject.split.doublets[[6]],
                                   dataObject.split.doublets[[7]], dataObject.split.doublets[[8]], dataObject.split.doublets[[9]], dataObject.split.doublets[[10]], dataObject.split.doublets[[11]],
                                   dataObject.split.doublets[[12]], dataObject.split.doublets[[13]], dataObject.split.doublets[[14]], dataObject.split.doublets[[15]], dataObject.split.doublets[[16]],
                                   dataObject.split.doublets[[17]], dataObject.split.doublets[[18]], dataObject.split.doublets[[19]], dataObject.split.doublets[[20]], dataObject.split.doublets[[21]], 
                                   dataObject.split.doublets[[22]], dataObject.split.doublets[[23]], dataObject.split.doublets[[24]], dataObject.split.doublets[[25]], dataObject.split.doublets[[26]],
                                   dataObject.split.doublets[[27]], dataObject.split.doublets[[28]], dataObject.split.doublets[[29]], dataObject.split.doublets[[30]], dataObject.split.doublets[[31]],
                                   dataObject.split.doublets[[32]], dataObject.split.doublets[[33]], dataObject.split.doublets[[34]], dataObject.split.doublets[[35]], dataObject.split.doublets[[36]],
                                   dataObject.split.doublets[[37]], dataObject.split.doublets[[38]], dataObject.split.doublets[[39]]),
                             project = paste0("cwow"))
# confirm number 
print(paste0(dim(dataObject.doublets)[2]," total doublets"))

## -------------Save----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
saveRDS(dataObject.singlets, paste0("../rObjects/",projectID,"_singlets.rds"), compress = FALSE)
saveRDS(dataObject.doublets, paste0("../rObjects/",projectID,"_doublets.rds"), compress = FALSE)

# Inspect
print(paste0("inspect singlets object"))
dataObject.singlets

## -------------Plot----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Set path
doublet_tables_dir <- "../results/DoubletFinder"

# --- 1. Read and Combine Data Tables  ---
# Get a list of all files matching your pattern
file_list <- list.files(path = doublet_tables_dir,
                        pattern = paste0(projectID, "_doubletFinder_table_.*"),
                        full.names = TRUE)

# Initialize an empty list to store data frames
all_samples_data <- list()

# Loop through each file, read it, extract sample name, and store in the list
for (file_path in file_list) {
  # Extract the base filename (e.g., "CWOW_cellbender_doubletFinder_table_LBD_ATS_F4.tsv")
  base_filename <- basename(file_path)
  
  # Extract the sample name using string manipulation
  # This assumes the sample name is always after "table_" and before the file extension
  sample_name_pattern <- paste0(projectID, "_doubletFinder_table_")
  sample_name <- gsub(sample_name_pattern, "", base_filename) # Remove the prefix
  sample_name <- sub("\\..*$", "", sample_name) # Remove file extension (e.g., .tsv)
  
  # Read the table
  df <- read_tsv(file_path, show_col_types = FALSE)
  
  # Add a 'Sample' column
  df$Sample <- sample_name
  
  # Add to our list
  all_samples_data[[sample_name]] <- df
}

# If you want to combine all data frames into a single data frame:
combined_doublet_data <- bind_rows(all_samples_data)
# Now match the sample order as in the dataObject
desired_sample_order <- unique(dataObject$Sample_ID)

# Convert the 'Sample' column in combined_doublet_data to a factor
# with levels set to the desired order
combined_doublet_data$Sample <- factor(combined_doublet_data$Sample,
                                       levels = desired_sample_order)

# --- 2. Data Preparation for Plotting ---
# Ensure 'DF_class' is a factor for consistent plotting order (e.g., Singlet then Doublet)
combined_doublet_data$DF_class <- factor(combined_doublet_data$DF_class, levels = c("Doublet", "Singlet"))

# Calculate the percentage of doublets for each sample
doublet_percentage_df <- combined_doublet_data %>%
  group_by(Sample) %>%
  dplyr::summarise(
    Total_Doublets = sum(Freq[DF_class == "Doublet"]),
    Total_Cells = sum(Freq),
    Percentage_Doublets = (Total_Doublets / Total_Cells) * 100
  ) %>%
  dplyr::select(Sample, Percentage_Doublets)

# Print the resulting table
print(doublet_percentage_df)

#---- Prepare data for labels ---
# Calculate the total height of each bar for label positioning
label_data <- combined_doublet_data %>%
  group_by(Sample) %>%
  dplyr::summarise(Freq = sum(Freq)) %>%
  # Join with the doublet percentage data
  left_join(doublet_percentage_df, by = "Sample")

# --- 3. Create the Stacked Bar Plot ---
stacked_bar_plot <- ggplot(combined_doublet_data, aes(x = Sample, y = Freq)) + # Removed fill from here
  geom_bar(stat = "identity", position = "stack", aes(fill = DF_class)) + # Added fill here
  scale_fill_manual(values = c("Singlet" = "#66C2A5", "Doublet" = "black")) +
  # Add the text labels
  geom_text(data = label_data,
            aes(x = Sample, y = Freq, # x is Sample, y is the total bar height
                label = paste0(round(Percentage_Doublets, 1), "%")), # Format the label
            vjust = -0.5, # Adjust vertical position slightly above the bar
            size = 3.5,   # Adjust text size as needed
            color = "black") + # Text color
  labs(title = "Percentage of doublets per sample",
       x = "Sample",
       y = "Nuclei Count",
       fill = "Cell Classification") +
  theme_classic() +
  scale_y_continuous(breaks = seq(0,30000, by = 2000), limits = c(0,30000)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), # Rotate x-axis labels for readability
        plot.title = element_text(hjust = 0.5)) # Center the plot title

# Print the plot
print(stacked_bar_plot)

# --- Optional: Save the plot ---
path <- paste0("../results/DoubletFinder/",projectID,
               "_stacked_barplot_all_samples")
pdf(paste0(path, ".pdf"), width = 13,height = 5)
print(stacked_bar_plot)
dev.off()
