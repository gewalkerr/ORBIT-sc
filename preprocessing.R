# ------------------ Seurat scRNAseq for ORBIT longitudinal study analysis ---------------------


# Load libraries
library(ggplot2)
library(tidyverse)
library(Seurat)
library(patchwork)
library(scales)

# ------------------ Generating Seurat Objects per Sample ---------------------

## Set file paths
setwd('~/orbit-sc')
data_path = "/home/modianoj/shared/RIS_Analysis/ORBIT_lymphoma_scRNAseq/Cellranger_Out/"

## Create Seurat objects


# Read in sample names
samples <- readLines("samples.txt")

# Initiate list to store objects
sample_objects <- list()

# Iterate over samples
for (sample_id in samples) {
  # Read in filtered data
  print(paste('Reading in ', sample_id, sep=""))
  # Get full dir
  directory = paste(data_path, sample_id,  "/outs/filtered_feature_bc_matrix/", sep="")
  data <- Read10X(data.dir = directory)
  
  # Generate Seurat object
  batch_key = substr(sample_id, 17, nchar(sample_id)-7) # Grab only unique ID
  seurat_obj <- CreateSeuratObject(counts = data, project = batch_key, min.cells = 3, min.features = 200)
  
  # Append curr object to list
  sample_objects[[batch_key]] <- seurat_obj

  # Save object
  print(paste("Saving ", sample_id, sep=""))
  saveRDS(seurat_obj, file = paste("objects/", batch_key, ".rds", sep=""))
}

# ---------------------------- Merge Samples to Compare ----------------------------

merged <- merge(sample_objects[[1]], y = sample_objects[2:length(sample_objects)], 
             add.cell.ids = names(sample_objects), project="Hank")

saveRDS(merged, "objects/samples_merged.rds")


# ---------------------- Initial Plots to Visualize Counts -------------------------

## Visualize number of cells per sample

# Preserve sample order
merged@meta.data$orig.ident <- factor(merged@meta.data$orig.ident, levels = unique(merged@meta.data$orig.ident))
# Plot
merged@meta.data %>% 
  ggplot(aes(x=orig.ident, fill=orig.ident)) + 
  geom_bar(color="black") +
  stat_count(geom = "text", colour = "black", size = 3.5, 
             aes(label = after_stat(count)),
             position=position_stack(vjust=1.05))+
  theme_classic() +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Number of Cells per Sample") +
  labs(x = "Sample", fill = "Sample")
ggsave("figures/preprocessing/cells_per_sample.png", width = 8, height = 5, dpi=400)

## Visualize nCount and nFeatures across samples

VlnPlot(merged, features = "nCount_RNA", layer="counts", group.by="orig.ident",
        alpha=0.2, pt.size = 0.1) + labs(x = "Sample", fill = "Sample")
ggsave("figures/preprocessing/raw_count_distribution.png", width = 8, height = 5, dpi=400)

VlnPlot(merged, features = "nFeature_RNA", layer="counts", group.by="orig.ident",
        alpha=0.2, pt.size = 0.1) + labs(x = "Sample", fill = "Sample")
ggsave("figures/preprocessing/raw_feat_distribution.png", width = 8, height = 5, dpi=400)

FeatureScatter(merged, "nCount_RNA", "nFeature_RNA", group.by = "orig.ident", pt.size = .3, raster=FALSE) +
  ggtitle("Raw gene count by feature count") +
  labs(fill = "Sample")
ggsave("figures/preprocessing/count_by_feature_scatter.png", width = 8, height = 5, dpi=400)

# Visualize MT counts across samples

merged[["percent.mt"]] <- PercentageFeatureSet(merged, pattern = "^MT-")
VlnPlot(merged, features = "percent.mt", layer="counts", alpha=0.2) +
  labs(x = "Sample", fill = "Sample") +
  geom_hline(yintercept = 10, color="red", linetype='dotted') +
  geom_hline(yintercept = 5, color="red", linetype='dotted')
ggsave("figures/preprocessing/mt_count_all_samples.png", width = 8, height = 5, dpi=400)

## Density plots
metadata <- merged@meta.data

# Features
metadata %>% 
  ggplot(aes(color=orig.ident, x=nFeature_RNA, fill= orig.ident)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 650,color="red",linetype="dotted") +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("nFeature Density Across Samples") +
  guides(color = guide_legend(title = "Sample"), 
         fill = guide_legend(title = "Sample"))
ggsave("figures/preprocessing/feature_density_all_samples.png", width = 8, height = 5, dpi=400)

# Counts
metadata %>% 
  ggplot(aes(color=orig.ident, x=nCount_RNA, fill= orig.ident)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 650,color="red",linetype="dotted") +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("nCount Density Across Samples") +
  guides(color = guide_legend(title = "Sample"), 
         fill = guide_legend(title = "Sample"))
ggsave("figures/preprocessing/count_density_all_samples.png", width = 8, height = 5, dpi=400)

# MT
metadata %>% 
  ggplot(aes(x=percent.mt,fill=orig.ident)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10(labels = label_comma()) + 
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("MT Count Density Across Samples") +
  guides(color = guide_legend(title = "Sample"), 
         fill = guide_legend(title = "Sample"))
ggsave("figures/preprocessing/mt_density_all_samples.png", width = 8, height = 5, dpi=400)


# ----------------------- Individual Samples Preprocessing -------------------------------
setwd('~/orbit-sc/objects/raw/')

# Read in single RDS
sample = 'pretx' 
object = readRDS('PreTx.rds')

# Calculate percent MT
object[["percent.mt"]] <- PercentageFeatureSet(object, pattern = "^MT-")

# Visualize
# Violin plots
VlnPlot(object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), layer="counts", ncol = 3)

# Scatter plots
plot1 <- FeatureScatter(object, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Set thresholds
mt = 5
feature_min = 300
feature_max = 5000
count_min = 600
count_max = 30000



# Density plots
p1 <- object@meta.data %>% 
  ggplot(aes(color=orig.ident, x=nFeature_RNA, fill= orig.ident)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() + 
  scale_x_log10() + 
  geom_vline(xintercept = feature_min,color="red",linetype="dotted") +
  geom_vline(xintercept = feature_max,color="red",linetype="dotted") +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  theme(legend.position = "none") +
  ggtitle("nFeature")
p2 <- object@meta.data %>% 
  ggplot(aes(color=orig.ident, x=nCount_RNA, fill= orig.ident)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() + 
  scale_x_log10() + 
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  theme(legend.position = "none") +
  ggtitle("nCount")
p3 <- object@meta.data %>% 
  ggplot(aes(x=percent.mt,fill=orig.ident)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10(labels = label_comma()) + 
  geom_vline(xintercept = mt,color="red",linetype="dotted") +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("MT Count")
p1 + p2 + p3
fig_path = paste("~/orbit-sc/figures/preprocessing/by_sample/", sample, "_densities.png", sep="")
ggsave(fig_path, width = 10, height = 4, dpi=400)

# Scatter (all metrics + thresholds)
ggplot(object@meta.data) +
  geom_point(aes(x=nCount_RNA,y=nFeature_RNA,fill=percent.mt > 5),shape=21,alpha=0.4) + 
  theme_classic() +
  scale_x_log10()+
  scale_y_log10()+
  #facet_grid(.~cond_tp) +
  geom_vline(xintercept = count_min,color="red",linetype="dotted")+
  geom_hline(yintercept=feature_min,color="red", linetype="dotted")+
  scale_fill_manual(values=c("FALSE"="lightblue", "TRUE"="purple"))  # Customize colors
fig_path = paste("~/orbit-sc/figures/preprocessing/by_sample/", sample, "_thresholds.png", sep="")
ggsave(fig_path, width = 6, height = 5, dpi=400)


## Filter out low quality cells
object <- subset(object, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)



# Visualize post filtering
VlnPlot(test, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), layer="counts", ncol = 3)

plot1 <- FeatureScatter(test, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(test, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2


# Save object
saveRDS(pre_tx, file = "......pre_tx.rds") # Edit file path










# Examine data
# A few key genes across certain # of columns
pbmc.data[c("CD3D", "TCL1A", "MS4A1"), 1:30]





