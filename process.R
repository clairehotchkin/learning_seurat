library(tidyverse)
library(Seurat)
library(patchwork)
library(umap)

# read in data
pbmc_data_filename = "filtered_gene_bc_matrices/hg19/"
pbmc_data <- Read10X(data.dir = pbmc_data_filename)

pbmc <- CreateSeuratObject(counts = pbmc_data, project = "pbmc3k", min.cells = 3, min.features = 200)

# subset to look at first 30 cells
pbmc_data[c("CD3D", "TCL1A", "MS4A1"), 1:30]

### FILTERING DATA ###

# use [[ to add columns to object metadata 
# use to store quality control stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# show QC metrics for first 5 cells 
head(pbmc@meta.data, 5)

# visualize QC metrics in order to filter cells 
    # filter for unique feature counts 200 < x > 2,5000 
    # filter for >5% mitochondrial counts 
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)

# use FeatureScatter to visualize feature-feature relationships
    # can also be used to visualize anything calculated by the object (ex. olumns in object metadata, PC scores)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

### NORMALIZING DATA ###

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
# normalized values stored in: pbmc[["RNA"]]@data

### FEATURE SELECTION ###

# calculate a subset of features that exhibit high cell-to-cell variation 
    # i.e. features that are highly expressed in some cells and lowly expressed in others 

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# identify the 10 most highly variable genes 
top10 <- head(VariableFeatures(pbmc), 10)

plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

### SCALING DATA ###
    # AKA apply linear transformation
    # ScaleData fnc: 
        # shifts expression of each gene so that the mean expression across cells is 0
        # scales expression of each gene so that variance across cells is 1 
    # results are stored in pbmc[["RNA"]]@scale.data

all_genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all_genes)


### LINEAR DIMENSION REDUCTION ###
# perform linear dimensional reduction
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# examine and visualize PCA results 
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")

DimPlot(pbmc, reduction = "pca")

# use dimHeatmap to explore sources of heterogeneity in a dataset 
    # useful when deciding which PCs to include for further analyses 
    # cells and features are ordered according to their PCA scores 
    # setting cells = number plots the 'extreme' cells on both ends of the spectrum
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)

# determine the dimensionality of the dataset 
    # the top PCs represent a robust compression of the dataset 
    # the question remains: how many components should we include? 10? 20? 100?
    # can implement a resampling test: 
        # randomly permute a subset of the data and rerun PCA
              # this constructs a 'null distribution' of feature scores 
        # repeat this procedure 
        # identify 'significant' PCs as those who have a strong enrichment of low p-value features
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)

# visualize using JackStraw Plot: compares distribution of p-values for each PC with a uniform distribution (dashed line)
    # significant PCs are those that have a "strong enrichment" of low p-value features
    # in the below plot, there is a sharp drop-off in significance after the first 10-12 PCs
JackStrawPlot(pbmc, dims = 1:15)

# elbow plots: another visualization tool
    # elbow plots rank PCs based on the percentage of variance explained by each one
    # in the below plot, there is an "elbow" around PC 9-10
    # this suggests that most of the "true signal" is captured in the first 10 PCs

ElbowPlot(pbmc)

### CLUSTERING CELLS ###
pbmc <- FindNeighbors(pbmc, dims = 1:10) # chooses first 10 PCs

# resolution parameter sets the 'granularity' of downstream clustering
    # increased values lead to more clusters 
    # set at 0.4-1.2 for datasets around 3k cells 
    # optimal resolution increases for larger datasets
pbmc <- FindClusters(pbmc, resolution = 0.5)

# look at cluster Ids of first 5 cells
# clusters can be found using Idents() fnc
head(Idents(pbmc), 5)

### NON-LINEAR DIMENSION REDUCTION ###
pbmc <- RunUMAP(pbmc, dims = 1:10)

# visualize. Can set 'label = TRUE' or use LabelClusters() to label individual clusters
DimPlot(pbmc, reduction = "umap")

# save object so that it can be easily loaded without the prior steps
saveRDS(pbmc, file = "../learning_seurat/pbmc_tutorial.rds")

### FINDING BIOMARKERS ###
# min.pct argument requires a feature to be detected at a minimum percentage in 
# either of the 2 groups of cells 
    # max.cells.per.ident can be set to downsample each identity class to have no more
    # cells than what its set to 
        # results in loss in power but speeds calculations 
cluster1.markers <- FindMarkers(pbmc, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 2, wt = avg_logFC) 

# one test for differential expression: ROC
    # ROC test returns the 'classification power' for any individual marker
cluster1.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

# visualize 
VlnPlot(pbmc, features = c("MS4A1", "CD79A"))

# plot raw counts
VlnPlot(pbmc, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)

FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", 
                               "CD8A"))

# heatmap plot
top10 <- pbmc.markers %>% 
  group_by(cluster) %>% 
  top_n(n = 10, wt = avg_logFC)

DoHeatmap(pbmc, features = top10$gene) + NoLegend()

### ASSIGNING CELL TYPE IDENTITY TO CLUSTERS ###
new.cluster.ids <- c("Naive CD4 T", "Memory CD4 T", "CD14+ Mono", "B", "CD8 T", "FCGR3A+ Mono", 
                    "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

saveRDS(pbmc, file = "../learning_seurat/pbmc3k_final.rds")
