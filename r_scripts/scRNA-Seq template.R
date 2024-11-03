install.packages("Matrix", type = "source")
install.packages("irlba", type = "source")

library(Seurat)
library(ggplot2)
library(sctransform)
library(Matrix)
library(irlba)
library(RSpectra)


pbmc.data <- Read10X(data.dir = "~/Downloads/filtered_gene_bc_matrices/hg19/")
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc

#count matrix 
#pbmc@assays$RNA@counts


#QC: The percentage of reads that map to the mitochondrial genome, MT can be ^Mt- or ^MT- depending on species of data.
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

#QC: nCount_RNA (The number of unique genes detected in each cell) and nFeature_RNA (the total number of molecules detected within a cell)
head(pbmc@meta.data)


#QC Metric visualisation
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#QC filtering
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)


#Normalisation. By default, we employ a global-scaling normalization method “LogNormalize” that normalizes the feature expression measurements for each cell by the total expression, multiplies this by a scale factor (10,000 by default), and log-transforms the result. 
#pbmc <- NormalizeData(pbmc)

#Alternative is The use of SCTransform (v2) replaces the need to run NormalizeData, FindVariableFeatures, or ScaleData (described below.)
#Transformed data will be available in the SCT assay, which is set as the default after running sctransform
#by default returns 3,000 variable features 
#remotes::install_version("matrixStats", version="1.1.0")
pbmc <- SCTransform(pbmc, vars.to.regress = "percent.mt", verbose = FALSE)


#normalised values - input for PCA
pbmc@assays$SCT@scale.data

#'corrected' UMI counts - for visualisation
pbmc@assays$SCT@counts

#log-normalised 'corrected' UMI counts - for visualisation
pbmc@assays$SCT@data

#Normalised+log-transformed values
#pbmc@assays$RNA@data

#identify high variable features. Default returns 2000 feature. No need to run if using SCTransform.
#pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)


#Scale data. Linear transformation prior to PCA.
#all.genes <- rownames(pbmc)
#pbmc <- ScaleData(pbmc, features = all.genes)

#Scaled values
#pbmc@assays$RNA@scale.data



#Dimensionality reduction and Clustering 
pbmc <- RunPCA(pbmc, verbose = FALSE)
pbmc <- RunUMAP(pbmc, dims = 1:30, verbose = FALSE)

#you can also determine dimension used for FindNeighbours (next step)
#ElbowPlot(pbmc)

#we should ideally also try different resolutions for clustering to see if cells are over-clustered or under-clustered!!
pbmc <- FindNeighbors(pbmc, dims = 1:30, verbose = FALSE)
#We find that setting this parameter between 0.4-1.2 typically returns good results for single-cell datasets of around 3K cells. Optimal resolution often increases for larger datasets.
pbmc <- FindClusters(pbmc, verbose = FALSE)

DimPlot(pbmc, reduction = 'umap', label = TRUE)


#you can do this in one step
#pbmc <- CreateSeuratObject(pbmc_data) %>%
#  PercentageFeatureSet(pattern = "^MT-", col.name = "percent.mt") %>%
#  SCTransform(vars.to.regress = "percent.mt") %>%
#  RunPCA() %>%
#  FindNeighbors(dims = 1:30) %>%
#  RunUMAP(dims = 1:30) %>%
#  FindClusters()


#VlnPlot(pbmc, features = c("CD8A", "GZMK", "CCL5", "S100A4", "ANXA1", "CCR7", "ISG15", "CD3D"),
#        pt.size = 0.2, ncol = 4)


# Identify gene markers - should specify we can specify thresholds for the minimum percentage of cells expressing the gene in either of the two groups of cells (min.pct) and minimum difference in expression between the two groups (min.dff.pct).
all_markers <-FindAllMarkers(pbmc, 
                             min.pct =  0.25, 
                             min.diff.pct = 0.25)


# Return top 10 markers for cluster specified 'x'
gen_marker_table <- function(x){
  all_markers[all_markers$cluster == x, ] %>%
    head(n=10)
}




# Create a data frame of results for clusters (here it's 11 but we want to be able to automatically specify)
library(purrr)
top10_markers <- map_dfr(0:11, gen_marker_table)

#write.csv(top10_markers, "results/top10_markers.csv", quote = F)



# we may not be able to distinguish between two clusters. in that case: we find differential markers between the two to help us distinguish.
# for example
markers_0vs1 <- FindMarkers(object = pbmc, ident.1 = 0, ident.2 = 1)


# name the clusters - need to change to automated
current_cluster_ids <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11)


#Here we can implement SingleR to determine cluster - not done yet.
#https://bioconductor.org/books/release/SingleRBook/introduction.html

#Could also use TransferData function to assign labels based on reference tissue labels: https://learn.gencore.bio.nyu.edu/seurat-integration-and-label-transfer/

#Then, rename the clusters to the new IDs - below is just an example I used EnrichR (not sure if correct). Should use output from SingleR or TransferData
new_cluster_ids <- c("CD4+ T cells", "Dendritic cells", "B cells", "CD8+ Mature T cells", "Macrophages", "Monocytes", "NK cells","Macrophages", "CD8+ Naive T cells", "Langerhans cells", "Endothelial cells", "Platelets")

#BiocManager::install("SingleR")
library(SingleR)


#BiocManager::install("celldex")

library(celldex)
#hpca.se <- celldex::HumanPrimaryCellAtlasData()
#hpca.se



# Changing IDs to cell type
pbmc@active.ident <- plyr::mapvalues(x = pbmc@active.ident, 
                                from = current_cluster_ids, 
                                to = new_cluster_ids)

# Add active idents to metadata
pbmc <- AddMetaData(pbmc, pbmc@active.ident, col.name = "Cluster")
pbmc <- AddMetaData(pbmc, pbmc$orig.ident, col.name = "Sample")

library(dittoSeq)

# Plot new DimPlot
dittoDimPlot(pbmc, "Cluster")
# Plot new DimPlot split by samples
dittoDimPlot(pbmc, "Cluster", split.by= "Sample")


# Plot expression of a specific_gene (use gene name, full capital except in some cases) per cluster different ways
dittoPlot(pbmc, "ENO1", group.by = "Cluster")
dittoRidgePlot(pbmc, "ENO1", group.by = "Cluster")
dittoBoxPlot(pbmc, "ENO1", group.by = "Cluster")

# Plot expression of a specific_gene (use gene name, full capital except in some cases) in Dimplot
dittoDimPlot(pbmc, "ENO1")

# Plot barplot showing percentage of each cell type in each sample
dittoBarPlot(pbmc, "Cluster", group.by = "Sample")



# Plot heatmap showing the expression of top 75 variable genes by cluster
dittoHeatmap(pbmc, VariableFeatures(pbmc)[1:75], annot.by = "Cluster",scaled.to.max = TRUE)
# OR also by sample and order by cluster
dittoHeatmap(pbmc, VariableFeatures(pbmc)[1:75],
             annot.by = c("Cluster", "Sample"), scaled.to.max = TRUE)
# OR also by sample and order by sample
dittoHeatmap(pbmc, VariableFeatures(pbmc)[1:75],
             annot.by = c("Cluster", "Sample"), order.by = "Sample", scaled.to.max = TRUE)


dittoDimPlot(pbmc, "IL2RA")

#Set similar genes to a group (e.g. genes identifying a cell type or differentially expressed genes)
Treg.genes <- c("CD3E", "CD4", "IL2RA")

#VlnPlot the expression of a group of genes together by cluster
multi_dittoPlot(pbmc, Treg.genes, group.by = "Cluster",vlnplot.lineweight = 0.2, jitter.size = 0.3)
#DimPlot the expression of a group of genes together by cluster
dittoDimPlot(pbmc, Treg.genes)
dittoDimPlot(pbmc, Treg.genes, split.by = "Sample")


#Plot modular score of the combined gene set by cluster, split by sample
dittoPlotVarsAcrossGroups(pbmc, Treg.genes, group.by = "Cluster", split.by = "Sample",
                          main = "Treg Markers", ylab = "Treg score")

#Add Modular score (in a different way to above)
pbmc <- AddModuleScore(pbmc, features = list(Treg.genes), name="Treg Score")

library(RColorBrewer)

#Plot modular score in Dimplot
FeaturePlot(pbmc,
            features = "Treg Score1", label = TRUE, repel = TRUE) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))



#Cell cell communication (Cellchat, Cellphonedb)
#Cell trajectories

