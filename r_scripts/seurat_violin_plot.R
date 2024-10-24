#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied: genome_pattern", call.=FALSE)
}

# install.packages("Matrix", type = "source")
# install.packages("irlba", type = "source")
# install.packages('Seurat')

library(Seurat)
library(ggplot2)
library(sctransform)
#may need to downgrade to Matrix 1.6-1
library(Matrix)
library(irlba)
library(RSpectra)


pbmc.data <- Read10X(data.dir = "/home/hananli/Documents/bioinformatics_platform/data/filtered_gene_bc_matrices/hg19/")
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
# pbmc

#count matrix 
# pbmc@assays$RNA@counts

genome_pattern <- args[1]
#QC: The percentage of reads that map to the mitochondrial genome, MT can be ^Mt- or ^MT- depending on species of data.
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = genome_pattern)

#QC: nCount_RNA (The number of unique genes detected in each cell) and nFeature_RNA (the total number of molecules detected within a cell)
# head(pbmc@meta.data)


#QC Metric visualisation
png(filename="/home/hananli/Documents/bioinformatics_platform/data/outfiles/out.png")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

