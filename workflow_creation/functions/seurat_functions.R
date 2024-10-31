# Seurat Functions for target.
read_data <- function(directory) {
    pbmc.data <- Read10X(data.dir = directory)
    pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
}

append_genome_pattern <- function(pbmc, genome_pattern) {
    pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = genome_pattern)
    pbmc
}

violin_plot <- function(pbmc) {
    VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
}

persist_plot <- function(plot, outfile) {
    ggsave(outfile, plot)
    outfile
}