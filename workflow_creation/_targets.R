library(targets)
tar_source("functions/seurat_functions.R")
tar_option_set(packages = c("Seurat", "Matrix", "irlba", "RSpectra", "ggplot2", "sctransform", "irlba"), error="null")

list(
    
    tar_target(read_seurat, read_data('/home/hananli/Documents/bioinformatics_platform/data/filtered_gene_bc_matrices/hg19/')), 
    
    tar_target(get_genome_pattern, append_genome_pattern(read_seurat, '^MT-')), 
    
    tar_target(v_plot, violin_plot(get_genome_pattern)), 
    
    tar_target(save_data, persist_plot(v_plot, '/home/hananli/Documents/bioinformatics_platform/data/outfiles/targets_out.png'), format = 'file')
    
)