setwd("~/Desktop/PhD/scRNA-Seq")
library(stats)
library(Seurat)
library(dplyr)
library(dplyr)
library(Matrix)
library(ggpubr)
library(ggplot2)
library(ggpubr)
library(gplots)
library(cowplot)
library(clustree)
library(ggraph)
library(patchwork)
library(dplyr)
library(data.table)
library(sctransform)
library(mixsqp)
library(devtools)

data_dir <- '~/Desktop/PhD/scRNA-Seq/D5D8D11_LN_Tumour_scRNASeq/Data/Compiled/'
list.files(data_dir)
Seurat_object <- Read10X(data.dir = data_dir)

metadata <- read.table("~/Desktop/PhD/scRNA-Seq/D5D8D11_LN_Tumour_scRNASeq/Data/Compiled/E-EHCA-2_new.sdrf", header=TRUE, sep="\t")
Seurat_object <- CreateSeuratObject(counts = Seurat_object, project = "E-EHCA-2.aggregated_filtered_normalised_counts", min.cells = 3, min.features = 200)
colOrd <- rownames(FetchData(Seurat_object,"ident"))
matched_metadata <- metadata %>% dplyr::slice(match(colOrd, Source.Name))
result <- matched_metadata[-1]
rownames(result) <- make.names(matched_metadata$Source.Name, unique=TRUE)
rownames(result)<-gsub("X","",rownames(result))
result

Seurat_object <- Read10X(data.dir = data_dir)
Seurat_object <- CreateSeuratObject(counts = Seurat_object, project = "E-EHCA-2.aggregated_filtered_counts", min.cells = 3, min.features = 200, meta.data = result)



#Seurat_object='matrix.mtx', rows='genes.tsv', cols='barcodes.tsv')
#https://bioinformatics.stackexchange.com/questions/5115/seurat-with-normalized-count-matrix
#https://www.biostars.org/p/460935/
#Seurat_object <- CreateSeuratObject(counts = data, project = "E-EHCA-2.aggregated_filtered_counts", meta.data = result, min.cells = 3, min.features = 200)
#Nd1','Nd2','Co1','Co2','Atp8','Atp6','Co3','Nd3','Nd4l','Nd4','Nd5','Nd6','Cytb','Rnr1','Rnr2'
#mtgenelist <- read.csv('~/Desktop/PhD/scRNA-Seq/D5D8D11_LN_Tumour_scRNASeq/Data/Compiled/mtgenelist.csv')[,1]
#Seurat_object <- PercentageFeatureSet(Seurat_object, features = c("ENSMUSG00000064370", "ENSMUSG00000064368", "ENSMUSG00000064367", "ENSMUSG00000064363", "ENSMUSG00000065947", "ENSMUSG00000064360", "ENSMUSG00000064358", "ENSMUSG00000064357", "ENSMUSG00000064356", "ENSMUSG00000064354", "ENSMUSG00000064351", "ENSMUSG00000064345", "ENSMUSG00000064341", "ENSMUSG00000064337", "ENSMUSG00000064339"), col.name = "percent.mt")

Seurat_object[["percent.mt"]] <- PercentageFeatureSet(object = Seurat_object, pattern = "^MT-")
remotes::install_version("matrixStats", version="1.1.0")
Seurat_object <- SCTransform(Seurat_object, method = "glmGampoi",vars.to.regress = "percent.mt", verbose = FALSE)


#mt.genes[!mt.genes %in% rownames(Seurat_object)]

Customgenelist <- read.csv('~/Desktop/PhD/scRNA-Seq/D5D8D11_LN_Tumour_scRNASeq/Data/Compiled/customgenelist.csv')[,1]

#Seurat_object <- FindVariableFeatures(object = Seurat_object)
#top10 <- head(VariableFeatures(Seurat_object), 10)
#plot1 <- VariableFeaturePlot(Seurat_object)
#plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
#plot2

DefaultAssay(Seurat_object) <- "SCT"
Seurat_object <- RunPCA(object = Seurat_object, features = Customgenelist, verbose = TRUE)
#Seurat_object <- RunUMAP(Seurat_object, dims = 1:28)



#Add custom cell identity to set Seurat cells: https://github.com/satijalab/seurat/issues/530
Seurat_object <- SetIdent(object = Seurat_object, value = "Factor.Value.inferred.cell.type...authors.labels.")


8482 -> saved.seed
set.seed(saved.seed)
RunTSNE(
  Seurat_object,
  dims=1:28,
  seed.use = saved.seed
) -> Seurat_object

DimPlot(Seurat_object,reduction = "tsne", pt.size = 1, label = TRUE, label.size = 3) + ggtitle("tSNE with default Perplexity (30)")
DimPlot(Seurat_object, reduction = "tsne", label=TRUE, label.size=2, split.by = "FactorValue.time.")
#TSNEPlot(object = Seurat_object, group.by = "Factor.Value.inferred.cell.type...authors.labels.")

#no. of reads in each cluster
VlnPlot(Seurat_object,features="nCount_RNA")

#no. of genes in each cluster
VlnPlot(Seurat_object,features="nFeature_RNA")



#Sub-clustering tumor Tregs
Seurat.tumourTcellcluster <- subset(Seurat_object, idents = 'tumor T cell')
Seurat.tumourTcellcluster <- FindNeighbors(Seurat.tumourTcellcluster, dims = 1:10, k.param = 5)
Seurat.tumourTcellcluster <- FindClusters(Seurat.tumourTcellcluster, resolution=0.9)
Seurat_object$sub_cluster <- as.character(Idents(Seurat_object))
Seurat_object$sub_cluster[Cells(Seurat.tumourTcellcluster)] <- paste("tumor",Idents(Seurat.tumourTcellcluster))
VlnPlot(Seurat_object, features = c("ENSMUSG00000039521", "ENSMUSG00000053977", "ENSMUSG00000023274", "ENSMUSG00000027276", "ENSMUSG00000026770"), group.by = "sub_cluster", pt.size = 0, combine = FALSE)
Seurat.tumourTcellcluster <- RenameIdents(Seurat.tumourTcellcluster, "0" = "CD8+ T cell", "1" = "Treg", "2" = "CD8+ T cell", "3" = "Teff", "4" = "CD8+ T cell", "5" = "Treg", "6" = "Treg", "7" = "Treg", "8" = "CD8+ T cell")
Seurat_object$sub_cluster[Cells(Seurat.tumourTcellcluster)] <- paste("tumor",Idents(Seurat.tumourTcellcluster))
DimPlot(Seurat_object, group.by = "sub_cluster")
Seurat_object <- SetIdent(object = Seurat_object, value = "sub_cluster")

DefaultAssay(Seurat_object) <- "RNA"

#expression of a specific gene (e.g. Jag1) in each original cluster or sub-cluster
VlnPlot(Seurat_object,features="ENSMUSG00000027276")
VlnPlot(Seurat_object,features="ENSMUSG00000027276", group.by = "FactorValue.time.")
VlnPlot(Seurat_object,features="ENSMUSG00000027276", group.by = 'sub_cluster')
#expression of a specific gene (e.g. Foxp3) in each cluster
VlnPlot(Seurat_object,features="ENSMUSG00000039521")
VlnPlot(Seurat_object,features="ENSMUSG00000039521", group.by = 'sub_cluster')
#expression of a specific gene (e.g. CD127; Tregs should be CD127-) in each cluster
VlnPlot(Seurat_object,features="ENSMUSG00000003882")

#Treg identity markers: Foxp3, Entpd1 (ENSMUSG00000048120; also expressed on CD8 T), Nt5e (ENSMUSG00000032420; also expressed on CD4 T), IL10 (ENSMUSG00000016529; also expressed on Th17 T)
#Treg markers (stability, activation, supression etc.)
#Ki67
VlnPlot(Seurat_object, features = c("ENSMUSG00000031004"), split.by = "FactorValue.time.", group.by = "sub_cluster",  pt.size = 0, combine = FALSE, idents = "tumor Treg", col=c("red", "green", "blue"))
#CD25/IL2RA
VlnPlot(Seurat_object, features = c("ENSMUSG00000026770"), split.by = "FactorValue.time.", group.by = "sub_cluster",  pt.size = 0, combine = FALSE, idents = "tumor Treg", col=c("red", "green", "blue"))
#CD27
VlnPlot(Seurat_object, features = c("ENSMUSG00000030336"), split.by = "FactorValue.time.", group.by = "sub_cluster",  pt.size = 0, combine = FALSE, idents = "tumor Treg", col=c("red", "green", "blue"))
#ICOS
VlnPlot(Seurat_object, features = c("ENSMUSG00000026009"), split.by = "FactorValue.time.", group.by = "sub_cluster",  pt.size = 0, combine = FALSE, idents = "tumor Treg", col=c("red", "green", "blue"))
#CTLA4 (suppression)
VlnPlot(Seurat_object, features = c("ENSMUSG00000026011"), split.by = "FactorValue.time.", group.by = "sub_cluster",  pt.size = 0, combine = FALSE, idents = "tumor Treg", col=c("red", "green", "blue"))
#FOXP3 (treg identity)
VlnPlot(Seurat_object, features = c("ENSMUSG00000039521"), split.by = "FactorValue.time.", group.by = "sub_cluster",  pt.size = 0, combine = FALSE, idents = "tumor Treg", col=c("red", "green", "blue"))
#IL10 (treg identity)
VlnPlot(Seurat_object, features = c("ENSMUSG00000016529"), split.by = "FactorValue.time.", group.by = "sub_cluster",  pt.size = 0, combine = FALSE, idents = "tumor Treg", col=c("red", "green", "blue"))
#Tigit (suppression)
VlnPlot(Seurat_object, features = c("ENSMUSG00000071552"), split.by = "FactorValue.time.", group.by = "sub_cluster",  pt.size = 0, combine = FALSE, idents = "tumor Treg", col=c("red", "green", "blue"))
#PD-1/PDCD1 (suppression)
VlnPlot(Seurat_object, features = c("ENSMUSG00000026285"), split.by = "FactorValue.time.", group.by = "sub_cluster",  pt.size = 0, combine = FALSE, idents = "tumor Treg", col=c("red", "green", "blue"))
#TGFb1
VlnPlot(Seurat_object, features = c("ENSMUSG00000002603"), split.by = "FactorValue.time.", group.by = "sub_cluster",  pt.size = 0, combine = FALSE, idents = "tumor Treg", col=c("red", "green", "blue"))
#CD39/Entpd1 (treg identity)
VlnPlot(Seurat_object, features = c("ENSMUSG00000048120"), split.by = "FactorValue.time.", group.by = "sub_cluster",  pt.size = 0, combine = FALSE, idents = "tumor Treg", col=c("red", "green", "blue"))
#JAG1
VlnPlot(Seurat_object, features = c("ENSMUSG00000027276"), split.by = "FactorValue.time.", group.by = "sub_cluster",  pt.size = 0, combine = FALSE, idents = "tumor Treg", col=c("red", "green", "blue"))
#ITGb8
VlnPlot(Seurat_object, features = c("ENSMUSG00000025321"), split.by = "FactorValue.time.", group.by = "sub_cluster",  pt.size = 0, combine = FALSE, idents = "tumor Treg", col=c("red", "green", "blue"))
#BTLA (suppression)
VlnPlot(Seurat_object, features = c("ENSMUSG00000052013"), split.by = "FactorValue.time.", group.by = "sub_cluster",  pt.size = 0, combine = FALSE, idents = "tumor Treg", col=c("red", "green", "blue"))
#LAG3 (suppression)
VlnPlot(Seurat_object, features = c("ENSMUSG00000030124"), split.by = "FactorValue.time.", group.by = "sub_cluster",  pt.size = 0, combine = FALSE, idents = "tumor Treg", col=c("red", "green", "blue"))
#Havcr2 (suppression)
VlnPlot(Seurat_object, features = c("ENSMUSG00000020399"), split.by = "FactorValue.time.", group.by = "sub_cluster",  pt.size = 0, combine = FALSE, idents = "tumor Treg", col=c("red", "green", "blue"))
#VSIR (suppression; low expression in paper)
VlnPlot(Seurat_object, features = c("ENSMUSG00000020101"), split.by = "FactorValue.time.", group.by = "sub_cluster",  pt.size = 1, combine = FALSE, idents = "tumor Treg", col=c("red", "green", "blue"))
#Tnfrsf18 (suppression)
VlnPlot(Seurat_object, features = c("ENSMUSG00000041954"), split.by = "FactorValue.time.", group.by = "sub_cluster",  pt.size = 1, combine = FALSE, idents = "tumor Treg", col=c("red", "green", "blue"))
#CD28 (co-stimulation)
VlnPlot(Seurat_object, features = c("ENSMUSG00000026012"), split.by = "FactorValue.time.", group.by = "sub_cluster",  pt.size = 1, combine = FALSE, idents = "tumor Treg", col=c("red", "green", "blue"))
#Tnfrsf9 (co-stimulation)
VlnPlot(Seurat_object, features = c("ENSMUSG00000028965"), split.by = "FactorValue.time.", group.by = "sub_cluster",  pt.size = 1, combine = FALSE, idents = "tumor Treg", col=c("red", "green", "blue"))
#Tnfrsf4 (co-stimulation)
VlnPlot(Seurat_object, features = c("ENSMUSG00000029075"), split.by = "FactorValue.time.", group.by = "sub_cluster",  pt.size = 1, combine = FALSE, idents = "tumor Treg", col=c("red", "green", "blue"))
#Nt5e (treg identity)
VlnPlot(Seurat_object, features = c("ENSMUSG00000032420"), split.by = "FactorValue.time.", group.by = "sub_cluster",  pt.size = 1, combine = FALSE, idents = "tumor Treg", col=c("red", "green", "blue"))
#gzmb (cytotoxicity)
VlnPlot(Seurat_object, features = c("ENSMUSG00000015437"), split.by = "FactorValue.time.", group.by = "sub_cluster",  pt.size = 1, combine = FALSE, idents = "tumor Treg", col=c("red", "green", "blue"))


DefaultAssay(object = Seurat_object) <- "RNA"

#FeaturePlot(Seurat_object, features = c("ENSMUSG00000031004", "ENSMUSG00000026770", "ENSMUSG00000030336", "ENSMUSG00000026009", "ENSMUSG00000026011", "ENSMUSG00000026770"), split.by = "FactorValue.time.", cols = c("grey", "red"), pt.size = 0, combine = FALSE)

plot1 <- VlnPlot(Seurat_object, features = c("ENSMUSG00000031004", "ENSMUSG00000026770", "ENSMUSG00000030336", "ENSMUSG00000026009"), split.by = "FactorValue.time.", group.by = "sub_cluster", pt.size = 0, combine = FALSE, idents = "tumor Treg", col=c("red", "green", "blue"))
wrap_plots(plots = plot1, ncol = 1)

Treg_actisupp_markers <- c("ENSMUSG00000031004", "ENSMUSG00000026770", "ENSMUSG00000030336", "ENSMUSG00000026009", "ENSMUSG00000026011", "ENSMUSG00000039521", "ENSMUSG00000016529", "ENSMUSG00000071552", "ENSMUSG00000026285", "ENSMUSG00000002603", "ENSMUSG00000048120", "ENSMUSG00000027276", "ENSMUSG00000025321", "ENSMUSG00000052013", "ENSMUSG00000030124", "ENSMUSG00000020399", "ENSMUSG00000020101", "ENSMUSG00000041954", "ENSMUSG00000026012", "ENSMUSG00000028965", "ENSMUSG00000029075", "ENSMUSG00000032420", "ENSMUSG00000015437")
#Treg_actisupp_markers <- read.csv('Treg_activation_suppression_markers.csv')[,1]
Treg_actsup_markers <- fread('Treg_activation_suppression_markers.csv')[,1]
Treg_actsup_markers = as.data.table(Treg_actsup_markers)
Treg_actisupp_markers_tumor_Treg_dotplotlabel <- c("Ki67","CD25/IL2RA","CD27", "ICOS", "CTLA4", "Foxp3", "IL10", "Tigit", "PD-1/PDCD1", "TGFb1", "CD39/Entpd1", "Jag1", "ITGb8", "BTLA", "LAG3", "Havcr2", "VSIR", "Tnfrsf18", "CD28", "Tnfrsf9", "Tnfrsf4", "Nt5e", "gzmb")
DotPlot(Seurat_object, features = Treg_actisupp_markers, cols = c("blue"), dot.scale = 8, group.by = "FactorValue.time.", idents = "tumor Treg") + RotatedAxis()+ scale_x_discrete(labels= Treg_actisupp_markers_tumor_Treg_dotplotlabel)

tumor_Treg <- WhichCells(Seurat_object, idents = 'tumor Treg')
DoHeatmap(Seurat_object, features = Treg_actisupp_markers, cells = tumor_Treg, group.by = "FactorValue.time.")

Treg_DC_crosstalk_markers <- c("ENSMUSG00000045382", "ENSMUSG00000050232", "ENSMUSG00000024621", "ENSMUSG00000025804", "ENSMUSG00000049103", "ENSMUSG00000035448", "ENSMUSG00000026285", "ENSMUSG00000020399", "ENSMUSG00000033467")
Treg_DC_crosstalk_markers_dotplotlabel <- c("CXCR4", "CXCR3", "CSF1R", "CCR1", "CCR2", "CCR3", "PDCD1", "HAVCR2", "TSLPR")
DotPlot(Seurat_object, features = Treg_DC_crosstalk_markers, cols = c("blue"), dot.scale = 8, group.by = "FactorValue.time.", idents = "tumor Treg") + RotatedAxis()+ scale_x_discrete(labels= Treg_DC_crosstalk_markers_dotplotlabel)


#Treg_new <- c("ENSMUSG00000064356", "ENSMUSG00000083365", "ENSMUSG00000083558", "ENSMUSG00000064360", "ENSMUSG00000034892")
#DotPlot(Seurat_object, features = Treg_new, cols = c("blue"), dot.scale = 8, group.by = "FactorValue.time.", idents = "tumor Treg") + RotatedAxis()


#FindMarkers(Seurat_object, ident.1 = "tumor T cell", ident.2 = "lymph node T cell")

#https://www.bioinformatics.babraham.ac.uk/training/10XRNASeq/seurat_workflow.html


#Sub-clustering LN Tregs
#Seurat.lnTcellcluster <- subset(Seurat_object, idents = 'lymph node T cell')
#Seurat.lnTcellcluster <- FindNeighbors(Seurat.lnTcellcluster, dims = 1:10, k.param = 5)
#Seurat.lnTcellcluster <- FindClusters(Seurat.lnTcellcluster, resolution=0.9)
#Seurat_object$sub_cluster <- as.character(Idents(Seurat_object))
#Seurat_object$sub_cluster[Cells(Seurat.lnTcellcluster)] <- paste("lymph node",Idents(Seurat.lnTcellcluster))
#VlnPlot(Seurat_object, features = c("ENSMUSG00000039521", "ENSMUSG00000053977", "ENSMUSG00000023274", "ENSMUSG00000027276", "ENSMUSG00000026770"), group.by = "sub_cluster", pt.size = 0, combine = FALSE)
#DimPlot(Seurat_object, group.by = "sub_cluster")
#Seurat.lnTcellcluster <- RenameIdents(Seurat.lnTcellcluster, "0" = "CD4+ Treg", "1" = "CD8+ T cell", "2" = "CD8+ Treg", "3" = "CD8+ T cell", "4" = "Teff", "5" = "CD8+ Treg", "6" = "Teff", "7" = "CD8+ T cell", "8" = "CD8+ T cell", "9" = "CD8+ T cell", "10" = "CD8+ T cell", "11" = "CD8+ T cell", "12" = "CD4+ Treg")
#Seurat_object$sub_cluster[Cells(Seurat.lnTcellcluster)] <- paste("lymph node",Idents(Seurat.lnTcellcluster))
#Seurat_object <- SetIdent(object = Seurat_object, value = "sub_cluster")





DimPlot(Seurat_object, group.by = "sub_cluster", split.by = "FactorValue.time.")


Seurat_object <- SetIdent(object = Seurat_object, value = "sub_cluster")
#"ENSMUSG00000039521" = Foxp3, "ENSMUSG00000053977" = CD8,  "ENSMUSG00000023274" = CD4, "ENSMUSG00000027276" = Jag1
#VlnPlot(Seurat_object, features = c("ENSMUSG00000039521", "ENSMUSG00000053977", "ENSMUSG00000023274", "ENSMUSG00000027276"), split.by = "FactorValue.time.", group.by = "sub_cluster", 
#        +         pt.size = 0, combine = FALSE)
VlnPlot(Seurat_object, features = c("ENSMUSG00000039521", "ENSMUSG00000053977", "ENSMUSG00000023274", "ENSMUSG00000027276"), split.by = "FactorValue.time.", group.by = "sub_cluster", pt.size = 0, combine = FALSE, idents="tumor Treg", col=c("red", "green", "blue"))



saveRDS(Seurat_object, file = "D5D8D11_LN_TumourNEW.rds")


#DE markers between Jag1+ and Jag1- tumor Treg
Seurat_object$Jag1.cell <- "Jag1.pos"
Seurat_object$Jag1.cell[WhichCells(Seurat_object, expression = ENSMUSG00000027276 < 5)] <- "Jag1.neg"
VlnPlot(Seurat_object, features = c("ENSMUSG00000027276"), split.by = "Jag1.cell", group.by = "sub_cluster",  pt.size = 0, combine = FALSE, idents = "tumor Treg", col=c("red", "green", "blue"))
DE_Markers_Jag1tumorTreg <- FindMarkers(Seurat_object, subset.ident = "tumor Treg", group.by = "Jag1.cell",ident.1 = "Jag1.neg", ident.2 = "Jag1.pos",  verbose = FALSE, assay = "SCT", slot = "data")
head(DE_Markers_Jag1tumorTreg, n = 20)
write.table(DE_Markers_Jag1tumorTreg,"DE_Markers_Jag1tumorTreg.txt",sep="\t",row.names=TRUE)

VlnPlot(Seurat_object, features = c("ENSMUSG00000006587"), split.by = "Jag1.cell", group.by = "sub_cluster",  pt.size = 1, combine = FALSE, idents = "tumor Treg", col=c("red", "green", "blue"))
DE_Markers_Jag1tumorTreg_dotplotlabel <- c("Gm49683","Jag1","Gm21854","Atp7b","Snai3","1700001P01Rik","Clvs2", "Upp1","Tshz3","Gm5457","Sall3","Srms","Asphd2","Tdgf1","Prtg","Anpep","Sema3d","Sgtb","Tlr7","Neto1")
DotPlot(Seurat_object, features = rownames(head(DE_Markers_Jag1tumorTreg, n=20)), cols = c("RdYlBu"), dot.scale = 8, group.by = "Jag1.cell", idents = "tumor Treg") + RotatedAxis() + scale_x_discrete(labels= DE_Markers_Jag1tumorTreg_dotplotlabel)
DotPlot(Seurat_object, features = rownames(head(DE_Markers_Jag1tumorTreg, n=20)), cols = c("RdYlBu"), dot.scale = 8, group.by = "Jag1.cell", split.by = "FactorValue.time.", idents = "tumor Treg") + RotatedAxis() + scale_x_discrete(labels= DE_Markers_Jag1tumorTreg_dotplotlabel)

Treg_cytokine_markers <- c("ENSMUSG00000016529","ENSMUSG00000025746", "ENSMUSG00000029375", "ENSMUSG00000027720", "ENSMUSG00000020383", "ENSMUSG00000055170", "ENSMUSG00000024401", "ENSMUSG00000000869", "ENSMUSG00000036117", "ENSMUSG00000025929", "ENSMUSG00000074695", "ENSMUSG00000040329", "ENSMUSG00000004371", "ENSMUSG00000003206", "ENSMUSG00000027776", "ENSMUSG00000039217")
Treg_cytokine_markers_tumor_Treg_dotplotlabel <- c("IL10", "IL6", "IL8", "IL2", "IL13", "IFNg", "TNFa", "IL4", "IL5", "IL17a", "IL22", "IL7", "IL11", "Ebi3", "IL12a", "IL18")
DotPlot(Seurat_object, features = Treg_cytokine_markers, cols = c("RdYlBu"), dot.scale = 8, group.by = "Jag1.cell", idents = "tumor Treg",  assay = "RNA") + RotatedAxis()+ scale_x_discrete(labels= Treg_cytokine_markers_tumor_Treg_dotplotlabel)
DotPlot(Seurat_object, features = Treg_cytokine_markers, cols = c("RdYlBu"), dot.scale = 8, group.by = "Jag1.cell", idents = "tumor Treg", split.by = "FactorValue.time.", assay = "RNA") + RotatedAxis()+ scale_x_discrete(labels= Treg_cytokine_markers_tumor_Treg_dotplotlabel)


Tumor_Treg_subset <- subset(x = Seurat_object, idents = c('tumor Treg'))
matrix<-Tumor_Treg_subset@assays$RNA@data
matrix_mod<-as.matrix(matrix)
gene<-as.numeric(matrix_mod["ENSMUSG00000027276",])
correlations<-apply(matrix_mod,1,function(x){cor(gene,x)})
write.table(correlations,"Jag1_correlations.txt",sep="\t")



#DE markers between Hopx+ and Hopx- tumor Treg
Seurat_object$Hopx.cell <- "Hopx.pos"
Seurat_object$Hopx.cell[WhichCells(Seurat_object, expression = ENSMUSG00000059325 < 150)] <- "Hopx.neg"
VlnPlot(Seurat_object, features = c("ENSMUSG00000059325"), split.by = "Hopx.cell", group.by = "sub_cluster",  pt.size = 0, combine = FALSE, idents = "tumor Treg", col=c("red", "green", "blue"))
DE_Markers_HopxtumorTreg <- FindMarkers(Seurat_object, subset.ident = "tumor Treg", group.by = "Hopx.cell",ident.1 = "Hopx.neg", ident.2 = "Hopx.pos",  verbose = FALSE)
head(DE_Markers_HopxtumorTreg, n = 20)

DotPlot(Seurat_object, features = c("ENSMUSG00000029075", "ENSMUSG00000028965"), cols = c("RdYlBu"), dot.scale = 8, group.by = "FactorValue.time.", split.by = "Hopx.cell", idents = "tumor Treg") + RotatedAxis()


#DE markers between D5, D8 and D11 tumor Treg (FindMarkers functino can only specify 2 identities max!!! If specify 3, will only compare idents 1 and 2. So in the case that you want to compare 3 conditions within same cluster, just specify cluster identity and write split by condition)
DE_Markers_tumor_Treg_D5vD8vD11 <- FindMarkers(Seurat_object, ident.1 = "tumor Treg", split.by ="FactorValue.time.", verbose = FALSE)
head(DE_Markers_tumor_Treg_D5vD8vD11, n = 20)
write.table(DE_Markers_tumor_Treg_D5vD8vD11,"DE_Markers_tumor_Treg_D5vD8vD11.txt",sep="\t",row.names=TRUE)

#Plot DE top hits between early vs. late tumor Tregs to see exactly when they are most differentially expressed
#VlnPlot(Seurat_object, features = c("ENSMUSG00000025321"), split.by = "FactorValue.time.", group.by = "sub_cluster", 
#        +         pt.size = 0, combine = FALSE)
#VlnPlot(Seurat_object, features = c("ENSMUSG00000108090"), split.by = "FactorValue.time.", group.by = "sub_cluster", 
#        +         pt.size = 0, combine = FALSE)
#VlnPlot(Seurat_object, features = c("ENSMUSG00000059325"), split.by = "FactorValue.time.", group.by = "sub_cluster", 
#        +         pt.size = 0, combine = FALSE)

#Plotting expression of top DE markers for D5, D8 and D11 tumor Treg (individually or top 20 altogether)
VlnPlot(Seurat_object, features = c("ENSMUSG00000025321"), split.by = "FactorValue.time.", group.by = "sub_cluster", 
        +         pt.size = 0, combine = FALSE, idents = "tumor Treg", col=c("red", "green", "blue"))
VlnPlot(Seurat_object, features = c("ENSMUSG00000108090"), split.by = "FactorValue.time.", group.by = "sub_cluster", 
        +         pt.size = 0, combine = FALSE, idents = "tumor Treg", col=c("red", "green", "blue"))
VlnPlot(Seurat_object, features = c("ENSMUSG00000059325"), split.by = "FactorValue.time.", group.by = "sub_cluster", 
        +         pt.size = 0, combine = FALSE, idents = "tumor Treg", col=c("red", "green", "blue"))
DE_Markers_tumor_Treg_D5vD8vD11_dotplotlabel <- c("Itgb8","Gm6637","Hopx","Tchp","Olfr1310","Traj27","Gstp3", "Rln3","Olfr1388","Trbj2-7","Trbv3","Sv2c","Nrn1","Marcks","Trav6-7-dv9","Stk26","Olfr1434","Ifitm2","D130007C19Rik","Clspn")
DotPlot(Seurat_object, features = rownames(head(DE_Markers_tumor_Treg_D5vD8vD11, n=20)), cols = c("RdYlBu"), dot.scale = 8, group.by = "FactorValue.time.", idents = "tumor Treg") + RotatedAxis() + scale_x_discrete(labels= DE_Markers_tumor_Treg_D5vD8vD11_dotplotlabel)


#DE markers between D5 vs. D8 tumor Treg only
Seurat_object <- SetIdent(object = Seurat_object, value = "sub_cluster")
DE_Markers_tumor_Treg_D5vD8 <- FindMarkers(Seurat_object, ident.1 = "5", ident.2 = "8", group.by ="FactorValue.time.", subset.ident = "tumor Treg", verbose = FALSE)
head(DE_Markers_tumor_Treg_D5vD8, n = 20)
write.table(DE_Markers_tumor_Treg_D5vD8,"DE_Markers_tumor_Treg_D5vD8.txt",sep="\t",row.names=TRUE)

VlnPlot(Seurat_object, features = c("ENSMUSG00000064356"), split.by = "FactorValue.time.", group.by = "sub_cluster", 
        +         pt.size = 0, combine = FALSE, idents = "tumor Treg", col=c("red", "green", "blue"))
VlnPlot(Seurat_object, features = c("ENSMUSG00000083365"), split.by = "FactorValue.time.", group.by = "sub_cluster", 
        +         pt.size = 0, combine = FALSE, idents = "tumor Treg", col=c("red", "green", "blue"))
DotPlot(Seurat_object, features = rownames(head(DE_Markers_tumor_Treg_D5vD8, n=20)), cols = c("RdYlBu"), dot.scale = 6, group.by = "FactorValue.time.", idents = "tumor Treg", cluster.idents = TRUE) + RotatedAxis()
DoHeatmap(Seurat_object, features = rownames(head(DE_Markers_tumor_Treg_D5vD8, n=20)), cells = tumor_Treg)


#DE markers between D5 vs. D11 tumor Treg only
Seurat_object <- SetIdent(object = Seurat_object, value = "sub_cluster")
DE_Markers_tumor_Treg_D5vD11 <- FindMarkers(Seurat_object, ident.1 = "5", ident.2 = "11", group.by ="FactorValue.time.", subset.ident = "tumor Treg", verbose = FALSE)
head(DE_Markers_tumor_Treg_D5vD11, n = 20)
write.table(DE_Markers_tumor_Treg_D5vD11,"DE_Markers_tumor_Treg_D5vD11.txt",sep="\t",row.names=TRUE)

VlnPlot(Seurat_object, features = c("ENSMUSG00000064356"), split.by = "FactorValue.time.", group.by = "sub_cluster", 
        +         pt.size = 0, combine = FALSE, idents = "tumor Treg", col=c("red", "green", "blue"))
VlnPlot(Seurat_object, features = c("ENSMUSG00000083365"), split.by = "FactorValue.time.", group.by = "sub_cluster", 
        +         pt.size = 0, combine = FALSE, idents = "tumor Treg", col=c("red", "green", "blue"))
DotPlot(Seurat_object, features = rownames(head(DE_Markers_tumor_Treg_D5vD11, n=20)), cols = c("RdYlBu"), dot.scale = 6, group.by = "FactorValue.time.", idents = "tumor Treg", cluster.idents = TRUE) + RotatedAxis()



#DE markers between D8 vs. D11 tumor Treg only (**no significant result)
Seurat_object <- SetIdent(object = Seurat_object, value = "sub_cluster")
DE_Markers_tumor_Treg_D8vD11 <- FindMarkers(Seurat_object, ident.1 = "8", ident.2 = "11", group.by ="FactorValue.time.", subset.ident = "tumor Treg", verbose = FALSE)
head(DE_Markers_tumor_Treg_D8vD11, n = 20)
write.table(DE_Markers_tumor_Treg_D8vD11,"DE_Markers_tumor_Treg_D8vD11.txt",sep="\t",row.names=TRUE)

#FOUND NO SIGNIFICANTLY DE MARKER HERE
#VlnPlot(Seurat_object, features = c("ENSMUSG00000023393"), split.by = "FactorValue.time.", group.by = "sub_cluster", 
#        +         pt.size = 0, combine = FALSE, col=c("red", "green", "blue"))
#DotPlot(Seurat_object, features = rownames(head(DE_Markers_tumor_Treg_D8vD11, n=20)), cols = c("RdYlBu"), dot.scale = 6, group.by = "FactorValue.time.", idents = "tumor Treg", cluster.idents = TRUE) + RotatedAxis()







saveRDS(Seurat_object, file = "D5D8D11_LN_Tumour_NEW.rds")




#Seurat_object Visualisation
DimPlot(Seurat_object, reduction = "umap", label = TRUE, split.by = "sel.k")
FeaturePlot(object = Seurat_object, features = "Il10", cols = c("grey", "blue", "red"), reduction = "umap", pt.size = 1, order = TRUE, label = FALSE)
FeaturePlot(object = Treg.combined, features = "Ebi3", cols = c("grey", "blue", "red"), reduction = "umap", pt.size = 1, order = TRUE, label = FALSE)
VlnPlot(Treg.combined, features = c("Il10","Ebi3",'Notch2','Jag1','Jag2','Notch1'), pt.size=0.2, ncol = 4)
DoHeatmap(Treg.combined, features = c('Il10', 'Ebi3','Notch2','Notch1','Jag1','Jag2'))

DoHeatmap(Treg.combined, features = c("Notch1", "Il10", "Jag1", "Ebi3", "Notch2", "Jag2", "Wnt4", "Acaa1a", "Acaa1b", "Acadl", "Acadm", "Acox1", "Acox2", "Acox3", "Acsbg1", "Acsbg2", "Acsbg3", "Acsl1", "Acsl3", "Acsl4", "Acsl5", "Acsl6", "Adipoq", "Angptl4", "Apoa1", "Apoa2", "Apoa5", "Apoc3", "Aqp7", "Cd36", "Cpt1a", "Cpt1b", "Cpt1c", "Cpt2", "Cyp27a1", "Cyp4a10", "Cyp4a12a", "Cyp4a12b", "Cyp4a14", "Cyp4a29", "Cyp4a30b", "Cyp4a31", "Cyp4a32", "Cyp7a1", "Cyp8b1", "Dbi", "Dbil5", "Ehhadh", "Fabp1", "Fabp2", "Fabp3", "Fabp4", "Fabp5", "Fabp6", "Fabp7", "Fads2", "Gk", "Gk2", "Gykl1", "Hmgcs1", "Hmgcs2", "Ilk", "Lpl", "Me1", "Me3", "Mmp1a", "Mmp1b", "Nr1h3", "Olr1", "Pck1", "Pck2", "Pdpk1", "Plin1", "Plin2", "Plin4", "Plin5", "Pltp", "Ppara", "Ppard", "Pparg", "Rxra", "Rxrb", "Rxrg", "Scd1", "Scd2", "Scd3", "Scd4", "Scp2", "Slc27a1", "Slc27a2", "Slc27a4", "Slc27a5", "Slc27a6", "Sorbs1", "Ubc", "Ucp1", "4930438A08Rik", "Adh1", "Adh4", "Adh5", "Adh7", "Aldh3a1", "Aldh3b1", "Aldh3b2", "Aldh3b3", "Aoc2", "Aoc3", "Aox1", "Aox2", "Aox3", "Aox4", "Comt", "Dbh", "Dct", "Ddc", "Fah", "Fahd1", "Got1", "Got1l1", "Got2", "Gstz1", "Hgd", "Hpd", "Il4i1", "Il4i1b", "Lao1", "Maoa", "Maob", "Mif", "Pnmt", "Tat", "Th", "Tomt", "Tpo", "Tyr", "Tyrp1","Actb", "Actg1", "Aifm1", "Akt1", "Akt2", "Akt3", "Apaf1", "Atf4", "Atm", "Bad", "Bak1", "Bax", "Bbc3", "Bcl2", "Bcl2a1a", "Bcl2a1b", "Bcl2a1c", "Bcl2a1d", "Bcl2l1", "Bcl2l11", "Bid", "Birc2", "Birc3", "Birc5", "Capn1", "Capn2", "Casp12", "Casp2", "Casp3", "Casp6", "Casp7", "Casp8", "Casp9", "Cflar", "Chuk", "Csf2rb", "Csf2rb2", "Ctsb", "Ctsc", "Ctsd", "Ctsf", "Ctsh", "Ctsk", "Ctsl", "Ctso", "Ctss", "Ctsw", "Ctsz", "Cycs", "Cyct", "Dab2ip", "Daxx", "Ddit3", "Dffa", "Dffb", "Diablo", "Eif2ak3", "Eif2s1", "Endog", "Ern1", "Fadd", "Fas", "Fasl", "Fos", "Gadd45a", "Gadd45b", "Gadd45g", "Gzmb", "Hras", "Hrk", "Htra2", "Ikbkb", "Ikbkg", "Il3", "Il3ra", "Itpr1", "Itpr2", "Itpr3", "Jun", "Kras", "Lmna", "Lmnb1", "Lmnb2", "Map2k1", "Map2k2", "Map3k14", "Map3k5", "Mapk1", "Mapk10", "Mapk3", "Mapk8", "Mapk9", "Mcl1", "Nfkb1", "Nfkbia", "Ngf", "Nras", "Ntrk1", "Parp1", "Parp2", "Parp3", "Parp4", "Pdpk1", "Pidd1", "Pik3ca", "Pik3cb", "Pik3cd", "Pik3r1", "Pik3r2", "Pik3r3", "Pmaip1", "Prf1", "Ptpn13","Raf1", "Rela", "Ripk1", "Septin4", "Spta1", "Sptan1", "Tnf", "Tnfrsf10b", "Tnfrsf1a", "Tnfsf10", "Tradd", "Traf1", "Traf2", "Trp53", "Tuba1a", "Tuba1b", "Tuba1c", "Tuba3a", "Tuba3b", "Tuba4a", "Tuba8", "Tubal3", "Xiap"))



DefaultAssay(Treg.combined) <- "RNA"
#find conserved markers between conditions (ctrl v.s. tumour) for a cluster
Cluster0.conserved.markers <- FindConservedMarkers(Treg.combined, ident.1 = 0, grouping.var = "condition")
Cluster1.conserved.markers <- FindConservedMarkers(Treg.combined, ident.1 = 1, grouping.var = "condition", verbose = FALSE)
Cluster2.conserved.markers <- FindConservedMarkers(Treg.combined, ident.1 = 2, grouping.var = "condition", verbose = FALSE)
Cluster3.conserved.markers <- FindConservedMarkers(Treg.combined, ident.1 = 3, grouping.var = "condition", verbose = FALSE)
Cluster4.conserved.markers <- FindConservedMarkers(Treg.combined, ident.1 = 4, grouping.var = "condition", verbose = FALSE)
Cluster5.conserved.markers <- FindConservedMarkers(Treg.combined, ident.1 = 5, grouping.var = "condition", verbose = FALSE)
Cluster6.conserved.markers <- FindConservedMarkers(Treg.combined, ident.1 = 6, grouping.var = "condition", verbose = FALSE)
head(Cluster0.conserved.markers)
head(Cluster1.conserved.markers)
head(Cluster2.conserved.markers)
head(Cluster3.conserved.markers)
head(Cluster4.conserved.markers)
head(Cluster5.conserved.markers)
head(Cluster6.conserved.markers)

#find enriched markers for every cluster compared to all remaining cells (cluster v.s. all other clusters)
cluster.markers <- FindAllMarkers(Treg.combined, only.pos = TRUE)
cluster.markers$pct.diff <- (cluster.markers$pct.1 - cluster.markers$pct.2)
#identify the top 20 markers
top20_cluster_markers <- cluster.markers %>% group_by(cluster) %>% top_n(n=20, wt=avg_log2FC)
write.csv(top20_cluster_markers, file='top20.cluster.markers.csv')


#identify the top 100 markers
top100_cluster_markers <- cluster.markers %>% group_by(cluster) %>% top_n(n=100, wt=avg_log2FC)
#Subsetting top 20 markers with adjusted p values lower than .05
top100pval <- subset(top100_cluster_markers, rowSums(top100_cluster_markers[5] < 0.05) > 0)

#ClusterProfiler
library(clusterProfiler)
library(org.Mm.eg.db)
library(AnnotationHub)

df <- top100pval[,7:6]
dfsample <- split(df$gene,df$cluster)
length(dfsample)
#The output of length(dfsample) returns how many clusters you have

dfsample$`0` = bitr(dfsample$`0`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
dfsample$`1` = bitr(dfsample$`1`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
dfsample$`2` = bitr(dfsample$`2`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
dfsample$`3` = bitr(dfsample$`3`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
dfsample$`4` = bitr(dfsample$`4`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
dfsample$`5` = bitr(dfsample$`5`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
dfsample$`6` = bitr(dfsample$`6`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")



#do the same here, a line like below for each cluster
genelist <- list("0" = dfsample$`0`$ENTREZID, 
                 "1" = dfsample$`1`$ENTREZID,
                 "2" = dfsample$`2`$ENTREZID,
                 "3" = dfsample$`3`$ENTREZID,
                 "4" = dfsample$`4`$ENTREZID,
                 "5" = dfsample$`5`$ENTREZID,
                 "6" = dfsample$`6`$ENTREZID,
                 
)

GOclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichGO", OrgDb = "org.Mm.eg.db")
dotplot(GOclusterplot)

KEGGclusterplot <- compareCluster(geneCluster = genelist, organism = "mmu", fun = "enrichKEGG", pvalueCutoff = 0.05)
dotplot(KEGGclusterplot)




theme_set(theme_cowplot())
treg.0 <- subset(Treg.combined, idents = 0)
treg.1 <- subset(Treg.combined, idents = 1)
treg.2 <- subset(Treg.combined, idents = 2)
treg.3 <- subset(Treg.combined, idents = 3)
treg.4 <- subset(Treg.combined, idents = 4)
treg.5 <- subset(Treg.combined, idents = 5)
treg.6 <- subset(Treg.combined, idents = 6)

Idents(treg.1) <- "condition"
avg.treg.1 <- as.Seurat_object.frame(log1p(AverageExpression(treg.1, verbose = FALSE)$RNA))
avg.treg.1$gene <- rownames(avg.treg.1)

genes.to.label = c("Notch1", "Il10", "Jag1", "Ebi3", "Rsp27", "Notch2", "Jag2", "Wnt4")
p1.1 <- ggplot(avg.treg.1, aes(CTRL, TUMOUR)) + geom_point() + ggtitle("Treg.1")
p1.1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE)
p1.1

Treg.combined$celltype.condition <- paste(Idents(Treg.combined), Treg.combined$condition, sep = "_")
Treg.combined$celltype <- Idents(Treg.combined)
Idents(Treg.combined) <- "celltype.condition"
#Identify differentially-expressed genes in cluster x between conditions
Differential_Markers.cluster0 <- FindMarkers(Treg.combined, ident.1 = "0_CTRL", ident.2 = "0_TUMOUR", verbose = FALSE)
head(Differential_Markers.cluster0, n = 15)
write.csv(Differential_Markers.cluster0, file="cluster0.DE.markers.csv")
Differential_Markers.cluster1 <- FindMarkers(Treg.combined, ident.1 = "1_CTRL", ident.2 = "1_TUMOUR", verbose = FALSE)
head(Differential_Markers.cluster1, n = 15)
write.csv(Differential_Markers.cluster1, file="cluster1.DE.markers.csv")
Differential_Markers.cluster2 <- FindMarkers(Treg.combined, ident.1 = "2_CTRL", ident.2 = "2_TUMOUR", verbose = FALSE)
head(Differential_Markers.cluster2, n = 15)
write.csv(Differential_Markers.cluster2, file="cluster2.DE.markers.csv")
Differential_Markers.cluster3 <- FindMarkers(Treg.combined, ident.1 = "3_CTRL", ident.2 = "3_TUMOUR", verbose = FALSE)
head(Differential_Markers.cluster3, n = 15)
write.csv(Differential_Markers.cluster3, file="cluster3.DE.markers.csv")
Differential_Markers.cluster4 <- FindMarkers(Treg.combined, ident.1 = "4_CTRL", ident.2 = "4_TUMOUR", verbose = FALSE)
head(Differential_Markers.cluster4, n = 15)
write.csv(Differential_Markers.cluster4, file="cluster4.DE.markers.csv")
Differential_Markers.cluster6 <- FindMarkers(Treg.combined, ident.1 = "6_CTRL", ident.2 = "6_TUMOUR", verbose = FALSE)
head(Differential_Markers.cluster6, n = 15)
write.csv(Differential_Markers.cluster6, file="cluster6.DE.markers.csv")

markers.to.plot <- c("Notch1", "Il10", "Jag1", "Ebi3", "Notch2", "Jag2", "Wnt4")
#Visualize how our markers of interest express in different clusters between conditions
DotPlot(Treg.combined, features = rev(markers.to.plot), cols = c("blue", "red"), dot.scale = 8) + RotatedAxis()
DimPlot(Treg.combined, label = TRUE)


#Combine all Treg subpopulations/clusters in one
#current.cluster.ids <- c(0, 1, 2, 3, 4, 5, 6)

#new.cluster.ids <- c('Tregs', 'Tregs',
#                     'Tregs', 'Tregs',
#                     'Tregs', 'Tregs', 'Tregs')

#Treg.combined@active.ident <- plyr::mapvalues(x = Treg.combined@active.ident, from = current.cluster.ids, to = new.cluster.ids)
#DimPlot(object = Treg.combined, reduction = "umap", label=TRUE, pt.size = 0.5)


Idents(Treg.combined) <- "celltype"
Treg.combined <- RenameIdents(Treg.combined, `0` = "Tregs", `1` = "Tregs", `2` = "Tregs", 
                              `3` = "Tregs", `4` = "Tregs", `5` = "Tregs", `6` = "Tregs")
DimPlot(Treg.combined, label = TRUE)

tregs <- subset(Treg.combined, idents = "Tregs")
Idents(tregs) <- "condition"
avg.tregs <- as.Seurat_object.frame(log1p(AverageExpression(tregs, verbose = FALSE)$RNA))
avg.tregs$gene <- rownames(avg.tregs)

genes.to.label = c("Notch1", "Il10", "Jag1", "Ebi3", "Rsp27", "Notch2", "Jag2", "Wnt4")
p2.1 <- ggplot(avg.tregs, aes(CTRL, TUMOUR)) + geom_point() + ggtitle("Tregs")
p2.1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE)
p2.1


Treg.combined <- subset(Treg.combined, idents = "Tregs")
Idents(Treg.combined) <- "condition"
#Identify differentially-expressed genes for all clusters combined between conditions
Differential_Markers.all.clusters <- FindMarkers(tregs, ident.1 = "CTRL", ident.2 = "TUMOUR", verbose = FALSE)
head(Differential_Markers.all.clusters, n = 100)


#Visualize how our markers of interest express between conditions (averaged across all cluster)
ppar.markers.to.plot <- c("Notch1", "Il10", "Jag1", "Ebi3", "Notch2", "Jag2", "Wnt4", "Acaa1a", "Acaa1b", "Acadl", "Acadm", "Acox1", "Acox2", "Acox3", "Acsbg1", "Acsbg2", "Acsbg3", "Acsl1", "Acsl3", "Acsl4", "Acsl5", "Acsl6", "Adipoq", "Angptl4", "Apoa1", "Apoa2", "Apoa5", "Apoc3", "Aqp7", "Cd36", "Cpt1a", "Cpt1b", "Cpt1c", "Cpt2", "Cyp27a1", "Cyp4a10", "Cyp4a12a", "Cyp4a12b", "Cyp4a14", "Cyp4a29", "Cyp4a30b", "Cyp4a31", "Cyp4a32", "Cyp7a1", "Cyp8b1", "Dbi", "Dbil5", "Ehhadh", "Fabp1", "Fabp2", "Fabp3", "Fabp4", "Fabp5", "Fabp6", "Fabp7", "Fads2", "Gk", "Gk2", "Gykl1", "Hmgcs1", 
                          "Hmgcs2", "Ilk", "Lpl", "Me1", "Me3", "Mmp1a", "Mmp1b", "Nr1h3", "Olr1", "Pck1", "Pck2", "Pdpk1", "Plin1", "Plin2", "Plin4", "Plin5", "Pltp", "Ppara", "Ppard", "Pparg", "Rxra", "Rxrb", "Rxrg", "Scd1", "Scd2", "Scd3", "Scd4", "Scp2", "Slc27a1", "Slc27a2", "Slc27a4", "Slc27a5", "Slc27a6", "Sorbs1", "Ubc", "Ucp1")
DotPlot(Treg.combined, features = rev(ppar.markers.to.plot), cols = c("blue", "red"), dot.scale = 8) + RotatedAxis()

plots <- VlnPlot(Treg.combined, features = c("Ctla4","Rpl37a","Tagln2","Vim","Ccl5","Rpl35"), split.by = "condition", group.by = "celltype", 
                 pt.size = 0, combine = FALSE)
plots

tyrosine.markers.to.plot <- c("4930438A08Rik", "Adh1", "Adh4", "Adh5", "Adh7", "Aldh3a1", "Aldh3b1", "Aldh3b2", "Aldh3b3", "Aoc2", "Aoc3", "Aox1", "Aox2", "Aox3", "Aox4", "Comt", "Dbh", "Dct", "Ddc", "Fah", "Fahd1", "Got1", "Got1l1", "Got2", "Gstz1", "Hgd", "Hpd", "Il4i1", "Il4i1b", "Lao1", "Maoa", "Maob", "Mif", "Pnmt", "Tat", "Th", "Tomt", "Tpo", "Tyr", "Tyrp1")

DotPlot(Treg.combined, features = rev(tyrosine.markers.to.plot), cols = c("blue", "red"), dot.scale = 8) + RotatedAxis()


noquote(toString(apoptosis_features))

apoptosis.markers.to.plot <- c("Actb", "Actg1", "Aifm1", "Akt1", "Akt2", "Akt3", "Apaf1", "Atf4", "Atm", "Bad", "Bak1", "Bax", "Bbc3", "Bcl2", "Bcl2a1a", "Bcl2a1b", "Bcl2a1c", "Bcl2a1d", "Bcl2l1", "Bcl2l11", "Bid", "Birc2", "Birc3", "Birc5", "Capn1", "Capn2", "Casp12", "Casp2", "Casp3", "Casp6", "Casp7", "Casp8", "Casp9", "Cflar", "Chuk", "Csf2rb", "Csf2rb2", "Ctsb", "Ctsc", "Ctsd", "Ctsf", "Ctsh", "Ctsk", "Ctsl", "Ctso", "Ctss", "Ctsw", "Ctsz", "Cycs", "Cyct", "Dab2ip", "Daxx", "Ddit3", "Dffa", "Dffb", "Diablo", "Eif2ak3", "Eif2s1", "Endog", "Ern1", "Fadd", "Fas", "Fasl", "Fos", "Gadd45a", "Gadd45b", "Gadd45g", "Gzmb", "Hras", "Hrk", "Htra2", "Ikbkb", "Ikbkg", "Il3", "Il3ra", "Itpr1", "Itpr2", "Itpr3", "Jun", "Kras", "Lmna", "Lmnb1", "Lmnb2", "Map2k1", "Map2k2", "Map3k14", "Map3k5", "Mapk1", "Mapk10", "Mapk3", "Mapk8", "Mapk9", "Mcl1", "Nfkb1", "Nfkbia", "Ngf", "Nras", "Ntrk1", "Parp1", "Parp2", "Parp3", "Parp4", "Pdpk1", "Pidd1", "Pik3ca", "Pik3cb", "Pik3cd", "Pik3r1", "Pik3r2", "Pik3r3", "Pmaip1", "Prf1", "Ptpn13","Raf1", "Rela", "Ripk1", "Septin4", "Spta1", "Sptan1", "Tnf", "Tnfrsf10b", "Tnfrsf1a", "Tnfsf10", "Tradd", "Traf1", "Traf2", "Trp53", "Tuba1a", "Tuba1b", "Tuba1c", "Tuba3a", "Tuba3b", "Tuba4a", "Tuba8", "Tubal3", "Xiap")

DotPlot(Treg.combined, features =rev(apoptosis.markers.to.plot), cols = c("blue", "red"), dot.scale = 8) + RotatedAxis()


plots2 <- VlnPlot(Treg.combined, features = c("Mif","Got2","Got1","Adh5","Actg1","Actb", "Ctsl", "Ctsk", "Ctsh", "Ctsf"), split.by = "condition", group.by = "celltype", 
                  pt.size = 0, combine = FALSE)
plots2
