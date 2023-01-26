.libPaths()
assign(".lib.loc", "E:/AG Scholz/R libraries2", envir = environment(.libPaths))
.libPaths()


library(Seurat)
library(patchwork)
library(dplyr)

#sessionInfo()
Matrix_1<- ReadMtx(  mtx = "T:/Dokumente/Hausarbeit/R_Datasets/Ovar/GSM6720925_ovarian_cancer2251.mtx.gz", 
                     cells = "T:/Dokumente/Hausarbeit/R_Datasets/Ovar/GSM6720925_ovarian_cancer_barcodes2251.tsv.gz",
                     features = "T:/Dokumente/Hausarbeit/R_Datasets/Ovar/GSM6720925_ovarian_cancer_genes2251.tsv.gz"
)

Matrix_2<- ReadMtx(  mtx = "T:/Dokumente/Hausarbeit/R_Datasets/Ovar/GSM6720926_2ovarian_cancer2267.mtx.gz", 
                     cells = "T:/Dokumente/Hausarbeit/R_Datasets/Ovar/GSM6720926_2ovarian_cancer_barcodes2267.tsv.gz",
                     features = "T:/Dokumente/Hausarbeit/R_Datasets/Ovar/GSM6720926_2ovarian_cancer_features2267.tsv.gz"
)
Matrix_3<- ReadMtx(  mtx = "T:/Dokumente/Hausarbeit/R_Datasets/Ovar/GSM6720927_3ovarian_cancer_2283.mtx.gz", 
                     cells = "T:/Dokumente/Hausarbeit/R_Datasets/Ovar/GSM6720927_3ovarian_cancer_barcodes_2283.tsv.gz",
                     features = "T:/Dokumente/Hausarbeit/R_Datasets/Ovar/GSM6720927_3ovarian_cancer_features_2283.tsv.gz"
)

Matrix_4<- ReadMtx(  mtx = "T:/Dokumente/Hausarbeit/R_Datasets/Ovar/GSM6720928_4ovarian_cancer_2293.mtx.gz", 
                     cells = "T:/Dokumente/Hausarbeit/R_Datasets/Ovar/GSM6720928_4ovarian_cancer_barcodes_2293.tsv.gz",
                     features = "T:/Dokumente/Hausarbeit/R_Datasets/Ovar/GSM6720928_4ovarian_cancer_features_2293.tsv.gz"
)

Matrix_5<- ReadMtx(
  mtx = "T:/Dokumente/Hausarbeit/R_Datasets/Ovar/GSM6720929_5ovarian_cancer_2380.mtx.gz", 
  cells = "T:/Dokumente/Hausarbeit/R_Datasets/Ovar/GSM6720929_5ovarian_cancer_barcodes_2380.tsv.gz",
  features = "T:/Dokumente/Hausarbeit/R_Datasets/Ovar/GSM6720929_5ovarian_cancer_features_2380.tsv.gz"
)
Matrix_6<- ReadMtx(  mtx = "T:/Dokumente/Hausarbeit/R_Datasets/Ovar/GSM6720930_6ovarian_cancer_2428.mtx.gz", 
                     cells = "T:/Dokumente/Hausarbeit/R_Datasets/Ovar/GSM6720930_6ovarian_cancer_barcodes_2428.tsv.gz",
                     features = "T:/Dokumente/Hausarbeit/R_Datasets/Ovar/GSM6720930_6ovarian_cancer_features_2428.tsv.gz"
)

Matrix_7<- ReadMtx(
  mtx = "T:/Dokumente/Hausarbeit/R_Datasets/Ovar/GSM6720931_7ovarian_cancer_matrix_2467.mtx.gz", 
  cells = "T:/Dokumente/Hausarbeit/R_Datasets/Ovar/GSM6720931_7ovarian_cancer_barcodes_2467.tsv.gz",
  features = "T:/Dokumente/Hausarbeit/R_Datasets/Ovar/GSM6720931_7ovarian_cancer_features_2467.tsv.gz"
)

Matrix_8<- ReadMtx(
  mtx = "T:/Dokumente/Hausarbeit/R_Datasets/Ovar/GSM6720932_8ovarian_cancer_matrix_2497.mtx.gz", 
  cells = "T:/Dokumente/Hausarbeit/R_Datasets/Ovar/GSM6720932_8ovarian_cancer_barcodes_2497.tsv.gz",
  features = "T:/Dokumente/Hausarbeit/R_Datasets/Ovar/GSM6720932_8ovarian_cancer_features_2497.tsv.gz"
)
Seurat_1<-CreateSeuratObject(counts = Matrix_1)
Seurat_1
Seurat_1[["origin"]]<-"OvarCancer1"
Seurat_1[["condition"]] <- "tumor"

Seurat_2<- CreateSeuratObject(counts = Matrix_2)
Seurat_2
Seurat_2[["origin"]]<-"OvarCancer2"
Seurat_2[["condition"]]<- "tumor"

Seurat_3<- CreateSeuratObject(counts = Matrix_3)
Seurat_3
Seurat_3[["origin"]]<-"OvarCancer3"
Seurat_3[["condition"]]<- "tumor"

Seurat_4<- CreateSeuratObject(counts = Matrix_4)
Seurat_4
Seurat_4[["origin"]]<- "OvarCancer4"
Seurat_4[["condition"]]<- "tumor"

Seurat_5<- CreateSeuratObject(counts = Matrix_5)
Seurat_5
Seurat_5[["origin"]] <- "OvarCancer5"
Seurat_5[["condition"]]<- "tumor"

Seurat_6<- CreateSeuratObject(counts = Matrix_6)
Seurat_6
Seurat_6[["origin"]]<-"OvarCancer6"
Seurat_6[["condition"]]<- "tumor"  

Seurat_7<- CreateSeuratObject(counts = Matrix_7)
Seurat_7
Seurat_7[["origin"]] <- "OvarCancer7"
Seurat_7[["condition"]]<- "tumor"

Seurat_8<- CreateSeuratObject(counts = Matrix_8)
Seurat_8
Seurat_8[["origin"]] <- "OvarCancer8"
Seurat_8[["condition"]]<- "tumor"



list_1<-list(Seurat_1,Seurat_2,Seurat_3, Seurat_4, Seurat_5, Seurat_6, Seurat_7, Seurat_8)

list_1 <- lapply(X = list_1, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures (x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = list_1)
list_1 <- lapply(X= list_1, FUN = function(x){
  x <- ScaleData (x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})


data.anchors <- FindIntegrationAnchors(object.list = list_1, anchor.features = features, reduction = "rpca")
data.combined <- IntegrateData(anchorset = data.anchors)

DefaultAssay(data.combined) <- "integrated"

saveRDS(data.combined, file ="C:/Users/Landstoh/Documents/seurat_ovary_cancer_merged.RDS", compress = FALSE)

data.combined <- readRDS(file="C:/Users/Landstoh/Documents/seurat_ovary_cancer_merged.RDS")

data.combined <- ScaleData(data.combined, verbose = TRUE)
data.combined <- RunPCA(data.combined, npcs = 15, verbose = TRUE)

DimPlot(data.combined, reduction = "pca")
ElbowPlot(data.combined)
data.combined <- RunUMAP(data.combined, reduction = "pca", dims = 1:15)
data.combined <- FindNeighbors(data.combined, reduction = "pca", dims = 1:15)
data.combined <- FindClusters(data.combined, resolution = 0.5)

p0 <- DimPlot(data.combined, reduction = "umap", label = TRUE, 
              repel = TRUE)
p0

p1 <- DimPlot(data.combined, reduction = "umap",label = TRUE,
              repel = TRUE, split.by = "origin")
p1

p2 <- DimPlot(data.combined, reduction = "umap", label = TRUE, 
              repel = TRUE, group.by = "origin")
p2
FeaturePlot(data.combined, features= c("WT1"),min.cutoff = "q10", max.cutoff = "q90",label = TRUE,
            repel = TRUE)

RidgePlot(data.combined, feature = c("WT1", "TP53"))
VlnPlot(data.combined, features = "WT1", pt.size = 0)

FeaturePlot(data.combined, features = c("WT1"))

marker_gene(data.combined)
marker_genes <- FindAllMarkers(data.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers.all <- tibble(marker_genes %>%
                        group_by(cluster) %>%
                        slice_max(n = 2, order_by = avg_log2FC))

saveRDS(markers.all, file = "C:/Users/Landstoh/Documents/integrated_tumor_data_all_markers", compress =FALSE )

new_cluster_names <- c("T Cells", "NK", "Fibroblasts", "Plasma cells", "Smooth muscle", "Granulocytes", "B cells", "Neurons", "Macrophages", "Endothelial cells", "Endothelial cells", "Granulocytes", "Schwann cells / Neurons", "Keratinocytes",
                       "neuron/endothelia", "plasma cells", "fibroblasts/myoepithelial cells","Erythroid cells", "activated B cells", "endometrial stroma cells", "Myeloid cells", "dendritic cells", "Immune cells", "granulocytes")
names(new_cluster_names) <- levels(data.combined)
data.combined <- RenameIdents(data.combined, new_cluster_names)
DimPlot(data.combined, reduction = "umap", label = TRUE, pt.size = 0.5)
VlnPlot(data.combined, features = "WT1", pt.size = 0)
FeaturePlot(data.combined, features= c("WT1"),min.cutoff = "q10", max.cutoff = "q90",label = TRUE,
            repel = TRUE)
markers <- readRDS(file = "C:/Users/Landstoh/Documents/integrated_tumor_data_all_markers")
Idents(data.combined) <- factor(Idents(data.combined), levels = c("T Cells", "NK", "Fibroblasts", "Plasma cells", "Smooth muscle", "Granulocytes", "B cells", "Neurons", "Macrophages", "Endothelial cells", "Endothelial cells", "Granulocytes", "Schwann cells / Neurons", "Keratinocytes",
                                                                  "neuron/endothelia", "plasma cells", "fibroblasts/myoepithelial cells","Erythroid cells", "activated B cells", "endometrial stroma cells", "Myeloid cells", "dendritic cells", "Immune cells", "granulocytes"))
markers_clusters_plot <- c("SLC38A1","IL7R", "GNLY", "CCL5", "COL1A1", "LUM", "IGHM", "IGKC","ACTA2", "RGS5", "WFDC2", "RHEX","CD19", "KCND2", "PRKG1", "C1QA", "SPP1", "ISG15", "IFI27", "CCL21", "CCL14", "AC024230.1", "FMN1", 
                           "NLRP3", "S100A9", "S100A8", "ZNF385D", "ADAMTS9", "IGLC2", "IGLC3", "CFD", "MGP", "HMGB2", "CENPF", "BANK1", "MS4A1", "MT1E", "ALDH1A3", "HLA-DPA1", "HLA-DPB1", "IRF4", "GZMB", "AL136456.1", "AREG", "TPSAB1", "TPSB2") 
DotPlot(data.combined, features = markers_clusters_plot, dot.scale = 8)
install.packages("hdf5r")

Matrix_9 <- Read10X_h5("T:/Dokumente/Hausarbeit/R_Datasets/Ovar/Ovar physiologisch/GSM3319041_sample_3-14_filtered_gene_bc_matrices_h5.h5")
Seurat_9<-CreateSeuratObject(counts = Matrix_9)
Seurat_9
Seurat_9[["origin"]]<-"OvarPhys1"
Seurat_9[["condition"]] <- "healthy tissue"

Matrix_10 <- Read10X_h5("T:/Dokumente/Hausarbeit/R_Datasets/Ovar/Ovar physiologisch/GSM3319042_sample_3-15_filtered_gene_bc_matrices_h5.h5")
Seurat_10<-CreateSeuratObject(counts = Matrix_10)
Seurat_10
Seurat_10[["origin"]]<-"OvarPhys2"
Seurat_10[["condition"]] <- "healthy tissue"

Matrix_11 <- Read10X_h5("T:/Dokumente/Hausarbeit/R_Datasets/Ovar/Ovar physiologisch/GSM3319046_sample_3-5_filtered_gene_bc_matrices_h5.h5")
Seurat_11<-CreateSeuratObject(counts = Matrix_11)
Seurat_11
Seurat_11[["origin"]]<-"OvarPhys3"
Seurat_11[["condition"]] <- "healthy tissue"

Matrix_12 <- Read10X_h5("T:/Dokumente/Hausarbeit/R_Datasets/Ovar/Ovar physiologisch/GSM3557968_sample5_B2_i10E_filtered_gene_bc_matrices_h5.h5")
Seurat_12<-CreateSeuratObject(counts = Matrix_12)
Seurat_12
Seurat_12[["origin"]]<-"OvarPhys4"
Seurat_12[["condition"]] <- "healthy tissue"

list_2<-list(Seurat_1,Seurat_2,Seurat_3, Seurat_4, Seurat_5, Seurat_6, Seurat_7, Seurat_8, Seurat_9, Seurat_10, Seurat_11, Seurat_12)

saveRDS(list_2, file = "C:/Users/Landstoh/Documents/list_tumor_and_phys", compress =FALSE )

list_2 <- lapply(X = list_2, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures (x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = list_2)
list_2 <- lapply(X= list_2, FUN = function(x){
  x <- ScaleData (x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})


data.anchors <- FindIntegrationAnchors(object.list = list_2, anchor.features = features, reduction = "rpca")
data.all.merged <- IntegrateData(anchorset = data.anchors)
saveRDS(data.anchors, file = "C:/Users/Landstoh/Documents/anchors_merged", compress = FALSE)


DefaultAssay(data.all.merged) <- "integrated"

saveRDS(data.all.merged, file ="C:/Users/Landstoh/Documents/tumor_phys_merged.RDS", compress = FALSE)
data.all.merged <- readRDS(file = "C:/Users/Landstoh/Documents/tumor_phys_merged.RDS")

data.all.merged <- ScaleData(data.all.merged, verbose = TRUE)
data.all.merged <- RunPCA(data.all.merged, npcs = 15, verbose = TRUE)

DimPlot(data.all.merged, reduction = "pca")
ElbowPlot(data.all.merged)
data.all.merged <- RunUMAP(data.all.merged, reduction = "pca", dims = 1:15)
data.all.merged <- FindNeighbors(data.all.merged, reduction = "pca", dims = 1:15)
data.all.merged <- FindClusters(data.all.merged, resolution = 0.5)

DefaultAssay(data.all.merged) <- "RNA"
p0 <- DimPlot(data.all.merged, reduction = "umap", label = TRUE, 
              repel = TRUE)
p0

p1 <- DimPlot(data.all.merged, reduction = "umap",
              repel = TRUE, group.by = "condition")
p1

p2 <- DimPlot(data.all.merged, reduction = "umap", label = TRUE, 
              repel = TRUE, split.by = "condition")
p2

p3 <- FeaturePlot(data.all.merged, features= c("WT1"),min.cutoff = "q10", max.cutoff = "q90",label = TRUE,
                  repel = TRUE, group.by = "condition")
p3

p4<- FeaturePlot(data.all.merged, features = c("WT1"),min.cutoff = "q10", max.cutoff = "q90",label = TRUE,
                 repel = TRUE, split.by = "condition")
p4

p5<- FeaturePlot(data.all.merged, features = c("EPCAM", "PAX8", "CD24"), min.cutoff = "q10", max.cutoff = "q90",label = TRUE,
                 repel = TRUE, split.by = "condition")

p5  #tumor cells

p6 <- FeaturePlot(data.all.merged, features = c("CD79A", "IGHG3", "IGKC"), min.cutoff = "q10", max.cutoff = "q90",label = TRUE,
                  repel = TRUE, split.by = "condition")
p6 #b cells

p7 <- FeaturePlot(data.all.merged, features = c("CD3D", "CD3E", "TRBC2"), min.cutoff = "q10", max.cutoff = "q90",label = TRUE,
                  repel = TRUE, split.by = "condition")
p7 #t cells

p8<- FeaturePlot(data.all.merged, features = c("AMH"),label = TRUE, min.cutoff = "q10", max.cutoff = "q90",
                 repel = TRUE)
p8
p9 <- DimPlot(data.all.merged, reduction = "umap", label = TRUE, 
              repel = TRUE, split.by = "condition", group.by ="origin")
p9

p10 <- FeaturePlot(data.all.merged, features = c("TP53"), label = TRUE, min.cutoff = "q10", max.cutoff = "q90", split.by = "condition")
p10


FeatureScatter(data.all.merged, feature1 = "WT1", feature2 = "EPCAM")

DotPlot


data.all.merged[["clusternumber"]] <- Idents(object = data.all.merged)
data.all.merged[["celltype"]]

new_cluster_names <- c("T Cells", "Fibroblasts", "unknown", "NK Cells", "Fibroblasts", "Epithelial Cancer", "Myeloid Cells", "Fibroblasts", "Myeloid Cells", "Endothelial Cells", "Ovar Stroma", "Epithelial Cancer", "B Cells", "Endothelial Cells",
                       "T Cells", "T Cells", "T Cells","Fibroblasts", "Fibroblasts", "B Cells", "Epithelial Cancer", "Fibroblasts", "Myeloid Cells", "Fibroblasts", "Myeloid Cells", "Myeloid Cells", "Granulocytes", "Endothelial Cells")
names(new_cluster_names) <- levels(data.all.merged)
data.all.merged <- RenameIdents(data.all.merged, new_cluster_names)
DimPlot(data.all.merged, reduction = "umap", label = TRUE, pt.size = 0.5, split.by = "condition")
VlnPlot(data.all.merged, features = "WT1", pt.size = 0, split.by = "condition")
VlnPlot(data.all.merged, features = "WT1", pt.size = 0)


data.all.merged <- readRDS(file = "C:/Users/Landstoh/Documents/tumor_phys_merged.RDS")

data.all.merged <- ScaleData(data.all.merged, verbose = TRUE)
data.all.merged <- RunPCA(data.all.merged, npcs = 15, verbose = TRUE)

DimPlot(data.all.merged, reduction = "pca")
ElbowPlot(data.all.merged)
data.all.merged <- RunUMAP(data.all.merged, reduction = "pca", dims = 1:15)
data.all.merged <- FindNeighbors(data.all.merged, reduction = "pca", dims = 1:15)
data.all.merged <- FindClusters(data.all.merged, resolution = 0.5)

DefaultAssay(data.all.merged) <- "RNA"
celltype_markers <- c("CD79A", "IGHG3", "CD3D", "TRBC2", "CD1D","CD1A", "CLDN5", "PECAM1", "COL1A1", "COL1A2", "STAR", "FOXL2", "EPCAM", "PAX8", "CD68", "LYZ", "NKG7", "KLRD1", "IL1RL1", "CPA3", "WT1")
DotPlot(data.all.merged, features = celltype_markers, cols = c("blue", "red"), dot.scale = 8) + RotatedAxis()

cluster2.markers <- FindMarkers(data.all.merged, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n= 5)

saveRDS(cluster2.markers, file ="C:/Users/Landstoh/Documents/cluster2_markers.RDS", compress = FALSE)

cluster3.markers <- FindMarkers(data.all.merged, ident.1 = 3, min.pct = 0.25)
head(cluster3.markers, n= 5)

saveRDS(cluster3.markers, file ="C:/Users/Landstoh/Documents/cluster3_markers.RDS", compress = FALSE)

cluster13.markers <- FindMarkers(data.all.merged, ident.1 = 13, min.pct = 0.25)
head(cluster13.markers, n= 5)

saveRDS(cluster13.markers, file ="C:/Users/Landstoh/Documents/cluster13_markers.RDS", compress = FALSE)

cluster15.markers <- FindMarkers(data.all.merged, ident.1 = 15, min.pct = 0.25)
head(cluster15.markers, n= 5)

saveRDS(cluster15.markers, file ="C:/Users/Landstoh/Documents/cluster15_markers.RDS", compress = FALSE)

cluster21.markers <- FindMarkers(data.all.merged, ident.1 = 21, min.pct = 0.25)
head(cluster21.markers, n= 5)

saveRDS(cluster21.markers, file ="C:/Users/Landstoh/Documents/cluster21_markers.RDS", compress = FALSE)

cluster23.markers <- FindMarkers(data.all.merged, ident.1 = 23, min.pct = 0.25)
head(cluster23.markers, n= 5)

saveRDS(cluster23.markers, file ="C:/Users/Landstoh/Documents/cluster2_markers.RDS", compress = FALSE)

cluster26.markers <- FindMarkers(data.all.merged, ident.1 = 26, min.pct = 0.25)
head(cluster26.markers, n= 5)

saveRDS(cluster26.markers, file ="C:/Users/Landstoh/Documents/cluster2_markers.RDS", compress = FALSE)

p11 <- FeaturePlot(data.all.merged, features = c("ISG15"), min.cutoff = "q10", max.cutoff = "q90", label = TRUE)
p11

p12 <- VlnPlot(data.all.merged, features = "WT1", split.by = "condition", pt.size = 0)
p12

p13 <- FeaturePlot (data.all.merged, features = c("WT1"), label = TRUE, repel = TRUE, split.by = "condition", min.cutoff = "q10", max.cutoff = "q90")
p13

p14 <- FeaturePlot(subset(data.all.merged, origin = "OvarCancer1"),  features = c("WT1"),label = TRUE, repel = TRUE, split.by  = "origin", min.cutoff = "q10", max.cutoff = "q90")
p14

p15<- FeaturePlot(data.all.merged, features = c("EGFR","WT1"), min.cutoff = "q10", max.cutoff = "q90")
p15

p16 <- FeaturePlot(data.all.merged, features = c("TP53"), label = TRUE, repel = TRUE, split.by = "condition")
p16

pdf("p14.pdf")
p14
dev.off()

ggsave("WT1_FeaturePlots_splitby_origin.pdf", plot = "p14", device = "pdf", path = "C:/Users/Landstoh/Documents/WT1_FeaturePlots_splitby_origin.pdf") 
ggsave("S:/AG/AG-Scholz-NGS/Daten/Helena/WT1_FeaturePlots_splitby_origin.pdf", plot = p14, device = "pdf")
sessionInfo()

library(scCustomize)
library(tidyverse)

p14 <- FeaturePlot(seurat_object = data.all.merged, features = c("WT1"),reduction = "umap", repel = TRUE,colors_use = viridis_inferno_dark_high, split.by = "origin", min.cutoff = "q10", max.cutoff = "q90", ncol = 2) + patchwork::plot_layout(ncol = 4, nrow = 3)
p14

p17 <- FeaturePlot_scCustom(seurat_object = data.all.merged, colors_use = viridis_inferno_dark_high, features = c("WT1"), split.by = "origin")+ patchwork::plot_layout(ncol = 4, nrow = 3)
p17  
