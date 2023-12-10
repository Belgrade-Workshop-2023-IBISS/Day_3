library(dplyr)
library(Seurat)
library(patchwork)
D59.data <- Read10X(data.dir = "D:/Users/Jovana-s26/Desktop/Emils-bioinfo-workshop/week3/D59_fetal_filtered_gene_bc_matrices/GRCh38")
D59 <- CreateSeuratObject(counts = D59.data, project = "retina", min.cells = 3, min.features = 200)

D59[["percent.mt"]] <- PercentageFeatureSet(D59, pattern = "^MT-")
D59[["percent.rb"]] <- PercentageFeatureSet(D59, pattern = "^RPS|^RPL|^MRPS|^MRPL", assay = 'RNA')
VlnPlot(D59, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 4)
plot1 <- FeatureScatter(D59, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(D59, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
D59 <- subset(D59, subset = nCount_RNA > 300 & nCount_RNA < 10000 & nFeature_RNA > 200 & nFeature_RNA < 3500 & percent.mt < 15 & percent.rb < 45)

D59 <- NormalizeData(D59, normalization.method = "LogNormalize", scale.factor = 10000)
D59 <- FindVariableFeatures(D59, selection.method = "vst", nfeatures = 2000)
D59 <- ScaleData(D59, features = rownames(D59))

s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes
D59 <- CellCycleScoring(D59, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
head(D59)
D59 <- ScaleData(D59, verbose = T, vars.to.regress = c('nCount_RNA', 'nFeature_RNA', 'percent.mt', "percent.rb","S.Score","G2M.Score"))

ProcessSeu <- function(Seurat){
  Seurat <- NormalizeData(Seurat)
  Seurat <- FindVariableFeatures(Seurat, selection.method = "vst", nfeatures = 2000)
  Seurat <- ScaleData(Seurat)
  Seurat <- RunPCA(Seurat)
  Seurat <- FindNeighbors(Seurat, dims = 1:20)
  Seurat <- FindClusters(Seurat, resolution = 0.5)
  Seurat <- RunUMAP(Seurat, dims = 1:20)
  Seurat <- RunTSNE(Seurat,  dims.use = 1:20)
  return (Seurat)
}  

D59 <- ProcessSeu(D59)
library(ggplot2)
DimPlot(D59, reduction = "umap", label = TRUE) + ggtitle("retina D59")

#remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
library(DoubletFinder)

RDoublet <- function(tmp){
  sweep.res.list <- paramSweep(tmp, PCs = 1:30, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  pKopt <- as.numeric(as.character(bcmvn$pK[bcmvn$BCmetric == max(bcmvn$BCmetric)]))
  pKopt <- pKopt[order(pKopt, decreasing = TRUE) ]
  pKopt <- pKopt[1]
  homotypic.prop <- modelHomotypic(tmp$seurat_clusters) 
  nExp_poi <- round(0.1*length(colnames(tmp)))
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  tmp <- doubletFinder(tmp, PCs = 1:30, pN = 0.25, pK = pKopt, nExp = nExp_poi, reuse.pANN = FALSE)
  tmp <- doubletFinder(tmp, PCs = 1:30, pN = 0.25, pK = pKopt, nExp = nExp_poi.adj, reuse.pANN = paste("pANN_0.25",pKopt,nExp_poi, sep="_"))
  return (tmp) 
}


D59 <- RDoublet(D59)
head(D59)
D59 <- subset(D59, subset = DF.classifications_0.25_0.2_328 == 'Singlet')
D59 <- subset(D59, subset = DF.classifications_0.25_0.2_289 == 'Singlet') 

D59$time.point <- "D59"
head(D59)
saveRDS(D59, file = "E:/Bioinformatics/D59first.rds")


D125C.data <- Read10X(data.dir = "D:/Users/Jovana-s26/Desktop/Emils-bioinfo-workshop/week3/D125Cfetal_filtered_gene_bc_matrices/GRCh38")
D125C <- CreateSeuratObject(counts = D125C.data, project = "retina", min.cells = 3, min.features = 200)

D125C[["percent.mt"]] <- PercentageFeatureSet(D125C, pattern = "^MT-")
D125C[["percent.rb"]] <- PercentageFeatureSet(D125C, pattern = "^RPS|^RPL|^MRPS|^MRPL", assay = 'RNA')
VlnPlot(D125C, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 4)
plot1 <- FeatureScatter(D125C, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(D125C, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
D125C <- subset(D125C, subset = nCount_RNA > 300 & nCount_RNA < 10000 & nFeature_RNA > 200 & nFeature_RNA < 2000 & percent.mt < 25 & percent.rb < 40)

D125C <- NormalizeData(D125C, normalization.method = "LogNormalize", scale.factor = 10000)
D125C <- FindVariableFeatures(D125C, selection.method = "vst", nfeatures = 2000)
D125C <- ScaleData(D125C, features = rownames(D125C))

D125C <- CellCycleScoring(D125C, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
head(D125C)
D125C <- ScaleData(D125C, verbose = T, vars.to.regress = c('nCount_RNA', 'nFeature_RNA', 'percent.mt', "percent.rb","S.Score","G2M.Score"))

D125C <- ProcessSeu(D125C)
DimPlot(D125C, reduction = "umap", label = TRUE) + ggtitle("retina D125C")

D125C <- RDoublet(D125C)
head(D125C)
D125C <- subset(D125C, subset = DF.classifications_0.25_0.005_1830 == 'Singlet')
D125C <- subset(D125C, subset = DF.classifications_0.25_0.005_1598 == 'Singlet') 

D125C$time.point <- "D125C"
head(D125C)
saveRDS(D125C, file = "E:/Bioinformatics/D125Cfirst.rds")



ProcessInt <- function(data.integrated){ 
  data.integrated <- ScaleData(data.integrated, verbose = FALSE)
  data.integrated <- RunPCA(data.integrated, npcs = 30, verbose = FALSE)
  data.integrated <- FindNeighbors(data.integrated, dims = 1:20)
  data.integrated <- FindClusters(data.integrated, resolution = 0.5)
  data.integrated <- RunUMAP(data.integrated, reduction = "pca", dims = 1:20)
  data.integrated <- RunTSNE(data.integrated,  dims.use = 1:20 )
}


integration_list <- list(D59, D125C)

features <- SelectIntegrationFeatures(object.list = integration_list)
data.anchors <- FindIntegrationAnchors(object.list = integration_list, anchor.features = features)

data.combined <- IntegrateData(anchorset = data.anchors)

retina.combined <- ProcessInt(data.combined)
DimPlot(retina.combined, reduction = "umap", label = TRUE) + ggtitle("retina combined")
DimPlot(retina.combined, reduction = "umap", split.by = "time.point", label = TRUE)


retina.markers <- FindAllMarkers(retina.combined, only.pos = TRUE)
library(tidyverse)
retina.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 3) %>%
  ungroup() -> top3

DoHeatmap(retina.combined, features = top3$gene) + NoLegend()

new.cluster.ids <- c("Rod", "1", "Microglia", "3", "4", "5", "6", "7", "8", "9", "10", "Astrocyte", "12", "13", "Microglia")  #we used https://doi.org/10.1016/j.xgen.2022.100164 and Data S1. Marker genes from scRNA-seq by cell type, related to Figure 1 for annotation
names(new.cluster.ids) <- levels(retina.combined)
retina.combined <- RenameIdents(retina.combined, new.cluster.ids)
DimPlot(retina.combined, reduction = "umap", label = TRUE, pt.size = 0.5) + ggtitle("clusters")

plot <- DimPlot(retina.combined, reduction = "umap", label = TRUE, label.size = 4.5) + xlab("UMAP_1") + ylab("UMAP_2") +
  theme(axis.title = element_text(size = 18), legend.text = element_text(size = 18)) + guides(colour = guide_legend(override.aes = list(size = 10)))
plot
FeaturePlot(retina.combined, features = c("RCVRN", "NRL", "TUBB4B"))
VlnPlot(retina.combined, features = c("RCVRN", "NRL", "TUBB4B"))


saveRDS(retina.combined, file = "D:/Users/Jovana-s26/Desktop/Emils-bioinfo-workshop/week3/retina.combined.rds")

