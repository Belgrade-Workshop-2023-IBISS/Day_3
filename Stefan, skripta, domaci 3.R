library(dplyr)
library(Seurat)
library(patchwork)

fetal59.data <- Read10X(data.dir = ("D:/Users/Stefan-s26/Desktop/domać/D59_fetal_filtered_gene_bc_matrices/GRCh38"))
fetal59 <- CreateSeuratObject(counts = fetal59.data, project = "pr1", min.cells = 3, min.features = 200)


fetal59[["percent.mt"]] <- PercentageFeatureSet(fetal59, pattern = "^MT-")
fetal59[["percent.rb"]] <- PercentageFeatureSet(fetal59, pattern = "^RPS|^RPL|^MRPS|^MRPL", assay = 'RNA')
VlnPlot(fetal59, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 4)
plot1 <- FeatureScatter(fetal59, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(fetal59, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
fetal59 <- subset(fetal59, subset = nCount_RNA > 300 & nCount_RNA < 10000 & nFeature_RNA > 200 & nFeature_RNA < 3500 & percent.mt < 15 & percent.rb < 45)


fetal59 <- NormalizeData(fetal59, normalization.method = "LogNormalize", scale.factor = 10000)
fetal59 <- FindVariableFeatures(fetal59, selection.method = "vst", nfeatures = 2000)
fetal59 <- ScaleData(fetal59, features = rownames(fetal59))


s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes
fetal59 <- CellCycleScoring(fetal59, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
head(fetal59)
fetal59 <- ScaleData(fetal59, verbose = T, vars.to.regress = c('nCount_RNA', 'nFeature_RNA', 'percent.mt', "percent.rb","S.Score","G2M.Score"))


ProcessSeu <- function(Seurat){
  Seurat <- NormalizeData(Seurat)
  Seurat <- FindVariableFeatures(Seurat, selection.method = "vst", nfeatures = 2000)
  Seurat <- ScaleData(Seurat)
  Seurat <- RunPCA(Seurat)
  Seurat <- FindNeighbors(Seurat, dims = 1:20)
  Seurat <- FindClusters(Seurat, resolution = 0.5)
  Seurat <- RunUMAP(Seurat, dims = 1:20)
  Seurat <- RunTSNE(Seurat,  dims.use = 1:20 )
  DimPlot(object = Seurat, reduction = "umap")
  return (Seurat)
} 


fetal59 <- ProcessSeu(fetal59)
DimPlot(fetal59, reduction = "umap") 


#remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
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




fetal59 <- RDoublet(fetal59)
head(fetal59)
fetal59 <- subset(fetal59, subset = DF.classifications_0.25_0.19_328 == 'Singlet')
fetal59 <- subset(fetal59, subset = DF.classifications_0.25_0.19_289 == 'Singlet')


saveRDS(fetal59, file = "D:/Users/Stefan-s26/Desktop/domać/D59_fetal_filtered_gene_bc_matrices/fetal59.rds")


fetal59$time.point <- "fetal59"
head(fetal59)

fetal125.data <- Read10X(data.dir = ("D:/Users/Stefan-s26/Desktop/domać/D125Cfetal_filtered_gene_bc_matrices/GRCh38"))
fetal125 <- CreateSeuratObject(counts = fetal125.data, project = "pr1", min.cells = 3, min.features = 200)


fetal125[["percent.mt"]] <- PercentageFeatureSet(fetal125, pattern = "^MT-")
fetal125[["percent.rb"]] <- PercentageFeatureSet(fetal125, pattern = "^RPS|^RPL|^MRPS|^MRPL", assay = 'RNA')
VlnPlot(fetal125, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 4)
plot1 <- FeatureScatter(fetal125, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(fetal125, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
fetal125 <- subset(fetal125, subset = nCount_RNA > 300 & nCount_RNA < 10000 & nFeature_RNA > 200 & nFeature_RNA < 3500 & percent.mt < 20 & percent.rb < 35)


fetal125 <- NormalizeData(fetal125, normalization.method = "LogNormalize", scale.factor = 10000)
fetal125 <- FindVariableFeatures(fetal125, selection.method = "vst", nfeatures = 2000)
fetal125 <- ScaleData(fetal125, features = rownames(fetal125))


s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes
fetal125 <- CellCycleScoring(fetal125, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
head(fetal125)
fetal125 <- ScaleData(fetal125, verbose = T, vars.to.regress = c('nCount_RNA', 'nFeature_RNA', 'percent.mt', "percent.rb","S.Score","G2M.Score"))


ProcessSeu <- function(Seurat){
  Seurat <- NormalizeData(Seurat)
  Seurat <- FindVariableFeatures(Seurat, selection.method = "vst", nfeatures = 2000)
  Seurat <- ScaleData(Seurat)
  Seurat <- RunPCA(Seurat)
  Seurat <- FindNeighbors(Seurat, dims = 1:20)
  Seurat <- FindClusters(Seurat, resolution = 0.5)
  Seurat <- RunUMAP(Seurat, dims = 1:20)
  Seurat <- RunTSNE(Seurat,  dims.use = 1:20 )
  DimPlot(object = Seurat, reduction = "umap")
  return (Seurat)
} 


fetal125 <- ProcessSeu(fetal125)
DimPlot(fetal125, reduction = "umap") 



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




fetal125 <- RDoublet(fetal125)
head(fetal125)
fetal125 <- subset(fetal125, subset = DF.classifications_0.25_0.005_1796 == 'Singlet')
fetal125 <- subset(fetal125, subset = DF.classifications_0.25_0.005_1561== 'Singlet')


saveRDS(fetal125, file = "D:/Users/Stefan-s26/Desktop/domać/D125Cfetal_filtered_gene_bc_matrices/GRCh38/fetal125.rds")
fetal125$time.point <- "fetal125"
head(fetal125)

#integration of data

#Choose the objects for integration

integration_list <- list(fetal59,fetal125)

features <- SelectIntegrationFeatures(object.list = integration_list)
data.anchors <- FindIntegrationAnchors(object.list = integration_list, anchor.features = features)

fetal.combined <- IntegrateData(anchorset = data.anchors)
ProcessInt <- function(data.integrated){
  data.integrated <- ScaleData(data.integrated, verbose = FALSE)
  data.integrated <- RunPCA(data.integrated, npcs = 30, verbose = FALSE)
  data.integrated <- FindNeighbors(data.integrated, dims = 1:20)
  data.integrated <- FindClusters(data.integrated, resolution = 0.5)
  data.integrated <- RunUMAP(data.integrated, reduction = "pca", dims = 1:20)
  data.integrated <- RunTSNE(data.integrated,  dims.use = 1:20 )
}

fetal.combined <- ProcessInt (fetal.combined)

DimPlot(fetal.combined, reduction = "umap") + ggtitle("fetal data combined")
DimPlot(fetal.combined, reduction = "umap", split.by = "time.point") + ggtitle("fetal data integrated")

retina.markers <- FindAllMarkers(fetal.combined, only.pos = TRUE)

library(tidyverse)

retina.markers %>%
  
  group_by(cluster) %>%
  
  dplyr::filter(avg_log2FC > 1) %>%
  
  slice_head(n = 6) %>%
  
  ungroup() -> top6



DoHeatmap(fetal.combined, features = top6$gene) 



new.cluster.ids <- c("Rod", "1", "Microglia", "3", "4", "5", "6", "7", "8", "9", "10", "Astrocyte", "12", "13", "Microglia", "15", "16")

names(new.cluster.ids) <- levels(fetal.combined)

fetal.combined <- RenameIdents(fetal.combined, new.cluster.ids)

DimPlot(fetal.combined, reduction = "umap", label = TRUE, pt.size = 0.5) + ggtitle("clusters")



plot <- DimPlot(fetal.combined, reduction = "umap", label = TRUE, label.size = 4.5) + xlab("UMAP_1") + ylab("UMAP_2") +
  
  theme(axis.title = element_text(size = 18), legend.text = element_text(size = 18)) + guides(colour = guide_legend(override.aes = list(size = 10)))

plot

FeaturePlot(retina.combined, features = c("RCVRN", "NRL", "TUBB4B"))

VlnPlot(retina.combined, features = c("RCVRN", "NRL", "TUBB4B"))





saveRDS(retina.combined, file = "D:/Users/Stefan-s26/Desktop/omać/fetal.combined.rds")
        