
fetal_59.data <- Read10X("D:/Users/Ana/Documents/Kursevi-i-obuke/Emils-bioinfo-workshop/week3/D59_fetal_filtered_gene_bc_matrices/GRCh38")
fetal_59 <- CreateSeuratObject(counts = fetal_59.data, project = "FetalRetina_Day59", min.cells = 3, min.features = 200)

fetal_125c.data <- Read10X("D:/Users/Ana/Documents/Kursevi-i-obuke/Emils-bioinfo-workshop/week3/D125Cfetal_filtered_gene_bc_matrices/GRCh38")
fetal_125c <- CreateSeuratObject(counts = fetal_125c.data, project = "FetalRetina_Day125c", min.cells = 3, min.features = 200)

####### Day 59 data (pre)processing

fetal_59[["percent.mt"]] <- PercentageFeatureSet(fetal_59, pattern = "^MT-")
fetal_59[["percent.rb"]] <- PercentageFeatureSet(fetal_59 , pattern = "^RPS|^RPL|^MRPS|^MRPL", assay = 'RNA')

VlnPlot(fetal_59, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 4)
FeatureScatter(fetal_59, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

fetal_59 <- NormalizeData(fetal_59)
fetal_59 <- FindVariableFeatures(fetal_59, selection.method = "vst")
fetal_59 <- ScaleData(fetal_59, features = rownames(fetal_59))
head(fetal_59)

s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes
fetal_59 <- CellCycleScoring(fetal_59, s.features = s.genes, g2m.features = g2m.genes, set.ident = T)
head(fetal_59)
fetal_59 <- ScaleData(fetal_59, verbose = T, vars.to.regress = c('nCount_RNA', 'nFeature_RNA', 'percent.mt', "percent.rb","S.Score","G2M.Score"))


fetal_59 <- subset(fetal_59, subset = nCount_RNA > 300 & nCount_RNA < 10000 & nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 15 & percent.rb < 30)
fetal_59

ProcessSeu <- function(Seurat){
  Seurat <- NormalizeData(Seurat)
  Seurat <- FindVariableFeatures(Seurat, selection.method = "vst", nfeatures = 2000)
  Seurat <- ScaleData(Seurat)
  Seurat <- RunPCA(Seurat)
  Seurat <- FindNeighbors(Seurat, dims = 1:20)
  Seurat <- FindClusters(Seurat, resolution = 0.5)
  Seurat <- RunUMAP(Seurat, dims = 1:20)
  Seurat <- RunTSNE(Seurat,  dims.use = 1:20 )
  return (Seurat)
}

fetal_59 <- ProcessSeu(fetal_59)
library(ggplot2)
DimPlot(fetal_59, reduction = "umap", label = TRUE) + ggtitle("Fetal Retina Day 59")


RDoublet <- function(tmp){
  sweep.res.list <- paramSweep(tmp, PCs = 1:30, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  pKopt <- as.numeric(as.character(bcmvn$pK[bcmvn$BCmetric == max(bcmvn$BCmetric)]))
  pKopt <- pKopt[order(pKopt, decreasing = TRUE) ]
  pKopt <- pKopt[1]
  homotypic.prop <- modelHomotypic(tmp$seurat_clusters) 
  nExp_poi <- round(0.1*length(colnames(tmp)))  ## Assuming 10% doublet formation rate 
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  tmp <- doubletFinder(tmp, PCs = 1:30, pN = 0.25, pK = pKopt, nExp = nExp_poi, reuse.pANN = FALSE)
  tmp <- doubletFinder(tmp, PCs = 1:30, pN = 0.25, pK = pKopt, nExp = nExp_poi.adj, reuse.pANN = paste("pANN_0.25",pKopt,nExp_poi, sep="_"))
  return (tmp) 
}

library(DoubletFinder)

fetal_59 <- RDoublet(fetal_59)
head(fetal_59)

fetal_59 <- subset(fetal_59, subset = DF.classifications_0.25_0.01_260 == 'Singlet') #remove the doublets
fetal_59 <- subset(fetal_59, subset = DF.classifications_0.25_0.01_231 == 'Singlet')
fetal_59$time_point <- "day59"

saveRDS(fetal_59, file = "D:/Users/Ana/Documents/Kursevi-i-obuke/Emils-bioinfo-workshop/week3/FetalRetina_Day59.rds")


###### Day 125c data (pre)processing

fetal_125c[["percent.mt"]] <- PercentageFeatureSet(fetal_125c, pattern = "^MT-")
fetal_125c[["percent.rb"]] <- PercentageFeatureSet(fetal_125c , pattern = "^RPS|^RPL|^MRPS|^MRPL", assay = 'RNA')

VlnPlot(fetal_125c, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 4)
FeatureScatter(fetal_125c, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

fetal_125c <- NormalizeData(fetal_125c)
fetal_125c <- FindVariableFeatures(fetal_125c, selection.method = "vst")
fetal_125c <- ScaleData(fetal_125c, features = rownames(fetal_59))
head(fetal_125c)

fetal_125c <- CellCycleScoring(fetal_125c, s.features = s.genes, g2m.features = g2m.genes, set.ident = T)
head(fetal_125c)
fetal_125c <- ScaleData(fetal_125c, verbose = T, vars.to.regress = c('nCount_RNA', 'nFeature_RNA', 'percent.mt', "percent.rb","S.Score","G2M.Score"))


fetal_125c <- subset(fetal_125c, subset = nCount_RNA > 300 & nCount_RNA < 10000 & nFeature_RNA > 200 & nFeature_RNA < 2000 & percent.mt < 25 & percent.rb < 30)
fetal_125c

fetal_125c <- ProcessSeu(fetal_125c)
DimPlot(fetal_125c, reduction = "umap", label = TRUE) + ggtitle("Fetal Retina Day 125C")

fetal_125c <- RDoublet(fetal_125c)
head(fetal_125c)

fetal_125c <- subset(fetal_125c, subset = DF.classifications_0.25_0.3_1628 == 'Singlet') #remove the doublets
fetal_125c <- subset(fetal_125c, subset = DF.classifications_0.25_0.3_1385 == 'Singlet')
fetal_125c$time_point <- "day125c"
head(fetal_125c)

saveRDS(fetal_125c, file = "D:/Users/Ana/Documents/Kursevi-i-obuke/Emils-bioinfo-workshop/week3/FetalRetina_Day125c.rds")


##### Data integration 

ProcessInt <- function(data.integrated){ 
  data.integrated <- ScaleData(data.integrated, verbose = FALSE)
  data.integrated <- RunPCA(data.integrated, npcs = 30, verbose = FALSE)
  data.integrated <- FindNeighbors(data.integrated, dims = 1:20)
  data.integrated <- FindClusters(data.integrated, resolution = 0.5)
  data.integrated <- RunUMAP(data.integrated, reduction = "pca", dims = 1:20)
  data.integrated <- RunTSNE(data.integrated,  dims.use = 1:20 )
}

integration_list <- list(fetal_59, fetal_125c)

features <- SelectIntegrationFeatures(object.list = integration_list)
data.anchors <- FindIntegrationAnchors(object.list = integration_list, anchor.features = features)

data.combined <- IntegrateData(anchorset = data.anchors)

fetal_combined <- ProcessInt(data.combined)

DimPlot(fetal_combined, reduction = "umap", label = TRUE) + ggtitle("Fetal integrated data")

DimPlot(fetal_combined, reduction = "umap", label = TRUE, split.by = "time_point") + ggtitle("Fetal integrated data")


####### Assigning cell type

fetal_markers <- FindAllMarkers(fetal_combined, only.pos = TRUE)

library(tidyverse)

fetal_markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(fetal_combined, features = top10$gene) + NoLegend()

write.csv(top10, "D:/Users/Ana/Documents/Kursevi-i-obuke/Emils-bioinfo-workshop/week3/top10-markers-in-every-cluster.csv", row.names = FALSE)

# Jovana found us a paper with retinal cell markers :) 

new.cluster.ids <- c("Rod", "1", "Microglia", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13")
names(new.cluster.ids) <- levels(fetal_combined)
fetal_combined <- RenameIdents(fetal_combined, new.cluster.ids)
DimPlot(fetal_combined, reduction = "umap", label = TRUE, pt.size = 0.5) + ggtitle("Retinal cell types")

# plot <- DimPlot(fetal_combined, reduction = "umap", label = TRUE, label.size = 4.5) + xlab("UMAP 1") + ylab("UMAP 2") +
#   theme(axis.title = element_text(size = 18), legend.text = element_text(size = 18)) + guides(colour = guide_legend(override.aes = list(size = 10)))
# plot
FeaturePlot(fetal_combined, features = c("RCVRN", "NRL", "TUBB4B")) # Rod markers
VlnPlot(fetal_combined, features = c("RCVRN", "NRL", "TUBB4B")) # Rod markers


saveRDS(fetal_combined, file = "D:/Users/Ana/Documents/Kursevi-i-obuke/Emils-bioinfo-workshop/week3/fetal_combined.rds")
