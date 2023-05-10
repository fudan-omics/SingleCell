##设置镜像
##安装Seurat
# install.packages("Seurat")

library(Matrix)
library(Seurat)
library(dplyr)
# setwd('D:/RNAseq/asthma Th2 scRNAseq_GSE131935/Rresult/') 
getwd()
matrix_dir = 'D:/repositories/SingleCell/Tasks/TestOnAPM/data/0422/' 
filepath=paste0(matrix_dir,"SS2_15.csv")
filepath
# raw_counts<-read.table(file=filepath, header = T, row.names=1, sep="," ,as.is=T)
raw_counts<-read.table(file=filepath, header = T, sep="," ,as.is=T)
raw_counts[1:6, 1:6]
raw_counts[,1]
rownames(raw_counts)=raw_counts[,1]

dim(raw_counts)
# write.csv(raw_counts,paste0(matrix_dir,"test.csv"),row.names=T)
pbmc <- CreateSeuratObject(counts = raw_counts, project = "m1_Th", min.cells = 3, min.features = 200)

pbmc
GetAssayData(object = pbmc, slot = 'counts')[9:15, 7:18]
min.cells = 3
min.features = 200

pbmc[["percent.mt"]] <- PercentageFeatureSet(object = pbmc, pattern = "^mt-") ##线粒体
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0.05)
pbmc <- subset(pbmc, subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & percent.mt < 8)
pbmc <- NormalizeData(pbmc)
dim(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
dim(pbmc)

pbmc <- ScaleData(pbmc)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap")
pbmc <- RunTSNE(pbmc, features = VariableFeatures(object = pbmc))
DimPlot(pbmc, reduction = "tsne")

##找到每个cluster的marker
markers <- FindMarkers(pbmc, ident.1 = 0,  only.pos = TRUE,min.pct = 0.25)
marker.path=paste0(matrix_dir,'markers_0.csv')
write.table(markers, file=marker.path, sep="\t", quote=FALSE, col.names=NA)

##直接根据marker画图

table <- read.table (marker.path, sep="\t",header = FALSE)
genes <- as.vector(table[1:4,1])
plots <- VlnPlot(pbmc, features = genes, group.by = "seurat_clusters", pt.size = 0, combine = FALSE)
a <- CombinePlots(plots = plots, ncol = 1)
ggsave("violin_test_1.png", a)

table <- read.table ( '/Users/apple/Desktop/Rwork/test1.csv', header = FALSE)
genes <- as.vector(table[8:14,1])
plots <- VlnPlot(pbmc, features = genes, group.by = "seurat_clusters", pt.size = 0, combine = FALSE)
a <- CombinePlots(plots = plots, ncol = 1)
ggsave("violin_test_2.png", a)

table <- read.table ( '/Users/apple/Desktop/Rwork/test1.csv', header = FALSE)
genes <- as.vector(table[15:21,1])
plots <- VlnPlot(pbmc, features = genes, group.by = "seurat_clusters", pt.size = 0, combine = FALSE)
a <- CombinePlots(plots = plots, ncol = 1)
ggsave("violin_test_3.png", a)
table <- read.table ( '/Users/apple/Desktop/Rwork/test1.csv', header = FALSE)
genes <- as.vector(table[22:28,1])
plots <- VlnPlot(pbmc, features = genes, group.by = "seurat_clusters", pt.size = 0, combine = FALSE)
a <- CombinePlots(plots = plots, ncol = 1)
ggsave("violin_test_4.png", a)

a <- FeaturePlot(pbmc, reduction = "tsne",features = "Cidec", cols = c("grey", "red"))
ggsave("genemap_1.png", a)

##标记
table <- read.table ( '/Users/apple/Desktop/Rwork/test1.csv', header = FALSE)
genes <- as.vector(table[1:8,1])
new.cluster.ids <- genes
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
a<-DimPlot(pbmc, reduction = "tsne", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave("tsne.png", a)









