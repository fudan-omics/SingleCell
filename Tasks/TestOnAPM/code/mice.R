##设置镜像
##安装Seurat
# install.packages("Seurat")

library(Matrix)
library(Seurat)
library(dplyr)
# setwd('D:/RNAseq/asthma Th2 scRNAseq_GSE131935/Rresult/') 
getwd()
matrix_dir = 'D:/repositories/SingleCell/Tasks/TestOnAPM/data/MICE/' 
filepath=paste0(matrix_dir,"umi-count-matrix.csv")
filepath
# raw_counts<-read.table(file=filepath, header = T, row.names=1, sep="," ,as.is=T)
raw_counts<-read.table(file=filepath, header = T, sep="," ,as.is=T)
colnames(raw_counts)[1:50]
raw_counts[1:20,1:2]
sample_df<-read.table(file='D:/repositories/SingleCell/Tasks/TestOnAPM/data/Mice/GSE155391_cell-level-metadata.csv', header = T, sep="," ,as.is=T)
raw_counts.list <- split(raw_counts, sample_df$treatment)
dim(raw_counts)
typeof(raw_counts)
str(raw_counts)
treatments=c(unique(sample_df$treatment))
treatments
dim(raw_counts[[treatments[1]]])
raw_counts[1:6, 1:6]
dim(raw_counts)
hdm_df=rbind(raw_counts.list$Control,raw_counts.list$HDM)

dim(raw_counts.list$HDM)
dim(raw_counts.list$Control)
# write.csv(raw_counts,paste0(matrix_dir,"test.csv"),row.names=T)
pbmc <- CreateSeuratObject(counts = raw_counts, project = "m1_Th", min.cells = 3, min.features = 200,names.field = 1)
table(Idents(object = pbmc))
pbmc.list <- SplitObject(pbmc, split.by = "orig.ident")
hdm.list=list(pbmc.list$X1,pbmc.list$X2)
names(hdm.list)=c('Control','HDM')
hdm.list
test.list=hdm.list
test.list <- lapply(X = test.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = test.list)
features
immune.anchors <- FindIntegrationAnchors(object.list = test.list, anchor.features = features)
immune.combined <- IntegrateData(anchorset = immune.anchors)
DefaultAssay(immune.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindClusters(immune.combined, resolution = 0.5)
DimPlot(immune.combined, reduction = "umap", split.by = "orig.ident")

dim(pbmc)
GetAssayData(object = pbmc, slot = 'counts')[9:15, 7:18]
min.cells = 3
min.features = 200

pbmc[["percent.mt"]] <- PercentageFeatureSet(object = pbmc, pattern = "^mt-") ##线粒体
pbmc
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0.05)
pbmc <- subset(pbmc, subset = nFeature_RNA > 300 & nFeature_RNA  < 4000 & percent.mt < 8)
pbmc <- NormalizeData(pbmc)
dim(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
dim(pbmc)
# Identify the 10 most highly variable genes
top10 <- head(x = VariableFeatures(object = pbmc), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(object = pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))

# process data
pbmc <- ScaleData(pbmc)
# linear reduce dimension
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
# non-linear reduce dimension
pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap")
pbmc <- RunTSNE(pbmc, features = VariableFeatures(object = pbmc))
DimPlot(pbmc, reduction = "tsne")

##找到每个cluster的marker
markers <- FindMarkers(pbmc, ident.1 = 0,  only.pos = TRUE,min.pct = 0.25)
rownames(markers)
raw_counts[rownames(markers),1][1:2]
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









