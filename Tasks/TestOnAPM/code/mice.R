##设置镜像
##安装Seurat
# install.packages("Seurat")
BiocManager::install("SingleR")
browseVignettes("SingleR")




library(Matrix)
library(Seurat)
library(dplyr)
library(ggplot2)
library(SingleR)

## Prepare Data ##
# # setwd('D:/RNAseq/asthma Th2 scRNAseq_GSE131935/Rresult/') 
# getwd()
# matrix_dir = 'D:/repositories/SingleCell/Tasks/TestOnAPM/data/MICE/' 
# filepath=paste0(matrix_dir,"umi-count-matrix.csv")
# filepath
# # raw_counts<-read.table(file=filepath, header = T, row.names=1, sep="," ,as.is=T)
# raw_counts<-read.table(file=filepath, header = T, , row.names=1,sep="," ,as.is=T)
# colnames(raw_counts)[1:5]
# rownames(raw_counts)[1:5]
# raw_counts[1:20,1:2]
# sample_df<-read.table(file='D:/repositories/SingleCell/Tasks/TestOnAPM/data/Mice/GSE155391_cell-level-metadata.csv', header = T, sep="," ,as.is=T)
# dim(raw_counts)
# typeof(raw_counts)
# # str(raw_counts)
# treatments=c(unique(sample_df$treatment))
# treatments
# dim(raw_counts[[treatments[1]]])
# raw_counts[1:6, 1:6]
# dim(raw_counts)

# write.csv(raw_counts,paste0(matrix_dir,"test.csv"),row.names=T)
# pbmc <- CreateSeuratObject(counts = raw_counts, project = "m1_Th", min.cells = 3, min.features = 200,names.field = 1)
# saveRDS(pbmc, file = "D:/repositories/SingleCell/Tasks/TestOnAPM/results/Mice/hdm_ctrl.rds")
## Prepare Data
pbmc=readRDS("D:/repositories/SingleCell/Tasks/TestOnAPM/results/Mice/hdm_ctrl.rds")


pbmc[["percent.mt"]] <- PercentageFeatureSet(object = pbmc, pattern = "^mt-") ##线粒体
pbmc
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),idents = c('Control','HDM'),ncol = 3, pt.size=0.05)


table(Idents(object = pbmc))
# Assign new labels to the cells in your Seurat object
Idents(pbmc) <- c("Control", "HDM",'APM','HMD_APM')[pbmc$orig.ident]
table(Idents(object = pbmc))

# Draw umap with updated labels
p1 <- DimPlot(pbmc, reduction = "umap", group.by = "orig.ident")
print(test_pbmc)
print(pbmc)
## Split Dataset
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
table(Idents(object = immune.combined))

DefaultAssay(immune.combined) <- "integrated"
immune.combined
# Run the standard workflow for visualization and clustering
table(Idents(object = immune.combined))

immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
immune.combined <- RunUMAP(immune.combined, dims = 1:30)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindClusters(immune.combined, resolution = 0.5) # orig.ident will be changed after this step
p1 <- DimPlot(immune.combined, reduction = "tsne", group.by = "orig.ident")
p2 <- DimPlot(immune.combined, reduction = "tsne", label = TRUE, repel = TRUE)
DimPlot(immune.combined, reduction = "umap", label = TRUE, repel = TRUE)
print(p1+p2)
print(p1)

levels(immune.combined)
cluster_list=list()
for (i in levels(immune.combined)){
    cluster_num=strtoi(i)
    cluster_list[cluster_num] <- list(FindMarkers(object = immune.combined, ident.1 = cluster_num, min.pct = 0.25))

}
cluster_list[[1]]
length(cluster_list)
cluster_list[[3]]
cluster1.markers <- list(FindMarkers(object = immune.combined, ident.1 = 1, min.pct = 0.25))
cluster1.markers
top_markers <- head(vector(cluster1.markers[order(cluster1.markers$p_val_adj), ]$gene))

cluster_list[1]=cluster1.markers
cluster_list[[1]]
immune.combined  <- RunTSNE(immune.combined, dims = 1:30)
DimPlot(immune.combined , reduction = "tsne")
# ggsave('D:/repositories/SingleCell/Tasks/TestOnAPM/results/Mice/p1.png',p1)
# ggsave('D:/repositories/SingleCell/Tasks/TestOnAPM/results/Mice/p2.png',p2)
Version(pbmc)
dim(pbmc)



##### Referenece ####
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


##### Reference #####
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









