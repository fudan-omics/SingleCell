##设置镜像
##安装Seurat
# install.packages("Seurat")

library(Matrix)
library(Seurat)
library(dplyr)
# setwd('D:/RNAseq/asthma Th2 scRNAseq_GSE131935/Rresult/') 
getwd()

##示例 https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE136916
matrix_dir = 'D:/repositories/SingleCell/Tasks/TestOnAPM/data/0422/' 
barcode.path<-paste0(matrix_dir,"barcodes.tsv")
genes.path<-paste0(matrix_dir,"features.tsv")
matrix.path<-paste0(matrix_dir,"matrix.mtx")
zebrafish.data <- readMM(file = matrix.path) ##mac上不能读压缩文件
gene.names = read.delim(genes.path,header = FALSE, stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path,header = FALSE,stringsAsFactors = FALSE)
colnames(zebrafish.data) = barcode.names$V1
rownames(zebrafish.data) = gene.names$V2 ##把示例中的V1改成V2
zebrafish.data[1:6, 1:6] ##check矩阵
dim(zebrafish.data) ##check矩阵

pbmc <- CreateSeuratObject(counts = zebrafish.data, project = "lung")
GetAssayData(object = pbmc, slot = 'counts')[9:15, 7:18]

pbmc[["percent.mt"]] <- PercentageFeatureSet(object = pbmc, pattern = "^mt-") ##线粒体
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0.05)
pbmc <- subset(pbmc, subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & percent.mt < 8)
pbmc <- NormalizeData(pbmc)

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
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
write.table(markers, file="/Users/apple/Desktop/Result/markers_0.tsv", sep="\t", quote=FALSE, col.names=NA)

##直接根据marker画图

table <- read.table ( '/Users/apple/Desktop/Rwork/test1.csv', header = FALSE)
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









