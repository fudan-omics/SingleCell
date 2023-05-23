# Inspiration
## Work Flow
### Prompt
#### Background
- 有一个分析的工作，我们想找跟哮喘发病相关的关键代谢基因，打算从单细胞测序数据入手，看看有没有表达显著改变的代谢基因
- 分析方法使用R语言的Seurat包
#### Question
1. 使用已有的data matrix完成该分析工作需要哪些关键步骤，请给出一个可用工作流。
2. Seurat包中有哪些部分是完成该工作必须学习的？
### Answer
1. 以下是一个可用的工作流：
- 数据质量控制和过滤（例如，检查细胞数目、基因数目和质量控制指标，并去除垃圾细胞和低质量细胞）
- 数据归一化和标准化（例如，使用normalize.total、log1p和ScaleData函数）
- 特征选择和降维（例如，使用FindVariableFeatures和PCA函数）
- 群集分析和细胞类型鉴定（例如，使用FindClusters和DimPlot函数）
- 群集标记和差异表达分析（例如，使用ClusterMarkers、VlnPlot和FindMarkers函数）

2. 以下是Seurat包中完成该工作必须学习的部分：
- 数据加载和预处理（例如，Read10X、CreateSeuratObject和FilterCells函数）
- 群集分析和细胞类型鉴定（例如，FindClusters和DimPlot函数）
- 差异表达分析和可视化（例如，FindMarkers、VlnPlot和DotPlot函数）
- 可视化工具（例如，FeaturePlot、UMAP和t-SNE函数） 

## Understand websites
### 是比较对照和实验
[https://satijalab.org/seurat/articles/integration_introduction.html]
The link you provided is to a tutorial on how to integrate multiple single-cell RNA sequencing (scRNA-seq) datasets using the Seurat R package. The tutorial begins by discussing the challenges of integrating scRNA-seq datasets, such as identifying cell populations that are present across multiple datasets. It then introduces the Seurat integration methods, which first identify cross-dataset pairs of cells that are in a matched biological state ('anchors'). These methods can be used both to correct for technical differences between datasets (i.e. batch effect correction), and to perform comparative scRNA-seq analysis of across experimental conditions. The tutorial then provides a step-by-step guide on how to use the Seurat integration methods, using a dataset of human immune cells (PBMC) in either a resting or interferon-stimulated state as an example. Finally, the tutorial discusses some of the limitations of the Seurat integration methods and provides some tips for troubleshooting common problems.

# Code
## Utils
### Data Split
#### Task
df1 has 2 columns named as [treatment,cell], column names of df2 is the same as values in cell in df1. There are four types of `treatment`, split df2 by treatment type.
## Seurat
### Change Label Names
#### Background
- Using Seurat package in R
- R version 4.2.2
- Seurat version 4.1.3
#### Task
- Using `p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "orig.ident")` in Seurat to draw umap 
- orig.ident of my Seurat Object named `pbmc` is 'X1' and 'X2'
- I want to convert 'X1' into 'Control' while 'X2' to 'HMC'
- How to change the orig.ident of `pbmc` or change labels in umap figures
#### Answer
```R
# Assign new labels to the cells in your Seurat object
Idents(pbmc) <- c("Control", "HMC")[pbmc$orig.ident]
# Draw umap with updated labels
p1 <- DimPlot(pbmc, reduction = "umap", group.by = "orig.ident")
```
### Cluster
#### Background
- Using Seurat package in R
- R version 4.2.2
- Seurat version 4.1.3
- seurat object named `immune.combined`
#### Task
- There is a seurat object `immune.combined` which is analysed by the standard flow and is clustered
- A list named `cluster_list` stores dataframes produced by `FindMarkers(object = immune.combined, ident.1 = cluster_num, min.pct = 0.25)` of each cluster
- I want every cluster in `immune.combined` is automatically annotated cell type by markers with `SingleR` package
#### Answer
```R
library(SingleR)
# 获取Seurat对象中的表达矩阵和基因名
expr_matrix <- GetAssayData(object = immune.combined, slot = "counts")
gene_names <- rownames(expr_matrix)
# 对每个聚类进行注释
for (cluster_num in unique(Idents(object = immune.combined))) {
  # 获取聚类标识为cluster_num的细胞的表达矩阵和基因名
  cluster_expr_matrix <- expr_matrix[, Idents(object = immune.combined) == cluster_num]
  cluster_gene_names <- gene_names
  
  # 获取该聚类的marker基因
  cluster_markers <- FindMarkers(object = immune.combined, ident.1 = cluster_num, min.pct = 0.25)
  top_markers <- head(cluster_markers[order(cluster_markers$p_val_adj), ], 10)$gene
  
  # 使用SingleR进行注释
  annotation <- SingleR(test = cluster_expr_matrix, labels = NULL, gene.names = cluster_gene_names, 
                        ref = "refdata-celltype", de.tags = top_markers)$labels
  # 将注释结果保存到聚类对应的数据框中
  cluster_list[[cluster_num]]$cell_type <- annotation
}

# 查看注释结果
cluster_list

```