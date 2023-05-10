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