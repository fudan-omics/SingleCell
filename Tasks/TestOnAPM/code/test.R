install.packages('devtools')
library(Seurat)
library(SeuratData)
library(patchwork)
library(devtools)
install_github('satijalab/seurat-data')
InstallData("ifnb")
LoadData("ifnb")

