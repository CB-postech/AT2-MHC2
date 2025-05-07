set.seed(12345)

BiocManager::install("DropletUtils")
install.packages("glmnet")
library(Matrix)
library(Seurat)
library(SeuratObject)
library(dplyr)
library(writexl)
library(readxl)
library(ggplot2)
library(readr)
library(Rcpp)
library(harmony)
library(glmnet) ## for makeX(): Convert data.frame to dgCMatrix
library(DropletUtils)
library(stringr)

file_path="D:/MJY/GSE155900/data"
setwd(file_path)

############### 01. Read raw data and create seurat object ##########################

# Create a list to store file paths
files <- list.files(pattern = "GSM.*.expression_matrix.txt.gz$", full.names = TRUE)
sample_ids <- sapply(files, function(x) str_extract(x, "P\\d+"))
seurat_list <-  lapply(seq_along(files), function(i) {
  data <- read.table(files[i], header = TRUE, row.names = 1)
  seurat_obj <- CreateSeuratObject(counts = data)
  return(seurat_obj)
})

merged_seurat <- merge(seurat_list[[1]], 
                       y = seurat_list[-1],
                       add.cell.ids = sample_ids,
                       project = "merged")

merged_seurat$sample <- ""
merged_seurat$sample <- substring(colnames(merged_seurat), 0,3)

############### 02. Pre-processing (Normalization, Scaling, Dimension reduction) ##########################
merged_seurat <- NormalizeData(merged_seurat)
merged_seurat <- FindVariableFeatures(merged_seurat)
merged_seurat <- ScaleData(merged_seurat)

print(merged_seurat)

###Dimension reduction
merged_seurat <- RunPCA(merged_seurat, features=VariableFeatures(object = merged_seurat))
ElbowPlot(merged_seurat, ndims=50)

n_pcs=20; cluster_resolution =0.8

merged_seurat <- FindNeighbors(merged_seurat, dims=1:n_pcs)
merged_seurat <- FindClusters(merged_seurat, resolution = cluster_resolution)
merged_seurat <- RunUMAP(merged_seurat, dims=1:n_pcs)

DimPlot(merged_seurat, label=T)
unique(merged_seurat@meta.data)
DimPlot(merged_seurat, group.by = "sample")

FeaturePlot(merged_seurat, features = "EPCAM")
FeaturePlot(merged_seurat, features = c("SFTPC", "SCGB3A2", "SOX9", "SOX2"))
FeaturePlot(merged_seurat, features = c("LYZ", "COL1A2", "ACTA2", "EPCAM", "CD3D", "CD79A", "PECAM1", "HBG1", "ID3", "UNC5C"))

merged_seurat$stage <- as.character(merged_seurat$sample)

Idents(merged_seurat) <- "stage"
merged_seurat$stage[WhichCells(merged_seurat, idents = c("P1_", "P2_", "P3_","P4_", "P6_", "P7_", "P8_", "P9_"))] <- '< 1yr'
merged_seurat$stage[WhichCells(merged_seurat, idents = c("P5_"))] <- '3-4yr'
merged_seurat$stage[WhichCells(merged_seurat, idents = c("P10"))] <- '2-3yr'

############### 03. Epithelial cell subsetting ##########################
Idents(merged_seurat) <- "seurat_clusters"
Epi_merged_seurat <- subset(merged_seurat, subset = seurat_clusters %in% c("24","7", "4", "3", "0", "32","9"))
DimPlot(Epi_merged_seurat, label = T)
DimPlot(Epi_merged_seurat, label = T, group.by = "stage")

### Normalization
Epi_merged_seurat <- NormalizeData(Epi_merged_seurat)

## feature selection; Find HVGs (Used 2000 genes with high cell-to-cell variation, which were calculated using the FindVariableFeatures function in seurat for further dimensionality reduction)
Epi_merged_seurat <- FindVariableFeatures(Epi_merged_seurat) # Default nFeatures = 2000

###Scaling; z-transformation (Effects of variable(percent.mt) were estimated and regressed out using a GLM (ScaleData function, model.use="linear"))
Epi_merged_seurat <- ScaleData(Epi_merged_seurat)

###Dimension reduction
Epi_merged_seurat <- RunPCA(Epi_merged_seurat, features=VariableFeatures(object = Epi_merged_seurat))
ElbowPlot(Epi_merged_seurat, ndims=50)

n_pcs=15; cluster_resolution =0.8

Epi_merged_seurat <- FindNeighbors(Epi_merged_seurat, dims=1:n_pcs)
Epi_merged_seurat <- FindClusters(Epi_merged_seurat, resolution = cluster_resolution)
Epi_merged_seurat <- RunUMAP(Epi_merged_seurat, dims=1:n_pcs)

DimPlot(Epi_merged_seurat, label=T)
DimPlot(Epi_merged_seurat, label=T, group.by = "stage")
FeaturePlot(Epi_merged_seurat, features = c('SOX2', 'SOX9', 'SFTPC', 'AGER', 'SCGB1A1', 'FAM183A'))
FeaturePlot(Epi_merged_seurat, features = c('LAMP3', 'ABCA3'))

saveRDS(Epi_merged_seurat, file="Epi_merged_seurat_normed.rds")

######### Batch correction (harmony)

harmony_Epi_merged_seurat = Epi_merged_seurat
harmony_Epi_merged_seurat <- RunHarmony(harmony_Epi_merged_seurat, "sample", max.iter = 30)
ElbowPlot(harmony_Epi_merged_seurat, reduction = "harmony", ndims = 30)

harmony.dims = 10

harmony_Epi_merged_seurat <- harmony_Epi_merged_seurat %>%
  RunUMAP(reduction = "harmony", dims = 1:harmony.dims) %>%
  FindNeighbors(reduction = "harmony", dims = 1:harmony.dims) %>%
  FindClusters()

DimPlot(harmony_Epi_merged_seurat, label = T)
DimPlot(harmony_Epi_merged_seurat, label = T, split.by = "stage")
FeaturePlot(harmony_Epi_merged_seurat, features = c('SOX2', 'SOX9', 'SFTPC', 'AGER', 'SCGB1A1', 'FAM183A'))
FeaturePlot(harmony_Epi_merged_seurat, features = c('EPCAM'))

######### 04. subset annotation (dividing timepoint for annotation is unnecessary) ######################

# DotPlot(harmony_Epi_merged_seurat, group.by = "seurat_clusters", features = c("C1orf194", "C20orf85", "CAPS", "FAM183A",
#                                                                               "MUC5B", "SCGB1A1", "SCGB3A1", "TFF3",
#                                                                               "AGBL4", "PPP1R9A", "SOX5",
#                                                                               "AGER", "CAV1", "CLIC3",
#                                                                               "PGC", "SFTPB", "SFTPC",
#                                                                               "HOXB5", "MFAP2", "MFAP4",
#                                                                               "KRT17", "MKI67", "TOP2A",
#                                                                               "S100A8")) +
#   theme(axis.text.x = element_text(angle=45, hjust=1))

# DotPlot(harmony_Epi_merged_seurat, group.by = "seurat_clusters", features = c("CD14","CD68","CST3", "LYZ", "PLAUR",
#                                                                               "COL14A1", "COL1A2", "LUM", "PDGFRA", "TCF21",
#                                                                               "ACTA2", "CNN1", "MYH11", "PDGFRB","TAGLN",
#                                                                               "CD3D", "CD3E", "KLRD1", "PTPRC",
#                                                                               "CD79A", "EBF1", "IGHM",
#                                                                               "CAVIN2", "IGFBP4", "PECAM1",
#                                                                               "HBA2", "HBB", "HBG1",
#                                                                               "ADGRB3", "CNTNAP2", "ROR1")) +
#   theme(axis.text.x = element_text(angle=45, hjust=1))

harmony_Epi_merged_seurat$annot <- as.character(harmony_Epi_merged_seurat$seurat_clusters)
Idents(harmony_Epi_merged_seurat) <- "annot"

harmony_Epi_merged_seurat$annot[WhichCells(harmony_Epi_merged_seurat, idents = c("0", "1", "2", "3", "4","6", "9"))] <- 'AT2'
harmony_Epi_merged_seurat$annot[WhichCells(harmony_Epi_merged_seurat, idents = c("5","8"))] <- 'AT1'
harmony_Epi_merged_seurat$annot[WhichCells(harmony_Epi_merged_seurat, idents = c("7","10", "11"))] <- 'Airway epithelial cells'

DimPlot(harmony_Epi_merged_seurat, label = T)
saveRDS(harmony_Epi_merged_seurat,"harmony_Epi_merged_annot_GSE155900.rds" )
