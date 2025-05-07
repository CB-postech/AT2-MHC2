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

########## GSE155900

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

######## GSE260769


file_path="C:/Users/owner/Desktop/MJY/GSE260769"
setwd(file_path)

######### 01. Read Raw data and create seurat object ######################
mydata_AL4 <- Read10X("data/AL4/")
colnames(mydata_AL4) <- paste0(substring(colnames(mydata_AL4),1,17), "AL4")
mydata_AL5 <- Read10X("data/AL5/")
colnames(mydata_AL5) <- paste0(substring(colnames(mydata_AL5),1,17), "AL5")
mydata_AL6 <- Read10X("data/AL6/")
colnames(mydata_AL6) <- paste0(substring(colnames(mydata_AL6),1,17), "AL6")
mydata_AL7 <- Read10X("data/AL7/")
colnames(mydata_AL7) <- paste0(substring(colnames(mydata_AL7),1,17), "AL7")
mydata_PT33 <- Read10X("data/PT33/")
colnames(mydata_PT33) <- paste0(substring(colnames(mydata_PT33),1,17), "PT33")
mydata_PT36 <- Read10X("data/PT36/")
colnames(mydata_PT36) <- paste0(substring(colnames(mydata_PT36),1,17), "PT36")
mydata_PT50 <- Read10X("data/PT50/")
colnames(mydata_PT50) <- paste0(substring(colnames(mydata_PT50),1,17), "PT50")
mydata_PT51 <- Read10X("data/PT51/")
colnames(mydata_PT51) <- paste0(substring(colnames(mydata_PT51),1,17), "PT51")
mydata_PT66 <- Read10X("data/PT66/")
colnames(mydata_PT66) <- paste0(substring(colnames(mydata_PT66),1,17), "PT66")
mydata_PT106 <- Read10X("data/PT106/")
colnames(mydata_PT106) <- paste0(substring(colnames(mydata_PT106),1,17), "PT106")

# samples <- list(mydata_AL4, mydata_AL5, mydata_AL6, mydata_AL7, 
#                 mydata_PT33, mydata_PT36, mydata_PT50, 
#                 mydata_PT51, mydata_PT66, mydata_PT106)

mydata <- cbind(mydata_AL4,mydata_AL5,mydata_AL6,mydata_AL7,mydata_PT33,mydata_PT36,mydata_PT50, mydata_PT51,mydata_PT66,mydata_PT106)

e.out <- emptyDrops(mydata)
myseurat <- CreateSeuratObject(counts = mydata)

myseurat$sample <- ""
myseurat$sample <- substring(colnames(myseurat), 18)

head(myseurat@assays)
head(myseurat)
tail(myseurat)
unique(myseurat$sample)

myseurat$patient_info <- as.character(myseurat$sample)
Idents(myseurat) <- "patient_info"
myseurat$patient_info[WhichCells(myseurat, idents='AL4', slot='counts')] <- 'A_47y'
myseurat$patient_info[WhichCells(myseurat, idents='AL5', slot='counts')] <- 'A_74y'
myseurat$patient_info[WhichCells(myseurat, idents='AL6', slot='counts')] <- 'A_54y'
myseurat$patient_info[WhichCells(myseurat, idents='AL7', slot='counts')] <- 'A_67y'
myseurat$patient_info[WhichCells(myseurat, idents='PT106', slot='counts')] <- 'P_10w'
myseurat$patient_info[WhichCells(myseurat, idents='PT33', slot='counts')] <- 'P_12w'
myseurat$patient_info[WhichCells(myseurat, idents='PT51', slot='counts')] <- 'P_17w'
myseurat$patient_info[WhichCells(myseurat, idents='PT66', slot='counts')] <- 'P_20w'
myseurat$patient_info[WhichCells(myseurat, idents='PT36', slot='counts')] <- 'P_25w'
myseurat$patient_info[WhichCells(myseurat, idents='PT50', slot='counts')] <- 'P_40w'

order <- c('P_10w', 'P_12w', 'P_17w', 'P_20w', 'P_25w', 'P_40w', 'A_47y', 'A_54y', 'A_67y', 'A_74y')
myseurat$patient_info <- factor(myseurat$patient_info, levels=order)

unique(myseurat$patient_info)

tail(myseurat$patient_info)

myseurat$percent.mt = PercentageFeatureSet(myseurat, pattern="^MT-")

######### 02. QC (filtered out nFeature_RNA<200, percent.mt > 10 by reference) ######################

myseurat_qc <- subset(myseurat, subset = nFeature_RNA > 200 & percent.mt < 10)

###Normalization (expression values were then scaled to 10000 transcripts per cell and log-transformed)
myseurat_qc <- NormalizeData(myseurat_qc, normalization.method = "LogNormalize", scale.factor = 10000)

### feature selection; Find HVGs (Used 2000 genes with high cell-to-cell variation, which were calculated using the FindVariableFeatures function in seurat for further dimensionality reduction)
myseurat_qc <- FindVariableFeatures(myseurat_qc) # Default nFeatures = 2000

### save RDS 
saveRDS(myseurat_qc, file = "GSE260769_qc.rds")


######### 03. Pre-processing (Normalization, Scaling, Dimension reduction) ######################
seurat1 <- readRDS("GSE260769_qc.rds")
head(seurat1@meta.data)

order <- c('P_10w', 'P_12w', 'P_17w', 'P_20w', 'P_25w', 'P_40w', 'A_47y', 'A_54y', 'A_67y', 'A_74y')
seurat1$patient_info <- factor(seurat1$patient_info, levels=order)

#VlnPlot(seurat1, features = c("nCount_RNA"), group.by = "patient_info", pt.size = 0.1)+
  theme_minimal()

#VlnPlot(seurat1, features = c("percent.mt"), group.by = "patient_info", pt.size = 0.1)+
  theme_minimal()

##Normalization (expression values were then scaled to 10000 transcripts per cell and log-transformed)
seurat1 <- NormalizeData(seurat1)

## feature selection; Find HVGs (Used 2000 genes with high cell-to-cell variation, which were calculated using the FindVariableFeatures function in seurat for further dimensionality reduction)
seurat1 <- FindVariableFeatures(seurat1) # Default nFeatures = 2000

###Scaling; z-transformation (Effects of variable(percent.mt) were estimated and regressed out using a GLM (ScaleData function, model.use="linear"))
seurat1 <- ScaleData(seurat1, vars.to.regress = "percent.mt", model.use = "linear")

###Dimension reduction
seurat1 <- RunPCA(seurat1, features=VariableFeatures(object = seurat1))
ElbowPlot(seurat1, ndims=50)

n_pcs=20; cluster_resolution =0.8

seurat1 <- FindNeighbors(seurat1, dims=1:n_pcs)
seurat1 <- FindClusters(seurat1, resolution = cluster_resolution)
seurat1 <- RunUMAP(seurat1, dims=1:n_pcs)

DimPlot(seurat1, label=T)

DimPlot(seurat1, label = T, group.by = "patient_info")
FeaturePlot(seurat1, features = c('EPCAM','SOX2', 'SOX9', 'SFTPC', 'SCGB1A1', 'AGER'))
FeaturePlot(seurat1, features = c('LYZ','COL1A2','ACTA2', 'CD3D','CD79A','PECAM1', 'ID3','HBG1'))

######### 04. Subsetting Epithelial cells ######################
head(seurat1@meta.data)

Epithelial_clusters <- c("13", "15", "35", "31", "28", "33", "21", "8", "18", "4", "23")
Idents(seurat1) <- "seurat_clusters"

Epi_seurat <- subset(seurat1, idents = Epithelial_clusters)
DimPlot(Epi_seurat, label = T, split.by = "patient_info")

saveRDS(Epi_seurat, file = "Epi_seurat_not normalized.rds")

### Normalization
Epi_seurat <- NormalizeData(Epi_seurat)

## feature selection; Find HVGs (Used 2000 genes with high cell-to-cell variation, which were calculated using the FindVariableFeatures function in seurat for further dimensionality reduction)
Epi_seurat <- FindVariableFeatures(Epi_seurat) # Default nFeatures = 2000

###Scaling; z-transformation (Effects of variable(percent.mt) were estimated and regressed out using a GLM (ScaleData function, model.use="linear"))
Epi_seurat <- ScaleData(Epi_seurat, vars.to.regress = "percent.mt", model.use = "linear")

###Dimension reduction
Epi_seurat <- RunPCA(Epi_seurat, features=VariableFeatures(object = Epi_seurat))
ElbowPlot(Epi_seurat, ndims=50)

n_pcs=20; cluster_resolution =0.4

Epi_seurat <- FindNeighbors(Epi_seurat, dims=1:n_pcs)
Epi_seurat <- FindClusters(Epi_seurat, resolution = cluster_resolution)
Epi_seurat <- RunUMAP(Epi_seurat, dims=1:n_pcs)

DimPlot(Epi_seurat, label=T)
DimPlot(Epi_seurat, label=T, split.by = "patient_info")

FeaturePlot(Epi_seurat, features = c('SOX2', 'SOX9', 'SFTPB', 'AGER', 'SCGB1A1', 'FOXJ1'))

saveRDS(Epi_seurat, file="Epi_seurat_normed.rds")

######### 05. subset annotation in each timepoint ######################
Idents(Epi_seurat) <- "patient_info"

pseudoglandular <- subset(Epi_seurat, ident = c("P_10w", "P_12w"))
canalicular <- subset(Epi_seurat, ident = c("P_17w", "P_20w", "P_25w"))
Saccular <- subset(Epi_seurat, ident = c("P_40w"))
Adult <- subset(Epi_seurat, ident = c("A_47y", "A_54y", "A_67y", "A_74y"))

###### pseudoglandular
pseudoglandular <- FindVariableFeatures(pseudoglandular, nfeatures = 2000)
pseudoglandular <- ScaleData(pseudoglandular)
pseudoglandular <- RunPCA(pseudoglandular, features = VariableFeatures(object = pseudoglandular))
ElbowPlot(pseudoglandular, ndims = 50)

n_pcs = 20; cluster_resolution=0.4
pseudoglandular <- FindNeighbors(pseudoglandular, dims = 1:n_pcs)
pseudoglandular <- FindClusters(pseudoglandular, resolution = cluster_resolution)
pseudoglandular <- RunUMAP(pseudoglandular, dims = 1:n_pcs)

DimPlot(pseudoglandular, label=T)
FeaturePlot(pseudoglandular, features = c('SOX2', 'SOX9', 'SFTPB', 'AGER', 'SCGB1A1', 'FOXJ1'))

FeaturePlot(pseudoglandular, features = c('SOX9', 'ID2', 'NKX2-1', 'FOXP2', 'HMGA2')) # undifferentiated lung progenitors
FeaturePlot(pseudoglandular, features = c('SOX2', 'SCGB3A2')) #proximal progenitors
FeaturePlot(pseudoglandular, features = c('SOX9', 'ID2', 'SFTPC')) #Distal progenitors
FeaturePlot(pseudoglandular, features = c('SFTPC', 'SFTPB', 'LAMP3', 'AGER')) #Early alveolar epithelial cells

DimPlot(pseudoglandular, label = T, group.by = "patient_info")

pseudoglandular$annot <- as.character(pseudoglandular$seurat_clusters)
Idents(pseudoglandular) <- "annot"
pseudoglandular$annot[WhichCells(pseudoglandular, idents = c("1", "7"))] <- 'Undifferentiated lung progenitors'
pseudoglandular$annot[WhichCells(pseudoglandular, idents = c("0", "5"))] <- 'Alveolar epithelial progenitors'
pseudoglandular$annot[WhichCells(pseudoglandular, idents = c("2", "9"))] <- 'Distal progenitors'
pseudoglandular$annot[WhichCells(pseudoglandular, idents = c("3", "4", "6"))] <- 'proximal progenitors'
pseudoglandular$annot[WhichCells(pseudoglandular, idents = c("8","10", "11"))] <- 'Unknown'

DimPlot(pseudoglandular, label = T)

###### Canalicular
canalicular <- FindVariableFeatures(canalicular, nfeatures = 2000)
canalicular <- ScaleData(canalicular)
canalicular <- RunPCA(canalicular, features = VariableFeatures(object = canalicular))
ElbowPlot(pseudoglandular, ndims = 50)

n_pcs = 20; cluster_resolution=0.4
canalicular <- FindNeighbors(canalicular, dims = 1:n_pcs)
canalicular <- FindClusters(canalicular, resolution = cluster_resolution)
canalicular <- RunUMAP(canalicular, dims = 1:n_pcs)

DimPlot(canalicular, label=T)
DimPlot(canalicular, label=T, group.by = "patient_info")

FeaturePlot(canalicular, features = c('SOX9', 'SFTPC', 'LAMP3')) # Alveolar epithelial progenitors
FeaturePlot(canalicular, features = c('AGER', 'PDPN', 'AQP5')) #Early AT1
FeaturePlot(canalicular, features = c('SFTPC', 'SFTPB', 'LAMP3', 'ABCA3')) #Early AT2
FeaturePlot(canalicular, features = c('SCGB3A2', 'FOXJ1', 'TP63')) #Airway epithelial cells

DimPlot(canalicular, label = T, group.by = "patient_info")

canalicular$annot <- as.character(canalicular$seurat_clusters)
Idents(canalicular) <- "annot"
canalicular$annot[WhichCells(canalicular, idents = c("3"))] <- 'Early AT1'
canalicular$annot[WhichCells(canalicular, idents = c("0"))] <- 'Early AT2'
canalicular$annot[WhichCells(canalicular, idents = c("2", "4"))] <- 'Airway epithelial cells'
canalicular$annot[WhichCells(canalicular, idents = c("6", "1", "5"))] <- 'Alveolar epithelial progenitors'

DimPlot(canalicular, label = T)

###### Saccular
Saccular <- FindVariableFeatures(Saccular, nfeatures = 2000)
Saccular <- ScaleData(Saccular)
Saccular <- RunPCA(Saccular, features = VariableFeatures(object = Saccular))
ElbowPlot(Saccular, ndims = 50)

n_pcs = 20; cluster_resolution=0.8
Saccular <- FindNeighbors(Saccular, dims = 1:n_pcs)
Saccular <- FindClusters(Saccular, resolution = cluster_resolution)
Saccular <- RunUMAP(Saccular, dims = 1:n_pcs)

DimPlot(Saccular, label=T)
DimPlot(Saccular, label=T, group.by = "patient_info")

FeaturePlot(Saccular, features = c('SFTPC', 'SFTPB', 'LAMP3', 'ABCA3')) #AT2
FeaturePlot(Saccular, features = c('SCGB3A2', 'FOXJ1', 'TP63', 'SCGB1A1', 'DNAI1')) #Airway epithelial cells

DimPlot(Saccular, label = T, group.by = "patient_info")

Saccular$annot <- as.character(Saccular$seurat_clusters)
Idents(Saccular) <- "annot"

Saccular$annot[WhichCells(Saccular, idents = c("2"))] <- 'AT2'
Saccular$annot[WhichCells(Saccular, idents = c("0", "1"))] <- 'Airway epithelial cells'

DimPlot(Saccular, label = T)


###### Adult
Adult <- FindVariableFeatures(Adult, nfeatures = 2000)
Adult <- ScaleData(Adult)
Adult <- RunPCA(Adult, features = VariableFeatures(object = Adult))
ElbowPlot(Adult, ndims = 50)

n_pcs = 20; cluster_resolution=0.4
Adult <- FindNeighbors(Adult, dims = 1:n_pcs)
Adult <- FindClusters(Adult, resolution = cluster_resolution)
Adult <- RunUMAP(Adult, dims = 1:n_pcs)

DimPlot(Adult, label=T)
DimPlot(Adult, label=T, group.by = "patient_info")

FeaturePlot(Adult, features = c('SFTPC', 'SFTPB', 'LAMP3', 'ABCA3')) #AT2
FeaturePlot(Adult, features = c('SCGB3A2', 'FOXJ1', 'TP63', 'SCGB1A1', 'DNAI1')) #Airway epithelial cells
FeaturePlot(Adult, features = c('AGER', 'PDPN', 'HOPX')) #AT1


Adult$annot <- as.character(Adult$seurat_clusters)
Idents(Adult) <- "annot"

Adult$annot[WhichCells(Adult, idents = c("6"))] <- 'AT2'
Adult$annot[WhichCells(Adult, idents = c("7"))] <- 'AT1'
Adult$annot[WhichCells(Adult, idents = c("0","1", "2", "3", "4", "5","8","9","10","11"))] <- 'Airway epithelial cells'

DimPlot(Adult, label = T)

######### 06. Merge every timepoints into one seurat object of GSE260769 ######################
merged_seurat_human <- merge(pseudoglandular, y=list(canalicular, Saccular, Adult))

merged_seurat_human$stage <- as.character(merged_seurat_human$patient_info)
Idents(merged_seurat_human) <- "stage"
merged_seurat_human$stage[WhichCells(merged_seurat_human, idents = c('P_10w', 'P_12w'))] <- '10-12pcw'
merged_seurat_human$stage[WhichCells(merged_seurat_human, idents = c('P_17w', 'P_20w', 'P_25w'))] <- '17-25pcw'
merged_seurat_human$stage[WhichCells(merged_seurat_human, idents = c('P_40w'))] <- '40pcw'
merged_seurat_human$stage[WhichCells(merged_seurat_human, idents = c('A_47y', 'A_54y', 'A_67y', 'A_74y'))] <- 'Adult'

merged_seurat_human <- FindVariableFeatures(merged_seurat_human, nfeatures = 2000)
merged_seurat_human <- ScaleData(merged_seurat_human)
merged_seurat_human <- RunPCA(merged_seurat_human, features = VariableFeatures(object = merged_seurat_human))
ElbowPlot(merged_seurat_human, ndims = 50)

n_pcs = 20; cluster_resolution=0.4
merged_seurat_human <- FindNeighbors(merged_seurat_human, dims = 1:n_pcs)
merged_seurat_human <- FindClusters(merged_seurat_human, resolution = cluster_resolution)
merged_seurat_human <- RunUMAP(merged_seurat_human, dims = 1:n_pcs)

DimPlot(merged_seurat_human, label=T)
DimPlot(merged_seurat_human, label = T, group.by = "annot")
DimPlot(merged_seurat_human, label = T, group.by = "patient_info") # Batch correction is not needed
DimPlot(merged_seurat_human, label = T, group.by = "stage")

saveRDS(merged_seurat_human, file = "Epi_subset_normed_annot_GSE260769.rds")
