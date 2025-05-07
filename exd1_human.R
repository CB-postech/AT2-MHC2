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


file_path="D:/MJY/Human MHCII_merged"
setwd(file_path)

############## 01. Read each pre-processed seurat object and create merged seurat object ####################

GSE155900_Epi <- readRDS("harmony_Epi_merged_annot_GSE155900.rds")
GSE260769_Epi <- readRDS("Epi_subset_normed_annot_GSE260769.rds")

colnames(GSE155900_Epi)
colnames(GSE260769_Epi)
DimPlot(GSE260769_Epi, label = T, group.by = "annot")
DimPlot(GSE155900_Epi, label = T, group.by = "annot")

merged_human <- merge(GSE260769_Epi, y=GSE155900_Epi)
head(merged_human@meta.data)
unique(merged_human$stage)

############## 02. Pre-processing (Normalization, Scaling, Dimension reduction) ####################

merged_human <- NormalizeData(merged_human)
merged_human <- FindVariableFeatures(merged_human)
merged_human <- ScaleData(merged_human)

###Dimension reduction
merged_human <- RunPCA(merged_human, features=VariableFeatures(object = merged_human))
ElbowPlot(merged_human, ndims=50)

n_pcs=20; cluster_resolution =0.8

merged_human <- FindNeighbors(merged_human, dims=1:n_pcs)
merged_human <- FindClusters(merged_human, resolution = cluster_resolution)
merged_human <- RunUMAP(merged_human, dims=1:n_pcs)

DimPlot(merged_human, label=T, group.by = "annot")
DimPlot(merged_human, label=T, group.by = "sample")
DimPlot(merged_human, label=T, group.by = "stage")

######### Batch correction (harmony)

harmony_merged_human = merged_human
harmony_merged_human <- RunHarmony(harmony_merged_human, "sample", max.iter = 30)
ElbowPlot(harmony_merged_human, reduction = "harmony", ndims = 30)

harmony.dims = 10

harmony_merged_human <- harmony_merged_human %>%
  RunUMAP(reduction = "harmony", dims = 1:harmony.dims) %>%
  FindNeighbors(reduction = "harmony", dims = 1:harmony.dims) %>%
  FindClusters()

DimPlot(harmony_merged_human, label = T, group.by = "annot")
DimPlot(harmony_merged_human, label = T, group.by = "stage")

ordered <- c('10-12pcw','17-25pcw','40pcw', '< 1yr', '2-3yr', '3-4yr', 'Adult')
harmony_merged_human$stage <- factor(harmony_merged_human$stage, levels = ordered)
DimPlot(harmony_merged_human, label = T, split.by = "stage")

saveRDS(harmony_merged_human, "harmony_merged_human.rds")

############## 03. AT2 subsetting ####################
Idents(harmony_merged_human) <- "annot"
AT2_seurat_human <- subset(harmony_merged_human, subset = annot %in% "AT2")

DimPlot(AT2_seurat_human, label = T)
print(table(AT2_seurat_human$stage))

############## 04. Calculate MHCII protein complex module score by using GOCC_MHCII Protein complex  ####################

ordered <- c('40pcw', '< 1yr', '2-3yr', '3-4yr', 'Adult')
AT2_seurat_human$stage <- factor(AT2_seurat_human$stage, levels=ordered)

gene_set_file <- "GOCC_MHC_CLASS_II_PROTEIN_COMPLEX.v2024.1.Hs.tsv"  
gene_set_data <- read.delim(gene_set_file)

gene_symbols_raw <- gene_set_data$GOCC_MHC_CLASS_II_PROTEIN_COMPLEX

gocc_mhc_class_ii_genes <- strsplit(gene_symbols_raw, split = ",")[[17]]

AT2_seurat_human <- JoinLayers(AT2_seurat_human)

AT2_seurat_human <- AddModuleScore(
  AT2_seurat_human, 
  features = list(gocc_mhc_class_ii_genes), 
  name = "GOCC_MHC_CLASS_II_Protein_Complex"
)

VlnPlot(AT2_seurat_human, 
        features = "GOCC_MHC_CLASS_II_Protein_Complex1", 
        group.by = "stage",
        pt.size = 0) + 
  ggtitle("GOCC MHC Class II Protein Complex Signature Score by Stage")

############## 05. AT2 marker selection (via marker gene selection) and calculate module score  ####################

unique(harmony_merged_human$annot)
Idents(harmony_merged_human) <- "annot"
AT2_markers <- FindMarkers(harmony_merged_human, ident.1 = "AT2", min.pct = 0.10)

AT2_markers_sig_FC_2.0 <- AT2_markers[AT2_markers$avg_log2FC > 2 & AT2_markers$p_val_adj < 0.05 & AT2_markers$pct.1 > 0.5,]
AT2_markers_sig_FC_2.0 <- rownames(AT2_markers_sig_FC_2.0)

### Trouble shooting
#AT2_markers_sig_FC_1.5 <- AT2_markers[AT2_markers$avg_log2FC > 1.5 & AT2_markers$p_val_adj < 0.05 & AT2_markers$pct.1 > 0.5,]
#AT2_markers_sig_FC_0.58 <- AT2_markers[AT2_markers$avg_log2FC > 0.58 & AT2_markers$p_val_adj < 0.05 & AT2_markers$pct.1 > 0.5,]
#AT2_markers_sig_FC_1.0 <- AT2_markers[AT2_markers$avg_log2FC > 1.0 & AT2_markers$p_val_adj < 0.05 & AT2_markers$pct.1 > 0.5,]

#AT2_markers_sig_FC_1.5 <- rownames(AT2_markers_sig_FC_1.5)
#AT2_markers_sig_FC_0.58<- rownames(AT2_markers_sig_FC_0.58)
#AT2_markers_sig_FC_1.0 <- rownames(AT2_markers_sig_FC_1.0)

ordered <- c('40pcw', '< 1yr', '2-3yr', '3-4yr', 'Adult')
AT2_seurat_human$stage <- factor(AT2_seurat_human$stage, levels=ordered)

AT2_seurat_human <- AddModuleScore(
  AT2_seurat_human, 
  features = list(AT2_markers_sig_FC_2.0), 
  name = "AT2_markers_sig_FC_2.0_sets"
)

VlnPlot(AT2_seurat_human, 
        features = "AT2_markers_sig_FC_2.0_sets1", 
        group.by = "stage",
        pt.size = 0) + 
  ggtitle("AT2_marker (log2FC > 2.0) Signature Score by Stage")
