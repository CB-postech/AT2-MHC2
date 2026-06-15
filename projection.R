### this script is for
### run Scran to normalize the data
### conda activate project_lung_exercise_R
# based on https://github.com/CB-postech/2022_KOGO_workshop/blob/main/KOGO_LCB_pipeline.md

library(Seurat)
library(scater)
library(SingleCellExperiment)
library(magrittr)
library(ggplot2)
library(BiocParallel)
library(future)
library(stringr)
library(dplyr)

source('/home/sjcho/yard/functions/R/seurat_count_to_normalization.R')
source('/home/sjcho/yard/functions/R/FeaturePlot_sjcho.R')
source('/home/sjcho/yard/functions/R/save_ggplot2_to_ppt.R')
source('/home/sjcho/projects/utils.R')

save_path = '/home/sjcho/projects/AT2_MHC2/code_for_publish/251208_revision/projection/outs/for_publication'
so.alv <- readRDS('/home/sjcho/projects/AT2_MHC2/20241113_after_cellbender/3.whole_cell_annotation/outs/3.1.epithelial_annotation/alv_annotation.rds')

so.alv$condition_dpi = factor(so.alv$condition_dpi, levels = c('naive_floxed', '7dpi_floxed', '14dpi_floxed', '30dpi_floxed', 
                                                                'naive_dAT2', '7dpi_dAT2', '14dpi_dAT2', '30dpi_dAT2'))

alv_cols = c(AT2 = "#8192ef", tAT2 = "#dfd13a", AT1 = "#f66b07", AT2.IFN = "#000000", AT2.prolif. = "#8bd690")
### 1. celltypist annotation
celltypist_result = fread('/home/sjcho/projects/AT2_MHC2/20241113_after_cellbender/manuscript_figures/figure4_and_exd6_20250225/alv_full_celltypist_result_overcluster1.csv')
celltypist_result = as.data.frame(celltypist_result)
rownames(celltypist_result) <- celltypist_result$V1
celltypist_result <- celltypist_result[, -1]

celltypist_result[celltypist_result[, 'predicted_labels'] == 'cycling.alv.epi', 'predicted_labels'] = 'AT2.prolif.'
celltypist_result[celltypist_result[, 'predicted_labels'] == 'IFN.response.alv.epi', 'predicted_labels'] = 'AT2.IFN'
celltypist_result[celltypist_result[, 'majority_voting'] == 'cycling.alv.epi', 'majority_voting'] = 'AT2.prolif.'
celltypist_result[celltypist_result[, 'majority_voting'] == 'IFN.response.alv.epi', 'majority_voting'] = 'AT2.IFN'

so.alv$predicted_labels = celltypist_result[Cells(so.alv), 'predicted_labels']
so.alv$predicted_labels = factor(so.alv$predicted_labels, levels = c('AT2', 'tAT2', 'AT1', 'AT2.prolif.', 'AT2.IFN'))
so.alv$majority_voting = celltypist_result[Cells(so.alv), 'majority_voting']
so.alv$majority_voting = factor(so.alv$majority_voting, levels = c('AT2', 'tAT2', 'AT1', 'AT2.prolif.', 'AT2.IFN')) 

p <- DimPlot(so.alv, group.by = 'majority_voting', cols = alv_cols)
ggsave(p, file = paste0(save_path, '/alv_dimplot.png'), width = 7, height = 6)

p <- fp_sjcho(subset(so.alv, condition == 'floxed'), features = c('Cd74', 'H2-Aa', 'H2-Ab1'), ncol = 3) & scale_colour_gradientn(colours  = rev(brewer.pal(11, "RdBu")))
ggsave(p, file = paste0(save_path, '/alv_Cd74_H2-Aa_H2-Ab1_floxed.png'), width = 9, height = 3)

################ Projection

###### A. MAER
###### 1. select hvg
query_obj <- readRDS('/home/sjcho/datas/reference_atlas/mouse_lung_alveolar/outs/1.3.a.alveolar_normalization_GSE2632927/merged_w_GSE262927_full.rds')
query_obj <- subset(query_obj, annotation == 'tAT2')
# so.262 <- readRDS('/home/sjcho/datas/public_data/GSE262927/GSE262927_rawX_symbol.rds')

hvg <- VariableFeatures(so.alv) 
hvg_intersected <- intersect(hvg, rownames(query_obj))
ref_mat <- as.matrix(so.alv[["RNA"]]$data[hvg_intersected, ])
query_mat <- as.matrix(query_obj[["RNA"]]$data[hvg_intersected, ])

###### 2. run knn to find nearest neighbors
knn_result <- KernelKnn::knn.index.dist(
  t(ref_mat), 
  t(query_mat), 
  k = 15, 
  method = "pearson_correlation"
)
knn_idx <- knn_result$test_knn_idx 

###### 3. calculate projection embeddings
ref_embeddings <- Embeddings(so.alv, reduction = "umap")

# 이웃들의 좌표 평균 계산
prj_coords <- t(apply(knn_idx, 1, function(x) colMeans(ref_embeddings[x, ])))
colnames(prj_coords) <- c("Prj_1", "Prj_2")
rownames(prj_coords) <- colnames(Cells(query_obj))

###### 4. visualization
df_ref <- as.data.frame(Embeddings(so.alv, reduction = "umap"))
colnames(df_ref) <- c("UMAP_1", "UMAP_2")
df_ref$Dataset <- "Reference"
df_ref$Group <- "Reference" 

# 2. Query 
df_query <- as.data.frame(prj_coords)
colnames(df_query) <- c("UMAP_1", "UMAP_2")
df_query$Dataset <- "Query"
df_query$Group <- "tAT2"

# 3. merge
df_total <- rbind(
  df_ref %>% mutate(Group = "Reference"),
  df_query
)

library(ggplot2)
library(MASS)
dens <- kde2d(df_query$UMAP_1, df_query$UMAP_2, n = 200)
df_query$density <- fields::interp.surface(
    list(x = dens$x, y = dens$y, z = dens$z),
    cbind(df_query$UMAP_1, df_query$UMAP_2)
)
df_query$density <- df_query$density / max(df_query$density, na.rm = TRUE)
df_query <- df_query[order(df_query$density), ]

p <- ggplot() +
    geom_point(data = df_ref,
               aes(x = UMAP_1, y = UMAP_2),
               color = "grey75", size = 0.4, alpha = 1) +
    geom_point(data = df_query,
               aes(x = UMAP_1, y = UMAP_2, color = density),
               size = 0.4, alpha = 1) +
    scale_color_distiller(palette = "Reds", direction = 1, name = "Density",
                          limits = c(0, NA)) +
    theme_bw() +
    labs(x = "UMAP1", y = "UMAP2") +
    theme(
        panel.grid = element_blank(),
        aspect.ratio = 1,
        axis.text = element_blank(),
        axis.ticks = element_blank()
    )
ggsave(p, file = paste0(save_path, '/projection_MAER_into_Inhouse_Reds.png'), width = 5, height = 4)
ggsave(p, file = paste0(save_path, '/projection_MAER_into_Inhouse_Reds.pdf'), width = 5, height = 4)

p <- ggplot() +
    geom_point(data = df_ref, 
                aes(x = UMAP_1, y = UMAP_2), 
               color = "grey75", size = 0.4, alpha = 1) +
    geom_point(data = df_query, 
                aes(x = UMAP_1, y = UMAP_2, color = "darkred"), 
                size = 0.4, alpha = 1) +
    theme_bw() +
    labs(x = "UMAP1", y = "UMAP2") +
    theme(
    panel.grid = element_blank(),
    aspect.ratio = 1,
    axis.text = element_blank(),
    axis.ticks = element_blank()
    ) + NoLegend()
ggsave(p, file = paste0(save_path, '/projection_of_query_onto_alveolar_reference_umap.png'), width = 5, height = 5)
ggsave(p, file = paste0(save_path, '/projection_of_query_onto_alveolar_reference_umap.pdf'), width = 5, height = 5)

###### B. GSE262927
###### 1. select hvg
query_obj <- readRDS('/home/sjcho/datas/reference_atlas/mouse_lung_alveolar/outs/1.3.a.alveolar_normalization_GSE2632927/merged_w_GSE262927_full.rds')
query_obj <- readRDS('/home/sjcho/datas/public_data/GSE262927/GSE262927_rawX_symbol.rds')
query_obj <- subset(query_obj, subset = author_cell_type == 'Alveolar_transitional')

hvg <- VariableFeatures(so.alv) 
hvg_intersected <- intersect(hvg, rownames(query_obj))
ref_mat <- as.matrix(so.alv[["RNA"]]$data[hvg_intersected, ])
query_mat <- as.matrix(query_obj[["RNA"]]$data[hvg_intersected, ])

###### 2. run knn to find nearest neighbors
knn_result <- KernelKnn::knn.index.dist(
  t(ref_mat), 
  t(query_mat), 
  k = 15, 
  method = "pearson_correlation"
)
knn_idx <- knn_result$test_knn_idx 

###### 3. calculate projection embeddings
ref_embeddings <- Embeddings(so.alv, reduction = "umap")

prj_coords <- t(apply(knn_idx, 1, function(x) colMeans(ref_embeddings[x, ])))
colnames(prj_coords) <- c("Prj_1", "Prj_2")
rownames(prj_coords) <- colnames(Cells(query_obj))

###### 4. visualization
# 1. Reference data frame
df_ref <- as.data.frame(Embeddings(so.alv, reduction = "umap"))
colnames(df_ref) <- c("UMAP_1", "UMAP_2")
df_ref$Dataset <- "Reference"
df_ref$Group <- "Reference" 

# 2. Query dataframe
df_query <- as.data.frame(prj_coords)
colnames(df_query) <- c("UMAP_1", "UMAP_2")
df_query$Dataset <- "Query"
df_query$Group <- "Alveolar_transitional"

# 3. merge
df_total <- rbind(
  df_ref %>% mutate(Group = "Reference"),
  df_query
)

library(ggplot2)
library(MASS)
dens <- kde2d(df_query$UMAP_1, df_query$UMAP_2, n = 200)
df_query$density <- fields::interp.surface(
    list(x = dens$x, y = dens$y, z = dens$z),
    cbind(df_query$UMAP_1, df_query$UMAP_2)
)
df_query$density <- df_query$density / max(df_query$density, na.rm = TRUE)
df_query <- df_query[order(df_query$density), ]

p <- ggplot() +
    geom_point(data = df_ref,
               aes(x = UMAP_1, y = UMAP_2),
               color = "grey75", size = 0.4, alpha = 1) +
    geom_point(data = df_query,
               aes(x = UMAP_1, y = UMAP_2, color = density),
               size = 0.4, alpha = 1) +
    scale_color_distiller(palette = "Reds", direction = 1, name = "Density",
                          limits = c(0, NA)) +
    theme_bw() +
    labs(x = "UMAP1", y = "UMAP2") +
    theme(
        panel.grid = element_blank(),
        aspect.ratio = 1,
        axis.text = element_blank(),
        axis.ticks = element_blank()
    )
ggsave(p, file = paste0(save_path, '/projection_GSE262927_into_Inhouse_Reds.png'), width = 5, height = 4)
ggsave(p, file = paste0(save_path, '/projection_GSE262927_into_Inhouse_Reds.pdf'), width = 5, height = 4)

p <- ggplot() +
    geom_point(data = df_ref, 
                aes(x = UMAP_1, y = UMAP_2), 
               color = "grey75", size = 0.4, alpha = 1) +
    geom_point(data = df_query, 
                aes(x = UMAP_1, y = UMAP_2, color = "darkred"), 
                size = 0.4, alpha = 1) +
    theme_bw() +
    labs(x = "UMAP1", y = "UMAP2") +
    theme(
    panel.grid = element_blank(),
    aspect.ratio = 1,
    axis.text = element_blank(),
    axis.ticks = element_blank()
    ) + NoLegend()
ggsave(p, file = paste0(save_path, '/projection_of_query_onto_GSE262927_umap.png'), width = 5, height = 5)
ggsave(p, file = paste0(save_path, '/projection_of_query_onto_GSE262927_umap.pdf'), width = 5, height = 5)

###### C. GSE141259
###### 1. select hvg
mtx <- Matrix::readMM("/home/sjcho/datas/public_data/GSE141259/GSE141259_WholeLung_rawcounts.mtx")
genes <- read.table('/home/sjcho/datas/public_data/GSE141259/GSE141259_WholeLung_genes.txt', stringsAsFactors = F)
rownames(mtx) <- genes$V1

cell_info <- read.csv('/home/sjcho/datas/public_data/GSE141259/GSE141259_WholeLung_cellinfo.csv', row.names = 1)
colnames(mtx) <- rownames(cell_info)

query_obj <- CreateSeuratObject(
  counts = mtx,
  meta.data = cell_info,
  min.cells = 0,
  min.features = 0
)
query_obj <- subset(query_obj, cells = Cells(query_obj)[query_obj$'cell.type' == 'Krt8 ADI'])
query_obj <- log_normalize(query_obj, save_path, data_name = 'GSE141259_ADI', pcs = 10, nfeatures = 2000)

hvg <- VariableFeatures(so.alv) 
hvg_intersected <- intersect(hvg, rownames(query_obj))
ref_mat <- as.matrix(so.alv[["RNA"]]$data[hvg_intersected, ])
query_mat <- as.matrix(query_obj[["RNA"]]$data[hvg_intersected, ])

###### 2. run knn to find nearest neighbors
knn_result <- KernelKnn::knn.index.dist(
  t(ref_mat), 
  t(query_mat), 
  k = 15, 
  method = "pearson_correlation"
)
knn_idx <- knn_result$test_knn_idx

###### 3. calculate projection embeddings
# Reference
ref_embeddings <- Embeddings(so.alv, reduction = "umap")

prj_coords <- t(apply(knn_idx, 1, function(x) colMeans(ref_embeddings[x, ])))
colnames(prj_coords) <- c("Prj_1", "Prj_2")
rownames(prj_coords) <- colnames(Cells(query_obj))

###### 4. visualization
# 1. Reference dataframe
df_ref <- as.data.frame(Embeddings(so.alv, reduction = "umap"))
colnames(df_ref) <- c("UMAP_1", "UMAP_2")
df_ref$Dataset <- "Reference"
df_ref$Group <- "Reference"

# 2. Query dataframe
df_query <- as.data.frame(prj_coords)
colnames(df_query) <- c("UMAP_1", "UMAP_2")
df_query$Dataset <- "Query"
df_query$Group <- "Alveolar_transitional"

# 3. data merge
df_total <- rbind(
  df_ref %>% mutate(Group = "Reference"),
  df_query
)

library(ggplot2)
library(MASS)
dens <- kde2d(df_query$UMAP_1, df_query$UMAP_2, n = 200)
df_query$density <- fields::interp.surface(
    list(x = dens$x, y = dens$y, z = dens$z),
    cbind(df_query$UMAP_1, df_query$UMAP_2)
)
df_query$density <- df_query$density / max(df_query$density, na.rm = TRUE)
df_query <- df_query[order(df_query$density), ]

p <- ggplot() +
    geom_point(data = df_ref,
               aes(x = UMAP_1, y = UMAP_2),
               color = "grey75", size = 0.4, alpha = 1) +
    geom_point(data = df_query,
               aes(x = UMAP_1, y = UMAP_2, color = density),
               size = 0.4, alpha = 1) +
    scale_color_distiller(palette = "Reds", direction = 1, name = "Density",
                          limits = c(0, NA)) +
    theme_bw() +
    labs(x = "UMAP1", y = "UMAP2") +
    theme(
        panel.grid = element_blank(),
        aspect.ratio = 1,
        axis.text = element_blank(),
        axis.ticks = element_blank()
    )
ggsave(p, file = paste0(save_path, '/projection_GSE141259_into_Inhouse_Reds.png'), width = 5, height = 4)
ggsave(p, file = paste0(save_path, '/projection_GSE141259_into_Inhouse_Reds.pdf'), width = 5, height = 4)

p <- ggplot() +
    geom_point(data = df_ref, 
                aes(x = UMAP_1, y = UMAP_2), 
               color = "grey75", size = 0.4, alpha = 1) +
    geom_point(data = df_query, 
                aes(x = UMAP_1, y = UMAP_2, color = "darkred"), 
                size = 0.4, alpha = 1) +
    theme_bw() +
    labs(x = "UMAP1", y = "UMAP2") +
    theme(
    panel.grid = element_blank(),
    aspect.ratio = 1,
    axis.text = element_blank(),
    axis.ticks = element_blank()
    ) + NoLegend()
ggsave(p, file = paste0(save_path, '/projection_of_query_onto_GSE141259_umap.png'), width = 5, height = 5)
ggsave(p, file = paste0(save_path, '/projection_of_query_onto_GSE141259_umap.pdf'), width = 5, height = 5)
