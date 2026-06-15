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

source('/home/sjcho/yard/functions/R/seurat_count_to_normalization.R')
source('/home/sjcho/yard/functions/R/FeaturePlot_sjcho.R')
source('/home/sjcho/yard/functions/R/save_ggplot2_to_ppt.R')

log_normalize <- function(so, save_path, data_name, pcs, nfeatures = 2000) {
    so <- NormalizeData(object = so)
    so <- FindVariableFeatures(object = so, nfeatures = nfeatures)
    so <- ScaleData(object = so)
    so <- RunPCA(object = so, npcs = 50)
    p <- ElbowPlot(so, ndims = 50, reduction = "pca")
    ggsave(p, file = paste0(save_path, '/', data_name, '_elbowplot.png'), width = 5, height = 5)
    so <- so %>%
        FindNeighbors(reduction = "pca") %>%
        FindClusters(resolution = 1) 
    so <- RunUMAP(so, dims = 1:pcs, reduction = 'pca')
    return(so)
}
plot_alveolar_epithelium_features <- function(so, save_path, data_name) {
    p <- fp_sjcho(so, features = c('Sftpc', 'Etv5', 'Hopx', 'Pdpn',
                                    'Aqp5', 'Sprr1a', 'Tnip3', 'Krt8', 'tdTomato', 'Foxj1', 'Muc5ac',
                                    'Scgb3a2', 'Krt5', 'Trp63', 'Cldn4', 
                                    'Mki67', 'Top2a', 'Nnat', 'Fn1', 'Plaur'), ncol = 5, order = T)
    ggsave(p, file = paste0(save_path, '/', data_name, '_alveolar_features.png'), width = 16, height = 12)
}

plot_lineage_features <- function(so, save_path, data_name) {
    p <- fp_sjcho(so, features = c('Sftpc', 'Etv5', 'Foxj1', 
                                    'Pdpn', 'Aqp5', 'Hopx',
                                    'Tnip3', 'Cldn4', 'Krt8',
                                    'Ptprc', 'Pecam1', 'Pdgfra'), ncol = 3)
    ggsave(p, file = paste0(save_path, '/', data_name, '_lineage_features.png'), width = 11, height = 12)
}
plot_metadata <- function(so, save_path, data_name) {
    p <- DimPlot(so, group.by = 'seurat_clusters', label = T, repel = T)
    ggsave(p, file = paste0(save_path, '/', data_name, '_clusters.png'), width = 7, height = 6)
    p <- DimPlot(so, group.by = 'seurat_clusters', split.by = 'dpi', label = T, ncol = 2)
    ggsave(p, file = paste0(save_path, '/', data_name, '_dpi.png'), width = 9, height = 8)
    p <- DimPlot(so, group.by = 'seurat_clusters', split.by = 'condition', label = T, ncol = 2)
    ggsave(p, file = paste0(save_path, '/', data_name, '_condition.png'), width = 9, height = 4)
    p <- DimPlot(so, group.by = 'seurat_clusters', split.by = 'condition_dpi', label = T, ncol = 4)
    ggsave(p, file = paste0(save_path, '/', data_name, '_condition_dpi.png'), width = 18, height = 8)
}

save_path = '/home/sjcho/projects/AT2_MHC2/20250714_revision/lineage_tracing/normalization_and_annotaiton/1.normalization/outs_lower250_mt5_umi2000'
data_name = 'lt_30dpi_'

sce_dAT2 <- readRDS('/home/sjcho/projects/AT2_MHC2/20250714_revision/lineage_tracing/QC/2.basic_QC/lower250/sce_dAT2_low_filter.rds')
sce_floxed <- readRDS('/home/sjcho/projects/AT2_MHC2/20250714_revision/lineage_tracing/QC/2.basic_QC/lower250/sce_floxed_low_filter.rds')

counts <- cbind(
  counts(sce_dAT2),
  counts(sce_floxed)
)

condition = c(rep('dAT2', ncol(counts(sce_dAT2))),
              rep('floxed', ncol(counts(sce_floxed))))

so <- CreateSeuratObject(counts, min.cells = 0, min.features = 0)
so$dpi <- '30dpi'
so$condition <- condition

so <- log_normalize(so, save_path, data_name, pcs = 15)

p <- DimPlot(so, group.by = 'condition', cols = c('gray50', 'darkred'))
ggsave(p, file = paste0(save_path, '/umap_condition.png'), width = 6.5, height = 6)
p <- DimPlot(so, group.by = 'condition', split.by = 'condition', cols = c('gray50', 'darkred'))
ggsave(p, file = paste0(save_path, '/umap_condition_split.png'), width = 13, height = 6)
p <- DimPlot(so, group.by = 'seurat_clusters', label = T, repel = T)
ggsave(p, file = paste0(save_path, '/', data_name, 'clusters.png'), width = 7, height = 6)

p <- fp_sjcho(so, features = c('Sox9', 'Col14a1', 'tdTomato'), ncol = 3, order = T)
ggsave(p, file = paste0(save_path, '/fp_sox9_col14a1_tdTomato.png'), width = 16, height = 5)

p <- fp_sjcho(so, features = c('Epcam', 'Pdgfra', 'Pdgfrb', 'Pecam1', 'Ptprc', 'tdTomato'), ncol = 3, order = T)
ggsave(p, file = paste0(save_path, '/fp_lineage_tdTomato.png'), width = 13, height = 8)

#### subset epithelial lineage
epithelial_cluster = c(6, 0, 1, 33, 29, 21, 8, 12, 18, 13, 7, 3, 4, 2, 17, 9)

p <- DimPlot(so, cells.highlight = Cells(so)[so$seurat_clusters %in% epithelial_cluster], sizes.highlight = 0.1, cols.highlight = 'darkred') + NoLegend()
ggsave(p, file = paste0(save_path, '/Dimplot_epithelial_clusters.png'), width = 6, height = 6)

save_path = '/home/sjcho/projects/AT2_MHC2/20250714_revision/lineage_tracing/normalization_and_annotaiton/1.normalization/outs_lower250_mt5_umi2000/epi_outs'

so.epi <- subset(so, subset = seurat_clusters %in% epithelial_cluster)
so.epi <- log_normalize(so.epi, save_path, paste0(data_name, 'epi'), pcs = 15)
so.epi$condition = factor(so.epi$condition, levels = c('floxed', 'dAT2'))

p <- DimPlot(so.epi, group.by = 'seurat_clusters', label = T, repel = T)
ggsave(p, file = paste0(save_path, '/', data_name, 'epi_clusters.png'), width = 7, height = 6)
p <- DimPlot(so.epi, group.by = 'condition', cols = c('gray50', 'darkred'))
ggsave(p, file = paste0(save_path, '/umap_condition.png'), width = 6.5, height = 6)
p <- DimPlot(so.epi, group.by = 'condition', split.by = 'condition', cols = c('gray50', 'darkred'))
ggsave(p, file = paste0(save_path, '/umap_condition_split.png'), width = 13, height = 6)

p <- fp_sjcho(so.epi, features = c('Sox9', 'Col14a1', 'tdTomato'), ncol = 3, order = T)
ggsave(p, file = paste0(save_path, '/fp_sox9_col14a1_tdTomato.png'), width = 16, height = 5)

p <- fp_sjcho(so.epi, features = c('Epcam', 'Pdgfra', 'Pdgfrb', 'Pecam1', 'Ptprc', 'tdTomato'), ncol = 3, order = T)
ggsave(p, file = paste0(save_path, '/fp_lineage_tdTomato.png'), width = 13, height = 8)

plot_alveolar_epithelium_features(so.epi, save_path, 'epi')

saveRDS(so.epi, file = paste0(save_path, '/so_epi_only.rds'))

library(sceasy)
library(reticulate)
use_condaenv('project_lung_exercise_R')
loompy <- reticulate::import('loompy')
so.epi[["RNA"]] <- as(so.epi[["RNA"]], "Assay")
sceasy::convertFormat(so.epi, from="seurat", to="anndata", drop_single_values=FALSE,
                       outFile= paste0('/home/sjcho/projects/AT2_MHC2/20250714_revision/lineage_tracing/normalization_and_annotaiton/1.normalization/outs_lower250_mt5_umi2000/epi_outs/epi.h5ad'))

############
alv_epi_clusters = c(12, 13, 11, 14, 17, 19)
save_path = '/home/sjcho/projects/AT2_MHC2/20250714_revision/lineage_tracing/normalization_and_annotaiton/1.normalization/outs_lower250_mt5_umi2000/alv_epi'

p <- DimPlot(so.epi, cells.highlight = Cells(so.epi)[so.epi$seurat_clusters %in% alv_epi_clusters], sizes.highlight = 0.1, cols.highlight = 'darkred') + NoLegend()
ggsave(p, file = paste0(save_path, '/Dimplot_alv_epi_clusters.png'), width = 6, height = 6)

so.alv_epi <- subset(so.epi, subset = seurat_clusters %in% alv_epi_clusters)

so.alv_epi <- log_normalize(so.alv_epi, save_path, paste0(data_name, 'alv_epi'), pcs = 10)
so.alv_epi$condition = factor(so.alv_epi$condition, levels = c('floxed', 'dAT2'))

p <- DimPlot(so.alv_epi, group.by = 'seurat_clusters', label = T, repel = T)
ggsave(p, file = paste0(save_path, '/', data_name, 'alv_epi_clusters.png'), width = 7, height = 6)
p <- DimPlot(so.alv_epi, group.by = 'condition', cols = c('gray50', 'darkred'))
ggsave(p, file = paste0(save_path, '/umap_condition.png'), width = 6.5, height = 6)
p <- DimPlot(so.alv_epi, group.by = 'condition', split.by = 'condition', cols = c('gray50', 'darkred'))
ggsave(p, file = paste0(save_path, '/umap_condition_split.png'), width = 13, height = 6)
p <- fp_sjcho(so.alv_epi, features = c('Sox9', 'Col14a1', 'tdTomato'), ncol = 3, order = T)
ggsave(p, file = paste0(save_path, '/fp_sox9_col14a1_tdTomato.png'), width = 16, height = 5)
p <- DimPlot(so.alv_epi, group.by = 'seurat_clusters', label = T, repel = T)
ggsave(p, file = paste0(save_path, '/alv_epi_clusters.png'), width = 7, height = 6)

p <- fp_sjcho(so.alv_epi, features = c('Epcam', 'Pdgfra', 'Pdgfrb', 'Ptprb', 'Pecam1', 'Ptprc', 'tdTomato'), ncol = 3, order = T)
ggsave(p, file = paste0(save_path, '/fp_lineage_tdTomato.png'), width = 13, height = 12)

plot_alveolar_epithelium_features(so.alv_epi, save_path, 'alv_epi')
plot_lineage_features(so.alv_epi, save_path, 'alv_epi')

####### v2
alv_epi_clusters_v2 = c(11, 12, 14, 9, 10, 15)

save_path = '/home/sjcho/projects/AT2_MHC2/20250714_revision/lineage_tracing/normalization_and_annotaiton/1.normalization/outs_lower250_mt5_umi2000/alv_epi_v2'

p <- DimPlot(so.alv_epi, cells.highlight = Cells(so.alv_epi)[!(so.alv_epi$seurat_clusters %in% alv_epi_clusters_v2)], sizes.highlight = 0.1, cols.highlight = 'darkred') + NoLegend()
ggsave(p, file = paste0(save_path, '/Dimplot_alv_epi_v2_clusters.png'), width = 6, height = 6)

so.alv_epi.v2 <- subset(so.alv_epi, subset = seurat_clusters %in% alv_epi_clusters_v2, invert = T)

so.alv_epi.v2 <- log_normalize(so.alv_epi.v2, save_path, paste0(data_name, 'alv_epi_v2'), pcs = 10)
so.alv_epi.v2$condition = factor(so.alv_epi.v2$condition, levels = c('floxed', 'dAT2'))

p <- DimPlot(so.alv_epi.v2, group.by = 'seurat_clusters', label = T, repel = T)
ggsave(p, file = paste0(save_path, '/', data_name, 'alv_epi_v2_clusters.png'), width = 7, height = 6)
p <- DimPlot(so.alv_epi.v2, group.by = 'condition', cols = c('gray50', 'darkred'))
ggsave(p, file = paste0(save_path, '/umap_condition.png'), width = 6.5, height = 6)
p <- DimPlot(so.alv_epi.v2, group.by = 'condition', split.by = 'condition', cols = c('gray50', 'darkred'))
ggsave(p, file = paste0(save_path, '/umap_condition_split.png'), width = 13, height = 6)
p <- fp_sjcho(so.alv_epi.v2, features = c('Sox9', 'Col14a1', 'tdTomato'), ncol = 3, order = T)
ggsave(p, file = paste0(save_path, '/fp_sox9_col14a1_tdTomato.png'), width = 16, height = 5)
p <- DimPlot(so.alv_epi.v2, group.by = 'seurat_clusters', label = T, repel = T)
ggsave(p, file = paste0(save_path, '/alv_epi_v2_clusters.png'), width = 7, height = 6)

p <- fp_sjcho(so.alv_epi.v2, features = c('Epcam', 'Pdgfra', 'Pdgfrb', 'Ptprb', 'Pecam1', 'Ptprc', 'tdTomato'), ncol = 3, order = T)
ggsave(p, file = paste0(save_path, '/fp_lineage_tdTomato.png'), width = 13, height = 12)

p <- VlnPlot(so.alv_epi.v2, features = 'tdTomato', pt.size = 0.05, group.by = 'seurat_clusters') + geom_boxplot(width = 0.2)
ggsave(p, file = paste0(save_path, '/tdTomato_vlnplot.png'), width = 8, height = 3)
p <- DotPlot(so.alv_epi.v2, features = 'tdTomato', group.by = 'seurat_clusters', scale = F) 
# x,y swap
p <- p + coord_flip()
ggsave(p, file = paste0(save_path, '/tdTomato_dotplot.png'), width = 8, height = 3)

plot_alveolar_epithelium_features(so.alv_epi.v2, save_path, 'alv_epi_v2')
plot_lineage_features(so.alv_epi.v2, save_path, 'alv_epi_v2')

dAT2_cells = Cells(so.alv_epi.v2)[so.alv_epi.v2$condition == 'floxed']
floxed_cells = Cells(so.alv_epi.v2)[so.alv_epi.v2$condition == 'dAT2']

so.alv_epi.v2@meta.data[dAT2_cells, 'condition'] = 'dAT2'
so.alv_epi.v2@meta.data[floxed_cells, 'condition'] = 'floxed'

p <- VlnPlot(subset(so.alv_epi.v2, seurat_clusters %in% c(1, 11, 2, 5, 3)), features = 'H2-Ab1', pt.size = 0.05, group.by = 'condition') + geom_boxplot(width = 0.2)
ggsave(p, file = paste0(save_path, '/H2-Ab1_vlnplot.png'), width = 4, height = 3)

saveRDS(so.alv_epi.v2, file = paste0(save_path, '/so_alv_epi.rds'))

###### for celltypist
library(sceasy)
library(reticulate)
use_condaenv('project_lung_exercise_R')
loompy <- reticulate::import('loompy')

so.alv_epi.v2[["RNA"]] <- as(so.alv_epi.v2[["RNA"]], "Assay")
sceasy::convertFormat(so.alv_epi.v2, 
                        from="seurat", 
                        to="anndata", 
                        drop_single_values=FALSE,
                        main_layer = 'data',
                        transfer_layers='counts',
                        outFile= paste0(save_path, 'so_alv_epi.h5ad')) 
