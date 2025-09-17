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
library(data.table)
library(RColorBrewer)

source('/home/sjcho/yard/functions/R/seurat_count_to_normalization.R')
source('/home/sjcho/yard/functions/R/FeaturePlot_sjcho.R')
source('/home/sjcho/yard/functions/R/save_ggplot2_to_ppt.R')

save_path = '/home/sjcho/projects/AT2_MHC2/code_for_publish/250724_tAT2/outs/2.test_are_there_pAT2'
so.alv <- readRDS('/home/sjcho/projects/AT2_MHC2/20241113_after_cellbender/3.whole_cell_annotation/outs/3.1.epithelial_annotation/alv_annotation.rds')

so.alv$condition_dpi = factor(so.alv$condition_dpi, levels = c('naive_floxed', '7dpi_floxed', '14dpi_floxed', '30dpi_floxed', 
                                                                'naive_dAT2', '7dpi_dAT2', '14dpi_dAT2', '30dpi_dAT2'))

alv_cols = c(AT2 = "#8192ef", tAT2 = "#dfd13a", AT1 = "#f66b07", AT2.INF = "#000000", AT2.prolif. = "#8bd690")

p <- DimPlot(so.alv, group.by = 'seurat_clusters')
ggsave(p, file = paste0(save_path, 'alv_dimplot.png'), width = 8, height = 6)

#################
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

################ pAT2 marker
pAT2_marker = c('Lcn2', 'Cxcl17', 'Lrg1', 'Itga7', 'Ptges', 'Glrx', 'Il33', 'Cd302')

so.alv <- AddModuleScore(so.alv, features = list(pAT2_marker), name = 'pAT2_score', assay = 'RNA')

p <- fp_sjcho(so.alv, features = 'pAT2_score1', pt.size = 0.1, order = T, ncol = 1)
p <- p + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
ggsave(p, file = paste0(save_path, '/pAT2_score.png'), width = 4, height = 3)

p <- FeaturePlot(so.alv, features = pAT2_marker, pt.size = 0.1, order = T, ncol = 4)
p <- p & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
ggsave(p, file = paste0(save_path, '/pAT2_score_featureplot.png'), width = 13, height = 6)