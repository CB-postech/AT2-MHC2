### conda activate project_lung_exercise_R

library(Seurat)
library(magrittr)
library(data.table)
library(ggplot2)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(ggsignif)

# load custom functions
source('/home/sjcho/yard/functions/R/seurat_count_to_normalization.R')
source('/home/sjcho/yard/functions/R/FeaturePlot_sjcho.R')
source('/home/sjcho/yard/functions/R/draw_proportion.R')
source('/home/sjcho/projects/AT2_MHC2/20241113_after_cellbender/utils.R')
source('/home/sjcho/yard/functions/R/draw_celltype_wise_dotplot.R')
source('/home/sjcho/yard/functions/R/gsea_sjcho.R')
source('/home/sjcho/yard/functions/R/save_ggplot2_to_ppt.R')

library(msigdbr)
h_gene_sets = msigdbr(species = "mouse", category = "H")
h_gene_sets$gs_name %>% unique

save_path = '/home/sjcho/projects/AT2_MHC2/code_for_publish/20250806_publish/inflammatory_signature_umap/outs'
so.alv <- readRDS('/home/sjcho/projects/AT2_MHC2/20241113_after_cellbender/3.whole_cell_annotation/outs/3.1.epithelial_annotation/alv_annotation.rds')

so.alv$condition_dpi = factor(so.alv$condition_dpi, levels = c('naive_floxed', '7dpi_floxed', '14dpi_floxed', '30dpi_floxed', 
                                                                'naive_dAT2', '7dpi_dAT2', '14dpi_dAT2', '30dpi_dAT2'))

alv_cols = c(AT2 = "#8192ef", tAT2 = "#dfd13a", AT1 = "#f66b07", AT2.INF = "#000000", AT2.prolif. = "#8bd690")

### 1. celltypist annotation
celltypist_result = fread('/home/sjcho/projects/AT2_MHC2/20241113_after_cellbender/manuscript_figures/figure4_and_exd6_20250225/alv_full_celltypist_result_overcluster1.csv')
celltypist_result = as.data.frame(celltypist_result)
rownames(celltypist_result) <- celltypist_result$V1
celltypist_result <- celltypist_result[, -1]

so.alv$predicted_labels = celltypist_result[Cells(so.alv), 'predicted_labels']
so.alv$predicted_labels[so.alv$predicted_labels == 'cycling.alv.epi'] <- 'AT2.prolif.'
so.alv$predicted_labels[so.alv$predicted_labels == 'IFN.response.alv.epi'] <- 'AT2.IFN'
so.alv$predicted_labels = factor(so.alv$predicted_labels, levels = c('AT2', 'tAT2', 'AT1', 'AT2.prolif.', 'AT2.IFN'))
so.alv$majority_voting = celltypist_result[Cells(so.alv), 'majority_voting']
so.alv$majority_voting[so.alv$majority_voting == 'cycling.alv.epi'] <- 'AT2.prolif.'
so.alv$majority_voting[so.alv$majority_voting == 'IFN.response.alv.epi'] <- 'AT2.IFN'
so.alv$majority_voting = factor(so.alv$majority_voting, levels = c('AT2', 'tAT2', 'AT1', 'AT2.prolif.', 'AT2.IFN'))

inflammatory_response = subset(h_gene_sets, gs_name == 'HALLMARK_INFLAMMATORY_RESPONSE')$gene_symbol

so.alv <- AddModuleScore(
    so.alv,
    features = list(inflammatory_response),
    name = 'HM_inflammatory')

p <- fp_sjcho(so.alv, features = 'HM_inflammatory1', ncol = 1, order = T)
p <- p + scale_colour_gradientn(colours = rev(c("#3d0000", "#8f0803", "#bb0c0c", "#f11b22", "#e76d71", "#ffdce1", "#fdf7f8", "#ffffff")))
ggsave(p, filename = file.path(save_path, '/HM_inflammatory.png'), width = 4.5, height = 4)
ggsave(p, filename = file.path(save_path, '/HM_inflammatory.pdf'), width = 4.5, height = 4)

p <- fp_sjcho(so.alv, features = 'HM_inflammatory1', ncol = 1, order = T)
ggsave(p, filename = file.path(save_path, '/HM_inflammatory.png'), width = 4.5, height = 4)
ggsave(p, filename = file.path(save_path, '/HM_inflammatory.pdf'), width = 4.5, height = 4)

p <- p + scale_color_gradientn(colors = rev(brewer.pal(11, "RdBu")))
ggsave(p, filename = file.path(save_path, '/HM_inflammatory_RdBu.png'), width = 4.5, height = 4)

p <- fp_sjcho(so.alv, features = 'HM_inflammatory1', split.by = 'condition_dpi') +
    patchwork::plot_layout(ncol = 4, nrow = 2)
p <- p & scale_color_gradientn(colors = rev(brewer.pal(11, "RdBu")))
ggsave(p, filename = file.path(save_path, '/HM_inflammatory_split_RdBu.png'), width = 12, height = 6)

p <- fp_sjcho(so.alv, features = 'HM_inflammatory1', split.by = 'condition_dpi') +
    patchwork::plot_layout(ncol = 4, nrow = 2)
ggsave(p, filename = file.path(save_path, '/HM_inflammatory_split_Reds.png'), width = 12, height = 6)

p <- fp_sjcho(subset(so.alv, majority_voting == 'AT2'), features = 'HM_inflammatory1', ncol = 1)
ggsave(p, filename = file.path(save_path, '/HM_inflammatory_AT2.png'), width = 4.5, height = 4)