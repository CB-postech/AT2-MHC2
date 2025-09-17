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
source('/home/sjcho/yard/functions/R/save_ggplot2_to_ppt.R')
source('/home/sjcho/yard/functions/R/draw_proportion.R')
source('/home/sjcho/projects/AT2_MHC2/20241113_after_cellbender/utils.R')
source('/home/sjcho/yard/functions/R/draw_celltype_wise_dotplot.R')
source('/home/sjcho/yard/functions/R/gsea_sjcho.R')

library(msigdbr)
h_gene_sets = msigdbr(species = "mouse", category = "H")
h_gene_sets$gs_name %>% unique

save_path = '/home/sjcho/projects/AT2_MHC2/code_for_publish/20250806_publish/anchor_AT2/outs'
so.alv <- readRDS('/home/sjcho/projects/AT2_MHC2/20241113_after_cellbender/3.whole_cell_annotation/outs/3.1.epithelial_annotation/alv_annotation.rds')

so.alv$condition_dpi = factor(so.alv$condition_dpi, levels = c('naive_floxed', 'naive_dAT2', '7dpi_floxed', '7dpi_dAT2', '14dpi_floxed', '14dpi_dAT2', '30dpi_floxed', '30dpi_dAT2'))

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

### 2. AT2 anchor RNAs at 30dpi
so.alv.AT2.30dpi <- subset(so.alv, majority_voting == 'AT2' & dpi == '30dpi')
genes <- c('Lcn2' ,'Il33', 'Lrg1', 'Dmkn')

Idents(so.alv.AT2.30dpi) <- 'condition'
p1 <- VlnPlot(so.alv.AT2.30dpi, features = 'Lcn2', group.by = 'condition', cols = c('darkgray', 'darkred'), pt.size = 0)
p1 <- p1 + geom_boxplot(width = 0.2, fill = 'white') + NoLegend()
p1 <- p1 + geom_signif(comparisons = list(c('floxed', 'dAT2')), map_signif_level = TRUE)
p1 <- p1 + ylim(so.alv.AT2.30dpi@assays$RNA$data['Lcn2', ] %>% expm1 %>% min, (so.alv.AT2.30dpi@assays$RNA$data['Lcn2', ] %>% max * 1.1))
p2 <- VlnPlot(so.alv.AT2.30dpi, features = 'Il33', group.by = 'condition', cols = c('darkgray', 'darkred'), pt.size = 0)
p2 <- p2 + geom_boxplot(width = 0.2, fill = 'white') + NoLegend()
p2 <- p2 + geom_signif(comparisons = list(c('floxed', 'dAT2')), map_signif_level = TRUE)
p2 <- p2 + ylim(so.alv.AT2.30dpi@assays$RNA$data['Il33', ] %>% expm1 %>% min, (so.alv.AT2.30dpi@assays$RNA$data['Il33', ] %>% max * 1.1))
p3 <- VlnPlot(so.alv.AT2.30dpi, features = 'Lrg1', group.by = 'condition', cols = c('darkgray', 'darkred'), pt.size = 0)
p3 <- p3 + geom_boxplot(width = 0.2, fill = 'white') + NoLegend()
p3 <- p3 + geom_signif(comparisons = list(c('floxed', 'dAT2')), map_signif_level = TRUE)
p3 <- p3 + ylim(so.alv.AT2.30dpi@assays$RNA$data['Lrg1', ] %>% expm1 %>% min, (so.alv.AT2.30dpi@assays$RNA$data['Lrg1', ] %>% max * 1.1))
p4 <- VlnPlot(so.alv.AT2.30dpi, features = 'Dmkn', group.by = 'condition', cols = c('darkgray', 'darkred'), pt.size = 0)
p4 <- p4 + geom_boxplot(width = 0.2, fill = 'white') + NoLegend()
p4 <- p4 + geom_signif(comparisons = list(c('floxed', 'dAT2')), map_signif_level = TRUE)
p4 <- p4 + ylim(so.alv.AT2.30dpi@assays$RNA$data['Dmkn', ] %>% expm1 %>% min, (so.alv.AT2.30dpi@assays$RNA$data['Dmkn', ] %>% max * 1.1))

p <- cowplot::plot_grid(p1, p2, p3, p4, ncol = 1)
ggsave(p, file = paste0(save_path, '/AT2_anchor_RNAs_30dpi.png'), width = 3, height = 12)
ggsave(p, file = paste0(save_path, '/AT2_anchor_RNAs_30dpi.pdf'), width = 3, height = 12)
ggplot2pptx(p, file = paste0(save_path, '/AT2_anchor_RNAs_30dpi.pptx'), width = 3, height = 12)

### 2. tAT2 anchor RNAs at 30dpi
so.alv.AT2.30dpi <- subset(so.alv, majority_voting == 'tAT2' & dpi == '30dpi')
genes <- c('Lcn2' ,'Il33', 'Lrg1', 'Dmkn')

Idents(so.alv.AT2.30dpi) <- 'condition'
p1 <- VlnPlot(so.alv.AT2.30dpi, features = 'Lcn2', group.by = 'condition', cols = c('darkgray', 'darkred'), pt.size = 0)
p1 <- p1 + geom_boxplot(width = 0.2, fill = 'white') + NoLegend()
p1 <- p1 + geom_signif(comparisons = list(c('floxed', 'dAT2')), map_signif_level = TRUE)
p1 <- p1 + ylim(so.alv.AT2.30dpi@assays$RNA$data['Lcn2', ] %>% expm1 %>% min, (so.alv.AT2.30dpi@assays$RNA$data['Lcn2', ] %>% max * 1.1))
p2 <- VlnPlot(so.alv.AT2.30dpi, features = 'Il33', group.by = 'condition', cols = c('darkgray', 'darkred'), pt.size = 0)
p2 <- p2 + geom_boxplot(width = 0.2, fill = 'white') + NoLegend()
p2 <- p2 + geom_signif(comparisons = list(c('floxed', 'dAT2')), map_signif_level = TRUE)
p2 <- p2 + ylim(so.alv.AT2.30dpi@assays$RNA$data['Il33', ] %>% expm1 %>% min, (so.alv.AT2.30dpi@assays$RNA$data['Il33', ] %>% max * 1.1))
p3 <- VlnPlot(so.alv.AT2.30dpi, features = 'Lrg1', group.by = 'condition', cols = c('darkgray', 'darkred'), pt.size = 0)
p3 <- p3 + geom_boxplot(width = 0.2, fill = 'white') + NoLegend()
p3 <- p3 + geom_signif(comparisons = list(c('floxed', 'dAT2')), map_signif_level = TRUE)
p3 <- p3 + ylim(so.alv.AT2.30dpi@assays$RNA$data['Lrg1', ] %>% expm1 %>% min, (so.alv.AT2.30dpi@assays$RNA$data['Lrg1', ] %>% max * 1.1))
p4 <- VlnPlot(so.alv.AT2.30dpi, features = 'Dmkn', group.by = 'condition', cols = c('darkgray', 'darkred'), pt.size = 0)
p4 <- p4 + geom_boxplot(width = 0.2, fill = 'white') + NoLegend()
p4 <- p4 + geom_signif(comparisons = list(c('floxed', 'dAT2')), map_signif_level = TRUE)
p4 <- p4 + ylim(so.alv.AT2.30dpi@assays$RNA$data['Dmkn', ] %>% expm1 %>% min, (so.alv.AT2.30dpi@assays$RNA$data['Dmkn', ] %>% max * 1.1))

p <- cowplot::plot_grid(p1, p2, p3, p4, ncol = 1)
ggsave(p, file = paste0(save_path, '/tAT2_anchor_RNAs_30dpi.pdf'), width = 3, height = 12)

### 2. AT1 anchor RNAs at 30dpi
so.alv.AT2.30dpi <- subset(so.alv, majority_voting %in% c('AT2', 'tAT2') & dpi == '30dpi')
genes <- c('Lcn2' ,'Il33', 'Lrg1', 'Dmkn')

Idents(so.alv.AT2.30dpi) <- 'condition'
p1 <- VlnPlot(so.alv.AT2.30dpi, features = 'Lcn2', group.by = 'condition', cols = c('darkgray', 'darkred'), pt.size = 0)
p1 <- p1 + geom_boxplot(width = 0.2, fill = 'white') + NoLegend()
p1 <- p1 + geom_signif(comparisons = list(c('floxed', 'dAT2')), map_signif_level = TRUE)
p1 <- p1 + ylim(so.alv.AT2.30dpi@assays$RNA$data['Lcn2', ] %>% expm1 %>% min, (so.alv.AT2.30dpi@assays$RNA$data['Lcn2', ] %>% max * 1.1))
p2 <- VlnPlot(so.alv.AT2.30dpi, features = 'Il33', group.by = 'condition', cols = c('darkgray', 'darkred'), pt.size = 0)
p2 <- p2 + geom_boxplot(width = 0.2, fill = 'white') + NoLegend()
p2 <- p2 + geom_signif(comparisons = list(c('floxed', 'dAT2')), map_signif_level = TRUE)
p2 <- p2 + ylim(so.alv.AT2.30dpi@assays$RNA$data['Il33', ] %>% expm1 %>% min, (so.alv.AT2.30dpi@assays$RNA$data['Il33', ] %>% max * 1.1))
p3 <- VlnPlot(so.alv.AT2.30dpi, features = 'Lrg1', group.by = 'condition', cols = c('darkgray', 'darkred'), pt.size = 0)
p3 <- p3 + geom_boxplot(width = 0.2, fill = 'white') + NoLegend()
p3 <- p3 + geom_signif(comparisons = list(c('floxed', 'dAT2')), map_signif_level = TRUE)
p3 <- p3 + ylim(so.alv.AT2.30dpi@assays$RNA$data['Lrg1', ] %>% expm1 %>% min, (so.alv.AT2.30dpi@assays$RNA$data['Lrg1', ] %>% max * 1.1))
p4 <- VlnPlot(so.alv.AT2.30dpi, features = 'Dmkn', group.by = 'condition', cols = c('darkgray', 'darkred'), pt.size = 0)
p4 <- p4 + geom_boxplot(width = 0.2, fill = 'white') + NoLegend()
p4 <- p4 + geom_signif(comparisons = list(c('floxed', 'dAT2')), map_signif_level = TRUE)
p4 <- p4 + ylim(so.alv.AT2.30dpi@assays$RNA$data['Dmkn', ] %>% expm1 %>% min, (so.alv.AT2.30dpi@assays$RNA$data['Dmkn', ] %>% max * 1.1))

p <- cowplot::plot_grid(p1, p2, p3, p4, ncol = 1)
ggsave(p, file = paste0(save_path, '/AT2_tAT2_anchor_RNAs_30dpi.pdf'), width = 3, height = 12)

### 2. total anchor RNAs at 30dpi
so.alv.AT2.30dpi <- subset(so.alv, dpi == '30dpi')
genes <- c('Lcn2' ,'Il33', 'Lrg1', 'Dmkn')

Idents(so.alv.AT2.30dpi) <- 'condition'
p1 <- VlnPlot(so.alv.AT2.30dpi, features = 'Lcn2', group.by = 'condition', cols = c('darkgray', 'darkred'), pt.size = 0)
p1 <- p1 + geom_boxplot(width = 0.2, fill = 'white') + NoLegend()
p1 <- p1 + geom_signif(comparisons = list(c('floxed', 'dAT2')), map_signif_level = TRUE)
p1 <- p1 + ylim(so.alv.AT2.30dpi@assays$RNA$data['Lcn2', ] %>% expm1 %>% min, (so.alv.AT2.30dpi@assays$RNA$data['Lcn2', ] %>% max * 1.1))
p2 <- VlnPlot(so.alv.AT2.30dpi, features = 'Il33', group.by = 'condition', cols = c('darkgray', 'darkred'), pt.size = 0)
p2 <- p2 + geom_boxplot(width = 0.2, fill = 'white') + NoLegend()
p2 <- p2 + geom_signif(comparisons = list(c('floxed', 'dAT2')), map_signif_level = TRUE)
p2 <- p2 + ylim(so.alv.AT2.30dpi@assays$RNA$data['Il33', ] %>% expm1 %>% min, (so.alv.AT2.30dpi@assays$RNA$data['Il33', ] %>% max * 1.1))
p3 <- VlnPlot(so.alv.AT2.30dpi, features = 'Lrg1', group.by = 'condition', cols = c('darkgray', 'darkred'), pt.size = 0)
p3 <- p3 + geom_boxplot(width = 0.2, fill = 'white') + NoLegend()
p3 <- p3 + geom_signif(comparisons = list(c('floxed', 'dAT2')), map_signif_level = TRUE)
p3 <- p3 + ylim(so.alv.AT2.30dpi@assays$RNA$data['Lrg1', ] %>% expm1 %>% min, (so.alv.AT2.30dpi@assays$RNA$data['Lrg1', ] %>% max * 1.1))
p4 <- VlnPlot(so.alv.AT2.30dpi, features = 'Dmkn', group.by = 'condition', cols = c('darkgray', 'darkred'), pt.size = 0)
p4 <- p4 + geom_boxplot(width = 0.2, fill = 'white') + NoLegend()
p4 <- p4 + geom_signif(comparisons = list(c('floxed', 'dAT2')), map_signif_level = TRUE)
p4 <- p4 + ylim(so.alv.AT2.30dpi@assays$RNA$data['Dmkn', ] %>% expm1 %>% min, (so.alv.AT2.30dpi@assays$RNA$data['Dmkn', ] %>% max * 1.1))

p <- cowplot::plot_grid(p1, p2, p3, p4, ncol = 1)
ggsave(p, file = paste0(save_path, '/total_anchor_RNAs_30dpi.pdf'), width = 3, height = 12)


### volcano plot
markers <- FindMarkers(so.alv.AT2.30dpi, ident.1 = 'dAT2', ident.2 = 'floxed', logfc.threshold = 0, min.pct = 0, verbose = FALSE)
markers <- markers[markers$p_val_adj != 1, ]