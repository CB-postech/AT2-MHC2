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

save_path = '/home/sjcho/projects/AT2_MHC2/code_for_publish/20250806_publish/cytokine_response/outs'
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

################# cytokine response
############### IL4, IL6, IL13, IFNa, IFNb, IFNg
library(msigdbr)
C5_gene_sets = msigdbr(species = "mouse", category = "C5")
# h_gene_sets = msigdbr(species = "mouse", category = "H")
h_gene_sets$gs_name %>% unique

GO_list = list()
GO_list$'HM_IL6_JAK_STAT3' = subset(h_gene_sets, gs_name == 'HALLMARK_IL6_JAK_STAT3_SIGNALING')$gene_symbol
GO_list$'HM_IL2_STAT5' = subset(h_gene_sets, gs_name == 'HALLMARK_IL2_STAT5_SIGNALING')$gene_symbol
GO_list$'HM_IFNa' = subset(h_gene_sets, gs_name == 'HALLMARK_INTERFERON_ALPHA_RESPONSE')$gene_symbol
GO_list$'HM_IFNg' = subset(h_gene_sets, gs_name == 'HALLMARK_INTERFERON_GAMMA_RESPONSE')$gene_symbol
GO_list$'HM_Tnfa_NFKB' = subset(h_gene_sets, gs_name == 'HALLMARK_TNFA_SIGNALING_VIA_NFKB')$gene_symbol

for (GO_term in grep('RESPONSE_TO_INTERLEUKIN', C5_gene_sets$gs_name %>% unique, value = T)) {
    GO_list[[GO_term]] = subset(C5_gene_sets, gs_name == GO_term)$gene_symbol
}
########## GOBP: IL-11 mediated signaling pathway gene set
GO_list[['Il-11_mediated_signaling_pathway']] = c('Il11ra3', 'Il11', 'Il11ra1', 'Il11ra2', 'Il6st', 'Jak1', 'Stat3')

for (GO_term in names(GO_list)) {
    so.alv <- AddModuleScore(so.alv, features = list(GO_list[[GO_term]]), name = GO_term)
}

so.alv.Intermediate <- subset(so.alv, majority_voting %in% c('AT2') & condition_dpi %in% c('14dpi_floxed', '14dpi_dAT2', '30dpi_floxed', '30dpi_dAT2'))
so.alv.Intermediate$condition_dpi <- factor(so.alv.Intermediate$condition_dpi, levels = c('14dpi_floxed', '14dpi_dAT2', '30dpi_floxed', '30dpi_dAT2'))

p <- VlnPlot(so.alv.Intermediate, features = 'HM_IFNa1', group.by = 'condition_dpi', pt.size = 0, cols = c('darkgray', 'darkred', 'darkgray', 'darkred')) + geom_boxplot(width = 0.2) + NoLegend()
p <- p + geom_signif(comparisons = list(c('14dpi_floxed', '14dpi_dAT2'), c('30dpi_floxed', '30dpi_dAT2')), map_signif_level = TRUE)
p <- p + ylim(so.alv.Intermediate$HM_IFNa1 %>% min, (so.alv.Intermediate$HM_IFNa1 %>% max * 1.1))
ggsave(p, file = paste0(save_path, '/AT2_IFNa_response.png'), width = 6, height = 5)
ggsave(p, file = paste0(save_path, '/AT2_IFNa_response.pdf'), width = 6, height = 5)
ggplot2pptx(p, 6, 5, paste0(save_path, '/AT2_IFNa_response.pptx'))

p <- VlnPlot(so.alv.Intermediate, features = 'HM_IFNg1', group.by = 'condition_dpi', pt.size = 0, cols = c('darkgray', 'darkred', 'darkgray', 'darkred')) + geom_boxplot(width = 0.2) + NoLegend()
p <- p + geom_signif(comparisons = list(c('14dpi_floxed', '14dpi_dAT2'), c('30dpi_floxed', '30dpi_dAT2')), map_signif_level = TRUE)
p <- p + ylim(so.alv.Intermediate$HM_IFNg1 %>% min, (so.alv.Intermediate$HM_IFNg1 %>% max * 1.1))
ggsave(p, file = paste0(save_path, '/AT2_IFNg_response.png'), width = 6, height = 5)
ggsave(p, file = paste0(save_path, '/AT2_IFNg_response.pdf'), width = 6, height = 5)
ggplot2pptx(p, 6, 5, paste0(save_path, '/AT2_IFNg_response.pptx'))

p <- VlnPlot(so.alv.Intermediate, features = 'BP_IL41', group.by = 'condition_dpi', pt.size = 0, cols = c('darkgray', 'darkred', 'darkgray', 'darkred')) + geom_boxplot(width = 0.2) + NoLegend()
p <- p + geom_signif(comparisons = list(c('14dpi_floxed', '14dpi_dAT2'), c('30dpi_floxed', '30dpi_dAT2')), map_signif_level = TRUE)
p <- p + ylim(so.alv.Intermediate$BP_IL41 %>% min, (so.alv.Intermediate$BP_IL41 %>% max * 1.1))
ggsave(p, file = paste0(save_path, '/AT2_IL4_response.png'), width = 6, height = 5)
ggsave(p, file = paste0(save_path, '/AT2_IL4_response.pdf'), width = 6, height = 5)
ggplot2pptx(p, 6, 5, paste0(save_path, '/AT2_IL4_response.pptx'))

p <- VlnPlot(so.alv.Intermediate, features = 'BP_IL61', group.by = 'condition_dpi', pt.size = 0, cols = c('darkgray', 'darkred', 'darkgray', 'darkred')) + geom_boxplot(width = 0.2) + NoLegend()
p <- p + geom_signif(comparisons = list(c('14dpi_floxed', '14dpi_dAT2'), c('30dpi_floxed', '30dpi_dAT2')), map_signif_level = TRUE)
p <- p + ylim(so.alv.Intermediate$BP_IL61 %>% min, (so.alv.Intermediate$BP_IL61 %>% max * 1.1))
ggsave(p, file = paste0(save_path, '/AT2_IL6_response.png'), width = 6, height = 5)
ggsave(p, file = paste0(save_path, '/AT2_IL6_response.pdf'), width = 6, height = 5)
ggplot2pptx(p, 6, 5, paste0(save_path, '/AT2_IL6_response.pptx'))

###################### TNFa, Il-1b, Il-11
### AT2
so.alv.Intermediate <- subset(so.alv, majority_voting %in% c('AT2') & condition_dpi %in% c('14dpi_floxed', '14dpi_dAT2', '30dpi_floxed', '30dpi_dAT2'))
so.alv.Intermediate$condition_dpi <- factor(so.alv.Intermediate$condition_dpi, levels = c('14dpi_floxed', '14dpi_dAT2', '30dpi_floxed', '30dpi_dAT2'))

p1 <- VlnPlot(so.alv.Intermediate, features = 'GOBP_RESPONSE_TO_INTERLEUKIN_11', group.by = 'condition_dpi', pt.size = 0, cols = c('darkgray', 'darkred', 'darkgray', 'darkred')) + geom_boxplot(width = 0.2) + NoLegend()
p1 <- p1 + geom_signif(comparisons = list(c('14dpi_floxed', '14dpi_dAT2'), c('30dpi_floxed', '30dpi_dAT2')), map_signif_level = TRUE)
p1 <- p1 + ylim(so.alv.Intermediate$GOBP_RESPONSE_TO_INTERLEUKIN_11 %>% min, (so.alv.Intermediate$GOBP_RESPONSE_TO_INTERLEUKIN_11 %>% max * 1.25))

p2 <- VlnPlot(so.alv.Intermediate, features = 'HM_Tnfa_NFKB1', group.by = 'condition_dpi', pt.size = 0, cols = c('darkgray', 'darkred', 'darkgray', 'darkred')) + geom_boxplot(width = 0.2) + NoLegend()
p2 <- p2 + geom_signif(comparisons = list(c('14dpi_floxed', '14dpi_dAT2'), c('30dpi_floxed', '30dpi_dAT2')), map_signif_level = TRUE)
p2 <- p2 + ylim(so.alv.Intermediate$HM_Tnfa_NFKB1 %>% min, (so.alv.Intermediate$HM_Tnfa_NFKB1 %>% max * 1.25))

p3 <- VlnPlot(so.alv.Intermediate, features = 'Il-11_mediated_signaling_pathway1', group.by = 'condition_dpi', pt.size = 0, cols = c('darkgray', 'darkred', 'darkgray', 'darkred')) + geom_boxplot(width = 0.2) + NoLegend()
p3 <- p3 + geom_signif(comparisons = list(c('14dpi_floxed', '14dpi_dAT2'), c('30dpi_floxed', '30dpi_dAT2')), map_signif_level = TRUE)
p3 <- p3 + ylim(so.alv.Intermediate$'Il-11_mediated_signaling_pathway1' %>% min, (so.alv.Intermediate$'Il-11_mediated_signaling_pathway1' %>% max * 1.25))

p <- cowplot::plot_grid(p1, p2, p3, ncol = 1)
ggsave(p, file = paste0(save_path, '/TNFa_IL1_IL11_AT2.pdf'), width = 4, height = 15)

### tAT2
so.alv.Intermediate <- subset(so.alv, majority_voting %in% c('tAT2') & condition_dpi %in% c('14dpi_floxed', '14dpi_dAT2', '30dpi_floxed', '30dpi_dAT2'))
so.alv.Intermediate$condition_dpi <- factor(so.alv.Intermediate$condition_dpi, levels = c('14dpi_floxed', '14dpi_dAT2', '30dpi_floxed', '30dpi_dAT2'))

p1 <- VlnPlot(so.alv.Intermediate, features = 'GOBP_RESPONSE_TO_INTERLEUKIN_11', group.by = 'condition_dpi', pt.size = 0, cols = c('darkgray', 'darkred', 'darkgray', 'darkred')) + geom_boxplot(width = 0.2) + NoLegend()
p1 <- p1 + geom_signif(comparisons = list(c('14dpi_floxed', '14dpi_dAT2'), c('30dpi_floxed', '30dpi_dAT2')), map_signif_level = TRUE)
p1 <- p1 + ylim(so.alv.Intermediate$GOBP_RESPONSE_TO_INTERLEUKIN_11 %>% min, (so.alv.Intermediate$GOBP_RESPONSE_TO_INTERLEUKIN_11 %>% max * 1.25))

p2 <- VlnPlot(so.alv.Intermediate, features = 'HM_Tnfa_NFKB1', group.by = 'condition_dpi', pt.size = 0, cols = c('darkgray', 'darkred', 'darkgray', 'darkred')) + geom_boxplot(width = 0.2) + NoLegend()
p2 <- p2 + geom_signif(comparisons = list(c('14dpi_floxed', '14dpi_dAT2'), c('30dpi_floxed', '30dpi_dAT2')), map_signif_level = TRUE)
p2 <- p2 + ylim(so.alv.Intermediate$HM_Tnfa_NFKB1 %>% min, (so.alv.Intermediate$HM_Tnfa_NFKB1 %>% max * 1.25))

p3 <- VlnPlot(so.alv.Intermediate, features = 'Il-11_mediated_signaling_pathway1', group.by = 'condition_dpi', pt.size = 0, cols = c('darkgray', 'darkred', 'darkgray', 'darkred')) + geom_boxplot(width = 0.2) + NoLegend()
p3 <- p3 + geom_signif(comparisons = list(c('14dpi_floxed', '14dpi_dAT2'), c('30dpi_floxed', '30dpi_dAT2')), map_signif_level = TRUE)
p3 <- p3 + ylim(so.alv.Intermediate$'Il-11_mediated_signaling_pathway1' %>% min, (so.alv.Intermediate$'Il-11_mediated_signaling_pathway1' %>% max * 1.25))

p <- cowplot::plot_grid(p1, p2, p3, ncol = 1)
ggsave(p, file = paste0(save_path, '/TNFa_IL1_IL11_tAT2.pdf'), width = 4, height = 15)

### AT1
so.alv.Intermediate <- subset(so.alv, majority_voting %in% c('AT1') & condition_dpi %in% c('14dpi_floxed', '14dpi_dAT2', '30dpi_floxed', '30dpi_dAT2'))
so.alv.Intermediate$condition_dpi <- factor(so.alv.Intermediate$condition_dpi, levels = c('14dpi_floxed', '14dpi_dAT2', '30dpi_floxed', '30dpi_dAT2'))

p1 <- VlnPlot(so.alv.Intermediate, features = 'GOBP_RESPONSE_TO_INTERLEUKIN_11', group.by = 'condition_dpi', pt.size = 0, cols = c('darkgray', 'darkred', 'darkgray', 'darkred')) + geom_boxplot(width = 0.2) + NoLegend()
p1 <- p1 + geom_signif(comparisons = list(c('14dpi_floxed', '14dpi_dAT2'), c('30dpi_floxed', '30dpi_dAT2')), map_signif_level = TRUE)
p1 <- p1 + ylim(so.alv.Intermediate$GOBP_RESPONSE_TO_INTERLEUKIN_11 %>% min, (so.alv.Intermediate$GOBP_RESPONSE_TO_INTERLEUKIN_11 %>% max * 1.25))

p2 <- VlnPlot(so.alv.Intermediate, features = 'HM_Tnfa_NFKB1', group.by = 'condition_dpi', pt.size = 0, cols = c('darkgray', 'darkred', 'darkgray', 'darkred')) + geom_boxplot(width = 0.2) + NoLegend()
p2 <- p2 + geom_signif(comparisons = list(c('14dpi_floxed', '14dpi_dAT2'), c('30dpi_floxed', '30dpi_dAT2')), map_signif_level = TRUE)
p2 <- p2 + ylim(so.alv.Intermediate$HM_Tnfa_NFKB1 %>% min, (so.alv.Intermediate$HM_Tnfa_NFKB1 %>% max * 1.25))

p3 <- VlnPlot(so.alv.Intermediate, features = 'Il-11_mediated_signaling_pathway1', group.by = 'condition_dpi', pt.size = 0, cols = c('darkgray', 'darkred', 'darkgray', 'darkred')) + geom_boxplot(width = 0.2) + NoLegend()
p3 <- p3 + geom_signif(comparisons = list(c('14dpi_floxed', '14dpi_dAT2'), c('30dpi_floxed', '30dpi_dAT2')), map_signif_level = TRUE)
p3 <- p3 + ylim(so.alv.Intermediate$'Il-11_mediated_signaling_pathway1' %>% min, (so.alv.Intermediate$'Il-11_mediated_signaling_pathway1' %>% max * 1.25))

p <- cowplot::plot_grid(p1, p2, p3, ncol = 1)
ggsave(p, file = paste0(save_path, '/TNFa_IL1_IL11_AT1.pdf'), width = 4, height = 15)

####################### Jun, Fos
### AT2
so.alv.Intermediate <- subset(so.alv, majority_voting %in% c('AT2') & condition_dpi %in% c('14dpi_floxed', '14dpi_dAT2', '30dpi_floxed', '30dpi_dAT2'))
p1 <- VlnPlot(so.alv.Intermediate, features = 'Jun', group.by = 'condition_dpi', pt.size = 0, cols = c('darkgray', 'darkred', 'darkgray', 'darkred')) + geom_boxplot(width = 0.2) + NoLegend()
p1 <- p1 + geom_signif(comparisons = list(c('14dpi_floxed', '14dpi_dAT2'), c('30dpi_floxed', '30dpi_dAT2')), map_signif_level = TRUE)
p1 <- p1 + ylim(so.alv.Intermediate@assays$RNA$data['Jun', ] %>% min, (so.alv.Intermediate@assays$RNA$data['Jun', ] %>% max * 1.1))

p2 <- VlnPlot(so.alv.Intermediate, features = 'Fos', group.by = 'condition_dpi', pt.size = 0, cols = c('darkgray', 'darkred', 'darkgray', 'darkred')) + geom_boxplot(width = 0.2) + NoLegend()
p2 <- p2 + geom_signif(comparisons = list(c('14dpi_floxed', '14dpi_dAT2'), c('30dpi_floxed', '30dpi_dAT2')), map_signif_level = TRUE)
p2 <- p2 + ylim(so.alv.Intermediate@assays$RNA$data['Fos', ] %>% min, (so.alv.Intermediate@assays$RNA$data['Fos', ] %>% max * 1.1))

p <- cowplot::plot_grid(p1, p2, ncol = 1)
ggsave(p, file = paste0(save_path, '/AT2_Jun_Fos.png'), width = 6, height = 10)
ggsave(p, file = paste0(save_path, '/AT2_Jun_Fos.pdf'), width = 6, height = 10)
ggplot2pptx(p, 6, 10, paste0(save_path, '/AT2_Jun_Fos.pptx'))

### tAT2
so.alv.Intermediate <- subset(so.alv, majority_voting %in% c('tAT2') & condition_dpi %in% c('14dpi_floxed', '14dpi_dAT2', '30dpi_floxed', '30dpi_dAT2'))
p1 <- VlnPlot(so.alv.Intermediate, features = 'Jun', group.by = 'condition_dpi', pt.size = 0, cols = c('darkgray', 'darkred', 'darkgray', 'darkred')) + geom_boxplot(width = 0.2) + NoLegend()
p1 <- p1 + geom_signif(comparisons = list(c('14dpi_floxed', '14dpi_dAT2'), c('30dpi_floxed', '30dpi_dAT2')), map_signif_level = TRUE)
p1 <- p1 + ylim(so.alv.Intermediate@assays$RNA$data['Jun', ] %>% min, (so.alv.Intermediate@assays$RNA$data['Jun', ] %>% max * 1.1))

p2 <- VlnPlot(so.alv.Intermediate, features = 'Fos', group.by = 'condition_dpi', pt.size = 0, cols = c('darkgray', 'darkred', 'darkgray', 'darkred')) + geom_boxplot(width = 0.2) + NoLegend()
p2 <- p2 + geom_signif(comparisons = list(c('14dpi_floxed', '14dpi_dAT2'), c('30dpi_floxed', '30dpi_dAT2')), map_signif_level = TRUE)
p2 <- p2 + ylim(so.alv.Intermediate@assays$RNA$data['Fos', ] %>% min, (so.alv.Intermediate@assays$RNA$data['Fos', ] %>% max * 1.1))

p <- cowplot::plot_grid(p1, p2, ncol = 1)
ggsave(p, file = paste0(save_path, '/tAT2_Jun_Fos.png'), width = 6, height = 10)
ggsave(p, file = paste0(save_path, '/tAT2_Jun_Fos.pdf'), width = 6, height = 10)
ggplot2pptx(p, 6, 10, paste0(save_path, '/tAT2_Jun_Fos.pptx'))

