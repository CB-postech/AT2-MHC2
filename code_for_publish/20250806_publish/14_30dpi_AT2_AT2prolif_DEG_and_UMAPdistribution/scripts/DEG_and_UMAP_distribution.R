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

save_path = '/home/sjcho/projects/AT2_MHC2/code_for_publish/20250806_publish/14_30dpi_AT2_AT2prolif_DEG_and_UMAPdistribution/outs'
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

p <- DimPlot(so.alv, cells.highlight = Cells(so.alv)[so.alv$majority_voting == 'AT1' & so.alv$dpi == '30dpi'], sizes.highlight = 0.5, cols.highlight = "#f66b07") + NoAxes() + NoLegend()
p <- p + ggtitle("")
p <- p + theme(
  panel.background = element_rect(fill = "transparent"), 
  plot.background = element_rect(fill = "transparent", color = NA),
  legend.background = element_rect(fill = "transparent"),
  legend.box.background = element_rect(fill = "transparent"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank()
)
ggsave(p, file = paste0(save_path, '/Dimplot_majority_voting_AT1_30dpi_nontraced.png'), dpi = 1000, width = 3, height = 3.25)

########### UMAP distribution of Pcan, Ki67, Cdkn1a / 3 cytokine responsiveness / HM inflammatory response

##### Pcan, Ki67, Cdkn1a
so.alv.14.30 <- subset(so.alv, dpi %in% c('14dpi', '30dpi'))

for (gene in c('Pcna', 'Mki67', 'Cdkn1a', 'Sftpc')) {
    p <- FeaturePlot(so.alv.14.30, features = gene, order = T, split = 'condition_dpi', pt.size = 1) + patchwork::plot_layout(ncol = 2, nrow = 2)
    p <- p & NoAxes()
    p <- p & scale_colour_gradientn(colours = rev(c("#3d0000", "#8f0803", "#bb0c0c", "#b9aeae", "#b8b8b8")), 
   # p <- p & scale_colour_gradientn(colours = rev(c("#3d0000", "#8f0803", "#bb0c0c", "#e76d71", "#ffdce1", "#fdf7f8", "#ffffff")), 
                                    limits = c(0, so.alv.14.30@assays$RNA$data[gene, ] %>% max %>% ceiling))
    p <- p & theme(legend.position = "right")
    ggsave(p, filename = paste0(save_path, '/splited_umap_', gene, '.png'), width = 10, height = 8)
    ggsave(p, filename = paste0(save_path, '/splited_umap_', gene, '.pdf'), width = 10, height = 8)
}

##### cytokine responsiveness
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

inflammatory_response = subset(h_gene_sets, gs_name == 'HALLMARK_INFLAMMATORY_RESPONSE')$gene_symbol
so.alv <- AddModuleScore(
    so.alv,
    features = list(inflammatory_response),
    name = 'HM_inflammatory')

########## GOBP: IL-11 mediated signaling pathway gene set
GO_list[['Il-11_mediated_signaling_pathway']] = c('Il11ra3', 'Il11', 'Il11ra1', 'Il11ra2', 'Il6st', 'Jak1', 'Stat3')

for (GO_term in names(GO_list)) {
    so.alv <- AddModuleScore(so.alv, features = list(GO_list[[GO_term]]), name = GO_term)
}

so.alv.14.30 <- subset(so.alv, dpi %in% c('14dpi', '30dpi'))
so.alv.14.30$'IL-1_response' = so.alv.14.30$GOBP_RESPONSE_TO_INTERLEUKIN_11
so.alv.14.30$'Tnfa_response' = so.alv.14.30$HM_Tnfa_NFKB1
so.alv.14.30$'Il-11_response' = so.alv.14.30$'Il-11_mediated_signaling_pathway1'
so.alv.14.30$'HM_inflammatory_response' = so.alv.14.30$'HM_inflammatory1'

for (feature in c('IL-1_response', 'Tnfa_response', 'Il-11_response', 'HM_inflammatory_response')) {
    p <- FeaturePlot(so.alv.14.30, features = feature, order = T, split = 'condition_dpi', pt.size = 1) + patchwork::plot_layout(ncol = 2, nrow = 2)
    p <- p & NoAxes()
    p <- p & scale_colour_gradientn(colours = rev(c("#3d0000", "#8f0803", "#bb0c0c", "#b9aeae", "#b8b8b8")), 
                                    limits = c(so.alv.14.30@meta.data[, feature] %>% min, so.alv.14.30@meta.data[, feature] %>% max))
    p <- p & theme(legend.position = "right")
    ggsave(p, filename = paste0(save_path, '/splited_umap_', feature, '.png'), width = 10, height = 8)
    ggsave(p, filename = paste0(save_path, '/splited_umap_', feature, '.pdf'), width = 10, height = 8)
}

###### 14dpi AT2 DEG / 14dpi AT2.prolif. DEG
markers.14.AT2 <- FindMarkers(so.alv, ident.1 = Cells(so.alv)[so.alv$dpi == '14dpi' & so.alv$majority_voting == 'AT2' & so.alv$condition == 'dAT2'],
                                ident.2 = Cells(so.alv)[so.alv$dpi == '14dpi' & so.alv$majority_voting == 'AT2' & so.alv$condition == 'floxed'],
                       logfc.threshold = 0.25, min.pct = 0.1)
write.csv(markers.14.AT2, file = paste0(save_path, '/markers_14dpi_AT2.csv'), row.names = T, quote = F, col.names = T)

markers.14.AT2.prolif. <- FindMarkers(so.alv, ident.1 = Cells(so.alv)[so.alv$dpi == '14dpi' & so.alv$majority_voting == 'AT2.prolif.' & so.alv$condition == 'dAT2'],
                                ident.2 = Cells(so.alv)[so.alv$dpi == '14dpi' & so.alv$majority_voting == 'AT2.prolif.' & so.alv$condition == 'floxed'],
                       logfc.threshold = 0.25, min.pct = 0.1)
write.csv(markers.14.AT2.prolif., file = paste0(save_path, '/markers_14dpi_AT2.prolif.csv'), row.names = T, quote = F, col.names = T)