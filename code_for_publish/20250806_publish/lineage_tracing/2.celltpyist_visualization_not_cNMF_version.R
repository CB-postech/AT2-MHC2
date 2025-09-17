### conda activate project_lung_exercise_R

library(Seurat)
library(magrittr)
library(data.table)
library(ggplot2)
library(stringr)
library(scater)
library(DropletUtils)
library(harmony)
library(pheatmap)
library(dplyr)
library(tidyr)
library(ggsignif)
library(EnhancedVolcano)
library(ggplotify)
library(reshape2)
library(RColorBrewer)

# load custom functions
source('/home/sjcho/yard/functions/R/seurat_count_to_normalization.R')
source('/home/sjcho/yard/functions/R/FeaturePlot_sjcho.R')
source('/home/sjcho/yard/functions/R/draw_proportion.R')
source('/home/sjcho/projects/AT2_MHC2/20241113_after_cellbender/utils.R')
source('/home/sjcho/yard/functions/R/draw_celltype_wise_dotplot.R')
source('/home/sjcho/yard/functions/R/gsea_sjcho.R')
source('/home/sjcho/yard/functions/R/save_ggplot2_to_ppt.R')

so <- readRDS('/home/sjcho/projects/AT2_MHC2/20250714_revision/lineage_tracing/normalization_and_annotaiton/1.normalization/outs_lower250_mt5_umi2000/alv_epi_v2/so_alv_epi.rds')
save_path = '/home/sjcho/projects/AT2_MHC2/code_for_publish/20250806_publish/lineage_tracing/outs'

celltypist_result <- read.csv('/home/sjcho/projects/AT2_MHC2/20250714_revision/lineage_tracing/normalization_and_annotaiton/1.normalization/outs_lower250_mt5_umi2000/alv_epi_v2/alv_full_celltypist_result_overcluster1.csv', 
                              row.names = 1, stringsAsFactors = FALSE)
so$predicted_labels = celltypist_result[Cells(so), 'predicted_labels']
so$predicted_labels[so$predicted_labels == 'cycling.alv.epi'] <- 'AT2.prolif.'
so$predicted_labels[so$predicted_labels == 'IFN.response.alv.epi'] <- 'AT2.IFN'
so$predicted_labels <- factor(so$predicted_labels, levels = c('AT1', 'AT2', 'tAT2', 'AT2.prolif.', 'AT2.IFN'))

so$majority_voting = celltypist_result[Cells(so), 'majority_voting']
so$majority_voting[so$majority_voting == 'cycling.alv.epi'] <- 'AT2.prolif.'
so$majority_voting[so$majority_voting == 'IFN.response.alv.epi'] <- 'AT2.IFN'
so$majority_voting <- factor(so$majority_voting, levels = c('AT1', 'AT2', 'tAT2', 'AT2.prolif.', 'AT2.IFN'))

alv_cols = c(AT2 = "#8192ef", tAT2 = "#dfd13a", AT1 = "#f66b07", AT2.IFN = "#000000", AT2.prolif. = "#8bd690")

result <- make.proportion.to.plot(so, col_meta = 'predicted_labels', row_meta = 'condition', colors = alv_cols, add_blackline = TRUE)
ggsave(result[[2]], file = paste0(save_path, '/alv_predicted_labels_proportion.png'), width = 4, height = 5)
result <- make.proportion.to.plot(so, col_meta = 'majority_voting', row_meta = 'condition', colors = alv_cols, add_blackline = TRUE)
ggsave(result[[2]], file = paste0(save_path, '/lineage_tracing_proportion_majority_voting.png'), width = 4, height = 5)
ggsave(result[[2]], file = paste0(save_path, '/lineage_tracing_proportion_majority_voting.pdf'), width = 4, height = 5)
ggplot2pptx(result[[2]], file = paste0(save_path, '/lineage_tracing_proportion_majority_voting.pptx'), width = 4, height = 5)

p <- DimPlot(so, group.by = 'majority_voting', cols = c(AT2 = "gray75", tAT2 = "gray75", AT1 = "#f66b07", AT2.IFN = "gray75", AT2.prolif. = "gray75")) + NoAxes() + NoLegend()
p <- p + ggtitle("")
p <- p + theme(
  panel.background = element_rect(fill = "transparent"), 
  plot.background = element_rect(fill = "transparent", color = NA),
  legend.background = element_rect(fill = "transparent"),
  legend.box.background = element_rect(fill = "transparent"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank()
)
ggsave(p, file = paste0(save_path, '/Dimplot_majority_voting_AT1.png'), dpi = 1000, width = 3, height = 3.25, bg = "transparent")
ggsave(p, file = paste0(save_path, '/Dimplot_majority_voting_AT1.pdf'), dpi = 1000, width = 3, height = 3.25, bg = "transparent")

# p <- DimPlot(so, group.by = 'predicted_labels', cols = cols) + NoAxes()
# ggsave(p, file = paste0(save_path, '/Dimplot_predicted_labels.png'), width = 4.5, height = 3)

# p <- DimPlot(so, cells.highlight = Cells(so)[so@assays$RNA$counts['tdTomato', ] > 0], cols.highlight = 'darkred', sizes.highlight = 0.5) + NoAxes() + NoLegend()
# ggsave(p, file = paste0(save_path, '/Dimplot_tdTomato.png'), width = 3, height = 3)

# result <- make.proportion.to.plot(subset(so, cells = Cells(so)[so@assays$RNA$counts['tdTomato', ] > 0]), col_meta = 'predicted_labels', row_meta = 'condition', colors = cols, add_blackline = TRUE)
# ggsave(result[[2]], file = paste0(save_path, '/alv_predicted_labels_proportion_tomato_positive.png'), width = 4, height = 5)
# result <- make.proportion.to.plot(subset(so, cells = Cells(so)[so@assays$RNA$counts['tdTomato', ] == 0]), col_meta = 'predicted_labels', row_meta = 'condition', colors = cols, add_blackline = TRUE)
# ggsave(result[[2]], file = paste0(save_path, '/alv_predicted_labels_proportion_tomato_negative.png'), width = 4, height = 5)

# so$tdTomato = ifelse(so@assays$RNA$counts['tdTomato', ] > 0, 'tdTomato_positive', 'tdTomato_negative')
# so$tdTomato = factor(so$tdTomato, levels = c('tdTomato_negative', 'tdTomato_positive'))

##################### contour plot
source('/home/sjcho/yard/functions/R/draw_contour.R')

p <- contour_plot_merge(so, name_of_dimension = 'umap',
            'condition',
            'majority_voting', 
            condition_color = c('floxed' = 'black', 'dAT2' = '#c00000'), 
            celltype_color = alv_cols, 
            nbin = 5, linewidth = 1)
ggsave(p, file = paste0(save_path, '/lineage_tracing_contour.pdf'), width = 6, height = 3)

p <- contour_plot_merge(so, name_of_dimension = 'umap',
            'condition',
            'majority_voting', 
            condition_color = c('floxed' = 'black', 'dAT2' = '#c00000'), 
            celltype_color = alv_cols, 
            nbin = 3.5, linewidth = 0.5)
ggsave(p, file = paste0(save_path, '/lineage_tracing_contour_v2.pdf'), width = 6, height = 3)

## palantir in python
library(sceasy)
library(reticulate)
use_condaenv('project_lung_exercise_R')
loompy <- reticulate::import('loompy')
so[["RNA"]] <- as(so[["RNA"]], "Assay")
sceasy::convertFormat(so, from="seurat", to="anndata", drop_single_values=FALSE,
                       outFile= paste0(save_path,'/full_so.h5ad'))

# Find the cell express Etv5 most
Cells(so)[max(so@assays$RNA$data['Etv5', ]) == so@assays$RNA$data['Etv5', ]]
# lt_30dpi_floxed_CAGCGAGGTCAGCGAG-1

p <- DimPlot(so, cells.highlight = Cells(so)[max(so@assays$RNA$data['Etv5', ]) == so@assays$RNA$data['Etv5', ]], 
             cols.highlight = 'darkred', sizes.highlight = 2) + NoAxes() + NoLegend()
ggsave(p, file = paste0(save_path, '/lineage_tracing_Dimplot_Etv5_max.png'), width = 3, height = 3)


### 2. AT2 anchor RNAs at 30dpi
so.AT2.30dpi <- subset(so, majority_voting == 'AT2' & dpi == '30dpi')
genes <- c('Lcn2' ,'Il33', 'Lrg1', 'Dmkn')

Idents(so.AT2.30dpi) <- 'condition'
p1 <- VlnPlot(so.AT2.30dpi, features = 'Lcn2', group.by = 'condition', cols = c('darkgray', 'darkred'), pt.size = 0)
p1 <- p1 + geom_boxplot(width = 0.2, fill = 'white') + NoLegend()
p1 <- p1 + geom_signif(comparisons = list(c('floxed', 'dAT2')), map_signif_level = TRUE)
p1 <- p1 + ylim(so.AT2.30dpi@assays$RNA$data['Lcn2', ] %>% expm1 %>% min, (so.AT2.30dpi@assays$RNA$data['Lcn2', ] %>% max * 1.1))
p2 <- VlnPlot(so.AT2.30dpi, features = 'Il33', group.by = 'condition', cols = c('darkgray', 'darkred'), pt.size = 0)
p2 <- p2 + geom_boxplot(width = 0.2, fill = 'white') + NoLegend()
p2 <- p2 + geom_signif(comparisons = list(c('floxed', 'dAT2')), map_signif_level = TRUE)
p2 <- p2 + ylim(so.AT2.30dpi@assays$RNA$data['Il33', ] %>% expm1 %>% min, (so.AT2.30dpi@assays$RNA$data['Il33', ] %>% max * 1.1))
p3 <- VlnPlot(so.AT2.30dpi, features = 'Lrg1', group.by = 'condition', cols = c('darkgray', 'darkred'), pt.size = 0)
p3 <- p3 + geom_boxplot(width = 0.2, fill = 'white') + NoLegend()
p3 <- p3 + geom_signif(comparisons = list(c('floxed', 'dAT2')), map_signif_level = TRUE)
p3 <- p3 + ylim(so.AT2.30dpi@assays$RNA$data['Lrg1', ] %>% expm1 %>% min, (so.AT2.30dpi@assays$RNA$data['Lrg1', ] %>% max * 1.1))
p4 <- VlnPlot(so.AT2.30dpi, features = 'Dmkn', group.by = 'condition', cols = c('darkgray', 'darkred'), pt.size = 0)
p4 <- p4 + geom_boxplot(width = 0.2, fill = 'white') + NoLegend()
p4 <- p4 + geom_signif(comparisons = list(c('floxed', 'dAT2')), map_signif_level = TRUE)
p4 <- p4 + ylim(so.AT2.30dpi@assays$RNA$data['Dmkn', ] %>% expm1 %>% min, (so.AT2.30dpi@assays$RNA$data['Dmkn', ] %>% max * 1.1))

p <- cowplot::plot_grid(p1, p2, p3, p4, ncol = 1)
ggsave(p, file = paste0(save_path, '/anchor_RNAs_30dpi_AT2.pdf'), width = 3, height = 12)

### 3. tAT2 anchor RNAs at 30dpi
so.AT2.30dpi <- subset(so, majority_voting == 'tAT2' & dpi == '30dpi')
genes <- c('Lcn2' ,'Il33', 'Lrg1', 'Dmkn')

Idents(so.AT2.30dpi) <- 'condition'
p1 <- VlnPlot(so.AT2.30dpi, features = 'Lcn2', group.by = 'condition', cols = c('darkgray', 'darkred'), pt.size = 0)
p1 <- p1 + geom_boxplot(width = 0.2, fill = 'white') + NoLegend()
p1 <- p1 + geom_signif(comparisons = list(c('floxed', 'dAT2')), map_signif_level = TRUE)
p1 <- p1 + ylim(so.AT2.30dpi@assays$RNA$data['Lcn2', ] %>% expm1 %>% min, (so.AT2.30dpi@assays$RNA$data['Lcn2', ] %>% max * 1.1))
p2 <- VlnPlot(so.AT2.30dpi, features = 'Il33', group.by = 'condition', cols = c('darkgray', 'darkred'), pt.size = 0)
p2 <- p2 + geom_boxplot(width = 0.2, fill = 'white') + NoLegend()
p2 <- p2 + geom_signif(comparisons = list(c('floxed', 'dAT2')), map_signif_level = TRUE)
p2 <- p2 + ylim(so.AT2.30dpi@assays$RNA$data['Il33', ] %>% expm1 %>% min, (so.AT2.30dpi@assays$RNA$data['Il33', ] %>% max * 1.1))
p3 <- VlnPlot(so.AT2.30dpi, features = 'Lrg1', group.by = 'condition', cols = c('darkgray', 'darkred'), pt.size = 0)
p3 <- p3 + geom_boxplot(width = 0.2, fill = 'white') + NoLegend()
p3 <- p3 + geom_signif(comparisons = list(c('floxed', 'dAT2')), map_signif_level = TRUE)
p3 <- p3 + ylim(so.AT2.30dpi@assays$RNA$data['Lrg1', ] %>% expm1 %>% min, (so.AT2.30dpi@assays$RNA$data['Lrg1', ] %>% max * 1.1))
p4 <- VlnPlot(so.AT2.30dpi, features = 'Dmkn', group.by = 'condition', cols = c('darkgray', 'darkred'), pt.size = 0)
p4 <- p4 + geom_boxplot(width = 0.2, fill = 'white') + NoLegend()
p4 <- p4 + geom_signif(comparisons = list(c('floxed', 'dAT2')), map_signif_level = TRUE)
p4 <- p4 + ylim(so.AT2.30dpi@assays$RNA$data['Dmkn', ] %>% expm1 %>% min, (so.AT2.30dpi@assays$RNA$data['Dmkn', ] %>% max * 1.1))

p <- cowplot::plot_grid(p1, p2, p3, p4, ncol = 1)
ggsave(p, file = paste0(save_path, '/anchor_RNAs_30dpi_tAT2.pdf'), width = 3, height = 12)

### 3. tAT2 & AT2 anchor RNAs at 30dpi
so.AT2.30dpi <- subset(so, majority_voting %in% c('tAT2', 'AT2') & dpi == '30dpi')
genes <- c('Lcn2' ,'Il33', 'Lrg1', 'Dmkn')

Idents(so.AT2.30dpi) <- 'condition'
p1 <- VlnPlot(so.AT2.30dpi, features = 'Lcn2', group.by = 'condition', cols = c('darkgray', 'darkred'), pt.size = 0)
p1 <- p1 + geom_boxplot(width = 0.2, fill = 'white') + NoLegend()
p1 <- p1 + geom_signif(comparisons = list(c('floxed', 'dAT2')), map_signif_level = TRUE)
p1 <- p1 + ylim(so.AT2.30dpi@assays$RNA$data['Lcn2', ] %>% expm1 %>% min, (so.AT2.30dpi@assays$RNA$data['Lcn2', ] %>% max * 1.1))
p2 <- VlnPlot(so.AT2.30dpi, features = 'Il33', group.by = 'condition', cols = c('darkgray', 'darkred'), pt.size = 0)
p2 <- p2 + geom_boxplot(width = 0.2, fill = 'white') + NoLegend()
p2 <- p2 + geom_signif(comparisons = list(c('floxed', 'dAT2')), map_signif_level = TRUE)
p2 <- p2 + ylim(so.AT2.30dpi@assays$RNA$data['Il33', ] %>% expm1 %>% min, (so.AT2.30dpi@assays$RNA$data['Il33', ] %>% max * 1.1))
p3 <- VlnPlot(so.AT2.30dpi, features = 'Lrg1', group.by = 'condition', cols = c('darkgray', 'darkred'), pt.size = 0)
p3 <- p3 + geom_boxplot(width = 0.2, fill = 'white') + NoLegend()
p3 <- p3 + geom_signif(comparisons = list(c('floxed', 'dAT2')), map_signif_level = TRUE)
p3 <- p3 + ylim(so.AT2.30dpi@assays$RNA$data['Lrg1', ] %>% expm1 %>% min, (so.AT2.30dpi@assays$RNA$data['Lrg1', ] %>% max * 1.1))
p4 <- VlnPlot(so.AT2.30dpi, features = 'Dmkn', group.by = 'condition', cols = c('darkgray', 'darkred'), pt.size = 0)
p4 <- p4 + geom_boxplot(width = 0.2, fill = 'white') + NoLegend()
p4 <- p4 + geom_signif(comparisons = list(c('floxed', 'dAT2')), map_signif_level = TRUE)
p4 <- p4 + ylim(so.AT2.30dpi@assays$RNA$data['Dmkn', ] %>% expm1 %>% min, (so.AT2.30dpi@assays$RNA$data['Dmkn', ] %>% max * 1.1))

p <- cowplot::plot_grid(p1, p2, p3, p4, ncol = 1)
ggsave(p, file = paste0(save_path, '/anchor_RNAs_30dpi_tAT2_AT2.pdf'), width = 3, height = 12)

### 4. alv anchor RNAs at 30dpi
so.AT2.30dpi <- subset(so, dpi == '30dpi')
genes <- c('Lcn2' ,'Il33', 'Lrg1', 'Dmkn')

Idents(so.AT2.30dpi) <- 'condition'
p1 <- VlnPlot(so.AT2.30dpi, features = 'Lcn2', group.by = 'condition', cols = c('darkgray', 'darkred'), pt.size = 0)
p1 <- p1 + geom_boxplot(width = 0.2, fill = 'white') + NoLegend()
p1 <- p1 + geom_signif(comparisons = list(c('floxed', 'dAT2')), map_signif_level = TRUE)
p1 <- p1 + ylim(so.AT2.30dpi@assays$RNA$data['Lcn2', ] %>% expm1 %>% min, (so.AT2.30dpi@assays$RNA$data['Lcn2', ] %>% max * 1.1))
p2 <- VlnPlot(so.AT2.30dpi, features = 'Il33', group.by = 'condition', cols = c('darkgray', 'darkred'), pt.size = 0)
p2 <- p2 + geom_boxplot(width = 0.2, fill = 'white') + NoLegend()
p2 <- p2 + geom_signif(comparisons = list(c('floxed', 'dAT2')), map_signif_level = TRUE)
p2 <- p2 + ylim(so.AT2.30dpi@assays$RNA$data['Il33', ] %>% expm1 %>% min, (so.AT2.30dpi@assays$RNA$data['Il33', ] %>% max * 1.1))
p3 <- VlnPlot(so.AT2.30dpi, features = 'Lrg1', group.by = 'condition', cols = c('darkgray', 'darkred'), pt.size = 0)
p3 <- p3 + geom_boxplot(width = 0.2, fill = 'white') + NoLegend()
p3 <- p3 + geom_signif(comparisons = list(c('floxed', 'dAT2')), map_signif_level = TRUE)
p3 <- p3 + ylim(so.AT2.30dpi@assays$RNA$data['Lrg1', ] %>% expm1 %>% min, (so.AT2.30dpi@assays$RNA$data['Lrg1', ] %>% max * 1.1))
p4 <- VlnPlot(so.AT2.30dpi, features = 'Dmkn', group.by = 'condition', cols = c('darkgray', 'darkred'), pt.size = 0)
p4 <- p4 + geom_boxplot(width = 0.2, fill = 'white') + NoLegend()
p4 <- p4 + geom_signif(comparisons = list(c('floxed', 'dAT2')), map_signif_level = TRUE)
p4 <- p4 + ylim(so.AT2.30dpi@assays$RNA$data['Dmkn', ] %>% expm1 %>% min, (so.AT2.30dpi@assays$RNA$data['Dmkn', ] %>% max * 1.1))

p <- cowplot::plot_grid(p1, p2, p3, p4, ncol = 1)
ggsave(p, file = paste0(save_path, '/anchor_RNAs_30dpi_alv.pdf'), width = 3, height = 12)

p <- fp_sjcho(so, features = c('Lcn2', 'Il33', 'Lrg1', 'Dmkn'), ncol =  2, order = T) & NoAxes()
ggsave(p, file = paste0(save_path, '/anchor_RNAs_30dpi_alv_featureplot.pdf'), width = 6, height = 5)

Idents(so) <- so$condition
markers.AT2 <- FindMarkers(so, ident.1 = 'dAT2', ident.2 = 'floxed', logfc.threshold = 0, min.pct = 0)

p <- fp_sjcho(so, features = c('log10UMI', 'log10nFeature', 'mt.pct'), ncol = 2)
ggsave(p, file = paste0(save_path, '/quality_control_featureplot.pdf'), width = 6, height = 5)

########### cytokine and Jun Fos
library(msigdbr)
C5_gene_sets = msigdbr(species = "mouse", category = "C5")
h_gene_sets = msigdbr(species = "mouse", category = "H")
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

for (GO_term in names(GO_list)) {
    so <- AddModuleScore(so, features = list(GO_list[[GO_term]]), name = GO_term)
}


###################### TNFa, Il-1b
### AT2
so.Intermediate <- subset(so, majority_voting %in% c('AT2'))
p1 <- VlnPlot(so.Intermediate, features = 'GOBP_RESPONSE_TO_INTERLEUKIN_11', group.by = 'condition', pt.size = 0, cols = c('darkgray', 'darkred', 'darkgray', 'darkred')) + geom_boxplot(width = 0.2) + NoLegend()
p1 <- p1 + geom_signif(comparisons = list(c('floxed', 'dAT2')), map_signif_level = TRUE)
p1 <- p1 + ylim(so.Intermediate$GOBP_RESPONSE_TO_INTERLEUKIN_11 %>% min, (so.Intermediate$GOBP_RESPONSE_TO_INTERLEUKIN_11 %>% max * 1.1))

p2 <- VlnPlot(so.Intermediate, features = 'HM_Tnfa_NFKB1', group.by = 'condition', pt.size = 0, cols = c('darkgray', 'darkred', 'darkgray', 'darkred')) + geom_boxplot(width = 0.2) + NoLegend()
p2 <- p2 + geom_signif(comparisons = list(c('floxed', 'dAT2')), map_signif_level = TRUE)
p2 <- p2 + ylim(so.Intermediate$HM_Tnfa_NFKB1 %>% min, (so.Intermediate$HM_Tnfa_NFKB1 %>% max * 1.1))

p <- cowplot::plot_grid(p1, p2, ncol = 1)
ggsave(p, file = paste0(save_path, '/AT2_TNFa_IL1.png'), width = 6, height = 10)
ggsave(p, file = paste0(save_path, '/AT2_TNFa_IL1.pdf'), width = 6, height = 10)
ggplot2pptx(p, 6, 10, paste0(save_path, '/AT2_TNFa_IL1.pptx'))

### tAT2
so.Intermediate <- subset(so, majority_voting %in% c('tAT2'))
p1 <- VlnPlot(so.Intermediate, features = 'GOBP_RESPONSE_TO_INTERLEUKIN_11', group.by = 'condition', pt.size = 0, cols = c('darkgray', 'darkred', 'darkgray', 'darkred')) + geom_boxplot(width = 0.2) + NoLegend()
p1 <- p1 + geom_signif(comparisons = list(c('floxed', 'dAT2')), map_signif_level = TRUE)
p1 <- p1 + ylim(so.Intermediate$GOBP_RESPONSE_TO_INTERLEUKIN_11 %>% min, (so.Intermediate$GOBP_RESPONSE_TO_INTERLEUKIN_11 %>% max * 1.1))

p2 <- VlnPlot(so.Intermediate, features = 'HM_Tnfa_NFKB1', group.by = 'condition', pt.size = 0, cols = c('darkgray', 'darkred', 'darkgray', 'darkred')) + geom_boxplot(width = 0.2) + NoLegend()
p2 <- p2 + geom_signif(comparisons = list(c('floxed', 'dAT2')), map_signif_level = TRUE)
p2 <- p2 + ylim(so.Intermediate$HM_Tnfa_NFKB1 %>% min, (so.Intermediate$HM_Tnfa_NFKB1 %>% max * 1.1))

p <- cowplot::plot_grid(p1, p2, ncol = 1)
ggsave(p, file = paste0(save_path, '/tAT2_TNFa_IL1.png'), width = 6, height = 10)
ggsave(p, file = paste0(save_path, '/tAT2_TNFa_IL1.pdf'), width = 6, height = 10)
ggplot2pptx(p, 6, 10, paste0(save_path, '/tAT2_TNFa_IL1.pptx'))

####################### Jun, Fos
### AT2
p1 <- VlnPlot(so.Intermediate, features = 'Jun', group.by = 'condition', pt.size = 0, cols = c('darkgray', 'darkred', 'darkgray', 'darkred')) + geom_boxplot(width = 0.2) + NoLegend()
p1 <- p1 + geom_signif(comparisons = list(c('floxed', 'dAT2')), map_signif_level = TRUE)
p1 <- p1 + ylim(so.Intermediate@assays$RNA$data['Jun', ] %>% min, (so.Intermediate@assays$RNA$data['Jun', ] %>% max * 1.1))

p2 <- VlnPlot(so.Intermediate, features = 'Fos', group.by = 'condition', pt.size = 0, cols = c('darkgray', 'darkred', 'darkgray', 'darkred')) + geom_boxplot(width = 0.2) + NoLegend()
p2 <- p2 + geom_signif(comparisons = list(c('floxed', 'dAT2')), map_signif_level = TRUE)
p2 <- p2 + ylim(so.Intermediate@assays$RNA$data['Fos', ] %>% min, (so.Intermediate@assays$RNA$data['Fos', ] %>% max * 1.1))

p <- cowplot::plot_grid(p1, p2, ncol = 1)
ggsave(p, file = paste0(save_path, '/AT2_Jun_Fos.png'), width = 6, height = 10)
ggsave(p, file = paste0(save_path, '/AT2_Jun_Fos.pdf'), width = 6, height = 10)
ggplot2pptx(p, 6, 10, paste0(save_path, '/AT2_Jun_Fos.pptx'))

p1 <- VlnPlot(so.Intermediate, features = 'H2-Ab1', group.by = 'condition', pt.size = 0, cols = c('darkgray', 'darkred', 'darkgray', 'darkred')) + geom_boxplot(width = 0.2) + NoLegend()
p1 <- p1 + geom_signif(comparisons = list(c('floxed', 'dAT2')), map_signif_level = TRUE)
p1 <- p1 + ylim(so.Intermediate@assays$RNA$data['H2-Ab1', ] %>% min, (so.Intermediate@assays$RNA$data['H2-Ab1', ] %>% max * 1.1))
ggsave(p1, file = paste0(save_path, '/AT2_H2-Ab1.png'), width = 6, height = 5)
ggsave(p1, file = paste0(save_path, '/AT2_H2-Ab1.pdf'), width = 6, height = 5)
ggplot2pptx(p1, 6, 5, paste0(save_path, '/AT2_H2-Ab1.pptx'))
