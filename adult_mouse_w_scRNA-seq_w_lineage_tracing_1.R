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

# load custom functions
source('/home/sjcho/yard/functions/R/seurat_count_to_normalization.R')
source('/home/sjcho/yard/functions/R/FeaturePlot_sjcho.R')
source('/home/sjcho/yard/functions/R/draw_proportion.R')
source('/home/sjcho/projects/AT2_MHC2/20241113_after_cellbender/utils.R')
source('/home/sjcho/yard/functions/R/draw_celltype_wise_dotplot.R')
source('/home/sjcho/yard/functions/R/gsea_sjcho.R')
source('/home/sjcho/yard/functions/R/save_ggplot2_to_ppt.R')

so <- readRDS('/home/sjcho/projects/AT2_MHC2/20250714_revision/lineage_tracing/normalization_and_annotaiton/1.normalization/outs_lower100/alv_epi_v3/so_alv_epi_v3.rds')
save_path = '/home/sjcho/projects/AT2_MHC2/20250714_revision/lineage_tracing/normalization_and_annotaiton/downstream/outs/2.celltypist/not_cNMF_annotation'

celltypist_result <- read.csv('/home/sjcho/projects/AT2_MHC2/20250714_revision/lineage_tracing/normalization_and_annotaiton/downstream/outs/2.celltypist/celltypist_alv_proportion_celltypist_result_not_cNMF.csv', 
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

condition_name = 'condition'
ref_condition = 'floxed'
exp_condition = 'dAT2'
celltype_name = 'majority_voting'
nbin = 10
type_color = c('floxed' = 'black', 'dAT2' = '#c00000')

p <- contour_plot_merge(so, name_of_dimension = 'umap',
            'condition',
            'majority_voting', 
            condition_color = c('floxed' = 'black', 'dAT2' = '#c00000'), 
            celltype_color = alv_cols, 
            nbin = 5, linewidth = 1)
ggsave(p, file = paste0(save_path, '/lineage_tracing_contour.png'), width = 6, height = 3)
ggsave(p, file = paste0(save_path, '/lineage_tracing_contour.pdf'), width = 6, height = 3)
ggplot2pptx(p, file = paste0(save_path, '/lineage_tracing_contour.pptx'), width = 6, height = 3)

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
# lt_30dpi_floxed_AAGTCACCAATGCCGA-1

p <- DimPlot(so, cells.highlight = Cells(so)[max(so@assays$RNA$data['Etv5', ]) == so@assays$RNA$data['Etv5', ]], 
             cols.highlight = 'darkred', sizes.highlight = 2) + NoAxes() + NoLegend()
ggsave(p, file = paste0(save_path, '/lineage_tracing_Dimplot_Etv5_max.png'), width = 3, height = 3)


### 2. AT2 anchor RNAs at 30dpi
so.alv.AT2.30dpi <- subset(so, majority_voting == 'AT2' & dpi == '30dpi')
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

############# cytokine and Jun Fos
