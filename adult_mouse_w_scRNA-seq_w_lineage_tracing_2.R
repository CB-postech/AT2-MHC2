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
library(ggpubr)

# load custom functions
source('/home/sjcho/yard/functions/R/seurat_count_to_normalization.R')
source('/home/sjcho/yard/functions/R/FeaturePlot_sjcho.R')
source('/home/sjcho/yard/functions/R/save_ggplot2_to_ppt.R')
source('/home/sjcho/yard/functions/R/draw_proportion.R')
source('/home/sjcho/projects/AT2_MHC2/20241113_after_cellbender/utils.R')
source('/home/sjcho/yard/functions/R/draw_celltype_wise_dotplot.R')
source('/home/sjcho/yard/functions/R/gsea_sjcho.R')

so <- readRDS('/home/sjcho/projects/AT2_MHC2/20250714_revision/lineage_tracing/normalization_and_annotaiton/1.normalization/outs_lower250_mt5_umi2000/alv_epi_v2/so_alv_epi.rds')
save_path = '/home/sjcho/projects/AT2_MHC2/code_for_publish/20250806_publish/lineage_tracing/outs/palantir_vis'

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

# Etv5 Vlion-plot
p <- VlnPlot(so, features = 'Etv5', group.by = 'seurat_clusters', pt.size = 0) + geom_boxplot(width = 0.2, fill = 'white')
ggsave(p, file = paste0(save_path, '/Etv5_VlnPlot_seurat_clusters.png'), width = 6, height = 4)
# cluster 1 has highest average-Etv5 expression

p <- DimPlot(so, group.by = 'seurat_clusters', label = T)
ggsave(p, file = paste0(save_path, '/seurat_clusters.png'), width = 6, height = 4)

# Find highest Etv5 expressed cell in cluster1
Etv5.cluster1 = so@assays$RNA$data['Etv5', Cells(so)[so$seurat_clusters == 1]]
Etv5.cluster1 = sort(Etv5.cluster1, decreasing = T)
# lt_30dpi_dAT2_CTCGATCCAAACCGTG-1

### load palantir result
branch <- fread('/home/sjcho/projects/AT2_MHC2/code_for_publish/20250806_publish/lineage_tracing/outs/palantir/branch_masks.csv', header = T)
branch <- as.data.frame(branch)
rownames(branch) <- branch$V1
branch <- branch[, -1]

pseudotime <- read.csv('/home/sjcho/projects/AT2_MHC2/code_for_publish/20250806_publish/lineage_tracing/outs/palantir/palantir_pseudotime.csv')
colnames(pseudotime) = c('cellname', 'pseudotime')
rownames(pseudotime) = pseudotime$cellname
pseudotime$condition = so$condition[rownames(pseudotime)]; pseudotime$condition = factor(pseudotime$condition, levels = c('dAT2', 'floxed'))
pseudotime$celltype = so$majority_voting[rownames(pseudotime)]; pseudotime$celltype = factor(pseudotime$celltype, levels = c('AT2', 'tAT2', 'AT1'))
pseudotime$dpi = so$dpi[rownames(pseudotime)]
so$pseudotime = pseudotime[colnames(so), 'pseudotime']

### path probability & pseudotime
probability <- read.csv('/home/sjcho/projects/AT2_MHC2/code_for_publish/20250806_publish/lineage_tracing/outs/palantir/palantir_fate_probabilities.csv')
rownames(probability) <- probability$X
so$path1_prob <- probability[colnames(so), 'lt_30dpi_dAT2_AATGGCTCAGTTGGCA.1']
so$path2_prob <- probability[colnames(so), 'lt_30dpi_dAT2_TCCACTTAGCTTAACG.1']

p <- FeaturePlot(so, features = 'path1_prob', pt.size = 0.75, order = T) + NoAxes() + NoLegend() + scale_colour_gradientn(colours = colorRampPalette(c("darkgray", 'green', "black"))(100))
p <- p + theme(plot.title = element_text(size = 0))
ggsave(p, file = paste0(save_path, '/path1_prob.png'), width = 4, height = 4)
ggsave(p, file = paste0(save_path, '/path1_prob.pdf'), width = 4, height = 4)
ggplot2pptx(p, 4, 4, paste0(save_path, '/path1_prob.pptx'))

p <- FeaturePlot(so, features = 'path2_prob', pt.size = 0.75, order = T) + NoAxes() + NoLegend() + scale_colour_gradientn(colours = colorRampPalette(c("darkgray", 'green', "black"))(100))
p <- p + theme(plot.title = element_text(size = 0))
ggsave(p, file = paste0(save_path, '/path2_prob.png'), width = 4, height = 4)
ggsave(p, file = paste0(save_path, '/path2_prob.pdf'), width = 4, height = 4)
ggplot2pptx(p, 4, 4, paste0(save_path, '/path2_prob.pptx'))

leg <- get_legend(FeaturePlot(so, features = 'pseudotime', pt.size = 0.75, order = T) + 
        scale_colour_gradientn(colours = colorRampPalette(c("darkgray", 'green', "black"))(100)))
ggplot2pptx(as.ggplot(leg), 1, 6, paste0(save_path, '/pseudotime_legend.pptx'))
ggsave(as.ggplot(leg), file = paste0(save_path, '/pseudotime_legend.png'), width = 1, height = 6)
ggsave(as.ggplot(leg), file = paste0(save_path, '/pseudotime_legend.pdf'), width = 1, height = 6)

p <- FeaturePlot(so, features = 'pseudotime', pt.size = 0.75, order = T) + NoAxes() + NoLegend() + scale_colour_gradientn(colours = colorRampPalette(c("darkgray", 'green', "black"))(100))
p <- p + geom_point(data = subset(so, cells = 'lt_30dpi_dAT2_CTCGATCCAAACCGTG-1')@reductions$umap@cell.embeddings, aes(x = umap_1, y = umap_2), color = 'darkred', size = 4)
p <- p + theme(plot.title = element_text(size = 0))
ggsave(p, file = paste0(save_path, '/pseudotime.png'), width = 4, height = 4)
ggsave(p, file = paste0(save_path, '/pseudotime.pdf'), width = 4, height = 4)
ggplot2pptx(p, 4, 4, paste0(save_path, '/pseudotime.pptx'))

# p1 <- FeaturePlot(subset(so, condition == 'floxed'), features = 'pseudotime', pt.size = 0.75, order = T) + NoAxes() + NoLegend() + scale_colour_gradientn(colours = colorRampPalette(c("darkgray", 'green', "black"))(100))
# p1 <- p1 + geom_point(data = subset(so, cells = Cells(so)[max(so@assays$RNA$data['Etv5', ]) == so@assays$RNA$data['Etv5', ]])@reductions$umap@cell.embeddings, aes(x = umap_1, y = umap_2), color = 'darkred', size = 4)
# p1 <- p1 + theme(plot.title = element_text(size = 0))
# p2 <- FeaturePlot(subset(so, condition == 'dAT2'), features = 'pseudotime', pt.size = 0.75, order = T) + NoAxes() + NoLegend() + scale_colour_gradientn(colours = colorRampPalette(c("darkgray", 'green', "black"))(100))
# p2 <- p2 + geom_point(data = subset(so, cells = Cells(so)[max(so@assays$RNA$data['Etv5', ]) == so@assays$RNA$data['Etv5', ]])@reductions$umap@cell.embeddings, aes(x = umap_1, y = umap_2), color = 'darkred', size = 4)
# p2 <- p2 + theme(plot.title = element_text(size = 0))
# ggsave(egg::ggarrange(p1, p2, ncol = 2), file = paste0(save_path, '/pseudotime_floxed_dAT2.png'), width = 8, height = 4)


### barcdoe & density plot
### check imported library
### egg, ggpubr shows different behavior (but, why?)
### use egg

### figure plot
# pseudotime barplot across celltype, dpi
# AT1 only
library(ggridges)
pseudotime$condition = factor(pseudotime$condition, levels = c('dAT2', 'floxed'))

p1 <- ggplot(pseudotime[pseudotime$celltype %in% c('AT1') & pseudotime$dpi == '30dpi', ], aes(x = pseudotime, y = condition, fill = condition)) +
    geom_density_ridges(scale = 3, from = 0, to = 1) + 
    labs(x = "Pseudotime", y = "Density") +
    scale_x_continuous(limits = c(0, 1), breaks = c(0.0, 0.25, 0.5, 0.75, 1)) + 
    scale_fill_manual(values = c('floxed' = 'darkgray', 'dAT2' = '#c00000'))
p3 <- ggplot(pseudotime[pseudotime$celltype == 'AT1' & pseudotime$dpi == '30dpi', ], aes(x = condition, y = pseudotime, fill = condition)) +
    geom_boxplot(alpha = 1, outlier.size = 0.25) + 
    theme_classic() +
    labs(x = "Condition", y = "Pseudotime") +
    scale_fill_manual(values = c('floxed' = 'darkgray', 'dAT2' = '#c00000')) + 
    scale_y_continuous(limits = c(0, 1), breaks = c(0.0, 0.25, 0.5, 0.75, 1)) +
    coord_flip() + NoLegend()
p1m <- p1 + theme_classic() + NoLegend() + theme(
      plot.title = element_text(size = 0),
      axis.text.x = element_text(size = 0),
      axis.title.x = element_text(size = 0), 
      axis.title.y = element_text(size = 0)
)
p3m <- p3 + NoLegend() + theme(
      plot.title = element_text(size = 0),
      axis.title.y = element_text(size = 0), 
      axis.title.x = element_text(size = 10),
      axis.text.x = element_text(size = 8)
)
p_AT1 <- as.ggplot(egg::ggarrange(p1m, p3m, nrow = 2, heights = c(2, 1)))

ggsave(p_AT1, file = paste0(save_path, '/AT1_pseudotime.png'), width = 4, height = 4)
ggsave(p_AT1, file = paste0(save_path, '/AT1_pseudotime.pdf'), width = 4, height = 4)
ggplot2pptx(p_AT1, 4, 4, paste0(save_path, '/AT1_pseudotime.pptx'))

t.test(pseudotime[pseudotime$celltype %in% c('AT1') & pseudotime$dpi == '30dpi' & pseudotime$condition == 'floxed', 'pseudotime'], 
        pseudotime[pseudotime$celltype %in% c('AT1') & pseudotime$dpi == '30dpi' & pseudotime$condition == 'dAT2', 'pseudotime'])
