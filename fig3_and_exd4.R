# run slignshot
### conda activate project_lung_exercise_R
# based on https://bioconductor.org/packages/devel/bioc/vignettes/slingshot/inst/doc/vignette.html#upstream-analysis

library(Seurat)
library(magrittr)
library(data.table)
library(stringr)
library(rstatix)
library(pheatmap)

library(ggplot2)
library(ggsignif)
library(dplyr)
library(RColorBrewer)
library(ggridges)
library(egg)
library(ggpubr)
library(viridis)
library(ggplotify)

# load custom functions
source('/home/sjcho/yard/functions/R/seurat_count_to_normalization.R')
source('/home/sjcho/yard/functions/R/FeaturePlot_sjcho.R')
source('/home/sjcho/yard/functions/R/save_ggplot2_to_ppt.R')
source('/home/sjcho/yard/functions/R/draw_proportion.R')
source('/home/sjcho/projects/AT2_MHC2/20241113_after_cellbender/utils.R')
source('/home/sjcho/yard/functions/R/draw_celltype_wise_dotplot.R')
source('/home/sjcho/yard/functions/R/draw_volcano.R')
source('/home/sjcho/projects/AT2_MHC2/20241113_after_cellbender/5.cellrank2/utils.cellrnak2.R')

save_path = '/home/sjcho/projects/AT2_MHC2/20241113_after_cellbender/manuscript_figures/figure4/figure4_more_visually_intuitive_outs/'
so.alv <- readRDS('/home/sjcho/projects/AT2_MHC2/20241113_after_cellbender/3.whole_cell_annotation/outs/3.1.epithelial_annotation/alv_annotation.rds')


so.alv$condition_dpi = factor(so.alv$condition_dpi, levels = c('naive_floxed', '7dpi_floxed', '14dpi_floxed', '30dpi_floxed', 
                                                                'naive_dAT2', '7dpi_dAT2', '14dpi_dAT2', '30dpi_dAT2'))

alv_cols = c(AT2 = "#8192ef", tAT2 = "#dfd13a", AT1 = "#f66b07", AT2.IFN = "#000000", AT2.prolif. = "#8bd690")

### celltypist prediction in python
# library(sceasy)
# library(reticulate)
# use_condaenv('project_lung_exercise_R')
# loompy <- reticulate::import('loompy')
# so.alv <- FindClusters(so.alv, resolution = 1.025)

# p <- DimPlot(so.alv, group.by = 'RNA_snn_res.1.025', split.by = 'condition_dpi', pt.size = 1.5, ncol = 4, label = TRUE) + NoLegend() & NoAxes()
# p <- p + theme(plot.title = element_text(size = 0)) & theme(strip.text.x = element_text(size = 0, face = "bold"))
# ggsave(p, file = paste0(save_path, 'test.png'), width = 13, height = 6)

# so.alv[["RNA"]] <- as(so.alv[["RNA"]], "Assay")
# sceasy::convertFormat(so.alv, from="seurat", to="anndata", drop_single_values=FALSE,
#                        outFile= paste0(save_path,'w_cycling_IFN.h5ad'))

### 1. celltypist annotation
saveRDS(so.alv, file = paste0('/home/sjcho/projects/AT2_MHC2/20241113_after_cellbender/manuscript_figures/', 'inhouse_adult_alveoalr_epithelium.rds'))

save_path = '/home/sjcho/projects/AT2_MHC2/20241113_after_cellbender/manuscript_figures/figure4/figure4_and_exd6_20250225/'

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

p <- DimPlot(so.alv, group.by = 'majority_voting', pt.size = 0.5, cols = alv_cols) + NoAxes() + NoLegend()
p <- p + theme(plot.title = element_text(size = 0))
ggsave(p, file = paste0(save_path, 'exd4.Cmajority_voting.png'), width = 4, height = 4)
ggsave(p, file = paste0(save_path, 'exd4.C.majority_voting.pdf'), width = 4, height = 4)
ggplot2pptx(p, 4, 4, paste0(save_path, 'exd4.C.majority_voting.pptx'))

### 2. pseudotime
so.alv.woInter_woPro <- subset(so.alv, cells = setdiff(Cells(so.alv), Cells(so.alv)[so.alv$majority_voting %in% c('AT2.IFN', 'AT2.prolif.')]))
so.alv.woInter_woPro <- log_normalize(so.alv.woInter_woPro, save_path, 'n_pcs_8', 8, 1000)

### annotation
p <- DimPlot(so.alv.woInter_woPro, group.by = 'majority_voting', cols = alv_cols, pt.size = 0.75) + NoAxes() + NoLegend()
p <- p + theme(plot.title = element_text(size = 0))
ggsave(p, file = paste0(save_path, 'figure3.i.majority_voting_wo_pro_IFN.png'), width = 4, height = 4)
ggsave(p, file = paste0(save_path, 'figure3.i.majority_voting_wo_pro_IFN.pdf'), width = 4, height = 4)
ggplot2pptx(p, 4, 4, paste0(save_path, 'figure3.i.majority_voting_wo_pro_IFN.pptx'))

leg <- get_legend(DimPlot(so.alv.woInter_woPro, group.by = 'majority_voting', cols = alv_cols))
ggplot2pptx(as_ggplot(leg), 1, 6, paste0(save_path, 'figure3.i.majority_voting_wo_pro_IFN_legend.pptx'))


# ## palantir in python
# library(sceasy)
# library(reticulate)
# use_condaenv('project_lung_exercise_R')
# loompy <- reticulate::import('loompy')
# so.alv.woInter_woPro[["RNA"]] <- as(so.alv.woInter_woPro[["RNA"]], "Assay")
# sceasy::convertFormat(so.alv.woInter_woPro, from="seurat", to="anndata", drop_single_values=FALSE,
#                        outFile= paste0(save_path,'wo_cycling_IFN.h5ad'))

# so.alv[["RNA"]] <- as(so.alv[["RNA"]], "Assay")
# sceasy::convertFormat(so.alv, from="seurat", to="anndata", drop_single_values=FALSE,
#                        outFile= paste0(save_path,'w_cycling_IFN.h5ad'))


### load palantir result
branch <- fread('/home/sjcho/projects/AT2_MHC2/20241113_after_cellbender/manuscript_figures/figure4/figure4_and_exd6_20250225/palantir/n_pcs8_branch_masks.csv', header = T)
branch <- as.data.frame(branch)
rownames(branch) <- branch$V1
branch <- branch[, -1]

pseudotime <- read.csv('/home/sjcho/projects/AT2_MHC2/20241113_after_cellbender/manuscript_figures/figure4/figure4_and_exd6_20250225/palantir/n_pcs8_palantir_pseudotime.csv')
colnames(pseudotime) = c('cellname', 'pseudotime')
rownames(pseudotime) = pseudotime$cellname
pseudotime$condition = so.alv.woInter_woPro$condition[rownames(pseudotime)]; pseudotime$condition = factor(pseudotime$condition, levels = c('dAT2', 'floxed'))
pseudotime$celltype = so.alv.woInter_woPro$majority_voting[rownames(pseudotime)]; pseudotime$celltype = factor(pseudotime$celltype, levels = c('AT2', 'tAT2', 'AT1'))
pseudotime$dpi = so.alv.woInter_woPro$dpi[rownames(pseudotime)]
pseudotime$branch_tAT2 = branch[rownames(pseudotime), '30dpi_dAT2_GAGGCCTGTCCACAGC-1']
pseudotime$branch_AT1 = branch[rownames(pseudotime), '30dpi_floxed_TCCCAGTTCCTACAAG-1']

so.alv.woInter_woPro$pseudotime = pseudotime[colnames(so.alv.woInter_woPro), 'pseudotime']


### pseudotime visualization (Fig 4.C)
start_cell = 'naive_dAT2_TAACGACCATGTCTAG-1'
so.start <- subset(so.alv.woInter_woPro, cells = start_cell)

### path probability & pseudotime
probability <- read.csv('/home/sjcho/projects/AT2_MHC2/20241113_after_cellbender/manuscript_figures/figure4/figure4_and_exd6_20250225/palantir/n_pcs8_palantir_fate_probabilities.csv')
rownames(probability) <- probability$X
so.alv.woInter_woPro$path1_prob <- probability[colnames(so.alv.woInter_woPro), 'X30dpi_dAT2_GAGGCCTGTCCACAGC.1']
so.alv.woInter_woPro$path2_prob <- probability[colnames(so.alv.woInter_woPro), 'X30dpi_floxed_TCCCAGTTCCTACAAG.1']

p <- FeaturePlot(so.alv.woInter_woPro, features = 'path1_prob', pt.size = 0.75, order = T) + NoAxes() + NoLegend() + scale_colour_gradientn(colours = colorRampPalette(c("darkgray", 'green', "black"))(100))
p <- p + theme(plot.title = element_text(size = 0))
ggsave(p, file = paste0(save_path, 'fig3.j.l.n_pcs_8_fp_path1_prob.png'), width = 4, height = 4)
ggsave(p, file = paste0(save_path, 'fig3.j.l.n_pcs_8_fp_path1_prob.pdf'), width = 4, height = 4)
ggplot2pptx(p, 4, 4, paste0(save_path, 'fig3.j.l.n_pcs_8_fp_path1_prob.pptx'))

p <- FeaturePlot(so.alv.woInter_woPro, features = 'path2_prob', pt.size = 0.75, order = T) + NoAxes() + NoLegend() + scale_colour_gradientn(colours = colorRampPalette(c("darkgray", 'green', "black"))(100))
p <- p + theme(plot.title = element_text(size = 0))
ggsave(p, file = paste0(save_path, 'fig3.j.r.n_pcs_8_fp_path2_prob.png'), width = 4, height = 4)
ggsave(p, file = paste0(save_path, 'fig3.j.r.n_pcs_8_fp_path2_prob.pdf'), width = 4, height = 4)
ggplot2pptx(p, 4, 4, paste0(save_path, 'fig3.j.r.n_pcs_8_fp_path2_prob.pptx'))

leg <- get_legend(FeaturePlot(so.alv.woInter_woPro, features = 'pseudotime', pt.size = 0.75, order = T) + 
        scale_colour_gradientn(colours = colorRampPalette(c("darkgray", 'green', "black"))(100)))
ggplot2pptx(as_ggplot(leg), 1, 6, paste0(save_path, 'fig3.j.r.n_pcs_8_fp_pseudotime_legend.pptx'))

p <- FeaturePlot(so.alv.woInter_woPro, features = 'pseudotime', pt.size = 0.75, order = T) + NoAxes() + NoLegend() + scale_colour_gradientn(colours = colorRampPalette(c("darkgray", 'green', "black"))(100))
p <- p + geom_point(data = so.start@reductions$umap@cell.embeddings, aes(x = umap_1, y = umap_2), color = 'darkred', size = 4)
p <- p + theme(plot.title = element_text(size = 0))
ggsave(p, file = paste0(save_path, 'fig3.j.t.n_pcs_8_fp_pseudotime.png'), width = 4, height = 4)
ggsave(p, file = paste0(save_path, 'fig3.j.t.n_pcs_8_fp_pseudotime.pdf'), width = 4, height = 4)
ggplot2pptx(p, 4, 4, paste0(save_path, 'fig3.j.t.n_pcs_8_fp_pseudotime.pptx'))


### barcdoe & density plot
### check imported library
### egg, ggpubr shows different behavior (but, why?)
### use egg

### figure plot
# pseudotime barplot across celltype, dpi (fig 4.D)
# AT1 only

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

ggsave(p_AT1, file = paste0(save_path, 'fig3.k.AT1_pseudotime.png'), width = 4, height = 4)
ggsave(p_AT1, file = paste0(save_path, 'fig3.k.AT1_pseudotime.pdf'), width = 4, height = 4)
ggplot2pptx(p_AT1, 4, 4, paste0(save_path, 'fig3.k.AT1_pseudotime.pptx'))

### contour plot (fig 3.i)
condition_name = 'condition'
ref_condition = 'floxed'
exp_condition = 'dAT2'
celltype_name = 'majority_voting'
nbin = 10
type_color = c('floxed' = 'black', 'dAT2' = '#c00000')

contour_plot(subset(so.alv.woInter_woPro, cells = Cells(so.alv.woInter_woPro)[so.alv.woInter_woPro$dpi == '30dpi']), name_of_dimension = 'umap',
            'condition', 'floxed', 'dAT2', 
            'majority_voting', 
            type_color = c('floxed' = 'black', 'dAT2' = '#c00000'), 
            celltype_color = alv_cols, 
            nbin = 5, linewidth = 0.5,
            save_path = save_path, data_name = 'fig3.i.dpi_30dpi', width = 5, height = 5)


### sox9
so.alv <- FindClusters(so.alv, resolution = 1.05)

cols = c('#f8766d', '#e68613', '#cd9600', '#aba300', '#7cae00', '#0cb702', '#00be67', '#00c19a', '#00bfc4', '#00b8e7', '#00a9ff', '#8494ff',
        '#c77cff', '#ed68ed', '#ff61cc', 'darkred')

p <- DimPlot(so.alv, group.by = 'RNA_snn_res.1.05', cols = cols, label = T, pt.size = 2) + NoLegend() + NoAxes()
p <- p + theme(plot.title = element_text(size = 0))
ggsave(p, file = paste0(save_path, 'exd.4.e.1.05.cluster_umap.png'), width = 6, height = 6)
ggsave(p, file = paste0(save_path, 'exd.4.e.1.05.cluster_umap.pdf'), width = 6, height = 6)
ggplot2pptx(p, 6, 6, paste0(save_path, 'exd.4.e.1.05.cluster_umap.pptx'))

leg <- get_legend(DimPlot(so.alv, group.by = 'RNA_snn_res.1.05', label = F) + NoAxes())
ggplot2pptx(as_ggplot(leg), 1, 6, paste0(save_path, 'exd.4.f.1.05.cluster_umap_legend.pptx'))

p <- DimPlot(so.alv, split.by = 'condition_dpi', cells.highlight = Cells(so.alv)[so.alv$RNA_snn_res.1.05 == '15'], cols.highlight = 'darkred', sizes.highlight = 1, pt.size = 1, ncol = 4) + NoLegend() & NoAxes()
p <- p + theme(strip.text.x = element_text(size = 0, face = "bold"))
ggsave(p, file = paste0(save_path, 'exd.4.F.cluster_umap_highlight15.png'), width = 13, height = 6)
ggsave(p, file = paste0(save_path, 'exd.4.F.cluster_umap_highlight15.pdf'), width = 13, height = 6)
ggplot2pptx(p, 13, 6, paste0(save_path, 'exd.4.F.cluster_umap_highlight15.pptx'))

### cluster-wise mean heatmap
### sox9, col14a1, club signature, Scgb1a1
# club cell siganture
load('/home/sjcho/projects/AT2_MHC2/20241113_after_cellbender/manuscript_figures/full_annotation/club_markers.RData')
# save as club_markers

club_markers.sig <- club_markers[club_markers$p_val_adj < 0.01 & club_markers$avg_log2FC > log2(2) & club_markers$pct.1 > 0.1, ]
write.csv(club_markers.sig, paste0(save_path, 'club_markers.sig.csv'))

so.alv <- AddModuleScore(so.alv, features = list(club_markers.sig %>% rownames), name = 'Club_cell_markers')

cluster_info <- so.alv$RNA_snn_res.1.05
sox9_expr <- expm1(so.alv@assays$RNA$data['Sox9', ])
col14a1_expr <- expm1(so.alv@assays$RNA$data['Col14a1', ])
club_signature <- so.alv@meta.data[, 'Club_cell_markers1']
# Scgb1a1_expr <- expm1(so.alv@assays$RNA$data['Cyp4f15', ])

cluster_means_sox9 <- tapply(sox9_expr, cluster_info, mean, na.rm = TRUE)
cluster_means_col14a1 <- tapply(col14a1_expr, cluster_info, mean, na.rm = TRUE)
cluster_means_club_signature <- tapply(club_signature, cluster_info, mean, na.rm = TRUE)
# cluster_means_Scgb1a1 <- tapply(Scgb1a1_expr, cluster_info, mean, na.rm = TRUE)

cluster_means_df <- data.frame(
  Sox9 = cluster_means_sox9,
  Col14a1 = cluster_means_col14a1,
  Club_signature = cluster_means_club_signature)
  # Scgb1a1 = cluster_means_Scgb1a1)

library(RColorBrewer)
my_palette <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(50))
pheatmap(cluster_means_df,
        cluster_rows = FALSE,
        cluster_cols = FALSE,
        show_rownames = TRUE,
        show_colnames = TRUE,
        fontsize = 15,
        display_numbers = FALSE,
        scale = "column",
        color = my_palette,
        breaks = seq(-3, 3, length.out = length(my_palette) + 1)
        # annotation_col = annotation_col,
        # annotation_colors = ann_colors
        # breaks = seq(-3, 3, length.out = 100),
        ) %>% as.ggplot -> p
ggsave(p, filename = paste0(save_path, 'exd4.G.sox9_col14a1_club_heatmap.png'), width = 3, height = 6)
ggsave(p, filename = paste0(save_path, 'exd4.G.sox9_col14a1_club_heatmap.pdf'), width = 3, height = 6)
ggplot2pptx(p, 3, 6, paste0(save_path, 'exd4.G.sox9_col14a1_club_heatmap.pptx'))


### whole cells
so <- readRDS('/home/sjcho/projects/AT2_MHC2/20241113_after_cellbender/manuscript_figures/full_annotation/full_annotation.rds')

so.full$annotation = factor(so.full$annotation, levels = c('T cell', 'Macrophage', 'Neutrophil', 'Erythrocyte', 
                                                            'AT2', 'tAT2', 'AT1', 'AT2.IFN', 'AT2.prolif', 
                                                            'Basal cell', 'Club cell', 'Ciliated cell', 'Serous cell', 'Basal.prolif',
                                                            'Stromal', 'Endothelial'))

alv_cols = c(AT2 = "#8192ef", tAT2 = "#dfd13a", AT1 = "#f66b07", AT2.IFN = "#000000", AT2.prolif = "#8bd690")
cols_full = c('#221677', '#d809e2', '#159b6c', 'red', 
              as.vector(alv_cols), 
              '#b2182b', '#9d9999', '#d6a877', '#128215', '#09359b', 
              '#513951', '#137c74')
names(cols_full) = levels(so.full$annotation)

p <- DimPlot(so.full, group.by = 'annotation', cols = cols_full) + NoAxes() + NoLegend() + theme(plot.title = element_text(size = 0))
ggsave(paste0(save_path, 'full_annotation_wo_legend.png'), p, width = 7, height = 7)
ggsave(paste0(save_path, 'full_annotation_wo_legend.pdf'), p, width = 7, height = 7)
ggplot2pptx(p, 7, 7, paste0(save_path, 'full_annotation_wo_legend.pptx'))

leg <- get_legend(DimPlot(so.full, group.by = 'annotation', cols = cols_full))
ggplot2pptx(as_ggplot(leg), 1, 7, paste0(save_path, 'full_annotation_legend.pptx'))

library(RColorBrewer)
markers <- c('Ptprc', 'Cd3e', 'Cd4', 'Cd68', 'C1qb', 'Fcgr3', 'S100a9', 'Hbb-bt', 'Hba-a2', 
            'Etv5', 'Abca3', 'Napsa', 'Cldn4', 'Tnip3', 'Krt8', 'Pdpn', 'Aqp5', 'Hopx', 
            'Krt5', 'Trp63', 'Scgb3a2', 'Cldn10', 'Foxj1', 'Ltf', 'Msln', 
            'Top2a', 'Mki67', 'Ifi44', 'Iigp1',
            'Pdgfra', 'Col1a1', 'Pecam1', 'Ptprb')
p <- DotPlot(so.full, features = markers, group.by = 'annotation', cols = cols_full) + NoAxes() + theme(plot.title = element_text(size = 0))
p <- p + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p <- p + scale_color_gradientn(colors = rev(brewer.pal(9, "RdBu")))
ggsave(paste0(save_path, 'full_annotation_markers.png'), p, width = 12, height = 6)
ggsave(paste0(save_path, 'full_annotation_markers.pdf'), p, width = 12, height = 6)
ggplot2pptx(p, 12, 6, paste0(save_path, 'full_annotation_markers.pptx'))
