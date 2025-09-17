### dot plot for 14dpi AT2
### the value is differential score between AT2, MHC2
set.seed(42)

library(Seurat)
library(scater)
library(SingleCellExperiment)
library(magrittr)
library(ggplot2)
library(stringr)
library(dplyr)
library(tidyr)
library(data.table)

save_path = '/home/sjcho/projects/AT2_MHC2/code_for_publish/20250806_publish/anti-inflammatory_Tcell/outs'

# load custom functions
source('/home/sjcho/yard/functions/R/seurat_count_to_normalization.R')
source('/home/sjcho/yard/functions/R/FeaturePlot_sjcho.R')
source('/home/sjcho/yard/functions/R/save_ggplot2_to_ppt.R')
source('/home/sjcho/yard/functions/R/draw_proportion.R')
source('/home/sjcho/projects/AT2_MHC2/20241113_after_cellbender/utils.R')
source('/home/sjcho/yard/functions/R/draw_celltype_wise_dotplot.R')
source('/home/sjcho/yard/functions/R/draw_proportion.R')
source('/home/sjcho/projects/AT2_MHC2/20241113_after_cellbender/manuscript_figures/chi-square_basal_exp/chi-sqaure_with_sampling.R')

so.cd4 <- readRDS('/home/sjcho/projects/AT2_MHC2/20241113_after_cellbender/4.Tcell_study/4.1.Tcell_subset/outs/4.5.full_CD4_transcriptomic_difference/woGD_cCd4.rds')

so.cd4$condition_dpi = paste0(so.cd4$dpi, '_', so.cd4$condition)
so.cd4$condition_dpi = factor(so.cd4$condition_dpi, levels = c('naive_floxed', '7dpi_floxed', '14dpi_floxed', '30dpi_floxed', 
                                                                'naive_dAT2', '7dpi_dAT2', '14dpi_dAT2', '30dpi_dAT2'))

so.cd4$annotation[so.cd4$annotation == 'Cxcr3+ Treg'] = 'Treg'
so.cd4$annotation[so.cd4$annotation == 'naive-like Treg'] = 'Treg'

so.cd4@meta.data[so.cd4$annotation == 'Cd4 Tcm', 'annotation'] = 'Cd4.T.cm'
so.cd4$annotation[so.cd4$annotation == 'Th1 Cd4 T cell'] = 'Cd4.T.Th1'
so.cd4$annotation[so.cd4$annotation == 'Th17 Cd4 T cell'] = 'Cd4.T.Th17'
so.cd4$annotation[so.cd4$annotation == 'cycling Cd4 T cell'] = 'Cd4.T.Th1.prolif.'
so.cd4$annotation[so.cd4$annotation == 'naive-like Cd4 T cell'] = 'Cd4.T.naive-like'
so.cd4$annotation[so.cd4$annotation == 'Proliferating Treg'] = 'Treg.prolif.'

Tcell_order = c( 'Treg', 'Treg.prolif.', 'Cd4.T.Th1', 'Cd4.T.Th1.prolif.', 'Cd4.T.naive-like', 'Cd4.T.Th17', 'Cd4.T.cm')
so.cd4$annotation = factor(so.cd4$annotation, levels = Tcell_order)

###### differential signature score
library(GO.db)
library(fgsea)
library(msigdbr)
library(org.Mm.eg.db)

# GOBPOFFSPRING : include all the childs
# GOBPCHILDREN : include all the lst level children (not childs's children)

GO_list <- list()
GO_terms <- list()
for (GO in c('GO:0050728', 'GO:0050729')) { # 
    child_terms <- GOBPCHILDREN[[GO]]

    genes <- AnnotationDbi::select(org.Mm.eg.db,
                                   keys = child_terms,
                                   columns = c("GO", "SYMBOL"),
                                   keytype = "GO")
    go_terms <- AnnotationDbi::select(GO.db,
                                      keys = unique(genes$GO),
                                      columns = c("GOID", "TERM"),
                                      keytype = "GOID")
    genes <- merge(genes, go_terms, by.x = "GO", by.y = "GOID", all.x = TRUE)
    genes <- genes[!is.na(genes$SYMBOL), ]
    if (is.null(GO_list)) {
        GO_list <- split(genes$SYMBOL, genes$TERM)
    } else {
        GO_list <- c(GO_list, split(genes$SYMBOL, genes$TERM))
    }
    GO_terms[[GO]] <- go_terms$TERM
}

### subset
so.cd4 <- AddModuleScore(so.cd4, features = GO_list, name = names(GO_list))
so.cd4.control <- subset(so.cd4, condition == 'floxed')
so.cd4.control$annotaiton_dpi = paste0(so.cd4.control$annotation, '_', so.cd4.control$dpi)
levels = c()
for (dpi in levels(so.cd4.control$dpi)) {
    for (annotation in levels(so.cd4.control$annotation)) {
        levels <- c(levels, paste0(annotation, '_', dpi))
    }
}
so.cd4.control$annotaiton_dpi <- factor(so.cd4.control$annotaiton_dpi, levels = levels)

########## proportion
Tcell_cols_order = c('Cd4.T.naive-like', 'Treg', 'Treg.prolif.', 'Cd4.T.Th1.prolif.', 'Cd4.T.Th1', 'Cd4.T.Th17', 'Cd4.T.cm')
Tcell_cols = c('#67cc7c', '#9d2c5c', '#016a01', '#caa1dd', '#53b1db', '#f90068', '#ffa82d'); names(Tcell_cols) = Tcell_cols_order

result <- make.proportion.to.plot(so.cd4, col_meta = 'annotation', row_meta = 'condition_dpi', colors = Tcell_cols, add_blackline = TRUE)
ggsave(result[[2]], file = paste0(save_path, '/CD4_proportion.png'), width = 8, height = 5)
ggsave(result[[2]], file = paste0(save_path, '/CD4_proportion.pdf'), width = 8, height = 5)
ggplot2pptx(result[[2]], file = paste0(save_path, '/CD4__proportion.pptx'), width = 8, height = 5)

########## signature score
p <- DotPlot(so.cd4.control, features = paste0(names(GO_list), seq(1, length(names(GO_list)))), group.by = 'dpi')
# RdBu
p <- p + scale_colour_gradientn(colours = rev(brewer.pal(11, "RdBu")))
p <- p + coord_flip()
ggsave(p, file = paste0(save_path, '/control_only_dpi_wise.png'), width = 12, height = 10)

df <- p$data
heatmap_plot <- ggplot(df, aes(x = id, y = features.plot, fill = avg.exp.scaled)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(colours = rev(brewer.pal(11, "RdBu"))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave(heatmap_plot, file = paste0(save_path, '/heatmap_control_only_dpi_wise.png'), width = 10, height = 5)

p <- DotPlot(so.cd4.control, features = paste0(names(GO_list), seq(1, length(names(GO_list)))), group.by = 'annotaiton_dpi')
# RdBu
p <- p + scale_colour_gradientn(colours = rev(brewer.pal(11, "RdBu")))
p <- p + coord_flip()
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave(p, file = paste0(save_path, '/control_only_dpi_and_annotation_wise.png'), width = 24, height = 10)

df <- p$data
heatmap_plot <- ggplot(df, aes(x = id, y = features.plot, fill = avg.exp.scaled)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(colours = rev(brewer.pal(11, "RdBu"))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

ggsave(heatmap_plot, file = paste0(save_path, '/heatmap_control_only_dpi_and_annotation_wise.png'), width = 16, height = 5)