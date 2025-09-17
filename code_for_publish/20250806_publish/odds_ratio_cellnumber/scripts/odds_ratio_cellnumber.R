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

alv_cols = c(AT2 = "#8192ef", tAT2 = "#dfd13a", AT1 = "#f66b07", AT2.IFN = "#000000", AT2.prolif. = "#8bd690")

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

odds_ratio_df = data.frame()
for (dpi in c('naive', '7dpi', '14dpi', '30dpi')) {
    for (celltype in levels(so.alv$majority_voting)) {
        control_cellnumber = length(Cells(so.alv)[so.alv$condition_dpi == paste0(dpi, '_floxed') & so.alv$majority_voting == celltype])
        dAT2_cellnumber = length(Cells(so.alv)[so.alv$condition_dpi == paste0(dpi, '_dAT2') & so.alv$majority_voting == celltype])
        control_total_cellnumber = length(Cells(so.alv)[so.alv$condition_dpi == paste0(dpi, '_floxed')])
        dAT2_total_cellnumber = length(Cells(so.alv)[so.alv$condition_dpi == paste0(dpi, '_dAT2')])

        control_cellnumber / control_total_cellnumber -> control_proportion
        dAT2_cellnumber / dAT2_total_cellnumber -> dAT2_proportion
        dAT2_proportion / control_proportion -> odds_ratio

        odds_ratio_df = rbind(odds_ratio_df, data.frame(
            dpi = dpi,
            majority_voting = celltype,
            control_proportion = control_proportion,
            dAT2_proportion = dAT2_proportion,
            odds_ratio = odds_ratio
        ))
    }
}

fold_change_df = data.frame()
for (dpi in c('naive', '7dpi', '14dpi', '30dpi')) {
    for (celltype in levels(so.alv$majority_voting)) {
        control_cellnumber = length(Cells(so.alv)[so.alv$condition_dpi == paste0(dpi, '_floxed') & so.alv$majority_voting == celltype])
        dAT2_cellnumber = length(Cells(so.alv)[so.alv$condition_dpi == paste0(dpi, '_dAT2') & so.alv$majority_voting == celltype])

        control_total = length(Cells(so.alv)[so.alv$condition_dpi == paste0(dpi, '_floxed')])
        dAT2_total = length(Cells(so.alv)[so.alv$condition_dpi == paste0(dpi, '_dAT2')])

        fold_change_df = rbind(fold_change_df, data.frame(
            dpi = dpi,
            majority_voting = celltype,
            fold_change = dAT2_cellnumber / control_cellnumber,
            log2_fold_change = log2(dAT2_cellnumber / control_cellnumber + 1e-10),
            odds_ratio = (dAT2_cellnumber / control_cellnumber) / (dAT2_total / control_total)
        ))
    }
}

result <- make.proportion.to.plot(so.alv, row_meta = 'condition_dpi', col_meta = 'majority_voting', colors = alv_cols, add_blackline = TRUE)
ggsave(result[[2]], filename = paste0(save_path, '/alv_majority_voting_condition_dpi_proportion.pdf'), width = 6, height = 4)