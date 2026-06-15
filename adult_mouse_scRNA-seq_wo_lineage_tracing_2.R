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
source('/home/sjcho/yard/functions/R/save_ggplot2_to_ppt.R')
source('/home/sjcho/yard/functions/R/draw_proportion.R')
source('/home/sjcho/projects/AT2_MHC2/20241113_after_cellbender/utils.R')
source('/home/sjcho/yard/functions/R/draw_celltype_wise_dotplot.R')
source('/home/sjcho/yard/functions/R/gsea_sjcho.R')

library(msigdbr)
h_gene_sets = msigdbr(species = "mouse", category = "H")
h_gene_sets$gs_name %>% unique

save_path = '/home/sjcho/projects/AT2_MHC2/20241113_after_cellbender/manuscript_figures/figure5_and_exd7_20250221/outs/'
so.alv <- readRDS('/home/sjcho/projects/AT2_MHC2/20241113_after_cellbender/3.whole_cell_annotation/outs/3.1.epithelial_annotation/alv_annotation.rds')

so.alv$condition_dpi = factor(so.alv$condition_dpi, levels = c('naive_floxed', '7dpi_floxed', '14dpi_floxed', '30dpi_floxed', 
                                                                'naive_dAT2', '7dpi_dAT2', '14dpi_dAT2', '30dpi_dAT2'))

alv_cols = c(AT2 = "#8192ef", tAT2 = "#dfd13a", AT1 = "#f66b07", AT2.INF = "#000000", AT2.prolif. = "#8bd690")

### 1. celltypist annotation
celltypist_result = fread('/home/sjcho/projects/AT2_MHC2/20241113_after_cellbender/manuscript_figures/figure4/figure4_and_exd6_20250225/alv_full_celltypist_result_overcluster1.csv')
celltypist_result = as.data.frame(celltypist_result)
rownames(celltypist_result) <- celltypist_result$V1
celltypist_result <- celltypist_result[, -1]

so.alv$predicted_labels = celltypist_result[Cells(so.alv), 'predicted_labels']
so.alv$predicted_labels = factor(so.alv$predicted_labels, levels = c('AT2', 'tAT2', 'AT1', 'AT2.prolif.', 'AT2.IFN'))
so.alv$majority_voting = celltypist_result[Cells(so.alv), 'majority_voting']
so.alv$majority_voting = factor(so.alv$majority_voting, levels = c('AT2', 'tAT2', 'AT1', 'AT2.prolif.', 'AT2.IFN'))


###################################################### fig4.b Augur for Alveolar epi
so.alv$annotation_dpi = paste0(so.alv$majority_voting, '_', so.alv$dpi)

library(Augur)
so.alv[['cell_type']] <- so.alv[['annotation_dpi']][[1]]
so.alv[['label']] <- so.alv[['condition']][[1]]
augur_default = calculate_auc(so.alv, n_threads = 48, subsample_size = 15, min_cells = 15)
saveRDS(augur_default, file = paste0(save_path, 'exd7.X_alv_Augur.rds'))

df <- data.frame(augur_default$AUC)
df <- df %>%
  mutate(dpi = gsub(".*_(\\d+dpi|naive)$", "\\1", cell_type),
         cell_type = gsub("_(\\d+dpi|naive)$", "", cell_type))

df_wide <- df %>%
  pivot_wider(names_from = cell_type, values_from = auc)

df_long <- df %>%
  mutate(dpi = factor(dpi, levels = unique(dpi))) %>% 
  pivot_wider(names_from = cell_type, values_from = auc) %>%
  pivot_longer(-dpi, names_to = "cell_type", values_to = "AUC")
df_long$cell_type <- factor(df_long$cell_type, levels = c("AT2", "tAT2", "AT1", "cycling.alv.epi", "IFN.response.alv.epi"))
df_long$dpi <- factor(df_long$dpi, levels = c("30dpi", "14dpi", "7dpi", "naive"))

p <- ggplot(df_long, aes(x = cell_type, y = dpi, fill = AUC)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(colors = c('gray95', '#f18e83', 'darkred'), na.value = "darkgray") +
  theme_minimal() +
  labs(title = "AUC Heatmap", x = "Cell Type", y = "DPI") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(p, file = paste0(save_path, 'fig4.b._alv_Augur_heatmap.png'), width = 5, height = 5)
ggsave(p, file = paste0(save_path, 'fig4.b._alv_Augur_heatmap.pdf'), width = 5, height = 5)
ggplot2pptx(p, 5, 5, paste0(save_path, 'fig4.b._alv_Augur_heatmap.pptx'))

so.14.AT2 <- subset(so.alv, cells = Cells(so.alv)[so.alv$dpi == '14dpi' & so.alv$majority_voting == 'AT2'])

########################################################## fig 4.c. 14dpi AT2 inflammatory response
library(GO.db)
library(fgsea)
library(org.Mm.eg.db)
library(msigdbr)
h_gene_sets = msigdbr(species = "mouse", category = "H")
child.terms <- GOBPCHILDREN[['GO:0006954']] # inflammatory response

genes <- AnnotationDbi::select(org.Mm.eg.db,
                                keys = child.terms,
                                columns = c("GO", "SYMBOL"),
                                keytype = "GO")
go_terms <- AnnotationDbi::select(GO.db,
                                  keys = unique(genes$GO),
                                  columns = c("GOID", "TERM"),
                                  keytype = "GOID")
genes <- merge(genes, go_terms, by.x = "GO", by.y = "GOID", all.x = TRUE)
genes <- genes[!is.na(genes$SYMBOL), ]
GO_list <- split(genes$SYMBOL, genes$TERM)
GO_list$'HM_inflammatory_response' = subset(h_gene_sets, gs_name == 'HALLMARK_INFLAMMATORY_RESPONSE')$gene_symbol

library(parallel)
min.cells = 1
p.value.list <- list()
diff.mean.list <- list()
mean_MS_list <- list()
mean_MS_anno_dpi_list <- list()

moduleScores_list <- mclapply(names(GO_list), function(term) {
  if (intersect(rownames(so.alv), GO_list[[term]]) %>% length < 2) {
    return(NULL)
  }
  temp_obj <- AddModuleScore(so.14.AT2, features = list(GO_list[[term]]), name = term)
  temp_obj@meta.data[, paste0(term, '1'), drop = FALSE]
}, mc.cores = 48)

moduleScores_list <- moduleScores_list[!sapply(moduleScores_list, is.null)]
moduleScores_df <- do.call(cbind, moduleScores_list)

control_cell = Cells(so.14.AT2)[so.14.AT2$condition == 'floxed']
dAT2_cell = Cells(so.14.AT2)[so.14.AT2$condition == 'dAT2']

module_results <- mclapply(colnames(moduleScores_df), function(module) {
    p_val <- wilcox.test(moduleScores_df[control_cell, module], 
                        moduleScores_df[dAT2_cell, module])$p.value
    diff_mean <- (mean(moduleScores_df[dAT2_cell, module]) - mean(moduleScores_df[control_cell, module]))
    mean_MS <- mean(moduleScores_df[, module])
    mean_MS_anno_dpi <- mean(moduleScores_df[c(dAT2_cell, control_cell), module])
    list(p.value = p_val, diff.mean = diff_mean, mean_MS = mean_MS, mean_MS_anno_dpi = mean_MS_anno_dpi)
}, mc.cores = 48)

p.value <- sapply(module_results, function(x) x$p.value)
diff.mean <- sapply(module_results, function(x) x$diff.mean)

df.p.diff.mean <- data.frame(p.value = p.value, diff.mean = diff.mean, row.names = colnames(moduleScores_df))
rownames(df.p.diff.mean) <- sub("1$", "", rownames(df.p.diff.mean))
df.p.diff.mean <- df.p.diff.mean[order(df.p.diff.mean$diff.mean), ]

## visualization
library(ggplot2)
library(dplyr)
library(tibble)

df.p.diff.mean <- df.p.diff.mean %>% rownames_to_column("Term")
df.p.diff.mean <- df.p.diff.mean %>% mutate(Term = factor(Term, levels = Term[order(diff.mean)]))

df.p.diff.mean <- df.p.diff.mean %>%
  mutate(logp = -log10(p.value),
         Significance = case_when(
           p.value > 0.05 ~ "not_significant",
           diff.mean > 0 ~ "dAT2 - floxed > 0",
           TRUE ~ "dAT2 - floxed < 0"
         ))

p <- ggplot(df.p.diff.mean, aes(x = Term, y = diff.mean)) +
  geom_point(aes(size = logp, color = Significance)) +
  scale_color_manual(values = c("dAT2 - floxed > 0" = "red", "dAT2 - floxed < 0" = "blue", "not_significant" = "gray50")) +
  coord_flip() + 
  labs(x = "", y = "differential score", size = "-log10(p-value)") + 
  geom_hline(yintercept = 0, color = 'black', linetype = 'dashed') +
  theme_classic() + 
  theme(axis.text.y = element_text(size = 6),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 6)
  ) +
 scale_y_continuous(breaks = c(-0.01, 0.00, 0.02, 0.04, 0.06)) +
 scale_size_continuous(range = c(1, 4))

ggsave(p, file = paste0(save_path, 'fig4.c.14dpi_AT2_inflammatory_response_childterm.png'), width = 6, height = 3)
ggsave(p, file = paste0(save_path, 'fig4.c.14dpi_AT2_inflammatory_response_childterm.pdf'), width = 6, height = 3)
ggplot2pptx(p, file = paste0(save_path, 'fig4.c.14dpi_AT2_inflammatory_response_childterm.pptx'), width = 6, height = 3)

#################################################### fig 4.d.
save_path = '/home/sjcho/projects/AT2_MHC2/20241113_after_cellbender/manuscript_figures/figure5_and_exd7_20250221/fig5.test.negative_regulation_of_inflammatory_response/'
so.cd4 <- readRDS('/home/sjcho/projects/AT2_MHC2/20241113_after_cellbender/4.Tcell_study/4.1.Tcell_subset/outs/4.5.full_CD4_transcriptomic_difference/woGD_cCd4.rds')

so.cd4$condition_dpi = paste0(so.cd4$dpi, '_', so.cd4$condition)
so.cd4$condition_dpi = factor(so.cd4$condition_dpi, levels = c('naive_floxed', '7dpi_floxed', '14dpi_floxed', '30dpi_floxed', 
                                                                'naive_dAT2', '7dpi_dAT2', '14dpi_dAT2', '30dpi_dAT2'))

Tcell_cols_order = c('Cd4.T.naive-like', 'Treg', 'Treg.prolif.', 'Cd4.T.Th1.prolif.', 'Cd4.T.Th1', 'Cd4.T.Th17', 'Cd4.T.cm')
Tcell_cols = c('#67cc7c', '#9d2c5c', '#016a01', '#caa1dd', '#53b1db', '#f90068', '#ffa82d'); names(Tcell_cols) = Tcell_cols_order
###### differential signature score
library(GO.db)
library(fgsea)
library(msigdbr)
library(org.Mm.eg.db)

# GOBPOFFSPRING : include all the childs
# GOBPCHILDREN : include all the lst level children (not childs's children)

child.terms3 <- GOBPCHILDREN[['GO:0050728']] # negative regulation of inflammatory response

genes <- AnnotationDbi::select(org.Mm.eg.db,
                                keys = child.terms3,
                                columns = c("GO", "SYMBOL"),
                                keytype = "GO")
go_terms <- AnnotationDbi::select(GO.db,
                                  keys = unique(genes$GO),
                                  columns = c("GOID", "TERM"),
                                  keytype = "GOID")
genes <- merge(genes, go_terms, by.x = "GO", by.y = "GOID", all.x = TRUE)

genes <- genes[!is.na(genes$SYMBOL), ]
GO_list <- split(genes$SYMBOL, genes$TERM)
GO_list[['negative regulation of inflammatory response']] = AnnotationDbi::select(org.Mm.eg.db,
                                                            keys = 'GO:0050728',
                                                            columns = c("GO", "SYMBOL"),
                                                            keytype = "GO")$SYMBOL


so.cd4.14 <- subset(so.cd4, cells = Cells(so.cd4)[so.cd4$dpi == '14dpi'])

moduleScores_list <- mclapply(names(GO_list), function(term) {
  if (intersect(rownames(so.cd4.14), GO_list[[term]]) %>% length < 2) {
    return(NULL)
  }
  temp_obj <- AddModuleScore(so.cd4.14, features = list(GO_list[[term]]), name = term)
  temp_obj@meta.data[, paste0(term, '1'), drop = FALSE]
}, mc.cores = 48)

#### condition wise difference
moduleScores_list <- moduleScores_list[!sapply(moduleScores_list, is.null)]
moduleScores_df <- do.call(cbind, moduleScores_list)

library(parallel)
min.cells = 1
p.value.list <- list()
t.statistic.list <- list()
diff.mean.list <- list()
mean_MS_list <- list()
mean_MS_anno_dpi_list <- list()
for (celltype in unique(so.cd4.14$annotation)) {
    control_cell = Cells(so.cd4.14)[so.cd4.14$annotation == celltype & so.cd4.14$condition == 'floxed']
    dAT2_cell = Cells(so.cd4.14)[so.cd4.14$annotation == celltype & so.cd4.14$condition == 'dAT2']

    if (control_cell %>% length < 15 | dAT2_cell %>% length < 15) {
        p.value.list[[celltype]] <- NA
        diff.mean.list[[celltype]] <- NA
    }
    else {
        module_results <- mclapply(colnames(moduleScores_df), function(module) {
          t.test.result <- t.test(moduleScores_df[dAT2_cell, module], 
                              moduleScores_df[control_cell, module])
          p_val <- t.test.result$p.value
          t.statistic <- t.test.result$statistic
          diff_mean <- (mean(moduleScores_df[control_cell, module]) - mean(moduleScores_df[dAT2_cell, module]))
          mean_MS <- mean(moduleScores_df[, module])
          mean_MS_anno_dpi <- mean(moduleScores_df[c(control_cell, dAT2_cell), module])
          list(p.value = p_val, t.statistic = t.statistic, diff.mean = diff_mean, mean_MS = mean_MS, mean_MS_anno_dpi = mean_MS_anno_dpi)
        }, mc.cores = 48)

        p.value <- sapply(module_results, function(x) x$p.value)
        t.statistic <- sapply(module_results, function(x) x$t.statistic)
        diff.mean <- sapply(module_results, function(x) x$diff.mean)
        mean_MS <- sapply(module_results, function(x) x$mean_MS)
        mean_MS_anno_dpi <- sapply(module_results, function(x) x$mean_MS_anno_dpi)

        p.value.list[[celltype]] <- p.value
        t.statistic.list[[celltype]] <- t.statistic
        diff.mean.list[[celltype]] <- diff.mean
        mean_MS_list[[celltype]] <- mean_MS
        mean_MS_anno_dpi_list[[celltype]] <- mean_MS_anno_dpi
    }
}

p.value.df <- as.data.frame(p.value.list)
rownames(p.value.df) <- colnames(moduleScores_df)
# p.value.df <- p.value.df[, order_level]
is.sig <- p.value.df < 0.05

t.statistic.df <- as.data.frame(t.statistic.list)
rownames(t.statistic.df) <- colnames(moduleScores_df)
# diff.mean.df <- diff.mean.df[, order_level]

t.statistic.sig <- t.statistic.df
t.statistic.sig[!is.sig] <- NA # only significant values

# find rows with NA only
rows_with_na_only <- which(apply(t.statistic.sig, 1, function(x) all(is.na(x))))

t.stastics.order <- t.statistic.sig
rownames(t.stastics.order) <- sub("1$", "", rownames(t.stastics.order))
ordered_row <- rownames(t.stastics.order)
ordered_row <- c(ordered_row[length(ordered_row)], ordered_row[1:length(ordered_row) - 1])
t.stastics.order <- t.stastics.order[ordered_row, ]

### dotplot
library(ggplot2)
library(reshape2)

p.value.df <- rownames_to_column(p.value.df, "Pathway")
t.statistic.df <- rownames_to_column(t.statistic.df, "Pathway")
p_long <- melt(p.value.df[c(1,2,3,4,7), ], id.vars = "Pathway", variable.name = "celltype", value.name = "pvalue")
t_long <- melt(t.statistic.df[c(1,2,3,4,7), ], id.vars = "Pathway", variable.name = "celltype", value.name = "tstat")

df <- merge(p_long, t_long, by = c("Pathway", "celltype"))
df$negLogP <- -log10(df$pvalue)

df$celltype <- factor(df$celltype, levels = c('Treg', 'Cd4.T.cm', 'Cd4.T.Th1', 'Cd4.T.naive.like', 'Cd4.T.Th17', 'Treg.prolif.', 'Cd4.T.Th1.prolif.'))

p <- ggplot(df, aes(x = celltype, y = Pathway)) +
    geom_point(data = subset(df, pvalue > 0.05),
                aes(x = celltype, y = Pathway, size = negLogP),
                color = "gray75") +
    geom_point(data = subset(df, pvalue <= 0.05),
                aes(x = celltype, y = Pathway, size = negLogP, color = tstat)) +
    scale_color_gradient2(
        low = "#0009f8", 
        mid = "white",   
        high = "#ff0000", 
        midpoint = 0,
        na.value = "gray50",
        limits = c(-1 * df$tstat %>% abs %>% max(na.rm = T), df$tstat %>% abs %>% max(na.rm = T))) +
    scale_size_continuous(range = c(0.1, 10)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
        size = "-log10(p-value)", 
        color = "t-statistic")
ggsave(p, file = paste0(save_path, 'fig4.d.negative.regulation.of.inflammatory.png'), width = 10, height = 8)
ggsave(p, file = paste0(save_path, 'fig4.d.negative.regulation.of.inflammatory.pdf'), width = 10, height = 8)
ggplot2pptx(p, 10, 8, paste0(save_path, 'fig4.d.negative.regulation.of.inflammatory.pptx'))

######################################### fig 4.g
### Code for fig 4.g is described in the chi-square test file

########################################## fig 4.h
library(readxl)
IL10_markers <- read_excel('/home/sjcho/projects/AT2_MHC2/20241113_after_cellbender/intestine organoid with IL-10.xlsx')
IL10_markers <- data.frame(IL10_markers)
IL10_genes = subset(IL10_markers, FDR < 0.05 & Log2.Fold.change.of.means > 0.25)$Gene.symbol

IL10_markers <- IL10_markers %>%
  arrange(desc(Log2.Fold.change.of.means)) %>%
  mutate(rank = row_number(),
         size_by_fdr = -log10(FDR))
p <- ggplot(IL10_markers, aes(x = rank, y = Log2.Fold.change.of.means)) +
  geom_point(aes(size = size_by_fdr), alpha = 0.7) +
  scale_size_continuous(name = "-log10(FDR)") +
  labs(x = "Rank", y = "Log2 Fold Change") +
  geom_hline(yintercept = log2(1.25), linetype = "dashed", color = "red") +
  theme_classic()
ggsave(p, file = paste0(save_path, 'fig4.h.IL10.response.png'), width = 15, height = 4)
ggsave(p, file = paste0(save_path, 'fig4.h.IL10.response.pdf'), width = 15, height = 4)
ggplot2pptx(p, 8, 4, paste0(save_path, 'fig4.h.IL10.response.pptx'))
leg <- get_legend(p)
ggplot2pptx(as.ggplot(leg), 1, 4, paste0(save_path, 'fig4.h.IL10.response.legend.pptx'))

##################################################### exd 6.a
library(future)
library(future.apply)
library(pbapply)
library(tidyverse)

plan(multisession, workers = 48)
options(future.globals.maxSize = 300 * 1024^3)

dpi_list <- c('naive', '7dpi', '14dpi', '30dpi')
celltype_list <- c('AT2', 'tAT2', 'AT1', 'AT2.prolif.', 'AT2.IFN')
combinations <- expand.grid(dpi = dpi_list, celltype = celltype_list, stringsAsFactors = FALSE)
log2FC_threshold <- log2(1.5)

results <- pblapply(1:nrow(combinations), function(i) {
  dpi <- combinations$dpi[i]
  celltype <- combinations$celltype[i]
  control_cell <- Cells(so.alv)[so.alv$condition == 'floxed' & so.alv$dpi == dpi & so.alv$majority_voting == celltype]
  dAT2_cell <- Cells(so.alv)[so.alv$condition == 'dAT2' & so.alv$dpi == dpi & so.alv$majority_voting == celltype]
  
  result_key <- paste0(celltype, '_', dpi)
  is_null <- length(control_cell) <= 15 | length(dAT2_cell) <= 15
  
  if (!is_null) {
    markers <- FindMarkers(so.alv, ident.1 = control_cell, ident.2 = dAT2_cell, min.pct = 0, logfc.threshold = 0)
    sig_genes <- sum((markers$avg_log2FC > log2FC_threshold | markers$avg_log2FC < -log2FC_threshold) & markers$p_val_adj < 0.05)
    cell_count <- length(control_cell) + length(dAT2_cell)
    deg_per_cell <- sig_genes / cell_count
  } else {
    deg_per_cell <- 0
    sig_genes <- 0
    cell_count <- 0
  }
  
  return(list(
    celltype = celltype,
    dpi = dpi,
    is_null = is_null,
    deg_per_cell = deg_per_cell,
    sig_genes = sig_genes,
    cell_count = cell_count
  ))
})

deg_long <- do.call(rbind, lapply(results, function(r) {
  data.frame(
    cell_type = r$celltype,
    dpi = r$dpi,
    DEG_per_cell = r$deg_per_cell,
    is_null = r$is_null,
    sig_genes = r$sig_genes,
    cell_count = r$cell_count
  )
}))

deg_long$DEG_per_cell[deg_long$is_null] <- NA

deg_long$dpi <- factor(deg_long$dpi, levels = c("30dpi", "14dpi", "7dpi", "naive"))
deg_long$cell_type <- factor(deg_long$cell_type, levels = c("AT2", "tAT2", "AT1", "AT2.prolif.", "AT2.IFN"))

p <- ggplot(deg_long, aes(x = cell_type, y = dpi, fill = DEG_per_cell)) +
  geom_tile(color = "black") +
  scale_fill_gradientn(colors = c('gray95', '#83a1f1', '#0d349e'), na.value = "darkgray") +
  theme_minimal() +
  labs(x = "Cell Type", y = "DPI", fill = "DEG per Cell") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(p, file = paste0(save_path, "exd6.a.DEG_per_cell_heatmap.pdf"), width = 5, height = 5)
ggsave(p, file = paste0(save_path, "exd6.a.DEG_per_cell_heatmap.png"), width = 5, height = 5)
ggplot2pptx(p, 5, 5, paste0(save_path, "exd6.a.DEG_per_cell_heatmap.pptx"))

########################### exd 6.b
Tcell_marker_order = c('naive-like Cd4 T cell', 'Cd4 Tem', 'Th1 Cd4 T cell', 'Th17 Cd4 T cell', 'effector Treg', 'naive-like Treg', 'Treg.prolif.', 'Th1 Cd4 T.prolif.')
markers <- c('Sell', 'Ccr7', 'Cd44', 'Itgb1', 'Tbx21', 'Ifng', 'Il17a', 'Ccr6', 'Foxp3', 'Cxcr3', 'Mki67', 'Top2a')

library(RColorBrewer)
so.cd4$marker_order = factor(so.cd4$annotation, levels = Tcell_marker_order)
p <- DotPlot(so.cd4, features = markers, group.by = 'marker_order', cols = Tcell_cols) + NoAxes() + theme(plot.title = element_text(size = 0))
p <- p + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p <- p + scale_color_gradientn(colors = rev(brewer.pal(9, "RdBu")))
ggsave(paste0(save_path, 'exd6.b.cd4_annotation_markers.png'), p, width = 6, height = 4)
ggsave(paste0(save_path, 'exd6.b.cd4_full_annotation_markers.pdf'), p, width = 6, height = 4)
ggplot2pptx(p, 6, 4, paste0(save_path, 'exd6.b.cd4_full_annotation_markers.pptx'))

############################ exd 6.c
p <- DimPlot(so.cd4, group.by = 'annotation', split.by = 'condition_dpi', cols = Tcell_cols, pt.size = 1.5, ncol = 4) + NoLegend()
p <- p + theme(plot.title = element_text(size = 0)) & theme(strip.text.x = element_text(size = 0, face = "bold"))
p <- p & NoAxes()
ggsave(p, file = paste0(save_path, 'exd6.c.annotation_split.png'), width = 13, height = 6)
ggsave(p, file = paste0(save_path, 'exd6.c.annotation_split.pdf'), width = 13, height = 6)
ggplot2pptx(p, 13, 6, paste0(save_path, 'exd6.c.annotation_split.pptx'))
