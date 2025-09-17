### conda activate project_lung_exercise_R

library(Seurat)
library(magrittr)
library(data.table)
library(ggplot2)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(ggsignif)
library(dplyr)

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

save_path = '/home/sjcho/projects/AT2_MHC2/code_for_publish/20250806_publish/not_traced_14dpi_AT2_AT2.prolif._inflammatory_cycling_ttest/outs'
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

########################################################## fig 4.c. 14dpi AT2 inflammatory response
library(GO.db)
library(fgsea)
library(org.Mm.eg.db)
library(msigdbr)
h_gene_sets = msigdbr(species = "mouse", category = "H")
child.terms <- c(GOBPCHILDREN[['GO:0006954']], GOBPCHILDREN[['GO:0007346']], GOBPCHILDREN[['GO:0044770']]) 
# 0006954 : inflammatory response / 0007346 : regulation of mitotic cell cycle / 0044770 : cell cycle phase transition

genes <- AnnotationDbi::select(org.Mm.eg.db,
                                keys = child.terms,
                                columns = c("GO", "SYMBOL"),
                                keytype = "GO")
go_terms <- AnnotationDbi::select(GO.db,
                                  keys = unique(genes$GO),
                                  columns = c("GOID", "TERM"),
                                  keytype = "GOID")
inflammatory_terms = go_terms[1:12, 'TERM']
cellcycle_terms = go_terms[13:30, 'TERM']

genes <- merge(genes, go_terms, by.x = "GO", by.y = "GOID", all.x = TRUE)
genes <- genes[!is.na(genes$SYMBOL), ]
GO_list <- split(genes$SYMBOL, genes$TERM)
GO_list$'HM_inflammatory_response' = subset(h_gene_sets, gs_name == 'HALLMARK_INFLAMMATORY_RESPONSE')$gene_symbol
GO_list$'HM_G2M_checkpoinit' = subset(h_gene_sets, gs_name == 'HALLMARK_G2M_CHECKPOINT')$gene_symbol

### delete if GO term only contain one gene
GO_list <- GO_list[sapply(GO_list, length) > 1]
inflammatory_terms = intersect(inflammatory_terms, names(GO_list)); inflammatory_terms = c(inflammatory_terms, 'HM_inflammatory_response')
cellcycle_terms = intersect(cellcycle_terms, names(GO_list)); cellcycle_terms = c(cellcycle_terms, 'HM_G2M_checkpoinit')

for (term in names(GO_list)) {
    so.alv <- AddModuleScore(so.alv, features = list(unlist(GO_list[term])), name = term)
}

############# gene x celltype x condition (bsAb, Combi) dotplot
############# It shows the differential expression of genes across different cell types and conditions

t.test.result.df <- data.frame(
  celltype = character(),
  term = character(),
  p_value = numeric(),
  test_statistic = numeric(),
  dpi = character()
)
for (celltype in c('AT2', 'AT2.prolif.', 'tAT2')) {
    for (dpi in c('14dpi', '30dpi')) {
        control_cells = Cells(so.alv)[so.alv$condition == 'floxed' & so.alv$dpi == dpi & so.alv$majority_voting == celltype]
        condition_cells = Cells(so.alv)[so.alv$condition == 'dAT2' & so.alv$dpi == dpi & so.alv$majority_voting == celltype]
        for (term in names(GO_list)) {
            t.test.result <- t.test(so.alv@meta.data[condition_cells, paste0(term, '1')], so.alv@meta.data[control_cells, paste0(term, '1')])
            p_value <- t.test.result$p.value
            test_statistic <- t.test.result$statistic

            t.test.result.df <- rbind(t.test.result.df, data.frame(
                celltype = celltype,
                term = term,
                p_value = p_value,
                test_statistic = test_statistic,
                dpi = dpi
        ))
        }
    }
}

plot.df <- t.test.result.df %>%
  mutate(
    logp = -log10(p_value),
    # 유의한 경우 색상 구분
    color_group = ifelse(p_value < 0.05, "sig", "nonsig")
  )
plot.df$celltype = factor(plot.df$celltype, levels = c('AT2', 'AT2.prolif.', 'tAT2'))

plot.df.vis <- plot.df
at2_terms <- plot.df.vis %>%
  filter(celltype == "AT2" & dpi == '14dpi') %>%
  mutate(group_order = case_when(
    term %in% inflammatory_terms ~ 1,
    term %in% cellcycle_terms ~ 2,
    TRUE ~ 3
  )) %>%
  arrange(group_order, p_value) %>%
  pull(term)

p <- ggplot(plot.df.vis, aes(
  x = celltype,
  y = term,
  size = logp,
  color = test_statistic
)) +
  geom_point(data = subset(plot.df.vis, color_group == "sig")) +
  geom_point(data = subset(plot.df.vis, color_group == "nonsig"), color = "gray50") +
  scale_size_continuous(range = c(1, 12)) +
  scale_color_distiller(
    palette = "RdBu",
    direction = -1,
    limits = c(-8, 8),
    oob = scales::squish   # ⬅️ 범위 밖 값도 -2, 2에 맞게 squeeze
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    panel.grid = element_blank()
  ) + facet_wrap(~ dpi, nrow = 1) +
  scale_y_discrete(limits = rev(at2_terms)) + 
  labs(
    x = "Celltype",
    size = "-log10(p)",
    color = "t.statistic"
  ) 
ggsave(p, file = paste0(save_path, '/cellcycle_inflammatory_dotplot_ttest.pdf'), width = 8, height = 11)
ggsave(p, file = paste0(save_path, '/cellcycle_inflammatory_dotplot_ttest.png'), width = 8, height = 11)

############# GESA instead of t.test
fgsea_list <- list()

for (celltype in c('AT2', 'AT2.prolif.', 'tAT2')) {
  for (dpi in c('14dpi', '30dpi')) {
    control_cells <- Cells(so.alv)[so.alv$condition == 'floxed' & 
                                   so.alv$dpi == dpi & 
                                   so.alv$majority_voting == celltype]
    condition_cells <- Cells(so.alv)[so.alv$condition == 'dAT2' & 
                                     so.alv$dpi == dpi & 
                                     so.alv$majority_voting == celltype]
    DEGs <- FindMarkers(
      so.alv,
      ident.1 = condition_cells,
      ident.2 = control_cells,
      log2fc.threshold = 0,
      min.pct = 0,
      min.diff.pct = -Inf
    )
    # delete no expression genes
    so.alv@assays$RNA$counts[, c(control_cells, condition_cells)] %>%
      Matrix::rowSums() %>%
      {names(.)[. == 0]} -> no_expression_genes

    DEGs <- DEGs[!rownames(DEGs) %in% no_expression_genes, ]
    log2fc <- DEGs$avg_log2FC
    names(log2fc) <- rownames(DEGs)

    fgseaRes <- fgsea(GO_list, log2fc)
    # NA 제거
    fgseaRes <- fgseaRes[!is.na(fgseaRes$padj), ]
    # NES 순으로 정렬
    fgseaRes <- fgseaRes[order(fgseaRes$NES), ]
    # celltype, dpi 정보 추가
    fgseaRes$celltype <- celltype
    fgseaRes$dpi <- dpi
    # 리스트에 저장
    fgsea_list[[paste(celltype, dpi, sep = "_")]] <- fgseaRes
    print(paste("Processed:", celltype, dpi))
  }
}

fgsea_all <- do.call(rbind, fgsea_list)


# fgsea_all 정리
plot_df <- fgsea_all %>%
  filter(!is.na(padj)) %>%  # padj가 NA인 것은 제거
  mutate(
    group = paste(celltype, dpi, sep = "_"),       # celltype과 dpi 합쳐서 x축에 표시할 그룹 만들기
    log10p = -log10(pval),                      # 점 크기
    NES_color = ifelse(pval <= 0.05, NES, NA)      # p.adj > 0.05면 색을 NA (나중에 회색 지정)
  )

p <- ggplot(plot_df, aes(x = celltype, y = pathway)) +
  geom_point(aes(size = log10p,
                 color = NES_color)) +
  scale_color_distiller(
    palette = "RdBu",
    direction = -1,
    limits = c(-2, 2),
    oob = scales::squish
  ) + facet_wrap(~ dpi, nrow = 1) +
  scale_size_continuous(range = c(1, 6)) +
  labs(
    x = "Celltype",
    y = "Pathway",
    size = "-log10(p-value)",
    color = "NES"
  ) + theme_classic() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1)
  )
ggsave(p, file = paste0(save_path, '/GSEA_celltype_dpi_dotplot.pdf'), width = 8, height = 11)
ggsave(p, file = paste0(save_path, '/GSEA_celltype_dpi_dotplot.png'), width = 8, height = 11)