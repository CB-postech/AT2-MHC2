### conda activate project_lung_exercise_R

library(Seurat)
library(magrittr)
library(data.table)
library(ggplot2)
library(reshape2)
library(ggplot2)
library(RColorBrewer)

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

save_path = '/home/sjcho/projects/AT2_MHC2/code_for_publish/20250806_publish/tAT2_pathway/outs'
so.alv <- readRDS('/home/sjcho/projects/AT2_MHC2/20241113_after_cellbender/3.whole_cell_annotation/outs/3.1.epithelial_annotation/alv_annotation.rds')

so.alv$condition_dpi = factor(so.alv$condition_dpi, levels = c('naive_floxed', '7dpi_floxed', '14dpi_floxed', '30dpi_floxed', 
                                                                'naive_dAT2', '7dpi_dAT2', '14dpi_dAT2', '30dpi_dAT2'))

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

result <- make.proportion.to.plot(so.alv, col_meta = 'majority_voting', row_meta = 'condition', colors = alv_cols, add_blackline = TRUE)
ggsave(result[[2]], file = paste0(save_path, '/proportion_majority_voting.png'), width = 4, height = 5)
ggsave(result[[2]], file = paste0(save_path, '/roportion_majority_voting.pdf'), width = 4, height = 5)
ggplot2pptx(result[[2]], file = paste0(save_path, '/proportion_majority_voting.pptx'), width = 4, height = 5)

########################################################## fig 4.c. 14dpi AT2 inflammatory response
library(GO.db)
library(fgsea)
library(org.Mm.eg.db)
library(msigdbr)
h_gene_sets = msigdbr(species = "mouse", category = "H")
GO_inflammatory_response <- GOBPCHILDREN[['GO:0006954']] # inflammatory response
GO_response_to_hypoxia <- GOBPCHILDREN[['GO:0071456']] # response to hypoxia
GO_glycolytic_process <- GOBPCHILDREN[['GO:0006096']] # glycolytic process
GO_cellular_senescence <- GOBPCHILDREN[['GO:0090398']] # cellular senescence
GO_signal_transduction_by_p53 <- GOBPCHILDREN[['GO:0072331']] # signal transduction by p53 class mediator

GO_list <- list()
GO_terms <- list()
for (GO in c('GO:0006954', 'GO:0071456', 'GO:0006096', 'GO:0090398', 'GO:0072331')) { # 
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

GO_list$'HM_inflammatory_response' = subset(h_gene_sets, gs_name == 'HALLMARK_INFLAMMATORY_RESPONSE')$gene_symbol
GO_list$'HM_p53' = subset(h_gene_sets, gs_name == 'HALLMARK_P53_PATHWAY')$gene_symbol
GO_list$'HM_glycolysis' = subset(h_gene_sets, gs_name == 'HALLMARK_GLYCOLYSIS')$gene_symbol
GO_list$'HM_hypoxia' = subset(h_gene_sets, gs_name == 'HALLMARK_HYPOXIA')$gene_symbol

GO_terms_order = c(GO_terms[['GO:0006954']], 'HM_inflammatory_response',
                   GO_terms[['GO:0071456']], 'HM_hypoxia',
                   GO_terms[['GO:0006096']], 'HM_glycolysis',
                   GO_terms[['GO:0090398']], 
                   GO_terms[['GO:0072331']], 'HM_p53')

result_df = data.frame(GO_term = character(),
                          dpi = character(),
                          celltype = character(),
                          p.value = numeric(),
                          t.statistic = numeric(),
                          stringsAsFactors = FALSE)
for (GO_term in names(GO_list)) {
    so.alv <- AddModuleScore(so.alv, features = list(GO_list[[GO_term]]), name = GO_term)
    for (dpi in c('14dpi', '30dpi')) {
        for (celltype in c('AT2', 'tAT2')) {
            floxed.cells = Cells(so.alv)[so.alv$condition == 'floxed' & so.alv$dpi == dpi & so.alv$majority_voting == celltype]
            dAT2.cells = Cells(so.alv)[so.alv$condition == 'dAT2' & so.alv$dpi == dpi & so.alv$majority_voting == celltype]

            t.test.result <- t.test(so.alv@meta.data[dAT2.cells, paste0(GO_term, '1')], 
                                    so.alv@meta.data[floxed.cells, paste0(GO_term, '1')])
            p.value <- t.test.result$p.value
            t.statistic <- t.test.result$statistic

            result_df <- rbind(result_df, 
                          data.frame(GO_term = GO_term, 
                                     dpi = dpi, 
                                     celltype = celltype, 
                                     p.value = p.value, 
                                     t.statistic = t.statistic))
            print(paste0("Processed GO term: ", GO_term, ", dpi: ", dpi, ", celltype: ", celltype))
        }
    }
}
result_df$dpi <- factor(result_df$dpi, levels = c("14dpi", "30dpi"))
result_df$GO_term <- factor(result_df$GO_term, levels = GO_terms_order)
result_df$celltype <- factor(result_df$celltype, levels = c('AT2', 'tAT2'))

max_tstat <- max(abs(result_df$t.statistic), na.rm = TRUE)

# 3) ggplot
p <- ggplot(result_df,
            aes(x = dpi,
                y = GO_term,
                size = -log10(p.value))) +
  geom_point(aes(color = ifelse(p.value <= 0.05, t.statistic, NA))) +
  scale_color_distiller(palette   = "RdBu",
                        direction = -1,
                        name      = "t statistic",
                        na.value  = "gray50",
                        limits    = c(-max_tstat, max_tstat)) +
  scale_size_continuous(range = c(3, 8),
                        name  = expression(-log[10](p.value))) +
  facet_wrap(~ celltype, nrow = 1) +
  theme_classic() +
  theme(strip.text      = element_text(size = 14, face = "bold"),
        axis.text.x     = element_text(size = 12),
        axis.text.y     = element_text(size = 12),
        legend.position = "right") +
  labs(x = NULL, y = NULL)
p <- p + theme(axis.text.x= element_text(angle=45, hjust=0, vjust=0))
ggsave(p, file = paste0(save_path, '/tAT2_pathway_t-test.png'), width = 16, height = 11)
ggsave(p, file = paste0(save_path, '/tAT2_pathway_t-test.pdf'), width = 12, height = 11)
ggplot2pptx(p, file = paste0(save_path, '/tAT2_pathway_t-test.pptx'), width = 16, height = 11)


################ representative genes
Hypoxia_glycolysis = c('Hif1a', 'Pgam1', 'Eno1', 'Aldoa', 'Myc', 'Pkm', 'Slc2a1', 'Cdkn1a', 'Slc16a3', 'Ndrg1')
senescence_p53 = c('Trp53', 'Mdm2', 'Ccnd1', 'Gdf15')

so.alv.AT2.tAT2.14.30 = subset(so.alv, majority_voting %in% c('AT2', 'AT2.prolif.') & dpi %in% c('14dpi', '30dpi'))
so.alv.AT2.tAT2.14.30$anno_condition_dpi = paste0(so.alv.AT2.tAT2.14.30$majority_voting, '_', so.alv.AT2.tAT2.14.30$condition_dpi)
so.alv.AT2.tAT2.14.30$anno_condition_dpi = factor(so.alv.AT2.tAT2.14.30$anno_condition_dpi, levels = c(
    'AT2_14dpi_floxed', 'AT2_14dpi_dAT2', 'AT2_30dpi_floxed', 'AT2_30dpi_dAT2',  'AT2.prolif._14dpi_floxed', 'AT2.prolif._14dpi_dAT2','AT2.prolif._30dpi_floxed', 'AT2.prolif._30dpi_dAT2'
))

p <- DotPlot(so.alv.AT2.tAT2.14.30, features = c(Hypoxia_glycolysis, senescence_p53), group.by = 'anno_condition_dpi')
p <- p + scale_color_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
p <- p + theme(axis.text.x = element_text(angle = 90))
ggsave(p, file = paste0(save_path, '/dotplot_AT2_tAT2_representative_genes.png'), width = 10, height = 6)
ggsave(p, file = paste0(save_path, '/dotplot_AT2_tAT2_representative_genes.pdf'), width = 10, height = 6)
ggplot2pptx(p, file = paste0(save_path, '/dotplot_AT2_tAT2_representative_genes.pptx'), width = 10, height = 6)

result_df = data.frame(gene = character(),
                          dpi = character(),
                          celltype = character(),
                          p.value = numeric(),
                          t.statistic = numeric(),
                          stringsAsFactors = FALSE)
for (gene in c(Hypoxia_glycolysis, senescence_p53)) {   
    for (dpi in c('14dpi', '30dpi')) {
        for (celltype in c('AT2.prolif.', 'AT2', 'tAT2')) {
            floxed.cells = Cells(so.alv)[so.alv$condition == 'floxed' & so.alv$dpi == dpi & so.alv$majority_voting == celltype]
            dAT2.cells = Cells(so.alv)[so.alv$condition == 'dAT2' & so.alv$dpi == dpi & so.alv$majority_voting == celltype]

            t.test.result <- t.test((so.alv@assays$RNA$data[gene, dAT2.cells]), 
                                    (so.alv@assays$RNA$data[gene, floxed.cells]))
            p.value <- t.test.result$p.value
            t.statistic <- t.test.result$statistic

            # wilcox.result <- wilcox.test(expm1(so.alv@assays$RNA$data[gene, dAT2.cells]), 
            #                             expm1(so.alv@assays$RNA$data[gene, floxed.cells]))
            # p.value <- wilcox.result$p.value
            # avg_log2FC = log2(mean(expm1(so.alv@assays$RNA$data[gene, dAT2.cells])) / mean(expm1(so.alv@assays$RNA$data[gene, floxed.cells])))

            result_df <- rbind(result_df, 
                          data.frame(gene = gene, 
                                     dpi = dpi, 
                                     celltype = celltype, 
                                     p.value = p.value, 
                                     t.statistic = t.statistic))
            print(paste0("Processed gene : ", gene, ", dpi: ", dpi, ", celltype: ", celltype))
        }
    }
}
result_df$dpi <- factor(result_df$dpi, levels = c("14dpi", "30dpi"))
result_df$celltype <- factor(result_df$celltype, levels = c('AT2', 'AT2.prolif.', 'tAT2'))

# 3) ggplot
max_tstat <- 8

# 3) ggplot
p <- ggplot(result_df,
            aes(x = dpi,
                y = gene,
                size = -log10(p.value))) +
  geom_point(aes(color = ifelse(p.value <= 0.05, t.statistic, NA))) +
  scale_color_distiller(palette   = "RdBu",
                        direction = -1,
                        # name      = "avg log2FC",
                        name      = "t-statistics",
                        na.value  = "gray50",
                        limits    = c(-max_tstat, max_tstat), 
                        oob = scales::squish) +
  scale_size_continuous(range = c(3, 8),
                        name  = expression(-log[10](p.value))) +
  facet_wrap(~ celltype, nrow = 1) +
  theme_classic() +
  theme(strip.text      = element_text(size = 14, face = "bold"),
        axis.text.x     = element_text(size = 12),
        axis.text.y     = element_text(size = 12),
        legend.position = "right") +
  labs(x = NULL, y = NULL)
p <- p + theme(axis.text.x= element_text(angle=45, hjust=0, vjust=0))
ggsave(p, file = paste0(save_path, '/tAT2_pathway_representative_genes_t-test.png'), width = 12, height = 6)
ggsave(p, file = paste0(save_path, '/tAT2_pathway_representative_genes_t-test.pdf'), width = 8, height = 6)
ggplot2pptx(p, file = paste0(save_path, '/tAT2_pathway_representative_genes_t-test.pptx'), width = 12, height = 6)


############# Inflammatory signature only
GO_list_inflammatory = GO_list[GO_terms_order][1:13]
GO_list_inflammatory = GO_list_inflammatory[!(is.na(names(GO_list_inflammatory)))]

result_df = data.frame(GO_term = character(),
                          dpi = character(),
                          celltype = character(),
                          p.value = numeric(),
                          t.statistic = numeric(),
                          stringsAsFactors = FALSE)
for (GO_term in names(GO_list_inflammatory)) {
    so.alv <- AddModuleScore(so.alv, features = list(GO_list_inflammatory[[GO_term]]), name = GO_term)
    for (dpi in c('14dpi', '30dpi')) {
        for (celltype in c('AT2', 'tAT2')) {
            floxed.cells = Cells(so.alv)[so.alv$condition == 'floxed' & so.alv$dpi == dpi & so.alv$majority_voting == celltype]
            dAT2.cells = Cells(so.alv)[so.alv$condition == 'dAT2' & so.alv$dpi == dpi & so.alv$majority_voting == celltype]

            t.test.result <- t.test(so.alv@meta.data[dAT2.cells, paste0(GO_term, '1')], 
                                    so.alv@meta.data[floxed.cells, paste0(GO_term, '1')])
            p.value <- t.test.result$p.value
            t.statistic <- t.test.result$statistic

            result_df <- rbind(result_df, 
                          data.frame(GO_term = GO_term, 
                                     dpi = dpi, 
                                     celltype = celltype, 
                                     p.value = p.value, 
                                     t.statistic = t.statistic))
            print(paste0("Processed GO term: ", GO_term, ", dpi: ", dpi, ", celltype: ", celltype))
        }
    }
}
result_df$dpi <- factor(result_df$dpi, levels = c("14dpi", "30dpi"))
result_df$GO_term <- factor(result_df$GO_term, levels = GO_terms_order)
result_df$celltype <- factor(result_df$celltype, levels = c('AT2', 'tAT2'))

max_tstat <- max(abs(result_df$t.statistic), na.rm = TRUE)

term_order = subset(result_df, celltype == 'AT2' & dpi == '14dpi')$GO_term[order(subset(result_df, celltype == 'AT2' & dpi == '14dpi')$p.value, decreasing = TRUE)]

# 3) ggplot
p <- ggplot(result_df,
            aes(x = dpi,
                y = GO_term,
                size = -log10(p.value))) +
  geom_point(aes(color = ifelse(p.value <= 0.05, t.statistic, NA))) +
  scale_color_gradient2(
    low      = "#0625f0",   # t.statistic가 음수(작을)일 때 파랑
    mid      = "white",     # 0일 때 흰색
    high     = "#fd1e1b",   # t.statistic가 양수(클)일 때 빨강
    midpoint = 0,
    limits   = c(-max_tstat, max_tstat),
    name     = "t statistic",
    na.value = "gray50"
  ) +
  scale_size_continuous(range = c(3, 8),
                        name  = expression(-log[10](p.value))) +
  facet_wrap(~ celltype, nrow = 1) +
  scale_y_discrete(limits = term_order) +
  theme_classic() +
  theme(strip.text      = element_text(size = 14, face = "bold"),
        axis.text.x     = element_text(size = 12),
        axis.text.y     = element_text(size = 12),
        legend.position = "right") +
  labs(x = NULL, y = NULL)
p <- p + theme(axis.text.x= element_text(angle=45, hjust=0, vjust=0))
ggsave(p, file = paste0(save_path, '/inflammatory_t-test.pdf'), width = 8, height = 5)

############# AT2, tAT2
GO_list_inflammatory = GO_list[GO_terms_order][1:13]
GO_list_inflammatory = GO_list_inflammatory[!(is.na(names(GO_list_inflammatory)))]

result_df = data.frame(GO_term = character(),
                          dpi = character(),
                          celltype = character(),
                          p.value = numeric(),
                          t.statistic = numeric(),
                          stringsAsFactors = FALSE)
for (GO_term in names(GO_list_inflammatory)) {
    so.alv <- AddModuleScore(so.alv, features = list(GO_list_inflammatory[[GO_term]]), name = GO_term)
    for (dpi in c('14dpi', '30dpi')) {
        floxed.cells = Cells(so.alv)[so.alv$condition == 'floxed' & so.alv$dpi == dpi & so.alv$majority_voting %in% c('AT2', 'tAT2')]
        dAT2.cells = Cells(so.alv)[so.alv$condition == 'dAT2' & so.alv$dpi == dpi & so.alv$majority_voting %in% c('AT2', 'tAT2')]

        t.test.result <- t.test(so.alv@meta.data[dAT2.cells, paste0(GO_term, '1')], 
                                so.alv@meta.data[floxed.cells, paste0(GO_term, '1')])
        p.value <- t.test.result$p.value
        t.statistic <- t.test.result$statistic

        result_df <- rbind(result_df, 
                        data.frame(GO_term = GO_term, 
                                    dpi = dpi, 
                                    celltype = celltype, 
                                    p.value = p.value, 
                                    t.statistic = t.statistic))
        print(paste0("Processed GO term: ", GO_term, ", dpi: ", dpi, ", celltype: ", celltype))
    }
}
result_df$dpi <- factor(result_df$dpi, levels = c("14dpi", "30dpi"))
result_df$GO_term <- factor(result_df$GO_term, levels = GO_terms_order)
result_df$celltype <- factor(result_df$celltype, levels = c('AT2', 'tAT2'))

max_tstat <- max(abs(result_df$t.statistic), na.rm = TRUE)

# 3) ggplot
p <- ggplot(result_df,
            aes(x = dpi,
                y = GO_term,
                size = -log10(p.value))) +
  geom_point(aes(color = ifelse(p.value <= 0.05, t.statistic, NA))) +
  scale_color_gradient2(
    low      = "#0625f0",   # t.statistic가 음수(작을)일 때 파랑
    mid      = "white",     # 0일 때 흰색
    high     = "#fd1e1b",   # t.statistic가 양수(클)일 때 빨강
    midpoint = 0,
    limits   = c(-max_tstat, max_tstat),
    name     = "t statistic",
    na.value = "gray50"
  ) +
  scale_size_continuous(range = c(3, 8),
                        name  = expression(-log[10](p.value))) +
  scale_y_discrete(limits = term_order) +
  theme_classic() +
  theme(strip.text      = element_text(size = 14, face = "bold"),
        axis.text.x     = element_text(size = 12),
        axis.text.y     = element_text(size = 12),
        legend.position = "right") +
  labs(x = NULL, y = NULL)
p <- p + theme(axis.text.x= element_text(angle=45, hjust=0, vjust=0))
ggsave(p, file = paste0(save_path, '/inflammatory_t-test_total_AT2_tAT2.pdf'), width = 8, height = 5)


############# total alv
GO_list_inflammatory = GO_list[GO_terms_order][1:13]
GO_list_inflammatory = GO_list_inflammatory[!(is.na(names(GO_list_inflammatory)))]

result_df = data.frame(GO_term = character(),
                          dpi = character(),
                          celltype = character(),
                          p.value = numeric(),
                          t.statistic = numeric(),
                          stringsAsFactors = FALSE)
for (GO_term in names(GO_list_inflammatory)) {
    so.alv <- AddModuleScore(so.alv, features = list(GO_list_inflammatory[[GO_term]]), name = GO_term)
    for (dpi in c('14dpi', '30dpi')) {
        floxed.cells = Cells(so.alv)[so.alv$condition == 'floxed' & so.alv$dpi == dpi]
        dAT2.cells = Cells(so.alv)[so.alv$condition == 'dAT2' & so.alv$dpi == dpi]

        t.test.result <- t.test(so.alv@meta.data[dAT2.cells, paste0(GO_term, '1')], 
                                so.alv@meta.data[floxed.cells, paste0(GO_term, '1')])
        p.value <- t.test.result$p.value
        t.statistic <- t.test.result$statistic

        result_df <- rbind(result_df, 
                        data.frame(GO_term = GO_term, 
                                    dpi = dpi, 
                                    celltype = celltype, 
                                    p.value = p.value, 
                                    t.statistic = t.statistic))
        print(paste0("Processed GO term: ", GO_term, ", dpi: ", dpi, ", celltype: ", celltype))
    }
}
result_df$dpi <- factor(result_df$dpi, levels = c("14dpi", "30dpi"))
result_df$GO_term <- factor(result_df$GO_term, levels = GO_terms_order)
result_df$celltype <- factor(result_df$celltype, levels = c('AT2', 'tAT2'))

max_tstat <- max(abs(result_df$t.statistic), na.rm = TRUE)

# 3) ggplot
p <- ggplot(result_df,
            aes(x = dpi,
                y = GO_term,
                size = -log10(p.value))) +
  geom_point(aes(color = ifelse(p.value <= 0.05, t.statistic, NA))) +
  scale_color_gradient2(
    low      = "#0625f0",   # t.statistic가 음수(작을)일 때 파랑
    mid      = "white",     # 0일 때 흰색
    high     = "#fd1e1b",   # t.statistic가 양수(클)일 때 빨강
    midpoint = 0,
    limits   = c(-max_tstat, max_tstat),
    name     = "t statistic",
    na.value = "gray50"
  ) +
  scale_size_continuous(range = c(3, 8),
                        name  = expression(-log[10](p.value))) +
  scale_y_discrete(limits = term_order) +
  theme_classic() +
  theme(strip.text      = element_text(size = 14, face = "bold"),
        axis.text.x     = element_text(size = 12),
        axis.text.y     = element_text(size = 12),
        legend.position = "right") +
  labs(x = NULL, y = NULL)
p <- p + theme(axis.text.x= element_text(angle=45, hjust=0, vjust=0))
ggsave(p, file = paste0(save_path, '/inflammatory_t-test_total_alv.pdf'), width = 8, height = 5)


############# cytokine
#fd1e1b
#0625f0