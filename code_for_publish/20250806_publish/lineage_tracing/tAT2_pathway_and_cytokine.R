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
save_path = '/home/sjcho/projects/AT2_MHC2/code_for_publish/20250806_publish/lineage_tracing/outs/tAT2_pathway'

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


#################### representative genes
Hypoxia_glycolysis = c('Hif1a', 'Pgam1', 'Eno1', 'Myc', 'Pkm', 'Slc2a1', 'Cdkn1a', 'Slc16a3', 'Ndrg1')
senescence_p53 = c('Trp53', 'Mdm2', 'Ccnd1', 'Gdf15')

so.AT2.tAT2.14.30 = subset(so, majority_voting %in% c('AT2', 'tAT2') & dpi %in% c( '30dpi'))
so.AT2.tAT2.14.30$anno_condition = paste0(so.AT2.tAT2.14.30$majority_voting, '_', so.AT2.tAT2.14.30$condition)
so.AT2.tAT2.14.30$anno_condition = factor(so.AT2.tAT2.14.30$anno_condition, levels = c(
    'AT2_floxed', 'AT2_dAT2', 'tAT2_floxed', 'tAT2_dAT2'
))

source('/home/sjcho/yard/functions/R/draw_heatmap_gene_expression_by_cluster.R')
p <- plot_gene_heatmap(so.AT2.tAT2.14.30, features = c(Hypoxia_glycolysis, senescence_p53), 
                        group.by = 'anno_condition', 
                        cols = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100), scale = 'column', 
                        save_path = save_path, feature_name = 'tAT2_pathway_representative_genes', 
                        gaps_row = c(2), gaps_col = NULL, 
                        show_colnames = TRUE, cluster_cols = FALSE,
                        width = 3, height = 1.5)
                        
### t-test
result_df = data.frame(gene = character(),
                          dpi = character(),
                          celltype = character(),
                          p.value = numeric(),
                          t.statistic = numeric(),
                          stringsAsFactors = FALSE)
for (gene in c(Hypoxia_glycolysis, senescence_p53)) {   
    for (dpi in c('30dpi')) {
        for (celltype in c('tAT2', 'AT2')) {
            floxed.cells = Cells(so)[so$condition == 'floxed' & so$dpi == dpi & so$majority_voting == celltype]
            dAT2.cells = Cells(so)[so$condition == 'dAT2' & so$dpi == dpi & so$majority_voting == celltype]

            t.test.result <- t.test(expm1(so@assays$RNA$data[gene, dAT2.cells]), 
                                    expm1(so@assays$RNA$data[gene, floxed.cells]))
            p.value <- t.test.result$p.value
            t.statistic <- t.test.result$statistic

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
result_df$celltype <- factor(result_df$celltype, levels = c('AT2', 'tAT2'))

# 3) ggplot
max_tstat <- max(abs(result_df$t.statistic), na.rm = TRUE)

# 3) ggplot
p <- ggplot(result_df,
            aes(x = dpi,
                y = gene,
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
ggsave(p, file = paste0(save_path, '/tAT2_pathway_representative_genes_t-test.png'), width = 4, height = 6)
ggsave(p, file = paste0(save_path, '/tAT2_pathway_representative_genes_t-test.pdf'), width = 8, height = 6)
ggplot2pptx(p, file = paste0(save_path, '/tAT2_pathway_representative_genes_t-test.pptx'), width = 12, height = 6)

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
    so <- AddModuleScore(so, features = list(GO_list_inflammatory[[GO_term]]), name = GO_term)
    for (dpi in c('30dpi')) {
        for (celltype in c('AT2', 'tAT2')) {
            floxed.cells = Cells(so)[so$condition == 'floxed' & so$dpi == dpi & so$majority_voting == celltype]
            dAT2.cells = Cells(so)[so$condition == 'dAT2' & so$dpi == dpi & so$majority_voting == celltype]

            t.test.result <- t.test(so@meta.data[dAT2.cells, paste0(GO_term, '1')], 
                                    so@meta.data[floxed.cells, paste0(GO_term, '1')])
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
result_df$dpi <- factor(result_df$dpi, levels = c("30dpi"))
result_df$GO_term <- factor(result_df$GO_term, levels = GO_terms_order)
result_df$celltype <- factor(result_df$celltype, levels = c('AT2', 'tAT2'))

max_tstat <- max(abs(result_df$t.statistic), na.rm = TRUE)

term_order = subset(result_df, celltype == 'AT2' & dpi == '30dpi')$GO_term[order(subset(result_df, celltype == 'AT2' & dpi == '30dpi')$p.value, decreasing = TRUE)]

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
  scale_y_discrete(limits = term_order) +
  theme_classic() +
  theme(strip.text      = element_text(size = 14, face = "bold"),
        axis.text.x     = element_text(size = 12),
        axis.text.y     = element_text(size = 12),
        legend.position = "right") +
  labs(x = NULL, y = NULL)
p <- p + theme(axis.text.x= element_text(angle=45, hjust=0, vjust=0))
ggsave(p, file = paste0(save_path, '/inflammatory_t-test.pdf'), width = 8, height = 5)

############## cytokine
################# cytokine response
############### IL4, IL6, IL13, IFNa, IFNb, IFNg
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
########## GOBP: IL-11 mediated signaling pathway gene set
GO_list[['Il-11_mediated_signaling_pathway']] = c('Il11ra3', 'Il11', 'Il11ra1', 'Il11ra2', 'Il6st', 'Jak1', 'Stat3')

for (GO_term in names(GO_list)) {
    so <- AddModuleScore(so, features = list(GO_list[[GO_term]]), name = GO_term)
}

### AT2
so.Intermediate <- subset(so, majority_voting %in% c('AT2') & condition %in% c('floxed', 'dAT2'))
so.Intermediate$condition <- factor(so.Intermediate$condition, levels = c('floxed', 'dAT2'))

p1 <- VlnPlot(so.Intermediate, features = 'GOBP_RESPONSE_TO_INTERLEUKIN_11', group.by = 'condition', pt.size = 0, cols = c('darkgray', 'darkred', 'darkgray', 'darkred')) + geom_boxplot(width = 0.2) + NoLegend()
p1 <- p1 + geom_signif(comparisons = list(c('floxed', 'dAT2')), map_signif_level = TRUE)
p1 <- p1 + ylim(so.Intermediate$GOBP_RESPONSE_TO_INTERLEUKIN_11 %>% min, (so.Intermediate$GOBP_RESPONSE_TO_INTERLEUKIN_11 %>% max * 1.25))

p2 <- VlnPlot(so.Intermediate, features = 'HM_Tnfa_NFKB1', group.by = 'condition', pt.size = 0, cols = c('darkgray', 'darkred', 'darkgray', 'darkred')) + geom_boxplot(width = 0.2) + NoLegend()
p2 <- p2 + geom_signif(comparisons = list(c('floxed', 'dAT2')), map_signif_level = TRUE)
p2 <- p2 + ylim(so.Intermediate$HM_Tnfa_NFKB1 %>% min, (so.Intermediate$HM_Tnfa_NFKB1 %>% max * 1.25))

p3 <- VlnPlot(so.Intermediate, features = 'Il-11_mediated_signaling_pathway1', group.by = 'condition', pt.size = 0, cols = c('darkgray', 'darkred', 'darkgray', 'darkred')) + geom_boxplot(width = 0.2) + NoLegend()
p3 <- p3 + geom_signif(comparisons = list(c('floxed', 'dAT2')), map_signif_level = TRUE)
p3 <- p3 + ylim(so.Intermediate$'Il-11_mediated_signaling_pathway1' %>% min, (so.Intermediate$'Il-11_mediated_signaling_pathway1' %>% max * 1.25))

p <- cowplot::plot_grid(p1, p2, p3, ncol = 1)
ggsave(p, file = paste0(save_path, '/TNFa_IL1_IL11_AT2.pdf'), width = 4, height = 15)

### tAT2
so.Intermediate <- subset(so, majority_voting %in% c('tAT2') & condition %in% c('floxed', 'dAT2'))
so.Intermediate$condition <- factor(so.Intermediate$condition, levels = c('floxed', 'dAT2'))

p1 <- VlnPlot(so.Intermediate, features = 'GOBP_RESPONSE_TO_INTERLEUKIN_11', group.by = 'condition', pt.size = 0, cols = c('darkgray', 'darkred', 'darkgray', 'darkred')) + geom_boxplot(width = 0.2) + NoLegend()
p1 <- p1 + geom_signif(comparisons = list(c('floxed', 'dAT2')), map_signif_level = TRUE)
p1 <- p1 + ylim(so.Intermediate$GOBP_RESPONSE_TO_INTERLEUKIN_11 %>% min, (so.Intermediate$GOBP_RESPONSE_TO_INTERLEUKIN_11 %>% max * 1.25))

p2 <- VlnPlot(so.Intermediate, features = 'HM_Tnfa_NFKB1', group.by = 'condition', pt.size = 0, cols = c('darkgray', 'darkred', 'darkgray', 'darkred')) + geom_boxplot(width = 0.2) + NoLegend()
p2 <- p2 + geom_signif(comparisons = list(c('floxed', 'dAT2')), map_signif_level = TRUE)
p2 <- p2 + ylim(so.Intermediate$HM_Tnfa_NFKB1 %>% min, (so.Intermediate$HM_Tnfa_NFKB1 %>% max * 1.25))

p3 <- VlnPlot(so.Intermediate, features = 'Il-11_mediated_signaling_pathway1', group.by = 'condition', pt.size = 0, cols = c('darkgray', 'darkred', 'darkgray', 'darkred')) + geom_boxplot(width = 0.2) + NoLegend()
p3 <- p3 + geom_signif(comparisons = list(c('floxed', 'dAT2')), map_signif_level = TRUE)
p3 <- p3 + ylim(so.Intermediate$'Il-11_mediated_signaling_pathway1' %>% min, (so.Intermediate$'Il-11_mediated_signaling_pathway1' %>% max * 1.25))

p <- cowplot::plot_grid(p1, p2, p3, ncol = 1)
ggsave(p, file = paste0(save_path, '/TNFa_IL1_IL11_tAT2.pdf'), width = 4, height = 15)

### AT1
so.Intermediate <- subset(so, majority_voting %in% c('AT1') & condition %in% c('floxed', 'dAT2'))
so.Intermediate$condition <- factor(so.Intermediate$condition, levels = c('floxed', 'dAT2'))

p1 <- VlnPlot(so.Intermediate, features = 'GOBP_RESPONSE_TO_INTERLEUKIN_11', group.by = 'condition', pt.size = 0, cols = c('darkgray', 'darkred', 'darkgray', 'darkred')) + geom_boxplot(width = 0.2) + NoLegend()
p1 <- p1 + geom_signif(comparisons = list(c('floxed', 'dAT2')), map_signif_level = TRUE)
p1 <- p1 + ylim(so.Intermediate$GOBP_RESPONSE_TO_INTERLEUKIN_11 %>% min, (so.Intermediate$GOBP_RESPONSE_TO_INTERLEUKIN_11 %>% max * 1.25))

p2 <- VlnPlot(so.Intermediate, features = 'HM_Tnfa_NFKB1', group.by = 'condition', pt.size = 0, cols = c('darkgray', 'darkred', 'darkgray', 'darkred')) + geom_boxplot(width = 0.2) + NoLegend()
p2 <- p2 + geom_signif(comparisons = list(c('floxed', 'dAT2')), map_signif_level = TRUE)
p2 <- p2 + ylim(so.Intermediate$HM_Tnfa_NFKB1 %>% min, (so.Intermediate$HM_Tnfa_NFKB1 %>% max * 1.25))

p3 <- VlnPlot(so.Intermediate, features = 'Il-11_mediated_signaling_pathway1', group.by = 'condition', pt.size = 0, cols = c('darkgray', 'darkred', 'darkgray', 'darkred')) + geom_boxplot(width = 0.2) + NoLegend()
p3 <- p3 + geom_signif(comparisons = list(c('floxed', 'dAT2')), map_signif_level = TRUE)
p3 <- p3 + ylim(so.Intermediate$'Il-11_mediated_signaling_pathway1' %>% min, (so.Intermediate$'Il-11_mediated_signaling_pathway1' %>% max * 1.25))

p <- cowplot::plot_grid(p1, p2, p3, ncol = 1)
ggsave(p, file = paste0(save_path, '/TNFa_IL1_IL11_AT1.pdf'), width = 4, height = 15)