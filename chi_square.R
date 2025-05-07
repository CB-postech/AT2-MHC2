### annotate T cell
# conda activate project_lung_exercise_R

library(Seurat)
library(magrittr)
library(data.table)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(ggsignif)
library(gprofiler2)
library(future)
library(pbmcapply)

set.seed(42)

# load custom functions
source('/home/sjcho/yard/functions/R/seurat_count_to_normalization.R')
source('/home/sjcho/yard/functions/R/FeaturePlot_sjcho.R')
source('/home/sjcho/yard/functions/R/save_ggplot2_to_ppt.R')
source('/home/sjcho/yard/functions/R/draw_proportion.R')
source('/home/sjcho/projects/AT2_MHC2/20241113_after_cellbender/utils.R')
source('/home/sjcho/yard/functions/R/draw_celltype_wise_dotplot.R')
source('/home/sjcho/yard/functions/R/draw_proportion.R')
source('/home/sjcho/projects/AT2_MHC2/20241113_after_cellbender/manuscript_figures/chi-square_basal_exp/chi-sqaure_with_sampling.R')

### ccr geneset
library(KEGGREST)
library(org.Mm.eg.db)
library(EnhancedVolcano)

pathway_genes <- keggLink("mmu", "pathway:mmu04060")
entrez_ids = sapply(strsplit(pathway_genes, ":"), function(x) x[2])
kegg_ccr <- mapIds(org.Mm.eg.db, keys = entrez_ids %>% as.vector, column = "SYMBOL", keytype = "ENTREZID") # ccr means cytokine, chemokine, and their receptor
kegg_ccr <- setdiff(kegg_ccr, 'Cd4')

cytokines <- c(grep('Ifn', kegg_ccr, value = T), grep('Il', kegg_ccr, value = T), grep('Tgf', kegg_ccr, value = T), grep('Bmp', kegg_ccr, value = T), grep('Gdf', kegg_ccr, value = T))
cytokines <- c(cytokines, 'Cntf', 'Csf1', 'Csf2', 'Csf3', 'Ctf1', 'Ctf2', 'Epo', 'Inhba', 'Inhbb', 'Inhbc', 'Inhbe', 'Lep')
cytokines <- c(cytokines, 'Mstn', 'Ngf', 'Nodal', 'Osm', 'Pf4', 'Ppbp', 'Prl', 'Pr15a1', 'Pr16a1', 'Thpo', 'Tnf')
cytokines <- c(cytokines, grep('Tnfsf', kegg_ccr, value = T), 'Tslp', 'Xcl1')
# not cytokines : Eda, Klk1b4, Rell1, Rell12, Fasl, Cd40lg, Inha, Relt

# remove receptor
# ra, rb, rb1, r1b, r1a, ra1, ra2, rc, re, rap, rnm, rl1, rg, st, r1, r2
cytokines <- cytokines[!grepl('r$', cytokines)]
cytokines <- cytokines[!grepl('ra$', cytokines)]
cytokines <- cytokines[!grepl('rc$', cytokines)]
cytokines <- cytokines[!grepl('re$', cytokines)]
cytokines <- cytokines[!grepl('rg$', cytokines)]
cytokines <- cytokines[!grepl('rap$', cytokines)]
cytokines <- cytokines[!grepl('rb$', cytokines)]
cytokines <- cytokines[!grepl('rl1$', cytokines)]
cytokines <- cytokines[!grepl('r1$', cytokines)]
cytokines <- cytokines[!grepl('r2$', cytokines)]
cytokines <- cytokines[!grepl('rb1$', cytokines)]
cytokines <- cytokines[!grepl('rb2$', cytokines)]
cytokines <- cytokines[!grepl('r1b$', cytokines)]
cytokines <- cytokines[!grepl('r1a$', cytokines)]
cytokines <- cytokines[!grepl('ra1$', cytokines)]
cytokines <- cytokines[!grepl('ra2$', cytokines)]
cytokines <- cytokines[!grepl('rl2$', cytokines)]
cytokines <- cytokines[!grepl('Il6st', cytokines)]

save_path = '/home/sjcho/projects/AT2_MHC2/20241113_after_cellbender/manuscript_figures/chi-square_basal_exp/outs_20250318/'
so.cd4 <- readRDS('/home/sjcho/projects/AT2_MHC2/20241113_after_cellbender/4.Tcell_study/4.1.Tcell_subset/outs/4.5.full_CD4_transcriptomic_difference/woGD_cCd4.rds')

so.cd4$condition_dpi = paste0(so.cd4$dpi, '_', so.cd4$condition)
so.cd4$condition_dpi = factor(so.cd4$condition_dpi, levels = c('naive_floxed', '7dpi_floxed', '14dpi_floxed', '30dpi_floxed', 
                                                                'naive_dAT2', '7dpi_dAT2', '14dpi_dAT2', '30dpi_dAT2'))
so.cd4$annotation[so.cd4$annotation == 'Cxcr3+ Treg'] = 'Treg'
so.cd4$annotation[so.cd4$annotation == 'naive-like Treg'] = 'Treg'
so.cd4$annotation[so.cd4$annotation == 'Proliferating Treg'] = 'Treg.prolif.'
so.cd4$annotation[so.cd4$annotation == 'cycling Cd4 T cell'] = 'Th1 Cd4 T.prolif.'

Tcell_cols_order = c('naive-like Cd4 T cell', 'Treg', 'Treg.prolif.', 'Th1 Cd4 T.prolif.', 'Th1 Cd4 T cell', 'Th17 Cd4 T cell', 'Cd4 Tcm')

n_iter <- 100 # number of sampling
num_cores <- 32
options('future.globals.maxSize' = 400 * 1014*1024^2) # 350GB

levels(so.cd4@meta.data[, 'condition']) = c('floxed', 'dAT2')

save_path = '/home/sjcho/projects/AT2_MHC2/20241113_after_cellbender/manuscript_figures/chi-square_basal_exp/chi-square-basal-level-outs-allDPI_Celltype/'

for (celltype in c('Th1 Cd4 T cell', 'Th17 Cd4 T cell', 'Cd4 Tcm')) {
  for (dpi in c('naive', '7dpi', '14dpi', '30dpi')) {
    if (celltype == 'Th1 Cd4 T cell' & dpi == 'naive') next
    if (celltype == 'Th1 Cd4 T cell' & dpi == '7dpi') next
    if (celltype == 'Th1 Cd4 T cell' & dpi == '14dpi') next
    print(paste0(celltype, ' ', dpi, ' Start'))
    so.use = subset(so.cd4, cells = Cells(so.cd4)[so.cd4$annotation == celltype & so.cd4$dpi == dpi])
    levels(so.use@meta.data[, 'condition']) = c('floxed', 'dAT2')
    # result <- result_list[[name]]
    result <- run_ChiSqaure_with_basal_exp_correction_mean(so.use, 'condition', n_iter = 100, target_umi_depth = 10000, target_umi_basal = 500, num_cores = num_cores, use_batch = FALSE, batch_size = 30)
    save(result, file = paste0(save_path, 'chi_square_result_mean_', celltype, '_', dpi, '.RData'))
    ### visualization
    chi_summary <- result[['chi_summary']]
    chi_summary_visual <- chi_summary
    chi_summary_visual$p_BH = ifelse(chi_summary_visual$p_BH < 10^-10, 10^-10, chi_summary_visual$p_BH)
    chi_summary_visual$log2FC = ifelse(chi_summary_visual$log2FC > 5, 5, ifelse(chi_summary_visual$log2FC < -5, -5, chi_summary_visual$log2FC))

    chi_summary_visual <- chi_summary_visual %>%
        mutate(is_cytokines = rownames(chi_summary_visual) %in% cytokines)

    # data filter; delete NA and sort by differential.std.residual
    chi_summary_visual <- chi_summary_visual %>%
        filter(!is.na(statistic) & !is.na(p_BH))

    highlighted_genes <- c('Il10')

    chi_summary_visual <- chi_summary_visual[chi_summary_visual$floxed_observed > 1 & chi_summary_visual$dAT2_observed > 1, ]
    chi_summary_visual$log10_oberseved = log10(chi_summary_visual$floxed_observed + chi_summary_visual$dAT2_observed + 1)
    # highlighted_genes <- kegg_ccr
    p <- ggplot(chi_summary_visual, aes(x = log2FC, y = -log10(p_BH))) +
        geom_point(data = subset(chi_summary_visual, !is_cytokines), 
                    color = "gray80", size = 2.5) +
        # draw kegg_ccr lately to make it on top
        geom_point(data = subset(chi_summary_visual, is_cytokines), 
                    color = "black", size = 5.5) +
        geom_point(data = subset(chi_summary_visual, rownames(chi_summary_visual) %in% highlighted_genes), 
            color = "red", size = 5.5) +
        # geom_text_repel(data = subset(chi_summary_visual, rownames(chi_summary_visual) %in% highlighted_genes),
        #                 aes(label = rownames(chi_summary_visual)[rownames(chi_summary_visual) %in% highlighted_genes]),
        #                 size = 6,
        #                 box.padding = 0.5,
        #                 segment.color = "grey50") +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    #        geom_vline(xintercept = 0, linetype = "dashed") +
        theme_classic() +
        theme(legend.position = "none") +
        theme(
        legend.position = "none",
            axis.title = element_text(size = 20),
            axis.text = element_text(size = 16),
            plot.title = element_text(size = 24, face = "bold")
        ) +
        # xlim(-2.5, 2.5) + 
        # ylim(0, 15) +
        labs(x = "log2FC", y = "-log10(adjusted p-value)")
    # ggsave(p, filename = paste0(save_path, 'test_14dpi_whole_Treg_chi_5000sample_300iter_differential.avg.png'), width = 8, height = 6)
    ggsave(p, filename = paste0(save_path, 'fig5.E.chi_mean.',  celltype, '_', dpi, '.png'), width = 8, height = 6)
    ggsave(p, filename = paste0(save_path, 'fig5.E.chi_mean.',  celltype, '_', dpi, '.pdf'), width = 8, height = 6)
    ggplot2pptx(p, 8, 6, paste0(save_path, 'fig5.E.chi_mean.',  celltype, '_', dpi, '.pptx'))
  }
}

############################ 14dpi Treg cytokines

### anti-inflammatory dotplot visualization

load('/home/sjcho/projects/AT2_MHC2/20241113_after_cellbender/manuscript_figures/chi-square_basal_exp/chi-square-basal-level-outs-allDPI_Celltype/chi_square_result_mean_Treg_14dpi.RData')

anti.inflammatory.cytokines = c('Il10', 'Il13', 'Il5', 'Il22', 'Il4', 'Tgfb1', 'Areg')
anti_inflammatory <- result$chi_summary[anti.inflammatory.cytokines, ]
anti_inflammatory$OR <- ((anti_inflammatory$floxed_expected + 1) / (anti_inflammatory$dAT2_expected + 1)) /
                                ((anti_inflammatory$floxed_observed + 1) / (anti_inflammatory$dAT2_observed + 1))
anti_inflammatory$log2OR <- log2(anti_inflammatory$OR)

######## further visualization requries

library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)

p_long <- p_BH_df %>% 
  tibble::rownames_to_column(var = "gene") %>% 
  pivot_longer(-gene, names_to = "celltype", values_to = "p_value")

log2OR_long <- log2OR_df %>% 
  tibble::rownames_to_column(var = "gene") %>% 
  pivot_longer(-gene, names_to = "celltype", values_to = "log2OR")

df_long <- left_join(p_long, log2OR_long, by = c("gene", "celltype"))

p <- ggplot(df_long, aes(x = celltype, y = gene)) +
  geom_point(aes(size = -log10(p_value), color = log2OR)) +
  scale_color_gradient2(
    low = "#0000ff", 
    mid = "white",   
    high = "#ff001e", 
    midpoint = 0,       
    na.value = "gray25"
  ) +
  labs(x = "Cell Type", 
       y = "Gene", 
       size = "-log10(p-value)", 
       color = "log2OR") +
  theme_minimal() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(p, file = paste0(save_path, 'fig4.g.14dpi.Treg.antiinfalmmatory.cytokines.png'))
