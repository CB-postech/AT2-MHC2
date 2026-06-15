### conda activate project_lung_exercise_R
### https://satijalab.org/seurat/articles/cell_cycle_vignette.html

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
source('/home/sjcho/yard/functions/R/draw_proportion.R')
source('/home/sjcho/projects/AT2_MHC2/20241113_after_cellbender/utils.R')
source('/home/sjcho/yard/functions/R/draw_celltype_wise_dotplot.R')
source('/home/sjcho/yard/functions/R/draw_volcano.R')
source('/home/sjcho/projects/AT2_MHC2/20241113_after_cellbender/5.cellrank2/utils.cellrnak2.R')
source('/home/sjcho/yard/functions/R/save_ggplot2_to_ppt.R')

save_path = '/home/sjcho/projects/AT2_MHC2/code_for_publish/251208_revision/cellcycle/outs'
so.alv <- readRDS('/home/sjcho/projects/AT2_MHC2/20241113_after_cellbender/3.whole_cell_annotation/outs/3.1.epithelial_annotation/alv_annotation.rds')

############# cell cycle genes from reactome #############
library(msigdbr)
C2_genes = msigdbr(species = "mouse", category = "C2")
s.genes = C2_genes %>% filter(gs_name == "REACTOME_S_PHASE") %>% pull(gene_symbol) %>% unique()
g2.genes = C2_genes %>% filter(gs_name == "REACTOME_G2_PHASE") %>% pull(gene_symbol) %>% unique()
m.genes = C2_genes %>% filter(gs_name == "REACTOME_M_PHASE") %>% pull(gene_symbol) %>% unique()

so.alv <- CellCycleScoring(so.alv, s.features = s.genes, g2m.features = c(g2.genes, m.genes), set.ident = TRUE)
so.alv$CC.prediction <- Idents(so.alv)

############ seurat CellCycleScoring
library(homologene)
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

s.genes.mm <- human2mouse(cc.genes.updated.2019$s.genes)$mouseGene
g2m.genes.mm <- human2mouse(cc.genes.updated.2019$g2m.genes)$mouseGene

celltypist_result = fread('/home/sjcho/projects/AT2_MHC2/20241113_after_cellbender/manuscript_figures/figure4_and_exd6_20250225/alv_full_celltypist_result_overcluster1.csv')
celltypist_result = as.data.frame(celltypist_result)
rownames(celltypist_result) <- celltypist_result$V1
celltypist_result <- celltypist_result[, -1]
so.alv$majority_voting <- celltypist_result$majority_voting[match(colnames(so.alv), rownames(celltypist_result))]

library(msigdbr)
h_gene_sets = msigdbr(species = "mouse", category = "H")
so.alv <- AddModuleScore(so.alv, features = list(subset(h_gene_sets, gs_name == 'HALLMARK_INFLAMMATORY_RESPONSE')$gene_symbol), name = "HM_inflammatory_response")
so.alv.14.AT2 <- subset(so.alv, cells = Cells(so.alv)[so.alv$majority_voting == "AT2" & so.alv$dpi == "14dpi"])
so.alv.30.AT2 <- subset(so.alv, cells = Cells(so.alv)[so.alv$majority_voting == "AT2" & so.alv$dpi == "30dpi"])

threshold_otsu <- function(x, nbins = 100) {
  x_valid <- x[!is.na(x)]
    
  h <- hist(x_valid, breaks = nbins, plot = FALSE)
  hist_counts <- as.numeric(h$counts)
  bin_centers <- h$mids
  
  weight1 <- cumsum(hist_counts)
  weight2 <- rev(cumsum(rev(hist_counts)))
  
  weighted_centers <- hist_counts * bin_centers
  mean1 <- cumsum(weighted_centers) / weight1
  mean2 <- rev(cumsum(rev(weighted_centers)) / cumsum(rev(hist_counts)))
  
  len_bins <- length(bin_centers)
  w1_sliced <- weight1[1:(len_bins - 1)]
  w2_sliced <- weight2[2:len_bins]
  m1_sliced <- mean1[1:(len_bins - 1)]
  m2_sliced <- mean2[2:len_bins]
  
  variance12 <- w1_sliced * w2_sliced * (m1_sliced - m2_sliced)^2
  idx <- which.max(variance12)
  threshold <- bin_centers[1:(len_bins - 1)][idx]
  
  return(threshold)
}

########## 14dpi AT2 high inflammatory response
so.alv.14.AT2.high_inflmmatory <- subset(so.alv.14.AT2, HM_inflammatory_response1 > threshold_otsu(so.alv.14.AT2$HM_inflammatory_response1))
so.alv.14.AT2.high_inflmmatory <- CellCycleScoring(so.alv.14.AT2.high_inflmmatory, s.features = s.genes.mm, g2m.features = g2m.genes.mm)

md <- data.frame(
  condition_dpi = unlist(so.alv.14.AT2.high_inflmmatory$condition_dpi),
  Phase = unlist(so.alv.14.AT2.high_inflmmatory$Phase)
)
df <- md %>%
  filter(condition_dpi %in% c("14dpi_floxed", "14dpi_dAT2")) %>%
  group_by(condition_dpi, Phase) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(condition_dpi) %>%
  mutate(prop = n / sum(n))

p <- ggplot(df, aes(x = condition_dpi, y = prop, fill = Phase)) +
  geom_col(position = "fill", width = 0.6) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("G1" = "#E76F51", "S" = "#2A9D8F", "G2M" = "#264653")) +
  labs(x = NULL, y = "Proportion", fill = "Phase") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(paste0(save_path, "/cellcycle_proportion_14dpi_AT2_inflammatory_high.pdf"), width = 3, height = 5)

########## 14dpi AT2 low inflammatory response
so.alv.14.AT2.low_inflammatory <- subset(so.alv.14.AT2, HM_inflammatory_response1 <= threshold_otsu(so.alv.14.AT2$HM_inflammatory_response1))
so.alv.14.AT2.low_inflammatory <- CellCycleScoring(so.alv.14.AT2.low_inflammatory, s.features = s.genes.mm, g2m.features = g2m.genes.mm)

md <- data.frame(
  condition_dpi = unlist(so.alv.14.AT2.low_inflammatory$condition_dpi),
  Phase = unlist(so.alv.14.AT2.low_inflammatory$Phase)
)
df <- md %>%
  filter(condition_dpi %in% c("14dpi_floxed", "14dpi_dAT2")) %>%
  group_by(condition_dpi, Phase) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(condition_dpi) %>%
  mutate(prop = n / sum(n))

p <- ggplot(df, aes(x = condition_dpi, y = prop, fill = Phase)) +
  geom_col(position = "fill", width = 0.6) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("G1" = "#E76F51", "S" = "#2A9D8F", "G2M" = "#264653")) +
  labs(x = NULL, y = "Proportion", fill = "Phase") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(paste0(save_path, "/cellcycle_proportion_14dpi_AT2_inflammatory_low.pdf"), width = 3, height = 5)

########## 30dpi AT2 high inflammatory response
so.alv.30.AT2.high_inflammatory <- subset(so.alv.30.AT2, HM_inflammatory_response1 > threshold_otsu(so.alv.30.AT2$HM_inflammatory_response1))
so.alv.30.AT2.high_inflammatory <- CellCycleScoring(so.alv.30.AT2.high_inflammatory, s.features = s.genes.mm, g2m.features = g2m.genes.mm)

md <- data.frame(
  condition_dpi = unlist(so.alv.30.AT2.high_inflammatory$condition_dpi),
  Phase = unlist(so.alv.30.AT2.high_inflammatory$Phase)
)
df <- md %>%
  filter(condition_dpi %in% c("30dpi_floxed", "30dpi_dAT2")) %>%
  group_by(condition_dpi, Phase) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(condition_dpi) %>%
  mutate(prop = n / sum(n))

p <- ggplot(df, aes(x = condition_dpi, y = prop, fill = Phase)) +
  geom_col(position = "fill", width = 0.6) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("G1" = "#E76F51", "S" = "#2A9D8F", "G2M" = "#264653")) +
  labs(x = NULL, y = "Proportion", fill = "Phase") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(paste0(save_path, "/cellcycle_proportion_30dpi_AT2_inflammatory_high.pdf"), width = 3, height = 5)

########## 30dpi AT2 low inflammatory response
so.alv.30.AT2.low_inflammatory <- subset(so.alv.30.AT2, HM_inflammatory_response1 <= threshold_otsu(so.alv.30.AT2$HM_inflammatory_response1))
so.alv.30.AT2.low_inflammatory <- CellCycleScoring(so.alv.30.AT2.low_inflammatory, s.features = s.genes.mm, g2m.features = g2m.genes.mm)

md <- data.frame(
  condition_dpi = unlist(so.alv.30.AT2.low_inflammatory$condition_dpi),
  Phase = unlist(so.alv.30.AT2.low_inflammatory$Phase)
)
df <- md %>%
  filter(condition_dpi %in% c("30dpi_floxed", "30dpi_dAT2")) %>%
  group_by(condition_dpi, Phase) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(condition_dpi) %>%
  mutate(prop = n / sum(n))

p <- ggplot(df, aes(x = condition_dpi, y = prop, fill = Phase)) +
  geom_col(position = "fill", width = 0.6) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("G1" = "#E76F51", "S" = "#2A9D8F", "G2M" = "#264653")) +
  labs(x = NULL, y = "Proportion", fill = "Phase") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(paste0(save_path, "/cellcycle_proportion_30dpi_AT2_inflammatory_low.pdf"), width = 3, height = 5)

AT2_14dpi_low_inflammatory_cells <- Cells(so.alv.14.AT2.low_inflammatory)
AT2_14dpi_high_inflammatory_cells <- Cells(so.alv.14.AT2.high_inflmmatory)
AT2_30dpi_low_inflammatory_cells <- Cells(so.alv.30.AT2.low_inflammatory)
AT2_30dpi_high_inflammatory_cells <- Cells(so.alv.30.AT2.high_inflammatory)

so.alv.AT2.14.30.floxed = subset(so.alv, cells = Cells(so.alv)[so.alv$condition == 'floxed' & so.alv$dpi %in% c('14dpi', '30dpi')])
so.alv.AT2.14.30.dAT2 = subset(so.alv, cells = Cells(so.alv)[so.alv$condition == 'dAT2' & so.alv$dpi %in% c('14dpi', '30dpi')])
p <- DimPlot(so.alv.AT2.14.30.floxed, cells.highlight = intersect(Cells(so.alv.AT2.14.30.floxed), AT2_14dpi_low_inflammatory_cells), cols.highlight = "blue", cols = "grey80", sizes.highlight = 1) + NoLegend()
p <- p + set_UMAP
ggsave(paste0(save_path, "/14dpi_AT2_low_inflammatory_cells_floxed.png"), width = 4, height = 4)
p <- DimPlot(so.alv.AT2.14.30.floxed, cells.highlight = intersect(Cells(so.alv.AT2.14.30.floxed), AT2_14dpi_high_inflammatory_cells), cols.highlight = "red", cols = "grey80", sizes.highlight = 1) + NoLegend()
p <- p + set_UMAP
ggsave(paste0(save_path, "/14dpi_AT2_high_inflammatory_cells_floxed.png"), width = 4, height = 4)

p <- DimPlot(so.alv.AT2.14.30.dAT2, cells.highlight = intersect(Cells(so.alv.AT2.14.30.dAT2), AT2_14dpi_low_inflammatory_cells), cols.highlight = "blue", cols = "grey80", sizes.highlight = 1) + NoLegend()
p <- p + set_UMAP
ggsave(paste0(save_path, "/14dpi_AT2_low_inflammatory_cells_dAT2.png"), width = 4, height = 4)
p <- DimPlot(so.alv.AT2.14.30.dAT2, cells.highlight = intersect(Cells(so.alv.AT2.14.30.dAT2), AT2_14dpi_high_inflammatory_cells), cols.highlight = "red", cols = "grey80", sizes.highlight = 1) + NoLegend()
p <- p + set_UMAP
ggsave(paste0(save_path, "/14dpi_AT2_high_inflammatory_cells_dAT2.png"), width = 4, height = 4)

p <- DimPlot(so.alv.AT2.14.30.floxed, cells.highlight = intersect(Cells(so.alv.AT2.14.30.floxed), AT2_30dpi_low_inflammatory_cells), cols.highlight = "blue", cols = "grey80", sizes.highlight = 1) + NoLegend()
p <- p + set_UMAP
ggsave(paste0(save_path, "/30dpi_AT2_low_inflammatory_cells_floxed.png"), width = 4, height = 4)
p <- DimPlot(so.alv.AT2.14.30.floxed, cells.highlight = intersect(Cells(so.alv.AT2.14.30.floxed), AT2_30dpi_high_inflammatory_cells), cols.highlight = "red", cols = "grey80", sizes.highlight = 1) + NoLegend()
p <- p + set_UMAP
ggsave(paste0(save_path, "/30dpi_AT2_high_inflammatory_cells_floxed.png"), width = 4, height = 4)

p <- DimPlot(so.alv.AT2.14.30.dAT2, cells.highlight = intersect(Cells(so.alv.AT2.14.30.dAT2), AT2_30dpi_low_inflammatory_cells), cols.highlight = "blue", cols = "grey80", sizes.highlight = 1) + NoLegend()
p <- p + set_UMAP
ggsave(paste0(save_path, "/30dpi_AT2_low_inflammatory_cells_dAT2.png"), width = 4, height = 4)
p <- DimPlot(so.alv.AT2.14.30.dAT2, cells.highlight = intersect(Cells(so.alv.AT2.14.30.dAT2), AT2_30dpi_high_inflammatory_cells), cols.highlight = "red", cols = "grey80", sizes.highlight = 1) + NoLegend()
p <- p + set_UMAP
ggsave(paste0(save_path, "/30dpi_AT2_high_inflammatory_cells_dAT2.png"), width = 4, height = 4)

############ merge all

make_md <- function(obj, dpi, inflam_level) {
  data.frame(
    condition_dpi = unlist(obj$condition_dpi),
    Phase = unlist(obj$Phase),
    dpi = dpi,
    inflammatory = inflam_level
  )
}

md_all <- bind_rows(
  make_md(so.alv.14.AT2.high_inflmmatory, "14dpi", "High"),
  make_md(so.alv.14.AT2.low_inflammatory, "14dpi", "Low"),
  make_md(so.alv.30.AT2.high_inflammatory, "30dpi", "High"),
  make_md(so.alv.30.AT2.low_inflammatory, "30dpi", "Low")
)

md_all$condition <- gsub("^\\d+dpi_", "", md_all$condition_dpi)
md_all$condition <- factor(md_all$condition, levels = c("floxed", "dAT2"))
md_all$inflammatory <- factor(md_all$inflammatory, levels = c("Low", "High"))
md_all$dpi <- factor(md_all$dpi, levels = c("14dpi", "30dpi"))

df_all <- md_all %>%
  group_by(dpi, inflammatory, condition, Phase) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(dpi, inflammatory, condition) %>%
  mutate(prop = n / sum(n))

p <- ggplot(df_all, aes(x = condition, y = prop, fill = Phase)) +
  geom_col(position = "fill", width = 0.6) +
  facet_grid(~ dpi + inflammatory, scales = "free_x") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("G1" = "#E76F51", "S" = "#2A9D8F", "G2M" = "#264653")) +
  labs(x = NULL, y = "Proportion", fill = "Phase") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(face = "bold"))
ggsave(p, paste0(save_path, "/cellcycle_proportion_combined.pdf"), width = 8, height = 5)
ggsave(p, paste0(save_path, "/cellcycle_proportion_combined.png"), width = 8, height = 5)
ggplot2pptx(p, filename = paste0(save_path, "/cellcycle_proportion_combined.pptx"), width = 8, height = 5)
