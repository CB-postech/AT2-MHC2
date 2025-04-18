library(dplyr)
library(Seurat)
library(tidyverse)
library(data.table)
library(biomaRt)

# load custom functions
source('/home/sjcho/yard/functions/R/seurat_count_to_normalization.R')
source('/home/sjcho/yard/functions/R/FeaturePlot_sjcho.R')
source('/home/sjcho/yard/functions/R/draw_proportion.R')
source('/home/sjcho/projects/AT2_MHC2/20241113_after_cellbender/utils.R')
source('/home/sjcho/yard/functions/R/draw_celltype_wise_dotplot.R')
source('/home/sjcho/yard/functions/R/gsea_sjcho.R')
source('/home/sjcho/yard/functions/R/save_ggplot2_to_ppt.R')

save_path = '/home/sjcho/projects/AT2_MHC2/20241113_after_cellbender/PND_Rcode_to_sjcho/figure1_sjcho_result/'

library(msigdbr)
BP <- msigdbr(species = "mouse", category = "C5")
MHC2_protein_complex_genes = subset(BP, gs_name == 'GOCC_MHC_CLASS_II_PROTEIN_COMPLEX')$gene_symbol

# no extra QC. used count matrix in corrresponding reference
so <- readRDS('/home/sjcho/projects/AT2_MHC2/20241113_after_cellbender/PND_Rcode_to_sjcho/figure1_sjcho_result/so_normalized.rds')

p <- DimPlot(so, group.by = 'cell_type')
ggsave(p, file = paste0(save_path, 'Fig1.a_celltype.png'), width = 9, height = 5)

### AT2 markers and MHC2 expression through developmental stages ###
load('/home/sjcho/datas/reference_atlas/mouse_lung_alveolar/outs/20250217_public_figures/markers/AT2_uninjury_marker.RData')
# from mouse alvoelar cell atlas
# save as AT2_uninjuryed_markers
marker.sig <- subset(AT2_uninjuryed_markers, avg_log2FC > log2(2) & p_val_adj < 0.05 & pct.1 > 0.5)

# E16.5-E18 / P1-P3 / p7-p14 / P14-p28/ adult
so$var_time = factor(so$var_time, levels = c('E12.5', 'E15.5', 'E17.5', 'P3', 'P7', 'P15', 'Adult'))
so.AT2 <- subset(so, cells = Cells(so)[so$cell_type == 'pulmonary alveolar type 2 cell'])
so.AT2 <- AddModuleScore(so.AT2, list(rownames(marker.sig)), name = 'AT2_uninjuryed_markers')
p1 <- VlnPlot(so.AT2, features = 'AT2_uninjuryed_markers1', group.by = 'var_time', cols = rep('gray30', length(unique(so$var_time))), pt.size = 0) + geom_boxplot(fill = 'white', width = 0.2)
p1 <- p1 + NoLegend()

so.AT2 <- AddModuleScore(so.AT2, list(MHC2_protein_complex_genes), name = 'GOCC_MHC_CLASS_II_PROTEIN_COMPLEX')
p2 <- VlnPlot(so.AT2, features = 'GOCC_MHC_CLASS_II_PROTEIN_COMPLEX1', group.by = 'var_time', cols = rep('gray60', length(unique(so$var_time))), pt.size = 0) + geom_boxplot(fill = 'white', width = 0.2)
p2 <- p2 + NoLegend()
p <- as.ggplot(egg::ggarrange(p1, p2, ncol = 2))
ggsave(p, file = paste0(save_path, 'exd1.b_AT2_GOCC_MHC_CLASS_II_PROTEIN_COMPLEX.png'), width = 8, height = 4)
ggsave(p, file = paste0(save_path, 'exd1.b_AT2_GOCC_MHC_CLASS_II_PROTEIN_COMPLEX.pdf'), width = 8, height = 4)
ggplot2pptx(p, file = paste0(save_path, 'exd1.b_AT2_GOCC_MHC_CLASS_II_PROTEIN_COMPLEX.pptx'), width = 8, height = 4)

p <- VlnPlot(so.AT2, features = 'H2-Ab1', group.by = 'var_time', cols = rep('gray10', length(unique(so$var_time))), pt.size = 0) + geom_boxplot(fill = 'gray95', width = 0.2)
p <- p + NoLegend()
ggsave(p, file = paste0(save_path, 'exd1.b_AT2_H2-Ab1.png'), width = 9, height = 5)
ggsave(p, file = paste0(save_path, 'exd1.b_AT2_H2-Ab1.pdf'), width = 9, height = 5)
ggplot2pptx(p, file = paste0(save_path, 'exd1.b_AT2_H2-Ab1.pptx'), width = 9, height = 5)

### exd 1.e. AT2 marker signature score & MHC2 expression correlation
## Acquire expression matrix
so.AT2 <- NormalizeData(so.AT2)
expr_matrix <- as.matrix(expm1(so.AT2@assays$RNA$data))

expr_matrix <- rbind(expr_matrix, as.vector(so.AT2$GOCC_MHC_CLASS_II_PROTEIN_COMPLEX1))
rownames(expr_matrix)[nrow(expr_matrix)] <- "GOCC_MHC_CLASS_II_PROTEIN_COMPLEX1"

##### correlation with GOCC_MHC_CLASS_II_PROTEIN_COMPLEX1
library(progressr)
handlers(global = TRUE)
handlers("txtprogressbar")  

p_value = numeric(nrow(expr_matrix))   
correlation = numeric(nrow(expr_matrix))

with_progress({
  p <- progressor(along = seq_len(nrow(expr_matrix)))

  results <- apply(expr_matrix, 1, function(x) {
    res <- cor.test(x, so.AT2$GOCC_MHC_CLASS_II_PROTEIN_COMPLEX1, method = 'pearson')
    p()
    return(c(res$estimate, res$p.value))
  })
})

correlation <- results[1, ]
p_value <- results[2, ]

df.cor = data.frame(correlation = correlation, p_value = p_value, row.names = rownames(expr_matrix))
df.cor = df.cor[df.cor$p_value %>% is.na() %>% `!`, ]

## gsea
df.cor.correlation <- df.cor$correlation
names(df.cor.correlation) <- rownames(df.cor)
GOCC = subset(BP, gs_subcat == 'GO:CC')
msigdbr_list = split(x = GOCC$gene_symbol, f = GOCC$gs_name) # load pathway

library(fgsea)
fgseaRes <- fgsea(msigdbr_list, df.cor.correlation)
# order by adjusted p.value
fgseaRes <- fgseaRes[order(fgseaRes$padj), ]
library(data.table)
fwrite(fgseaRes, paste0(save_path, 'exd1.e.GSEA_with_MHC2.csv'))
fgseaRes[fgseaRes$NES > 0, ]$pathway %>% head
fgseaRes[fgseaRes$pathway == 'GOCC_ALVEOLAR_LAMELLAR_BODY', ]
fgseaRes[fgseaRes$pathway == 'GOCC_LAMELLAR_BODY', ]

### top barplot
fgseaRes[fgseaRes$NES > 0, ]$pathway -> t
# GOCC_MULTIVESICULAR_BODY
# GOCC_SECRETORY_GRANULE_MEMBRANE 

NES.pos <- fgseaRes[fgseaRes$NES > 0 & fgseaRes$padj < 0.05, ]
NES.pos <- NES.pos[order(NES.pos$NES, decreasing = T), ]

top30 <- NES.pos %>%
  arrange(desc(NES)) %>%
  slice_head(n = 30)

top30$pathway <- factor(top30$pathway, levels = rev(top30$pathway))

p <- ggplot(top30, aes(x = pathway, y = NES)) +
  geom_bar(stat = "identity", fill = "darkred") +
  coord_flip() +
  theme_classic() +
  labs(
    title = "",
    x = "",
    y = "Normalized Enrichment Score"
  )
ggsave(p, file = paste0(save_path, 'exd1.c.top30_GSEA.png'), width = 10, height = 5)
ggsave(p, file = paste0(save_path, 'exd1.c.top30_GSEA.pdf'), width = 10, height = 5)
ggplot2pptx(p, file = paste0(save_path, 'exd1.c.top30_GSEA.pptx'), width = 10, height = 5)
