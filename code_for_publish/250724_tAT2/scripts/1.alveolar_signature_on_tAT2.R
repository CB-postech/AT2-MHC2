### this script is for
### run Scran to normalize the data
### conda activate project_lung_exercise_R
# based on https://github.com/CB-postech/2022_KOGO_workshop/blob/main/KOGO_LCB_pipeline.md

library(Seurat)
library(scater)
library(SingleCellExperiment)
library(magrittr)
library(ggplot2)
library(BiocParallel)
library(future)
library(stringr)

source('/home/sjcho/yard/functions/R/seurat_count_to_normalization.R')
source('/home/sjcho/yard/functions/R/FeaturePlot_sjcho.R')
source('/home/sjcho/yard/functions/R/save_ggplot2_to_ppt.R')

save_path = '/home/sjcho/projects/AT2_MHC2/code_for_publish/250724_tAT2/outs/1.alveolar_signature_on_tAT2'
so.alv <- readRDS('/home/sjcho/projects/AT2_MHC2/20241113_after_cellbender/3.whole_cell_annotation/outs/3.1.epithelial_annotation/alv_annotation.rds')

so.alv$condition_dpi = factor(so.alv$condition_dpi, levels = c('naive_floxed', '7dpi_floxed', '14dpi_floxed', '30dpi_floxed', 
                                                                'naive_dAT2', '7dpi_dAT2', '14dpi_dAT2', '30dpi_dAT2'))

alv_cols = c(AT2 = "#8192ef", tAT2 = "#dfd13a", AT1 = "#f66b07", AT2.INF = "#000000", AT2.prolif. = "#8bd690")

p <- DimPlot(so.alv, group.by = 'seurat_clusters')
ggsave(p, file = paste0(save_path, 'alv_dimplot.png'), width = 8, height = 6)

### 1. celltypist annotation
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

### pseudotime
pseudotime <- read.csv('/home/sjcho/projects/AT2_MHC2/20241113_after_cellbender/manuscript_figures/figure4_and_exd6_20250225/palantir/n_pcs8_palantir_pseudotime.csv')
colnames(pseudotime) = c('cellname', 'pseudotime')
rownames(pseudotime) = pseudotime$cellname
pseudotime$condition = so.alv$condition[rownames(pseudotime)]; pseudotime$condition = factor(pseudotime$condition, levels = c('dAT2', 'floxed'))
pseudotime$celltype = so.alv$majority_voting[rownames(pseudotime)]; pseudotime$celltype = factor(pseudotime$celltype, levels = c('AT2', 'tAT2', 'AT1'))
pseudotime$dpi = so.alv$dpi[rownames(pseudotime)]

so.alv$pseudotime = pseudotime[Cells(so.alv), 'pseudotime']

### load markers
load('/home/sjcho/datas/reference_atlas/mouse_lung_alveolar/outs/1.5.markers/MAST/markers_level0_AT2.Rdata')
AT2.markers <- markers %>% subset(avg_log2FC > log2(2) & p_val_adj < 0.05 & pct.1 > 0.3)

load('/home/sjcho/datas/reference_atlas/mouse_lung_alveolar/outs/1.5.markers/MAST/markers_level0_AT1.Rdata')
AT1.markers <- markers %>% subset(avg_log2FC > log2(2) & p_val_adj < 0.05 & pct.1 > 0.3)

load('/home/sjcho/datas/reference_atlas/mouse_lung_alveolar/outs/1.5.markers/MAST/markers_level0_Transitory.Rdata')
tAT2.markers <- markers %>% subset(avg_log2FC > log2(2) & p_val_adj < 0.05 & pct.1 > 0.3)

#### marker expression & pseudotime at 14dpi
so.14.tAT2 <- subset(so.alv, dpi == '14dpi' & majority_voting == 'tAT2')

so.14.tAT2 <- AddModuleScore(so.14.tAT2, features = list(rownames(tAT2.markers)), name = 'tAT2_signature')
so.14.tAT2 <- AddModuleScore(so.14.tAT2, features = list(rownames(AT1.markers)), name = 'AT1_signature')
so.14.tAT2 <- AddModuleScore(so.14.tAT2, features = list(rownames(AT2.markers)), name = 'AT2_signature')

t.test(so.14.tAT2@meta.data[Cells(so.14.tAT2)[so.14.tAT2$condition == 'dAT2'], 'AT1_signature1'], so.14.tAT2@meta.data[Cells(so.14.tAT2)[so.14.tAT2$condition == 'floxed'], 'AT1_signature1'])
t.test(so.14.tAT2@meta.data[Cells(so.14.tAT2)[so.14.tAT2$condition == 'dAT2'], 'AT2_signature1'], so.14.tAT2@meta.data[Cells(so.14.tAT2)[so.14.tAT2$condition == 'floxed'], 'AT2_signature1'])
t.test(so.14.tAT2@meta.data[Cells(so.14.tAT2)[so.14.tAT2$condition == 'dAT2'], 'tAT2_signature1'], so.14.tAT2@meta.data[Cells(so.14.tAT2)[so.14.tAT2$condition == 'floxed'], 'tAT2_signature1'])
t.test(so.14.tAT2@meta.data[Cells(so.14.tAT2)[so.14.tAT2$condition == 'dAT2'], 'pseudotime'], so.14.tAT2@meta.data[Cells(so.14.tAT2)[so.14.tAT2$condition == 'floxed'], 'pseudotime'])

### marker expression & pseudotime at 30dpi
so.30.tAT2 <- subset(so.alv, dpi == '30dpi' & majority_voting == 'tAT2')

so.30.tAT2 <- AddModuleScore(so.30.tAT2, features = list(rownames(tAT2.markers)), name = 'tAT2_signature')
so.30.tAT2 <- AddModuleScore(so.30.tAT2, features = list(rownames(AT1.markers)), name = 'AT1_signature')
so.30.tAT2 <- AddModuleScore(so.30.tAT2, features = list(rownames(AT2.markers)), name = 'AT2_signature')

t.test(so.30.tAT2@meta.data[Cells(so.30.tAT2)[so.30.tAT2$condition == 'dAT2'], 'AT1_signature1'], so.30.tAT2@meta.data[Cells(so.30.tAT2)[so.30.tAT2$condition == 'floxed'], 'AT1_signature1'])
t.test(so.30.tAT2@meta.data[Cells(so.30.tAT2)[so.30.tAT2$condition == 'dAT2'], 'AT2_signature1'], so.30.tAT2@meta.data[Cells(so.30.tAT2)[so.30.tAT2$condition == 'floxed'], 'AT2_signature1'])
t.test(so.30.tAT2@meta.data[Cells(so.30.tAT2)[so.30.tAT2$condition == 'dAT2'], 'tAT2_signature1'], so.30.tAT2@meta.data[Cells(so.30.tAT2)[so.30.tAT2$condition == 'floxed'], 'tAT2_signature1'])
t.test(so.30.tAT2@meta.data[Cells(so.30.tAT2)[so.30.tAT2$condition == 'dAT2'], 'pseudotime'], so.30.tAT2@meta.data[Cells(so.30.tAT2)[so.30.tAT2$condition == 'floxed'], 'pseudotime'])


### inflammatory signature at 14, 30dpi AT2, tAT2, AT1
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

result_df = data.frame(GO_term = character(),
                          dpi = character(),
                          celltype = character(),
                          p.value = numeric(),
                          t.statistic = numeric(),
                          stringsAsFactors = FALSE)
for (GO_term in names(GO_list)) {
    so.alv <- AddModuleScore(so.alv, features = list(GO_list[[GO_term]]), name = GO_term)
    
    for (dpi in c('14dpi', '30dpi')) {
        for (celltype in c('AT2', 'tAT2', 'AT1')) {
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

# result_df에 x축 순서를 지정한 factor 변수 추가
result_df$group <- factor(
  paste(result_df$dpi, result_df$celltype),
  levels = c(
    "14dpi AT2",  "30dpi AT2",
    "14dpi tAT2", "30dpi tAT2",
    "14dpi AT1",  "30dpi AT1"
  )
)

result_df$dpi <- factor(result_df$dpi, levels = c("14dpi", "30dpi"))

library(ggplot2)
library(RColorBrewer)

# 3) ggplot
p <- ggplot() +
  # (a) p.value >= 0.05: t.statistic 색상
  geom_point(
    data = subset(result_df, p.value <= 0.05),
    aes(
      x     = dpi,
      y     = GO_term,
      size  = -log10(p.value),
      color = t.statistic
    ),
    position = 'identity'
  ) +
  # (b) p.value < 0.05: gray50 고정
  geom_point(
    data = subset(result_df, p.value > 0.05),
    aes(
      x     = dpi,
      y     = GO_term,
      size  = -log10(p.value)
    ),
    color    = "gray50",
    position = 'identity'
  ) +
  # 색상·크기 스케일
  scale_color_distiller(
    palette   = "RdBu",
    direction = -1,
    name      = "t statistic"
  ) +
  scale_size_continuous(
    range = c(2, 8),
    name  = expression(-log[10](p.value))
  ) +
  # 패널 분리
  facet_wrap(~ celltype, nrow = 1) +
  # 레이블 및 테마
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_classic() +
  theme(
    strip.text      = element_text(size = 14, face = "bold"),
    axis.text.x     = element_text(size = 12),
    axis.text.y     = element_text(size = 12),
    legend.position = "right"
  )
ggsave(p, file = paste0(save_path, '/inflammatory_response_t-test.png'), width = 16, height = 8)

############# DEG test 14dpi, 30dpi tAT2
floxed.tAT2.14dpi.cells = Cells(so.alv)[so.alv$condition == 'floxed' & so.alv$dpi == '14dpi' & so.alv$majority_voting == 'tAT2']
dAT2.tAT2.14dpi.cells = Cells(so.alv)[so.alv$condition == 'dAT2' & so.alv$dpi == '14dpi' & so.alv$majority_voting == 'tAT2']

DEG.14.tAT2 = FindMarkers(so.alv, ident.1 = floxed.tAT2.14dpi.cells, ident.2 = dAT2.tAT2.14dpi.cells, logfc.threshold = 0, min.pct = 0, verbose = FALSE)

floxed.tAT2.30dpi.cells = Cells(so.alv)[so.alv$condition == 'floxed' & so.alv$dpi == '30dpi' & so.alv$majority_voting == 'tAT2']
dAT2.tAT2.30dpi.cells = Cells(so.alv)[so.alv$condition == 'dAT2' & so.alv$dpi == '30dpi' & so.alv$majority_voting == 'tAT2']

DEG.30.tAT2 = FindMarkers(so.alv, ident.1 = floxed.tAT2.30dpi.cells, ident.2 = dAT2.tAT2.30dpi.cells, logfc.threshold = 0, min.pct = 0, verbose = FALSE)