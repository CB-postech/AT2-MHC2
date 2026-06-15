library(Matrix)
library(ggplot2)
library(cowplot)
library(ggrepel)
library(pheatmap)
library(flashier)
library(gbcd)
library(Seurat)
library(magrittr)
library(readxl)
library(org.Hs.eg.db)
library(dplyr)
library(tidyr)
library(data.table)
library(stringr)

source('/home/sjcho/yard/functions/R/FeaturePlot_sjcho.R')
source('/home/sjcho/projects/AT2_MHC2/utils.R')
source('/home/sjcho/projects/utils.R')
source('/home/sjcho/yard/functions/R/save_ggplot2_to_ppt.R')

save_path = '/home/sjcho/projects/AT2_MHC2/code_for_publish/251208_revision/human_relevance/GSE171524'

# library(sceasy)
# library(reticulate)
# use_condaenv('project_lung_exercise_R')
# loompy <- reticulate::import('loompy')

# ## Anndata to Seurat
sceasy::convertFormat(
  "/home/sjcho/datas/public_data/GSE171524/GSE171524_raw.h5ad",
  from = "anndata",
  to = "seurat",
  outFile = "/home/sjcho/datas/public_data/GSE171524/GSE171524_raw.rds"
)

so <- readRDS('/home/sjcho/datas/public_data/GSE171524/GSE171524_raw.rds')
var_info <- read.csv('/home/sjcho/datas/public_data/GSE171524/GSE171524_var_info.csv', row.names = 1)
# make a new column which contains unique gene symbols
var_info <- var_info %>%
  tibble::rownames_to_column("gene_id") %>% 
  add_count(feature_name, name = "n_count") %>% 
  mutate(new_feature_name = if_else(n_count > 1, gene_id, feature_name)) %>% 
  select(-n_count) %>%
  tibble::column_to_rownames("gene_id")

count <- so@assays$RNA$counts
rownames(count) <- var_info$new_feature_name
so <- CreateSeuratObject(counts = count, meta.data =so@meta.data, min.cells = 0, min.features = 0)  
so$cell_type_main %>% table

################# alv.epi. subset
so.alv.epi <- subset(so, cell_type_intermediate %in% c('AT1', 'AT2'))
so.alv.epi <- log_normalize_harmony(so.alv.epi, save_path = save_path, data_name = 'GSE171524_alv_epi', pcs = 20, batch_name = 'donor_id')
so.alv.epi[["percent.mt"]] <- PercentageFeatureSet(so.alv.epi, pattern = "^MT-")
so.alv.epi[["percent.ribo"]] <- PercentageFeatureSet(so.alv.epi, pattern = "^RP[SL]")
so.alv.epi[['log10UMI']] <- log10(so.alv.epi$nCount_RNA + 1)

p <- fp_sjcho(so.alv.epi, features = c('percent.mt', 'log10UMI', 'nFeature_RNA', 'percent.ribo', 'KRT8', 'CLDN4', 'TNIP3'), ncol = 4)
ggsave(p, file = paste0(save_path, '/alv_epi_qc.png'), width = 13, height = 6)
p <- DimPlot(so.alv.epi, cells.highlight = Cells(so.alv.epi)[so.alv.epi$log10UMI < 3], sizes.highlight = 0.5) + NoLegend()
ggsave(p, file = paste0(save_path, '/alv_epi_low_umi.png'), width = 5, height = 5)
p <- DimPlot(so.alv.epi, cells.highlight = Cells(so.alv.epi)[so.alv.epi$percent.mt > 5], sizes.highlight = 0.1) + NoLegend()
ggsave(p, file = paste0(save_path, '/alv_epi_high_mt.png'), width = 5, height = 5)

so.alv.epi.v2 <- subset(so.alv.epi, log10UMI >= 3 & percent.mt <= 5)
so.alv.epi.v2 <- log_normalize_harmony(so.alv.epi.v2, save_path = save_path, data_name = 'GSE171524_alv_epi', pcs = 20, batch_name = 'donor_id')
so.alv.epi.v2 <- FindClusters(so.alv.epi.v2, resolution = 1.5)

p <- DimPlot(so.alv.epi.v2, group.by = 'cell_type_fine', label = TRUE, repel = TRUE) + NoLegend()
ggsave(p, file = paste0(save_path, '/v2_alv_epi_celltype.png'), width = 5, height = 5)
p <- DimPlot(so.alv.epi.v2, group.by = 'donor_id') + NoLegend()
ggsave(p, file = paste0(save_path, '/v2_alv_epi_donor.png'), width = 5, height = 5)
p <- fp_sjcho(so.alv.epi.v2, features = c('MKI67', 'TOP2A', 'PDGFRA', 'PDGFRB', 'PTPRC', 'PECAM1'), ncol = 3, order = T)
ggsave(p, file = paste0(save_path, '/v2_alv_epi_lineage_markers.png'), width = 10, height = 5)
p <- fp_sjcho(so.alv.epi.v2, features = c('percent.mt', 'log10UMI', 'nFeature_RNA', 'percent.ribo', 'KRT8', 'CLDN4', 'TNIP3'), ncol = 4)
ggsave(p, file = paste0(save_path, '/v2_alv_epi_qc.png'), width = 13, height = 6)
p <- DimPlot(so.alv.epi.v2, group.by = 'seurat_clusters', label = T)
ggsave(p, file = paste0(save_path, '/v2_alv_epi_cluster.png'), width = 6, height = 5)

##### delete cluster17 : immune doublet
so.alv.epi.v3 <- subset(so.alv.epi.v2, seurat_clusters != 17)
so.alv.epi.v3 <- log_normalize_harmony(so.alv.epi.v3, save_path = save_path, data_name = 'GSE171524_alv_epi_v3', pcs = 15, batch_name = 'donor_id')

p <- fp_sjcho(so.alv.epi.v3, features = c('percent.mt', 'log10UMI', 'nFeature_RNA', 'percent.ribo', 'KRT8', 'CLDN4', 'TNIP3'), ncol = 4)
ggsave(p, file = paste0(save_path, '/v3_alv_epi_qc.png'), width = 13, height = 6)
p <- DimPlot(so.alv.epi.v3, group.by = 'cell_type_fine', label = TRUE, repel = TRUE) + NoLegend()
ggsave(p, file = paste0(save_path, '/v3_alv_epi_celltype.png'), width = 5, height = 5)
p <- DimPlot(so.alv.epi.v3, group.by = 'donor_id') + NoLegend()
ggsave(p, file = paste0(save_path, '/v3_alv_epi_donor.png'), width = 5, height = 5)
p <- fp_sjcho(so.alv.epi.v3, features = c('MKI67', 'TOP2A', 'PDGFRA', 'PDGFRB', 'PTPRC', 'PECAM1', 'KRT8', 'CLDN4', 'CDKN1A', 'THBS1', 'CDH2', 'PCDH7'), ncol = 3, order = T)
ggsave(p, file = paste0(save_path, '/v3_alv_epi_lineage_markers.png'), width = 10, height = 13)
p <- DimPlot(so.alv.epi.v3, group.by = 'seurat_clusters', label = T)
ggsave(p, file = paste0(save_path, '/v3_alv_epi_cluster.png'), width = 6, height = 5)

markers.8 <- FindMarkers(so.alv.epi.v3, ident.1 = 8, only.pos = T)
markers.9 <- FindMarkers(so.alv.epi.v3, ident.1 = 9, only.pos = T)

####### delete cluster8 : stromal doublet (TBHS1, CDH2, COL1A1, PCDH7)
so.alv.epi.v4 <- subset(so.alv.epi.v3, seurat_clusters != 8)
so.alv.epi.v4 <- log_normalize_harmony(so.alv.epi.v4, save_path = save_path, data_name = 'GSE171524_alv_epi_v4', pcs = 20, batch_name = 'donor_id')
so.alv.epi.v4 <- FindClusters(so.alv.epi.v4, resolution = 2)

p <- fp_sjcho(so.alv.epi.v4, features = c('percent.mt', 'log10UMI', 'nFeature_RNA', 'percent.ribo', 'KRT8', 'CLDN4', 'CXCL17', 'LCN2', 'SFTPC', 'NAPSA', 'ABCA3', 'PDPN', 'AGER'), ncol = 4, order = T)
ggsave(p, file = paste0(save_path, '/v4_alv_epi_qc.png'), width = 13, height = 9)
p <- DimPlot(so.alv.epi.v4, group.by = 'cell_type_fine', label = TRUE, repel = TRUE) + NoLegend()
ggsave(p, file = paste0(save_path, '/v4_alv_epi_celltype.png'), width = 5, height = 5)
p <- DimPlot(so.alv.epi.v4, group.by = 'donor_id') + NoLegend()
ggsave(p, file = paste0(save_path, '/v4_alv_epi_donor.png'), width = 5, height = 5)
p <- fp_sjcho(so.alv.epi.v4, features = c('MKI67', 'TOP2A', 'PDGFRA', 'PDGFRB', 'PTPRC', 'PECAM1', 'KRT8', 'CLDN4', 'CDKN1A'), ncol = 3, order = T)
ggsave(p, file = paste0(save_path, '/v4_alv_epi_lineage_markers.png'), width = 10, height = 13)
p <- DimPlot(so.alv.epi.v4, group.by = 'seurat_clusters', label = T)
ggsave(p, file = paste0(save_path, '/v4_alv_epi_cluster.png'), width = 6, height = 5)
p <- VlnPlot(so.alv.epi.v4, features = c('KRT8', 'CLDN4', 'NDRG1', 'LCN2', 'TNIP3', 'CXCL17'), group.by = 'seurat_clusters', pt.size = 0) & geom_boxplot(width = 0.5, outlier.shape = NA)
ggsave(p, file = paste0(save_path, '/v4_alv_epi_markers_vln.png'), width = 10, height = 5)

# annotation
# cluster 17 : tAT2
# cluster 14, 6 : pAT2
# cluster 2, 9, 5, 4, 15, 11, 19, 16, 22 : AT1
# cluster 8, 7, 13, 12, 1, 21, 10, 20, 3, 18, 0 : AT2
so.alv.epi.v4$anno <- as.vector(so.alv.epi.v4$seurat_clusters)
so.alv.epi.v4$anno[so.alv.epi.v4$seurat_clusters %in% c(17, 14, 6)] <- 'tAT2'
so.alv.epi.v4$anno[so.alv.epi.v4$seurat_clusters %in% c(2, 9, 5, 4, 15, 11, 19, 16, 22)] <- 'AT1'
so.alv.epi.v4$anno[so.alv.epi.v4$seurat_clusters %in% c(8, 7, 13, 12, 1, 21, 10, 20, 3, 18, 0)] <- 'AT2'
so.alv.epi.v4$anno <- factor(so.alv.epi.v4$anno, levels = c('AT2', 'tAT2', 'AT1'))

p <- DotPlot(so.alv.epi.v4, features = c('SFTPC', 'ETV5', 'ABCA3', 'CDKN1A', 'KRT8', 'CLDN4', 'AGER', 'TIMP3', 'SPOCK2'), group.by = 'anno') + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p <- p + set_color_RdBu
ggsave(p, file = paste0(save_path, '/figure_v4_alv_epi_anno_dot.png'), width = 8, height = 5)
ggsave(p, file = paste0(save_path, '/figure_v4_alv_epi_anno_dot.pdf'), width = 8, height = 5)
saveRDS(p, file = paste0(save_path, '/figure_v4_alv_epi_anno_dot.rds'))

cols.use = c('AT2' = '#7687e8', 'tAT2' = '#dbcb45', 'AT1' = '#e3783d')
p <- DimPlot(so.alv.epi.v4, group.by = 'anno', cols = cols.use) + NoLegend() + set_UMAP
ggsave(p, file = paste0(save_path, '/figure_v4_alv_epi_anno_dimplot_color.png'), width = 5, height = 5)
ggsave(p, file = paste0(save_path, '/figure_v4_alv_epi_anno_dimplot_color.pdf'), width = 5, height = 5)

saveRDS(so.alv.epi.v4, file = paste0(save_path, '/so_alv_epi_v4.rds'))

so.alv.epi.v4[["RNA"]] <- as(so.alv.epi.v4[["RNA"]], "Assay")
sceasy::convertFormat(
  so.alv.epi.v4,
  from = "seurat",
  to = "anndata",
  main_layer = "data",        
  transfer_layers = c("counts"),
  outFile= paste0(save_path, '/so_alv_epi_v4.h5ad')
)

p <- VlnPlot(so.alv.epi.v4, features = c('ETV5'), group.by = 'seurat_clusters') + NoLegend() + geom_boxplot(width = 0.3)
ggsave(p, file = paste0(save_path, '/v4_alv_ETV5_cluster.png'), width = 6, height = 5)
cells_cluster3 = Cells(so.alv.epi.v4)[which(so.alv.epi.v4$seurat_clusters == 3)]
cell_ETV5_highest_in_cluster3 = names(which(so.alv.epi.v4@assays$RNA@data['ETV5', cells_cluster3] == max(so.alv.epi.v4@assays$RNA@data['ETV5', cells_cluster3])))
p <- DimPlot(so.alv.epi.v4, cells.highlight = cell_ETV5_highest_in_cluster3, sizes.highlight = 3) + NoLegend()
ggsave(p, file = paste0(save_path, '/v4_ETV5_highest_cell.png'), width = 5, height = 5)

so.alv.epi.v4 <- readRDS('/home/sjcho/projects/AT2_MHC2/code_for_publish/251208_revision/human_relevance/GSE171524/so_alv_epi_v4.rds')

################# Palantir visualization
### load palantir result
branch <- fread('/home/sjcho/projects/AT2_MHC2/code_for_publish/251208_revision/human_relevance/GSE171524/branch_masks.csv', header = T)
branch <- as.data.frame(branch)
rownames(branch) <- branch$V1
branch <- branch[, -1]

pseudotime <- read.csv('/home/sjcho/projects/AT2_MHC2/code_for_publish/251208_revision/human_relevance/GSE171524/palantir_pseudotime.csv')
colnames(pseudotime) = c('cellname', 'pseudotime')
rownames(pseudotime) = pseudotime$cellname
so.alv.epi.v4$pseudotime = pseudotime[colnames(so.alv.epi.v4), 'pseudotime']

cell_ETV5_highest_in_cluster0 = 'GGGAGTAAGAAACTGT-1_7'
so.start <- subset(so.alv.epi.v4, cells = cell_ETV5_highest_in_cluster0)

### path probability & pseudotime
probability <- read.csv('/home/sjcho/projects/AT2_MHC2/code_for_publish/251208_revision/human_relevance/GSE171524/palantir_fate_probabilities.csv')
rownames(probability) <- probability$X
so.alv.epi.v4$path1_prob <- probability[colnames(so.alv.epi.v4), 'ACTATGGGTCTGCGCA.1_1']
so.alv.epi.v4$path2_prob <- probability[colnames(so.alv.epi.v4), 'TAATTCCAGATTGACA.1_18']

p <- FeaturePlot(so.alv.epi.v4, features = 'path1_prob', pt.size = 0.75, order = T) + NoAxes() + NoLegend() + scale_colour_gradientn(colours = colorRampPalette(c("darkgray", 'green', "black"))(100))
p <- p + theme(plot.title = element_text(size = 0))
ggsave(p, file = paste0(save_path, '/fig.path1.png'), width = 4, height = 4)
ggsave(p, file = paste0(save_path, '/fig.path1.pdf'), width = 4, height = 4)
ggplot2pptx(p, 4, 4, paste0(save_path, '/fig.path1.pptx'))

p <- FeaturePlot(so.alv.epi.v4, features = 'path2_prob', pt.size = 0.75, order = T) + NoAxes() + NoLegend() + scale_colour_gradientn(colours = colorRampPalette(c("darkgray", 'green', "black"))(100))
p <- p + theme(plot.title = element_text(size = 0))
ggsave(p, file = paste0(save_path, '/fig.path2.png'), width = 4, height = 4)
ggsave(p, file = paste0(save_path, '/fig.path2.pdf'), width = 4, height = 4)
ggplot2pptx(p, 4, 4, paste0(save_path, '/fig.path2.pptx'))
leg <- get_legend(FeaturePlot(so.alv.epi.v4, features = 'pseudotime', pt.size = 0.75, order = T) + 
        scale_colour_gradientn(colours = colorRampPalette(c("darkgray", 'green', "black"))(100)))
ggplot2pptx(ggpubr::as_ggplot(leg), 1, 6, paste0(save_path, '/fig_pseudotime_legend.pptx'))

p <- FeaturePlot(so.alv.epi.v4, features = 'pseudotime', pt.size = 0.75, order = T) + NoAxes() + NoLegend() + scale_colour_gradientn(colours = colorRampPalette(c("darkgray", 'green', "black"))(100))
p <- p + geom_point(data = so.start@reductions$umap@cell.embeddings, aes(x = umap_1, y = umap_2), color = 'darkred', size = 4)
p <- p + theme(plot.title = element_text(size = 0))
ggsave(p, file = paste0(save_path, '/fig.pseudotiem.png'), width = 4, height = 4)
ggsave(p, file = paste0(save_path, '/fig.pseudotiem.pdf'), width = 4, height = 4)
ggplot2pptx(p, 4, 4, paste0(save_path, '/fig.pseudotiem.pptx'))

############ signature score
library(GO.db)
library(fgsea)
library(org.Hs.eg.db)
library(msigdbr)
h_gene_sets = msigdbr(species = "human", category = "H")
GOBP_sets = msigdbr(species = "human", category = "C5", subcategory = "BP")
child.terms <- GOBPCHILDREN[['GO:0006954']] # inflammatory response

genes <- AnnotationDbi::select(org.Hs.eg.db,
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
GO_list_use = GO_list[c('acute inflammatory response', 'HM_inflammatory_response', 'inflammatory response to antigenic stimulus')]

so.alv.epi.v4 <- AddModuleScore(so.alv.epi.v4, features = GO_list_use, name = names(GO_list_use))
p <- FeaturePlot(so.alv.epi.v4, features = paste0(names(GO_list_use), seq(1, length(GO_list_use))), ncol = 4, order = T)
p <- p & set_color_RdBu
ggsave(p, file = paste0(save_path, '/test.signature.whole.epi.png'), width = 16, height = 16)

so.tAT2 <- subset(so.alv.epi.v4, subset = anno %in% c('tAT2'))
plot_df <- so.tAT2@meta.data %>%
  dplyr::select(pseudotime, paste0(names(GO_list_use), seq(1, length(GO_list_use)))) %>% 
  pivot_longer(cols = -pseudotime, names_to = "Term", values_to = "Score")

df_cleaned <- plot_df %>%
  mutate(Term = Term %>%
    str_replace_all("_", " ") %>%
    str_remove("\\d+$") %>%
    str_replace_all("(to|of|in) ", "\\1\n")
  )

p <- ggplot(df_cleaned, aes(x = pseudotime, y = Score)) +
  geom_point(alpha = 1, size = 0.1, color = 'black') +  
  geom_smooth(method = "loess", color = "red", size = 1) + 
  facet_wrap(~Term, scales = "free_y") +
  labs(x = "Pseudotime", y = "Score") +
  theme_classic() +
  theme(legend.position = "none",
        strip.text = element_text(size = 8, face = "bold"))
ggsave(p, file = paste0(save_path, '/fig.pseudotime.signature.png'), width = 10, height = 3)
ggsave(p, file = paste0(save_path, '/fig.pseudotime.signature.pdf'), width = 10, height = 3)
saveRDS(p, file = paste0(save_path, '/fig.pseudotime.signature.rds'))

### 
plot_df <- so.alv.epi.v4@meta.data %>%
  dplyr::select(pseudotime, paste0(names(GO_list_use), seq(1, length(GO_list_use)))) %>% 
  pivot_longer(cols = -pseudotime, names_to = "Term", values_to = "Score")

df_cleaned <- plot_df %>%
  mutate(Term = Term %>%
    str_replace_all("_", " ") %>%
    str_remove("\\d+$") %>%
    str_replace_all("(to|of|in) ", "\\1\n")
  )

p <- ggplot(df_cleaned, aes(x = pseudotime, y = Score)) +
  geom_point(alpha = 1, size = 0.1, color = 'black') +  
  geom_smooth(method = "loess", color = "red", size = 1) + 
  facet_wrap(~Term, scales = "free_y") +
  labs(x = "Pseudotime", y = "Score") +
  theme_classic() +
  theme(legend.position = "none",
        strip.text = element_text(size = 8, face = "bold"))
ggsave(p, file = paste0(save_path, '/fig.pseudotime.signature.total.png'), width = 10, height = 3)
ggsave(p, file = paste0(save_path, '/fig.pseudotime.signature.total.pdf'), width = 10, height = 3)
saveRDS(p, file = paste0(save_path, '/fig.pseudotime.signature.rds'))
