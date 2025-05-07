######## QC scRNA-seq datasets to build Mouse Alveolar Epithelium Reference
######## GSE261793, GSE145031, GSE204787
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

source('/home/sjcho/datas/reference_atlas/mouse_lung_alveolar/scripts/utils.R')
source('/home/sjcho/yard/functions/R/seurat_count_to_normalization.R')
source('/home/sjcho/yard/functions/R/FeaturePlot_sjcho.R')
source('/home/sjcho/yard/functions/R/save_ggplot2_to_ppt.R')
save_path = '/home/sjcho/datas/reference_atlas/mouse_lung_alveolar/outs/1.alveolar_regeneration_datasets_only/'

# run GSEA
library(fgsea)
library(org.Mm.eg.db)
library(msigdbr)
h_gene_sets = msigdbr(species = "mouse", category = 'H')
bp_gene_sets = msigdbr(species = "mouse", category = 'C5', subcategory = 'BP')

### GSE261793
so_list_GSE261793 <- list(
    'annotattaed' = readRDS('/home/sjcho/datas/public_data/mouse_lung_alveolar/GSE261793/GSE261793_processed_seurat_annotated.rds'),
    'antibody' = readRDS('/home/sjcho/datas/public_data/mouse_lung_alveolar/GSE261793/GSE261793_processed_seurat_antibody.rds')
)

# have to change cell barcode to contain sample information
count_list_GSE261793 <- list(
    'annotattaed' = as.matrix(so_list_GSE261793[['annotattaed']]@assays$RNA$counts),
    'antibody' = as.matrix(so_list_GSE261793[['antibody']]@assays$RNA$counts)
)

meta_list_GSE261793 <- list(
    'annotattaed' = cbind(so_list_GSE261793[['annotattaed']]@meta.data[, c('gender', 'sample')], 'injury_type' = 'bleomycin', 'dpi' = 12, 'dataset' = 'GSE261793', 'condition' = 'NULL', 'sort' = 'NULL'),
    'antibody' = cbind('gender' = 'NULL', 'sample' = so_list_GSE261793[['antibody']]@meta.data[, 'sample'], 'injury_type' = 'bleomycin', 'dpi' = 12, 'dataset' = 'GSE261793', 'condition' = 'NULL', 'sort' = 'NULL')
)

### GSE145031
subdir_name = c('GSM4304609', 'GSM4304610', 'GSM4304611', 'GSM4304612', 'GSM4304613', 'GSM4304614')
dpi = c(0, 0, 14, 14, 28, 28); names(dpi) = subdir_name
sort = c('Tomato', 'nonTomato', 'Tomato', 'nonTomato', 'Tomato', 'nonTomato'); names(sort) = subdir_name
sce_list_GSE145031 <- list()
meta_list_GSE145031 <- list()
for (subdir in subdir_name) {
    dir_path = paste0('/home/sjcho/datas/public_data/mouse_lung_alveolar/GSE145031/', subdir)
    files = list.files(dir_path, full.names = T)
    lapply(files, rename_files)

    # sce_list_GSE145031[[subdir]] <- run_dropletUtils(dir_path, save_path, paste0('GSE145031_', subdir))
    cellnumbers = sce_list_GSE145031[[subdir]] %>% colnames %>% length
    meta_list_GSE145031[[subdir]] <- data.frame('sample' = rep(subdir, cellnumbers), 
                                                'dpi' = rep(dpi[subdir], cellnumbers), 
                                                'sort' = rep(sort[subdir], cellnumbers), 
                                                'injury_type' = rep('bleomycin', cellnumbers), 
                                                'dataset' = rep('GSE145031', cellnumbers),
                                                'condition' = rep('NULL', cellnumbers),
                                                'gender' = rep('NULL', cellnumbers))
}
count_list_GSE145031 <- lapply(sce_list_GSE145031, function(x) counts(x))

# now we get rid of empty droplets

### GSE204787
# its only offer filtered h5
count_list_GSE204787 <- list()
meta_list_GSE204787 <- list()

tmp = Read10X_h5('/home/sjcho/datas/public_data/mouse_lung_alveolar/GSE204787/GSM6191345_Tfcp2l1_KO_SPCCRE_YFPpos.h5')
colnames(tmp) = paste0('Tfcp2l1_KO_SPCCRE_YFPneg_', colnames(tmp))
cellnumbers = tmp %>% colnames %>% length
count_list_GSE204787[['Tfcp2l1_KO_SPCCRE_YFPneg']] = tmp
meta_list_GSE204787[['Tfcp2l1_KO_SPCCRE_YFPneg']] = data.frame('sample' = rep('Tfcp2l1_KO_SPCCRE_YFPneg', cellnumbers),
                                                                'dpi' = rep(14, cellnumbers),
                                                                'injury_type' = rep('H1N1_PR8', cellnumbers),
                                                                'dataset' = rep('GSE204787', cellnumbers),
                                                                'condition' = rep('Tfcp2l1_KO', cellnumbers),
                                                                'sort' = rep('NULL', cellnumbers),
                                                                'gender' = rep('NULL', cellnumbers))

tmp = Read10X_h5('/home/sjcho/datas/public_data/mouse_lung_alveolar/GSE204787/GSM6191346_WT_SPCCRE_YFPpos_14dpi.h5')
colnames(tmp) = paste0('WT_SPCCRE_YFPpos_14dpi_', colnames(tmp))
cellnumbers = tmp %>% colnames %>% length
count_list_GSE204787[['WT_SPCCRE_YFPpos_14dpi']] = tmp
meta_list_GSE204787[['WT_SPCCRE_YFPpos_14dpi']] = data.frame('sample' = rep('WT_SPCCRE_YFPpos_14dpi', cellnumbers),
                                                             'dpi' = rep(14, cellnumbers),
                                                             'injury_type' = rep('H1N1_PR8', cellnumbers),
                                                             'dataset' = rep('GSE204787', cellnumbers),
                                                             'condition' = rep('WT', cellnumbers),
                                                             'sort' = rep('NULL', cellnumbers),
                                                             'gender' = rep('NULL', cellnumbers))

### Merge it all
count_list <- c(count_list_GSE261793, count_list_GSE145031, count_list_GSE204787)
common_genes <- Reduce(intersect, lapply(count_list, rownames))
count_list <- lapply(count_list, function(x) x[common_genes, ])
# merge it
count <- do.call(cbind, count_list)
meta <- do.call(rbind, c(meta_list_GSE261793, meta_list_GSE145031, meta_list_GSE204787))

ref_cols <- names(meta_list_GSE261793[[1]])
meta_list_GSE261793 <- lapply(meta_list_GSE261793, function(x) x[, ref_cols])
meta_list_GSE145031 <- lapply(meta_list_GSE145031, function(x) x[, ref_cols])
meta_list_GSE204787 <- lapply(meta_list_GSE204787, function(x) x[, ref_cols])

meta <- do.call(rbind, c(meta_list_GSE261793, meta_list_GSE145031, meta_list_GSE204787))

so <- CreateSeuratObject(counts = count, meta.data = meta)
so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^mt-")

so.raw <- CreateSeuratObject(counts = count, meta.data = meta)
so.raw[["percent.mt"]] <- PercentageFeatureSet(so.raw, pattern = "^mt-")
so.raw[['log10_UMI']] <- log10(so.raw[['nCount_RNA']][[1]])

p <- ggplot(so.raw[['percent.mt']]) + geom_histogram(aes(x = percent.mt), bins = 250) + theme_classic()
p <- p + ggtitle('MT percent') + 
    theme(title = element_text(size = 20),
          axis.text = element_text(size = 15))
p <- p + geom_vline(xintercept = 7.5, linetype = 'dashed', color = 'red')
ggsave(p, file = paste0(save_path, 'histogram_percent_mt.png'), width = 12, height = 6)
p <- ggplot(so.raw[['log10_UMI']]) + geom_histogram(aes(x = log10_UMI), bins = 250) + theme_classic()
p <- p + ggtitle('log10_UMI') + 
    theme(title = element_text(size = 20),
          axis.text = element_text(size = 15))
p <- p + geom_vline(xintercept = 3, linetype = 'dashed', color = 'red')
ggsave(p, file = paste0(save_path, 'histogram_log10_UMIt.png'), width = 12, height = 6)

so <- subset(so.raw, nCount_RNA > 1000 & percent.mt < 15) 
# from 106438 
# 15 : 80574
# 13.5 : 75082
save_path_mt = '/home/sjcho/datas/public_data/mouse_lung_alveolar/outs/mt15/'
# so <- seurat_count_to_normalization(so, metas = ref_cols, 'alveoalr_merged', pcs = 20, fig_save = save_path, batch = 'dataset')

so <- NormalizeData(object = so)
so <- FindVariableFeatures(object = so)
so <- ScaleData(object = so)
so <- RunPCA(object = so, npcs = 20)
so <- RunHarmony(so, 'dataset')
so <- so %>%
    FindNeighbors(reduction = "harmony") %>%
    FindClusters(resolution = 1) 
so <- RunUMAP(so, dims = 1:20, reduction = 'harmony')

so <- FindClusters(so, resolution = 0.5)
p <- fp_sjcho(so, features = c('Epcam', 'Pecam1', 'Ptprc', 'Pdgfra'))
ggsave(p, file = paste0(save_path_mt, 'alveolar_merged_featureplot.png'), width = 10, height = 10)
p <- DimPlot(so, group.by = 'dataset', label = T, repel = T)
ggsave(p, file = paste0(save_path_mt, 'alveolar_merged_dataset.png'), width = 7, height = 6)
p <- DimPlot(so, group.by = 'seurat_clusters', label = T, repel = T)
ggsave(p, file = paste0(save_path_mt, 'alveolar_merged_seurat_clusters.png'), width = 7, height = 6)
p <- DimPlot(so, group.by = 'sample', label = T, repel = T)
ggsave(p, file = paste0(save_path_mt, 'alveolar_merged_seurat_samples.png'), width = 7, height = 6)
p <- fp_sjcho(so, features = c('Sftpc', 'Etv5', 'Abca3', 'Pdpn', 'Aqp5', 'Hopx', 'Foxj1', 'Scgb1a1', 'Scgb3a2', 'Krt5', 'Trp63', 'Cldn4', 'Sprr1a', 'Tnip3', 'Krt8'), ncol = 4)
ggsave(p, file = paste0(save_path_mt, 'alveolar_merged_epi_features.png'), width = 13, height = 15)
p <- fp_sjcho(so, features = c('Sftpc', 'Etv5', 'Foxj1', 
                                'Pdpn', 'Aqp5', 'Hopx',
                                'Tnip3', 'Cldn4', 'Krt8',
                                'Ptprc', 'Pecam1', 'Pdgfra'), ncol = 3)
ggsave(p, file = paste0(save_path_mt, 'alveolar_merged_epi_features2.png'), width = 11, height = 12)

so[['log10_UMI']] <- log10(so[['nCount_RNA']][[1]])
p <- fp_sjcho(so, features = c('log10_UMI', 'percent.mt', 'nFeature_RNA'), ncol = 2, order = T)
ggsave(p, file = paste0(save_path_mt, 'alveolar_merged_tech.png'), width = 10, height = 10)

p <- DimPlot(so, cells.highlight = Cells(so)[so[['percent.mt']][[1]] > 10], sizes.highlight = 0.1)
ggsave(p, file = paste0(save_path_mt, 'alveolar_merged_percent_mt10.png'), width = 7, height = 6)

alveolar_clusters <- c(0, 1, 2, 6, 7, 8, 12)
p <- DimPlot(so, cells.highlight = Cells(so)[so[['seurat_clusters']][[1]] %in% alveolar_clusters], sizes.highlight = 0.05, cols.highlight = 'darkred') + NoLegend()
ggsave(p, file = paste0(save_path_mt, 'alveolar_merged_alveolar_clusters.png'), width = 6, height = 6)

##### subset so.alv
so.alv <- CreateSeuratObject(so@assays$RNA$counts[, Cells(so)[so[['seurat_clusters']][[1]] %in% alveolar_clusters]], min.cells = 0, min.features = 0)
so.alv@meta.data[, ref_cols] <- so@meta.data[Cells(so)[so[['seurat_clusters']][[1]] %in% alveolar_clusters], ref_cols]
so.alv <- NormalizeData(object = so.alv)
so.alv <- FindVariableFeatures(object = so.alv)
so.alv <- ScaleData(object = so.alv)
so.alv <- RunPCA(object = so.alv, npcs = 10)
so.alv <- RunHarmony(so.alv, 'dataset')
so.alv <- so.alv %>%
    FindNeighbors(reduction = "harmony") %>%
    FindClusters(resolution = 1) 
so.alv <- RunUMAP(so.alv, dims = 1:10, reduction = 'harmony')
so.alv[["percent.mt"]] <- PercentageFeatureSet(so.alv, pattern = "^mt-")

so.alv <- FindClusters(so.alv, resolution = 1.5)

save_path_mt = '/home/sjcho/datas/public_data/mouse_lung_alveolar/outs/mt15/alveolar_only_'
p <- fp_sjcho(so.alv, features = c('Epcam', 'Pecam1', 'Ptprc', 'Pdgfra'))
ggsave(p, file = paste0(save_path_mt, 'alveolar_merged_featureplot.png'), width = 10, height = 10)
p <- DimPlot(so.alv, group.by = 'dataset')
ggsave(p, file = paste0(save_path_mt, 'alveolar_merged_dataset.png'), width = 7, height = 6)
p <- DimPlot(so.alv, group.by = 'RNA_snn_res.1.5', label = T, repel = T)
ggsave(p, file = paste0(save_path_mt, 'seurat_clusters_1.5.png'), width = 7, height = 6)
p <- DimPlot(so.alv, group.by = 'sample', label = T, repel = T)
ggsave(p, file = paste0(save_path_mt, 'alveolar_merged_seurat_samples.png'), width = 7, height = 6)
p <- fp_sjcho(so.alv, features = c('Sftpc', 'Etv5', 'Abca3', 'Pdpn', 'Aqp5', 'Hopx', 'Foxj1', 'Scgb1a1', 'Scgb3a2', 'Krt5', 'Trp63', 'Cldn4', 'Sprr1a', 'Tnip3', 'Krt8'), ncol = 4)
ggsave(p, file = paste0(save_path_mt, 'alveolar_merged_epi_features.png'), width = 13, height = 15)
p <- fp_sjcho(so.alv, features = c('Sftpc', 'Etv5', 'Abca3', 
                                'Pdpn', 'Aqp5', 'Hopx',
                                'Tnip3', 'Cldn4', 'Sprr1a',
                                'Krt8', 'Cdkn1a', 'Ly6a'), ncol = 3)
ggsave(p, file = paste0(save_path_mt, 'alveolar_merged_epi_features2.png'), width = 11, height = 12)
p <- fp_sjcho(so.alv, features = c('Top2a', 'Mki67'), ncol = 2)
ggsave(p, file = paste0(save_path_mt, 'alveolar_merged_proliferation.png'), width = 12, height = 6)
p <- fp_sjcho(so.alv, features = c('Lcn2', 'Cxcl17', 'Lrg1', 'Itga7', 'Ptges',' Glrx'), order = T, ncol = 3)
ggsave(p, file = paste0(save_path_mt, 'alveolar_merged_activatedAT2.png'), width = 13, height = 9)

so.alv[['log10_UMI']] <- log10(so.alv[['nCount_RNA']][[1]])
p <- fp_sjcho(so.alv, features = c('log10_UMI', 'percent.mt', 'nFeature_RNA'), ncol = 2, order = T)
ggsave(p, file = paste0(save_path_mt, 'alveolar_merged_tech.png'), width = 10, height = 10)

### HM gene set
# p53, EMT, Krt8
h_gene_sets[, 'gs_name'] %>% table %>% names
so.alv = AddModuleScore(so.alv, features = subset(h_gene_sets, gs_name == 'HALLMARK_P53_PATHWAY')[, 'gene_symbol'], name = 'HM_p53')
so.alv = AddModuleScore(so.alv, features = subset(h_gene_sets, gs_name == 'HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION')[, 'gene_symbol'], name = 'HM_EMT')

HM_ms_df = data.frame('tmp' = rep(0, length(colnames(so.alv))))
rownames(HM_ms_df) = colnames(so.alv)
for (hm_gene_set in h_gene_sets[, 'gs_name'] %>% table %>% names) {
    so.alv = AddModuleScore(so.alv, features = subset(h_gene_sets, gs_name == hm_gene_set)[, 'gene_symbol'], name = 'tmp')
    HM_ms_df[[hm_gene_set]] = so.alv[['tmp1']][[1]]
}
HM_ms_df[, 'tmp'] = NULL
save(HM_ms_df, file = paste0(save_path, 'alveolar_HM_ms_df.rdata'))

p <- fp_sjcho(so.alv, features = c('HM_p531', 'HM_EMT1', 'Krt8'), ncol = 2, order = T)
ggsave(p, file = paste0(save_path_mt, 'p53_EMT_Krt8.png'), width = 10, height = 10)

### based on above HM gene set, we can annotate the cell type

so.alv[['annotation']] = 'AT2'
so.alv[['annotation']][Cells(so.alv)[so.alv[['RNA_snn_res.1']][[1]] %in% c(6)], ] = 'Krt8_trasitory_cells'
so.alv[['annotation']][Cells(so.alv)[so.alv[['RNA_snn_res.1']][[1]] %in% c(8, 15, 16)], ] = 'AT1'
so.alv[['annotation']][Cells(so.alv)[so.alv[['RNA_snn_res.1']][[1]] %in% c(11, 12)], ] = 'proliferating_AT2/Krt8'

so.alv[['annotation']] = factor(so.alv[['annotation']][[1]], levels = c('AT2', 'AT1', 'Krt8_trasitory_cells', 'proliferating_AT2/Krt8'))
p <- DimPlot(so.alv, group.by = 'annotation', label = F, cols = c('#670404', '#f38b2b', 'red', '#ffa3a3'))
ggsave(p, file = paste0(save_path_mt, 'annotation.png'), width = 8.5, height = 6)

saveRDS(so.alv, file = paste0(save_path_mt, 'seurat.rds'))
saveRDS(so, file = paste0(save_path, 'full_merged_seurat_mt15.rds'))

so.alv[['annotation_cluster']] <- paste0(so.alv[['annotation']][[1]], '_', so.alv[['RNA_snn_res.1']][[1]])
annotation_cluster_list <- split(rownames(so.alv[['annotation_cluster']]), so.alv[['annotation_cluster']])
cluster_list <- split(rownames(so.alv[['RNA_snn_res.1.5']]), so.alv[['RNA_snn_res.1.5']])

HM_ms_df_mean_by_annotation_cluster <- lapply(cluster_list, function(rows) {
 colMeans(HM_ms_df[rows, , drop=FALSE])
})
HM_ms_df_mean_by_annotation_cluster %>% as.data.frame -> HM_ms_df_mean_by_annotation_cluster

variances_HM = apply(HM_ms_df, 2, var)
variance_df_HM = variances_HM %>% as.vector %>% as.data.frame
colnames(variance_df_HM) = 'variance'
variance_df_HM$rank = rank(-variances_HM)
p <- ggplot(variance_df_HM, aes(x = variance)) + geom_histogram(bins = 100) + theme_classic()
ggsave(p, file = paste0(save_path, 'HM_ms_variance.png'), width = 6, height = 6)
p <- ggplot(variance_df_HM, aes(x = rank, y = variance)) + geom_point() + theme_classic()
ggsave(p, file = paste0(save_path, 'HM_ms_variance_dotplot.png'), width = 6, height = 6)

HM_ms_df_mean_by_annotation_cluster_subset = HM_ms_df_mean_by_annotation_cluster[variance_df_HM$rank < 26, ]

my_palette <- colorRampPalette(c("darkblue", "white", "darkred"))(500)
png(paste0(save_path, 'HM_ms_pheatmap_test.png'),width = 3500, height = 3500, res = 300)
pheatmap_result <- pheatmap(HM_ms_df_mean_by_annotation_cluster_subset,
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        show_rownames = TRUE,
        show_colnames = TRUE,
        fontsize = 15,
        display_numbers = FALSE,
        scale = "row",
        color = my_palette,
        cutree_cols = 5
        # breaks = seq(-3, 3, length.out = 100),
        )
dev.off()

so.alv[['annotation']] <- 'AT2'
so.alv[['annotation']][Cells(so.alv)[so.alv[['RNA_snn_res.1.5']][[1]] %in% c(8, 13)], ] = 'p53_high_Krt8_transitory'
so.alv[['annotation']][Cells(so.alv)[so.alv[['RNA_snn_res.1.5']][[1]] %in% c(10, 15, 21)], ] = 'p53_low_Krt8_transitory'
so.alv[['annotation']][Cells(so.alv)[so.alv[['RNA_snn_res.1.5']][[1]] %in% c(2, 4, 7, 12)], ] = 'primed_AT2'
so.alv[['annotation']][Cells(so.alv)[so.alv[['RNA_snn_res.1.5']][[1]] %in% c(20, 22)], ] = 'AT1'
so.alv[['annotation']][Cells(so.alv)[so.alv[['RNA_snn_res.1.5']][[1]] %in% c(14, 18, 19)], ] = 'proliferating_Alveolar_epithelial'
so.alv[['annotation']] <- factor(so.alv[['annotation']][[1]], levels = c(
    'AT2', 'primed_AT2', 'p53_high_Krt8_transitory', 'p53_low_Krt8_transitory', 'AT1', 'proliferating_Alveolar_epithelial'
))
alv_cols = c('#c93a65', '#e29aad', '#440f0c', '#c18d81', '#ce1400', '#191412')
p <- DimPlot(so.alv, group.by = 'annotation', label = F, cols = alv_cols, raster = F)
ggsave(p, file = paste0(save_path_mt, 'heirachial_clustering.png'), width = 9, height = 6)
ggplot2pptx(p, 9, 6, paste0(save_path_mt, 'heirachial_clustering.pptx'))
### GO BP
bp_gene_sets[, 'gs_name'] %>% table %>% names %>% head

bp_ms_df = data.frame('tmp' = rep(0, length(colnames(so.alv))))
rownames(bp_ms_df) = colnames(so.alv)

bp_gene_sets_names = names(table(bp_gene_sets[, 'gs_name']))
bp_ms_df = sapply(seq_along(bp_gene_sets_names), function(i) {
  print(i)
  bp_gene_set = bp_gene_sets_names[i]
  genes = subset(bp_gene_sets, gs_name == bp_gene_set)$gene_symbol
  if(length(intersect(rownames(so.alv), genes)) >= 10) {
    so.alv = AddModuleScore(so.alv, features = list(genes), name = 'tmp')
    return(so.alv[['tmp1']][[1]])
  }
  return(rep(-100, ncol(so.alv)))
})
save(bp_ms_df, file = paste0(save_path, 'alveolar_bp_ms_df.rdata'))

rownames(bp_ms_df) = colnames(so.alv)
colnames(bp_ms_df) = bp_gene_sets_names
bp_ms_df = bp_ms_df[, bp_ms_df[1, ] != -100]
# find top300 variance
variances = apply(bp_ms_df, 2, var)
variance_df = variances %>% as.vector %>% as.data.frame
colnames(variance_df) = 'variance'
p <- ggplot(variance_df, aes(x = variance)) + geom_histogram(bins = 150) + theme_classic()
p <- p + xlim(0, 0.03)
ggsave(p, file = paste0(save_path, 'bp_ms_variance.png'), width = 6, height = 6)
bp_ms_df = bp_ms_df[, order(variances, decreasing = T)]
bp_ms_df_top500 = bp_ms_df[, 1:50]
bp_ms_df_over002 = bp_ms_df[, variances > 0.01]

so.alv[['annotation_1.5']] = paste0(so.alv[['annotation']][[1]], '_', so.alv[['RNA_snn_res.1.5']][[1]])
annotation_cluster_list <- split(colnames(so.alv), so.alv[['annotation_1.5']])
bp_ms_df_mean_by_annotation_cluster <- lapply(annotation_cluster_list, function(rows) {
 colMeans(bp_ms_df_over002[rows, ])
})
bp_ms_df_mean_by_annotation_cluster %>% as.data.frame -> bp_ms_df_mean_by_annotation_cluster

my_palette <- colorRampPalette(c("darkblue", "white", "darkred"))(500)
png(paste0(save_path, 'bp_ms_pheatmap.png'),width = 3500, height = 5000, res = 300)
pheatmap_result <- pheatmap(bp_ms_df_mean_by_annotation_cluster,
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        show_rownames = FALSE,
        show_colnames = TRUE,
        fontsize = 15,
        display_numbers = FALSE,
        scale = "row",
        color = my_palette
        # breaks = seq(-3, 3, length.out = 100),
        )
dev.off()

########## by HVG
hvgs <- VariableFeatures(so.alv)
hvgs_exp <- so.alv@assays$RNA$data[hvgs, ]
hvgs_df_mean_by_annotation_cluster <- lapply(annotation_cluster_list, function(rows) {
 colMeans(t(hvgs_exp)[rows, ])
})
hvgs_df_mean_by_annotation_cluster %>% as.data.frame -> hvgs_df_mean_by_annotation_cluster

my_palette <- colorRampPalette(c("darkblue", "white", "darkred"))(500)
png(paste0(save_path, 'hvg_pheatmap.png'),width = 3500, height = 5000, res = 300)
pheatmap_result <- pheatmap(hvgs_df_mean_by_annotation_cluster,
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        show_rownames = FALSE,
        show_colnames = TRUE,
        fontsize = 15,
        display_numbers = FALSE,
        scale = "row",
        color = my_palette
        # breaks = seq(-3, 3, length.out = 100),
        )
dev.off()


# GOBP_PEPTIDE_ANTIGEN_ASSEMBLY_WITH_MHC_CLASS_II_PROTEIN_COMPLEX
# GOBP_MHC_CLASS_II_BIOSYNTHETIC_PROCESS
MHC2_genesets = c('GOBP_PEPTIDE_ANTIGEN_ASSEMBLY_WITH_MHC_CLASS_II_PROTEIN_COMPLEX', 'GOBP_MHC_CLASS_II_BIOSYNTHETIC_PROCESS', 
'GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_EXOGENOUS_PEPTIDE_ANTIGEN_VIA_MHC_CLASS_II', 'GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_PEPTIDE_OR_POLYSACCHARIDE_ANTIGEN_VIA_MHC_CLASS_II')

for (geneset_name in MHC2_genesets) {
    geneset = subset(bp_gene_sets, gs_name == geneset_name)$gene_symbol
    so.alv = AddModuleScore(so.alv, features = list(geneset), name = geneset_name)
    # p <- VlnPlot(so.alv, features = paste0(geneset_name, 1), group.by = 'annotation', pt.size = 0, cols = alv_cols)
    # ggsave(p, file = paste0(save_path_mt, geneset_name, '_vilolin.png'), width = 10, height = 6)
}
p <- fp_sjcho(so.alv, features = paste0(MHC2_genesets, 1), ncol = 2, order = T)
ggsave(p, file = paste0(save_path_mt, 'MHC2_genesets.png'), width = 10, height = 10)
p <- fp_sjcho(so.alv, features = paste0(MHC2_genesets, 1), ncol = 2, order = T, min.cutoff = 'q50')
ggsave(p, file = paste0(save_path_mt, 'MHC2_genesets.png'), width = 10, height = 10)

# convert it to anndata
so.alv[['annotation_vector']] <- as.vector(so.alv[['annotation']][[1]])
library(sceasy)
library(reticulate)
use_condaenv('project_lung_exercise_R')
loompy <- reticulate::import('loompy')
so.alv[["RNA"]] <- as(so.alv[["RNA"]], "Assay")
sceasy::convertFormat(so.alv, from="seurat", to="anndata", drop_single_values=FALSE,
                       outFile= paste0(save_path_mt,'_seurat.h5ad'))



### Get airway & club signatures
so.ref <- readRDS('/home/sjcho/datas/reference_atlas/mouse_lung/refined_reference.rds')
Idents(so.ref) <- so.ref[['annotation.h']][[1]]
club_marker = FindMarkers(so.ref, ident.1 = c('Club cell'), ident.2 = c('AT2', 'AT1', 'Ciliated cell', 'PNEC', 'Mucous/Goblet cell', 'Basal cell'), only.pos = T, log2FC.threshold = 0.5)
airway_marker = FindMarkers(so.ref, ident.1 = c('Club cell', 'Ciliated cell', 'PNEC', 'Mucous/Goblet cell', 'Basal cell'), ident.2 = c('AT2', 'AT1'), only.pos = T, log2FC.threshold = 0.5)
save(club_marker, file = paste0('/home/sjcho/datas/reference_atlas/mouse_lung/', 'club_marker_log2FC05.rdata'))
save(airway_marker, file = paste0('/home/sjcho/datas/reference_atlas/mouse_lung/', 'airway_marker_log2FC05.rdata'))

### Get markers
Idents(so.alv) <- so.alv[['annotation']][[1]]
for (celltype in so.alv[['annotation']] %>% table %>% names) {
    markers = FindMarkers(so.alv, ident.1 = celltype, only.pos = T, log2FC.threshold = 0.5)
    save(markers, file = paste0(save_path_mt, celltype, '_markers_log2FC05.rdata'))
}

saveRDS(so, file = paste0(save_path, 'full_merged_seurat_mt15.rds'))
saveRDS(so.alv, file = paste0(save_path_mt, 'alv_only.rds'))

