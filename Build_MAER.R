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

so.merged <- so

save_path_mt10 = '/home/sjcho/datas/reference_atlas/mouse_lung_alveolar/outs/1.3.log_normalization_clean/mt10/'
so.merged_10 <- subset(so.merged, cells = Cells(so.merged)[so.merged[['percent.mt']] < 10])
so.merged_10 <- log_normalize(so.merged_10, save_path_mt10, 'merged_mt10', pcs = 20, nfeatures = 2000)
plot_alveolar_epithelium_features(so.merged_10, save_path_mt10, 'merged_mt10')
plot_lineage_features(so.merged_10, save_path_mt10, 'merged_mt10')
p <- fp_sjcho(so.merged_10, features = c('log10_UMI', 'percent.mt', 'nFeature_RNA'))
ggsave(p, file = paste0(save_path_mt10, 'log10_UMI_percent_mt_nFeature_RNA.png'), width = 10, height = 10)

# annotation
# 18 : stromal
# 19 : endothelial
# 5, 6, 12, 16, 20, 21 : immune
# 1, 4, 9, 17, 24 : ciliated
# 23, 14 : club
# others : alveolar epithelial
so.merged_10[['annotation']] = 'alveolar_epithelial'
so.merged_10[['annotation']][Cells(so.merged_10)[so.merged_10[['seurat_clusters']][[1]] %in% c(18)], ] = 'stromal'
so.merged_10[['annotation']][Cells(so.merged_10)[so.merged_10[['seurat_clusters']][[1]] %in% c(19)], ] = 'endothelial'
so.merged_10[['annotation']][Cells(so.merged_10)[so.merged_10[['seurat_clusters']][[1]] %in% c(5, 6, 12, 16, 20, 21)], ] = 'immune'
so.merged_10[['annotation']][Cells(so.merged_10)[so.merged_10[['seurat_clusters']][[1]] %in% c(1, 4, 9, 17, 24)], ] = 'ciliated'
so.merged_10[['annotation']][Cells(so.merged_10)[so.merged_10[['seurat_clusters']][[1]] %in% c(23, 14)], ] = 'club'

so.alv_mt10 <- subset(so.merged_10, cells = Cells(so.merged_10)[so.merged_10[['annotation']] == 'alveolar_epithelial'])
save_path_alv_mt10 = '/home/sjcho/datas/reference_atlas/mouse_lung_alveolar/outs/1.3.log_normalization_clean/alveolar_mt10/'
so.alv_mt10 <- log_normalize(so.alv_mt10, save_path_alv_mt10, 'alveolar_mt10', pcs = 10, nfeatures = 1000)
plot_alveolar_epithelium_features(so.alv_mt10, save_path_alv_mt10, 'alveolar_mt10')
plot_lineage_features(so.alv_mt10, save_path_alv_mt10, 'alveolar_mt10')
p <- fp_sjcho(so.alv_mt10, features = c('log10_UMI', 'percent.mt', 'nFeature_RNA'))
ggsave(p, file = paste0(save_path_alv_mt10, 'log10_UMI_percent_mt_nFeature_RNA.png'), width = 10, height = 10)

immune_cells = Cells(so.alv_mt10)[so.alv_mt10[['seurat_clusters']][[1]] %in% c('10')]

so.alv_mt10_del10 <- subset(so.alv_mt10, cells = Cells(so.alv_mt10)[so.alv_mt10[['seurat_clusters']] != 10])
save_path_alv_mt10_del10 = '/home/sjcho/datas/reference_atlas/mouse_lung_alveolar/outs/1.3.log_normalization_clean/alveolar_mt10_wo_immune/'
so.alv_mt10_del10 <- log_normalize(so.alv_mt10_del10, save_path_alv_mt10_del10, 'alveolar_mt10_wo_immune', pcs = 10, nfeatures = 1000)
plot_alveolar_epithelium_features(so.alv_mt10_del10, save_path_alv_mt10_del10, 'alveolar_mt10_wo_immune')
plot_lineage_features(so.alv_mt10_del10, save_path_alv_mt10_del10, 'alveolar_mt10_wo_immune')
p <- fp_sjcho(so.alv_mt10_del10, features = c('log10_UMI', 'percent.mt', 'nFeature_RNA'))
ggsave(p, file = paste0(save_path_alv_mt10_del10, 'log10_UMI_percent_mt_nFeature_RNA.png'), width = 10, height = 10)

immune_cells = Cells(so.alv_mt10_del10)[so.alv_mt10_del10[['seurat_clusters']][[1]] %in% c('17')]

so.alv_mt10_del1017 <- subset(so.alv_mt10_del10, cells = Cells(so.alv_mt10_del10)[so.alv_mt10_del10[['seurat_clusters']] != 17])
save_path_alv_mt10_del1017 = '/home/sjcho/datas/reference_atlas/mouse_lung_alveolar/outs/1.3.log_normalization_clean/alveolar_mt10_wo_immune_twice/'
so.alv_mt10_del1017 <- log_normalize(so.alv_mt10_del1017, save_path_alv_mt10_del1017, 'alveolar_mt10_wo_immune', pcs = 7, nfeatures = 1000)
plot_alveolar_epithelium_features(so.alv_mt10_del1017, save_path_alv_mt10_del1017, 'alveolar_mt10_wo_immune')
plot_lineage_features(so.alv_mt10_del1017, save_path_alv_mt10_del1017, 'alveolar_mt10_wo_immune')
p <- fp_sjcho(so.alv_mt10_del1017, features = c('log10_UMI', 'percent.mt', 'nFeature_RNA'))
ggsave(p, file = paste0(save_path_alv_mt10_del1017, 'log10_UMI_percent_mt_nFeature_RNA.png'), width = 10, height = 10)

p <- fp_sjcho(so.alv_mt10_del1017, features = c('H2-Ab1'), ncol = 1)
ggsave(p, file = paste0(save_path_alv_mt10_del1017, 'MHC2.png'), width = 3.5, height = 3)

so.merged.3 <- so.alv_mt10_del1017

meta.data <- fread('/home/sjcho/datas/reference_atlas/mouse_lung_alveolar/GSE262927/GSE262927_CellMetaData.csv')
meta.data <- as.data.frame(meta.data)
rownames(meta.data) <- paste0(meta.data$cb, '_', meta.data$orig.ident)
meta.data$orig.ident %>% table %>% names -> seq_to_ident

key <- gsub(".*?(\\d{3}).*", "\\1", seq_to_ident)
seq_to_ident <- setNames(seq_to_ident, key)

h5_files = '/home/sjcho/datas/reference_atlas/mouse_lung_alveolar/GSE262927/h5_files'
file_list = list.files(h5_files, full.names = T)

h5s <- list()
for (file in file_list) {
    h5 <- Read10X_h5(file)
    key <- sub(".*?(\\d{3})\\.h5$", "\\1", file)
    if (as.character(key) %in% names(seq_to_ident)) {
        ident <- seq_to_ident[[as.character(key)]]
        colnames(h5) <- paste0(colnames(h5), '_', ident)

    } else {
        ident <- NULL
        colnames(h5) <- paste0(colnames(h5), '_', key)
    }
    
    h5s[[key]] <- h5
}

merged_count <- do.call(cbind, h5s)
merged_count_used <- merged_count[, rownames(meta.data)] 
# 123189 cells left

so.tmp <- CreateSeuratObject(counts = merged_count_used, project = 'GSE262927', meta.data = meta.data, min.cells = 0, min.features = 0)
so.GSE262927 <- subset(so.tmp, cells = Cells(so)[so$celltype %in% c('AT1', 'AT1_AT2', 'AT2')])


common.genes = intersect(rownames(so.merged.3), rownames(so.GSE262927))
merged.count <- cbind(so.merged.3@assays$RNA$counts[common.genes,], so.GSE262927@assays$RNA$counts[common.genes,])

so.merged <- CreateSeuratObject(counts = merged.count, project = 'alveolar_merged', min.cells = 0, min.features = 0)
so.merged@meta.data[Cells(so.merged.3), c('sex', 'sample', 'injury_type', 'dpi', 'dataset', 'annotation')] = lapply(so.merged.3@meta.data[, c('gender', 'sample', 'injury_type', 'dpi', 'dataset', 'annotation')], as.character)
so.merged@meta.data[Cells(so.GSE262927), c('sex', 'sample', 'dpi', 'annotation')] = lapply(so.GSE262927@meta.data[, c('sex','experimental_group', 'sacrifice_day', 'subtype')], as.character)
# dataset
so.merged@meta.data[Cells(so.GSE262927), 'dataset'] = 'GSE262927'
so.merged@meta.data[Cells(so.GSE262927), 'injury_type'] = 'H1N1'
so.merged@meta.data[Cells(so.GSE262927), 'original_annotation'] = so.GSE262927@meta.data[, 'subtype']

plot_alveolar_epithelium_features <- function(so, save_path, data_name, reduction = 'umap') {
    p <- fp_sjcho(so, features = c('Epcam', 'Pecam1', 'Ptprc', 'Pdgfra'), reduction = reduction)
    ggsave(p, file = paste0(save_path, data_name, '_featureplot.png'), width = 10, height = 10)
    p <- DimPlot(so, group.by = 'dataset', label = T, repel = T, reduction = reduction)
    ggsave(p, file = paste0(save_path, data_name, '_dataset.png'), width = 7, height = 6)
    p <- DimPlot(so, group.by = 'seurat_clusters', label = T, repel = T, reduction = reduction)
    ggsave(p, file = paste0(save_path, data_name, '_seurat_clusters.png'), width = 7, height = 6)
    p <- DimPlot(so, group.by = 'sample', label = T, repel = T, reduction = reduction)
    ggsave(p, file = paste0(save_path, data_name, '_seurat_samples.png'), width = 7, height = 6)
    p <- fp_sjcho(so, features = c('Sftpc', 'Etv5', 'Abca3', 'Pdpn',
                                    'Aqp5', 'Foxj1', 'Scgb1a1',
                                    'Scgb3a2', 'Krt5', 'Trp63', 'Cldn4', 
                                    'Sprr1a', 'Tnip3', 'Krt8', 'Mki67', 'Top2a'), ncol = 4, reduction = reduction)
    ggsave(p, file = paste0(save_path, data_name, '_epi_features.png'), width = 13, height = 15)
    p <- fp_sjcho(so, features = c('Sftpc', 'Etv5', 'Foxj1', 
                                    'Pdpn', 'Aqp5', 'Hopx',
                                    'Tnip3', 'Cldn4', 'Krt8',
                                    'Ptprc', 'Pecam1', 'Pdgfra'), ncol = 3, reduction = reduction)
    ggsave(p, file = paste0(save_path, data_name, '_epi_features2.png'), width = 11, height = 12)
}
plot_lineage_features <- function(so, save_path, data_name, reduction = reduction) {
    p <- fp_sjcho(so, features = c('Sftpc', 'Etv5', 'Foxj1', 
                                    'Pdpn', 'Aqp5', 'Hopx',
                                    'Tnip3', 'Cldn4', 'Krt8',
                                    'Ptprc', 'Pecam1', 'Pdgfra'), ncol = 3, reduction = reduction)
    ggsave(p, file = paste0(save_path, data_name, '_lineage_features.png'), width = 11, height = 12)
}

so.merged <- log_normalize(so.merged, save_path, 'w_GSE262927', 6, nfeatures = 1000)

########### visulalize
 
study_cols = c('forestgreen', 'skyblue', 'brown', 'orange')

p <- DimPlot(so.merged, group.by = 'seurat_clusters', split.by = 'dataset', ncol = 4, pt.size = 0.75) & NoAxes()
ggsave(p, file = paste0(save_path, 'w_GSE262927_dataset.png'), width = 17, height = 6)
ggsave(p, file = paste0(save_path, 'w_GSE262927_dataset.pdf'), width = 17, height = 6)

plot_alveolar_epithelium_features(so.merged, save_path, 'w_GSE262927')
plot_lineage_features(so.merged, save_path, 'w_GSE262927')
p <- fp_sjcho(so.merged, features = c('Epcam', 'Cdh1', 'Pecam1', 'Ptprb', 'Pdgfra', 'Col1a1', 'Cdh5', 'Vim', 'Ptprc'), ncol = 3, order = T)
ggsave(p, file = paste0(save_path, 'w_GSE262927_other_lineage.png'), width = 10, height = 9)
p <- fp_sjcho(so.merged, features = c('Ifitm3', 'Ly6a', 'Kit', 'Cd34'), ncol = 2, order = T)
ggsave(p, file = paste0(save_path, 'w_GSE262927_doublet_specific.png'), width = 10, height = 9)

p <- DimPlot(so.merged, group.by = 'dpi', label = T, repel = T, split.by = 'dataset', ncol = 4) & NoAxes()
ggsave(p, file = paste0(save_path, 'w_GSE262927_dpi.png'), width = 13, height = 3)

so.merged@meta.data[Cells(so.merged)[so.merged$RNA_snn_res.1 %in% c('15', '11')], 'dpi']

# delete study specific doublet
so.merged.wo.d <- subset(so.merged, cells = Cells(so.merged)[so.merged$RNA_snn_res.1 %in% c('15', '11')], invert = T)
so.merged.wo.d <- log_normalize(so.merged.wo.d, save_path, 'w_GSE262927_wo_doublet', 6, nfeatures = 1000)

plot_alveolar_epithelium_features(so.merged.wo.d, save_path, 'w_GSE262927_wo_doublet')
plot_lineage_features(so.merged.wo.d, save_path, 'w_GSE262927_wo_doublet')

p <- fp_sjcho(so.merged.wo.d, features = c('Epcam', 'Cdh1', 'Pecam1', 'Ptprb', 'Pdgfra', 'Col1a1', 'Cdh5', 'Vim', 'Ptprc'), ncol = 3, order = T)
ggsave(p, file = paste0(save_path, 'w_GSE262927_wo_doublet_other_lineage.png'), width = 10, height = 9)

p <- fp_sjcho(so.merged.wo.d, features = c('Top2a', 'Mki67'), ncol = 2, order = T)
ggsave(p, file = paste0(save_path, 'w_GSE262927_wo_doublet_proliferating.png'), width = 6, height = 3)

resolution = 1
so.merged.wo.d <- FindClusters(so.merged.wo.d, resolution = resolution)

### annotate and visualize
### annotate based on canonical markers
so <- so.merged.wo.d

so$annotation = as.vector(so$seurat_clusters)
so$annotation[so$seurat_clusters %in% c(2, 11, 3, 4, 10, 1, 0, 5)] = 'AT2'
so$annotation[so$seurat_clusters %in% c(9)] = 'AT2.IFN'
so$annotation[so$seurat_clusters %in% c(12, 13, 15)] = 'AT2.prolif.'
so$annotation[so$seurat_clusters %in% c(6, 16)] = 'tAT2'
so$annotation[so$seurat_clusters %in% c(7, 8, 14, 17)] = 'AT1'
so$tmp1 = NULL
saveRDS(so, file = paste0(save_path, 'merged_w_GSE262927_full.rds'))

### AT2 : 2, 11, 3, 4, 10, 1, 0, 5
### interferon.sti.alv.epi : 9
### cycling.alv.epi : 12, 13, 15
### tAT2 : 6, 16
### AT1 : 7, 8, 14, 17

# AT2, interferon.sti.alv.epi, cycling.alv.epi, tAT2, AT1
so$clusters <- paste0('c', as.vector(so$seurat_clusters))
so$clusters <- factor(so$clusters, levels = paste0('c', c(c(2, 11, 3, 4, 0, 5, 10, 1), c(9), c(12, 13, 15), c(6, 16), c(7, 8, 14, 17))))

cols_clusters = c('#538dba', '#32869e', '#002c5b', '#846ded', '#22007a', '#0a2466', '#a4dafc', '#8b95ca', '#000000', '#bdfccf', '#163d2e', '#11c195', '#f9e79f', '#b9a708', '#eaa85d', '#b2460c', '#be9464', '#ff7118')
names(cols_clusters) = levels(so$clusters)

### extd.fig6.x.cluster & dataset

p <- DimPlot(so, group.by = 'clusters', cols = cols_clusters, pt.size = 1) + NoAxes() + NoLegend() + theme(plot.title = element_text(size = 0))
ggsave(p, filename = paste0(save_path, 'sup1.A.cluster.png'), width = 6, height = 6)
ggsave(p, filename = paste0(save_path, 'sup1.A.cluster.pdf'), width = 6, height = 6)
ggplot2pptx(p, 6, 6, paste0(save_path, 'sup1.A.cluster.pptx'))

leg <- get_legend(DimPlot(so, group.by = 'clusters', cols = cols_clusters))
ggplot2pptx(as_ggplot(leg), 1, 6, paste0(save_path, 'sup1.A.cluster.legend.pptx'))

p <- DimPlot(so, group.by = 'clusters', split.by = 'dataset', ncol = 4, cols = cols_clusters, pt.size = 0.75) + NoAxes() + NoLegend() + theme(plot.title = element_text(size = 0)) & theme(strip.text.x = element_text(size = 0, face = "bold"))
ggsave(p, filename = paste0(save_path, 'sup1.B.dataset_cluster.png'), width = 14, height = 3.5)
ggsave(p, filename = paste0(save_path, 'sup1.B.dataset_cluster.pdf'), width = 14, height = 3.5)
ggplot2pptx(p, 14, 3.5, paste0(save_path, 'sup1.B.dataset_cluster.pptx'))

p <- DimPlot(so, group.by = 'clusters', split.by = 'dataset', ncol = 4, cols = cols_clusters, pt.size = 0.75) + NoAxes() + NoLegend()
ggsave(p, filename = paste0(save_path, 'sup1.B.dataset_cluster_with_title.png'), width = 14, height = 3.5)

### extd.fig6.xx.heatmap

h_gene_sets = msigdbr(species = "mouse", category = 'H')

HM_ms_df = data.frame('tmp' = rep(0, length(colnames(so))))
rownames(HM_ms_df) = colnames(so)
for (hm_gene_set in h_gene_sets[, 'gs_name'] %>% table %>% names) {
    so = AddModuleScore(so, features = as.list(subset(h_gene_sets, gs_name == hm_gene_set)[, 'gene_symbol']), name = 'tmp')
    HM_ms_df[[hm_gene_set]] = so[['tmp1']][[1]]
    print(hm_gene_set)
}
HM_ms_df[, 'tmp'] = NULL

resolution = 1
cluster_list <- split(rownames(so[[paste0('RNA_snn_res.', resolution)]]), so[[paste0('RNA_snn_res.', resolution)]])
HM_ms_df_mean_by_annotation_cluster <- lapply(cluster_list, function(rows) {
    colMeans(HM_ms_df[rows, , drop=FALSE])
})

HM_ms_df_mean_by_annotation_cluster %>% as.data.frame -> HM_ms_df_mean_by_annotation_cluster
colnames(HM_ms_df_mean_by_annotation_cluster) = gsub('X', 'c', colnames(HM_ms_df_mean_by_annotation_cluster))

rownames(HM_ms_df_mean_by_annotation_cluster) = gsub('HALLMARK_', '', rownames(HM_ms_df_mean_by_annotation_cluster))

write.csv(HM_ms_df_mean_by_annotation_cluster, file = paste0(save_path, 'ext.6.e.HM_ms_df_mean_by_annotation_cluster.csv'))

annotation_col = data.frame(
    CellType = c(
        rep('AT2', 8),  # c2~c5
        rep('tAT2', 2),                    
        rep('AT1', 4),                 # c7~c8
        'IFN.response.alv.epi',       # c9
        rep('cycling.alv.epi',3)  # c12~c13
    )
)
rownames(annotation_col) = c(paste0('c', c(c(2, 11, 3, 4, 10, 1, 0, 5), c(6, 16), c(7, 8, 14, 17), c(9), c(12, 13, 15))))

ann_colors = list(
    CellType = c(AT2 = "#8192ef", tAT2 = "#fab76d", AT1 = "#e74c3c", IFN.response.alv.epi = "#000000", cycling.alv.epi = "#8bd690")
)

library(RColorBrewer)
my_palette <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(50))
pheatmap(HM_ms_df_mean_by_annotation_cluster,
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        show_rownames = TRUE,
        show_colnames = TRUE,
        fontsize = 15,
        display_numbers = FALSE,
        scale = "row",
        color = my_palette,
        breaks = seq(-3, 3, length.out = length(my_palette) + 1)
        # annotation_col = annotation_col,
        # annotation_colors = ann_colors
        # breaks = seq(-3, 3, length.out = 100),
        ) %>% as.ggplot -> p
ggsave(p, filename = paste0(save_path, 'sup1.C.HM_ms_pheatmap_1.', '_wo_HM.png'), width = 12, height = 11)
ggsave(p, filename = paste0(save_path, 'sup1.C.HM_ms_pheatmap_1.', '_wo_HM.pdf'), width = 12, height = 11)
ggplot2pptx(p, 12, 11, paste0(save_path, 'sup1.C.HM_ms_pheatmap_1.', '_wo_HM.pptx'))

### dotplot
library(RColorBrewer)
so$clusters <- factor(so$clusters, levels = paste0('c', c(c(0, 1, 2, 3, 4, 5, 10, 11), c(6, 16), c(7, 8, 14, 17), c(12, 13, 15), c(9))))
p <- DotPlot(so, features = c('Etv5', 'Abca3', 'Napsa', 'Cldn4', 'Tnip3', 'Krt8', 'Pdpn', 'Aqp5', 'Hopx', 'Top2a', 'Mki67', 'Ifi44', 'Iigp1'), group.by = 'clusters')
p <- p + theme_classic()
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p <- p + scale_color_gradientn(colors = rev(brewer.pal(9, "RdBu")))
ggsave(p, file = paste0(save_path, 'sup1.D.w_GSE262927_w_pro_clusters_dotplot.png'), width = 7, height = 5)
ggsave(p, file = paste0(save_path, 'sup1.D.w_GSE262927_w_pro_clusters_dotplot.pdf'), width = 7, height = 5)
ggplot2pptx(p, 6, 5, paste0(save_path, 'sup1.D.w_GSE262927_w_pro_clusters_dotplot.pptx'))

### celltypist
library(sceasy)
library(reticulate)
use_condaenv('project_lung_exercise_R')
loompy <- reticulate::import('loompy')
so[["RNA"]] <- as(so[["RNA"]], "Assay")
sceasy::convertFormat(so, from="seurat", to="anndata", drop_single_values=FALSE,
                       outFile= paste0('/home/sjcho/datas/reference_atlas/mouse_lung_alveolar/outs/20250217_public_figures/','merged_w_GSE262927_full.h5ad'))

### MHC2 expression level
gobp = msigdbr(species = "mouse", category = 'C5')

MHC2_geneset = subset(gobp, gs_name == 'GOCC_MHC_CLASS_II_PROTEIN_COMPLEX')$gene_symbol %>% as.vector
so <- AddModuleScore(so, features = list(MHC2_geneset), name = 'MHC2_geneset')

p <- FeaturePlot(so, features = c('H2-Ab1', 'MHC2_geneset1'), order = T, pt.size = 0.25) & NoAxes()
p <- p & scale_colour_gradientn(colours = rev(c("#2b0247", "#6800ad","gray90")))
p <- p & theme(plot.title = element_text(size = 0))
ggsave(p, filename = paste0(save_path, 'extd.fig6.X.MHC2_geneset.png'), width = 8, height = 4)

p <- StackedVlnPlot(so, features = c('H2-Ab1', 'MHC2_geneset1'), group.by = 'clusters', cols = cols_clusters, pt.size = 0) 
p <- p & NoLegend() & geom_boxplot(fill = 'white', width = 0.2, outlier.size = 0.5)
p <- p & theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(p, filename = paste0(save_path, 'sup1.E.MHC2_geneset_vlnplot.png'), width = 9, height = 6)
ggsave(p, filename = paste0(save_path, 'sup1.E.MHC2_geneset_vlnplot.pdf'), width = 9, height = 6)
ggplot2pptx(p, 9, 6, paste0(save_path, 'sup1.E.MHC2_geneset_vlnplot.pptx'))

### UMAP split.by damage type
so$injury_type[so$dpi == '0'] = 'No Damage'
so$injury_type[so$injury_type == 'H1N1_PR8'] = 'H1N1'
so$injury_type = factor(so$injury_type, levels = c('No Damage', 'H1N1', 'bleomycin'))

p <- DimPlot(so, group.by = 'clusters', split.by = 'injury_type', cols = cols_clusters, pt.size = 0.5, ncol = 3) & NoAxes() & NoLegend() & theme(plot.title = element_text(size = 0))
p <- p & theme(strip.text.x = element_text(size = 0, face = "bold"))
ggsave(p, filename = paste0(save_path, 'sup1.B.injury_type.png'), width = 9, height = 3)
ggsave(p, filename = paste0(save_path, 'sup1.B.injury_type.pdf'), width = 9, height = 3)
ggplot2pptx(p, 6, 6, paste0(save_path, 'sup1.B.injury_type.pptx'))

### markers
so <- readRDS('/home/sjcho/datas/reference_atlas/mouse_lung_alveolar/outs/20250217_public_figures/merged_w_GSE262927_full.rds')
Idents(so) <- so$annotation

for (celltype in unique(so$annotation)) {
    marker <- FindMarkers(so, ident.1 = celltype, min.pct = 0, logfc.threshold = 0, test.use = 'MAST', latet.var = 'dataset')
    save(marker, file = paste0('/home/sjcho/datas/reference_atlas/mouse_lung_alveolar/outs/20250217_public_figures/markers/', celltype, '_marker.RData'))
} 

AT2_markers <- FindMarkers(so, ident.1 = 'AT2', min.pct = 0, logfc.threshold = 0, test.use = 'MAST', latet.var = c('dataset', 'injury_type'))
save(AT2_markers, file = paste0('/home/sjcho/datas/reference_atlas/mouse_lung_alveolar/outs/20250217_public_figures/markers/AT2_marker.RData'))

AT2_uninjuryed_markers <- FindMarkers(so, ident.1 = Cells(so)[so$annotation == 'AT2' & so$injury_type == 'No Damage'], 
                                        ident.2 = Cells(so)[so$annotation %in% c('AT1', 'tAT2') & so$injury_type == 'No Damage'], 
                                        min.pct = 0, logfc.threshold = 0, test.use = 'MAST', latet.var = c('dataset', 'injury_type'))
save(AT2_uninjuryed_markers, file = paste0('/home/sjcho/datas/reference_atlas/mouse_lung_alveolar/outs/20250217_public_figures/markers/AT2_uninjury_marker.RData'))


so <- AddModuleScore(so, features = list(subset(AT2_markers, avg_log2FC > 1 & p_val_adj < 0.05 & pct.1 > 0.3) %>% rownames), name = 'AT2_marker')
p <- fp_sjcho(so, features = c('AT2_marker1'), order = T, pt.size = 0.25, ncol = 1) & NoAxes()
# p <- p & scale_colour_gradientn(colours = rev(c("#2b0247", "#6800ad","gray90")))
p <- p & theme(plot.title = element_text(size = 0))
ggsave(p, filename = paste0(save_path, 'extd.fig6.X.AT2_marker.png'), width = 5, height = 4)

so$annotation_injury = paste0(so$annotation, '_', so$injury_type)
p <- VlnPlot(so, features = 'Il33', group.by = 'annotation_injury', pt.size = 0, cols = rep('darkgray', length(unique(so$annotation_injury)))) + geom_boxplot(fill = 'white', width = 0.2)
p <- p + NoLegend()
ggsave(p, filename = paste0(save_path, 'extd.fig6.X.Il33_vlnplot.png'), width = 8, height = 4)

save_path = '/home/sjcho/datas/reference_atlas/mouse_lung_alveolar/outs/20250217_public_figures/AT2_marker_test/'
for (gene in (subset(AT2_markers, avg_log2FC > 1 & p_val_adj < 0.05 & pct.1 > 0.3) %>% rownames)) {
    p <- VlnPlot(so, features = gene, group.by = 'annotation_injury', pt.size = 0, cols = rep('darkgray', length(unique(so$annotation_injury)))) + geom_boxplot(fill = 'white', width = 0.2)
    p <- p + NoLegend()
    ggsave(p, filename = paste0(save_path, 'extd.fig6.X.', gene, '_vlnplot.png'), width = 8, height = 4)
}
