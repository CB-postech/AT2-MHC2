library(harmony)

log_normalize <- function(so, save_path, data_name, pcs, nfeatures = 2000) {
    so <- NormalizeData(object = so)
    so <- FindVariableFeatures(object = so, nfeatures = nfeatures)
    so <- ScaleData(object = so)
    so <- RunPCA(object = so, npcs = 50)
    p <- ElbowPlot(so, ndims = 50, reduction = "pca")
    ggsave(p, file = paste0(save_path, data_name, '_elbowplot.png'), width = 5, height = 5)
    so <- RunHarmony(so, 'dataset')
    so <- so %>%
        FindNeighbors(reduction = "harmony") %>%
        FindClusters(resolution = 1) 
    so <- RunUMAP(so, dims = 1:pcs, reduction = 'harmony')
    return(so)
}

rename_files <- function(file_path) {
  file_name <- basename(file_path)
  
  if (str_detect(file_name, "_barcodes.tsv.gz$")) {
    new_name <- "barcodes.tsv.gz"
  } else if (str_detect(file_name, "_gene.tsv.gz$")) {
    new_name <- "features.tsv.gz"
  } else if (str_detect(file_name, ".mtx.gz$")) {
    new_name <- "matrix.mtx.gz"
  } else {
    return(NULL)
  }
  new_path <- file.path(dirname(file_path), new_name)
  file.rename(file_path, new_path)
}

run_dropletUtils <- function(file_path, save_path, sample_info) {
    tmp_count <- Read10X(dir_path)
    colnames(tmp_count) <- paste0(sample_info, '_', colnames(tmp_count))

    rawsce <- SingleCellExperiment(assays = list(counts = tmp_count))
    br.out <- barcodeRanks(counts(rawsce))
    
    png(paste0(save_path, 'DropletUtils_', sample_info, '.png'))

    plot(br.out$rank, br.out$total, log = "xy", xlab = "Rank", ylab = "Total", main = sample_info)

    o <- order(br.out$rank)
    lines(br.out$rank[o], br.out$fitted[o], col = "red")

    e.out <- emptyDrops(counts(rawsce))  ## Cells that have UMI counts lower than 100 are empty cells.
    table(Sig=e.out$FDR <= 0.05, Limited=e.out$Limited)
    is.cell <- e.out$FDR <= 0.05

    print(sum(is.cell, na.rm=TRUE))
    print(table(br.out$rank == sum(is.cell, na.rm=TRUE)))

    abline(h=min(br.out$fitted[o], na.rm=TRUE), col="red", lty=2)
    abline(h=metadata(br.out)$knee, col="dodgerblue", lty=2)
    abline(h=metadata(br.out)$inflection, col="forestgreen", lty=2)
    legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen", "red"), legend=c("knee", "inflection", "FDR_0.05"))
    dev.off()

    # colnames(rawsce) = colData(rawsce)$Barcode
    rawsce <- rawsce[,which(e.out$FDR <= 0.05)]
    return(rawsce)
}

celltype_order = c('AT2', 'pAT2', 'Transitory_1', 'Transitory_2', 'AT1_1', 'AT1_2')
alv_cols = c('#c93a65', '#e29aad','#a50c16', '#6f1566', 'red', '#ffdbdb')
names(alv_cols) = celltype_order

add_damage_type <- function(so) {
  celltype_order = c('AT2', 'pAT2', 'Transitory_1', 'Transitory_2', 'AT1_1', 'AT1_2')
  sample_names = so[['sample']] %>% unique %>% unlist %>% as.vector
  damage_type = c('NoDamage', 'NoDamage', 'Bleomycin', 'Bleomycin', 'Bleomycin', 'Bleomycin', 'NoDamage', 'Bleomycin', 'Bleomycin', 'NoDamage', 'NoDamage', 'H1N1', 'H1N1', 'H1N1', 'H1N1', 'NoDamage', 'NoDamage')
  names(damage_type) = sample_names
  so[['damage_type']] <- damage_type[so[['sample']][[1]]] %>% as.vector
  so[['damage_type']] <- factor(so[['damage_type']][[1]], levels = c('NoDamage', 'Bleomycin', 'H1N1'))
  so[['damage_annotation']] <- paste0(so[['damage_type']][[1]], '_', so[['annotation']][[1]])
  level_order = c()
  for(celltype in celltype_order) {
      for(damage in c('NoDamage', 'Bleomycin', 'H1N1')) {
          level_order = c(level_order, paste0(damage, '_', celltype))
      }
  }
  so@meta.data[, 'damage_annotation'] <- factor(so[['damage_annotation']][[1]], 
                                      levels = level_order)
  return(so)
}

# stacked violin plot from https://divingintogeneticsandgenomics.rbind.io/post/stacked-violin-plot-for-visualizing-single-cell-data-in-seurat/
library(Seurat)
library(patchwork)
library(ggplot2)

## remove the x-axis text and tick
## plot.margin to adjust the white space between each plot.
## ... pass any arguments to VlnPlot in Seurat
modify_vlnplot<- function(obj, 
                          feature, 
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  p<- VlnPlot(obj, features = feature, pt.size = pt.size, ... )  + 
    xlab("") + ylab(feature) + ggtitle("") + 
    theme(legend.position = "none", 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.title.y = element_text(size = rel(1), angle = 0), 
          axis.text.y = element_text(size = rel(1)), 
          plot.margin = plot.margin ) 
  return(p)
}

## extract the max value of the y axis
extract_max<- function(p){
  ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}


## main function
StackedVlnPlot<- function(obj, features,
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  
  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  
  # Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(), axis.ticks.x = element_line())
  
  # change the y-axis tick to only max value 
  ymaxs<- purrr::map_dbl(plot_list, extract_max)
  plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x + 
                            scale_y_continuous(breaks = c(y)) + 
                            expand_limits(y = y))

  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}


Draw_HM_Heatmap <- function(df, save_path, filename, cutree_cols = 3) {
    
    # my_palette <- colorRampPalette(c("darkblue", "white", "darkred"))(500)
    # ## with hallmark
    # png(paste0(save_path, filename, '.png'),width = 4000, height = 3500, res = 300)
    # pheatmap(df,
    #         cluster_rows = TRUE,
    #         cluster_cols = TRUE,
    #         show_rownames = TRUE,
    #         show_colnames = TRUE,
    #         fontsize = 15,
    #         display_numbers = FALSE,
    #         scale = "row",
    #         color = my_palette,
    #         cutree_cols = cutree_cols
    #         # breaks = seq(-3, 3, length.out = 100),
    #         )
    # dev.off()

    # pdf(paste0(save_path, filename, '.pdf'),width = 13, height = 11)
    # pheatmap(df,
    #         cluster_rows = TRUE,
    #         cluster_cols = TRUE,
    #         show_rownames = TRUE,
    #         show_colnames = TRUE,
    #         fontsize = 15,
    #         display_numbers = FALSE,
    #         scale = "row",
    #         color = my_palette,
    #         cutree_cols = cutree_cols
    #         # breaks = seq(-3, 3, length.out = 100),
    #         )
    # dev.off()

    ## without hallmark
    rownames(df) = gsub('HALLMARK_', '', rownames(df))
    pheatmap(df,
            cluster_rows = TRUE,
            cluster_cols = TRUE,
            show_rownames = TRUE,
            show_colnames = TRUE,
            fontsize = 15,
            display_numbers = FALSE,
            scale = "row",
            color = my_palette,
            cutree_cols = cutree_cols
            # breaks = seq(-3, 3, length.out = 100),
            ) %>% as.ggplot -> p
    ggsave(p, filename = paste0(save_path, filename, '_wo_HM.png'), width = 11, height = 11)
    ggsave(p, filename = paste0(save_path, filename, '_wo_HM.pdf'), width = 11, height = 11)
    ggplot2pptx(p, 11, 11, paste0(save_path, filename, '_wo_HM.pptx'))

    pdf(paste0(save_path, filename, '_wo_HM.pdf'),width = 11, height = 11)
    pheatmap(df,
            cluster_rows = TRUE,
            cluster_cols = TRUE,
            show_rownames = TRUE,
            show_colnames = TRUE,
            fontsize = 15,
            display_numbers = FALSE,
            scale = "row",
            color = my_palette,
            cutree_cols = cutree_cols
            # breaks = seq(-3, 3, length.out = 100),
            )
    dev.off()
}

log_normalize <- function(so, save_path, data_name, pcs, nfeatures = 2000) {
    so <- NormalizeData(object = so)
    so <- FindVariableFeatures(object = so, nfeatures = nfeatures)
    so <- ScaleData(object = so)
    so <- RunPCA(object = so, npcs = 50)
    p <- ElbowPlot(so, ndims = 50, reduction = "pca")
    ggsave(p, file = paste0(save_path, data_name, '_elbowplot.png'), width = 5, height = 5)
    so <- RunHarmony(so, 'dataset')
    so <- so %>%
        FindNeighbors(reduction = "harmony") %>%
        FindClusters(resolution = 1) 
    so <- RunUMAP(so, dims = 1:pcs, reduction = 'harmony')
    return(so)
}
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
plot_lineage_features <- function(so, save_path, data_name) {
    p <- fp_sjcho(so, features = c('Sftpc', 'Etv5', 'Foxj1', 
                                    'Pdpn', 'Aqp5', 'Hopx',
                                    'Tnip3', 'Cldn4', 'Krt8',
                                    'Ptprc', 'Pecam1', 'Pdgfra'), ncol = 3, reduction = reduction)
    ggsave(p, file = paste0(save_path, data_name, '_lineage_features.png'), width = 11, height = 12)
}
