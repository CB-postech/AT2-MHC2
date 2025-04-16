### run dropletUtils to detect empty droplets
# based on https://github.com/CB-postech/2022_KOGO_workshop/blob/main/KOGO_QC1-DropletUtils.md

library(Seurat)
library(scater)
library(DropletUtils)
library(scran)

save_path = '/home/sjcho/projects/AT2_MHC2/20241113_after_cellbender/PND_Rcode_to_sjcho/inhouse_PND1_7_from_cellranger/preprocessing/outs/lower_300/'
cellranger_outputs <- "/home/sjcho/projects/AT2_MHC2/20241113_after_cellbender/PND_Rcode_to_sjcho/inhouse_PND1_7_from_cellranger/cellranger"

filenames <- list.dirs(path = cellranger_outputs, full.names = FALSE, recursive = FALSE)
print(filenames)

filtered_sce_list <- list()
for(dpi_condition in filenames) {
    rawsce <- read10xCounts(paste0(cellranger_outputs, '/', dpi_condition, '/outs/raw_feature_bc_matrix'), type = "sparse", compressed = TRUE)
    br.out <- barcodeRanks(counts(rawsce))
    
    png(paste0(save_path, 'DropletUtils_', dpi_condition, '.png'))

    plot(br.out$rank, br.out$total, log = "xy", xlab = "Rank", ylab = "Total", main = dpi_condition)

    o <- order(br.out$rank)
    lines(br.out$rank[o], br.out$fitted[o], col = "red")

    e.out <- emptyDrops(counts(rawsce), lower = 300)  ## Cells that have UMI counts lower than 100 are empty cells.
    table(Sig=e.out$FDR <= 0.05, Limited=e.out$Limited)
    is.cell <- e.out$FDR <= 0.05

    print(sum(is.cell, na.rm=TRUE))
    print(table(br.out$rank == sum(is.cell, na.rm=TRUE)))

    abline(h=min(br.out$fitted[o], na.rm=TRUE), col="red", lty=2)
    abline(h=metadata(br.out)$knee, col="dodgerblue", lty=2)
    abline(h=metadata(br.out)$inflection, col="forestgreen", lty=2)
    legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen", "red"), legend=c("knee", "inflection", "FDR_0.05"))
    dev.off()

    colnames(rawsce) = colData(rawsce)$Barcode
    rawsce <- rawsce[,which(e.out$FDR <= 0.05)]
    filtered_sce_list[[dpi_condition]] <- rawsce
}

saveRDS(filtered_sce_list, file = paste0(save_path, 'filtered_sce_list.rds'))

intersect(paste0('PND1___', Cells(filtered_sce_list[[1]])), cells)

### basic_QC
library(Seurat)
library(scater)
library(SingleCellExperiment)
library(magrittr)
library(ggplot2)
library(org.Mm.eg.db)

save_path = '/home/sjcho/projects/AT2_MHC2/20241113_after_cellbender/PND_Rcode_to_sjcho/inhouse_PND1_7_from_cellranger/preprocessing/outs/preprocessing/'
save_path = '/home/sjcho/projects/AT2_MHC2/20241113_after_cellbender/PND_Rcode_to_sjcho/inhouse_PND1_7_from_cellranger/preprocessing/outs/preprocessing/detected_filtering/'
cellranger_outputs <- "/home/sjcho/projects/AT2_MHC2/20241113_after_cellbender/PND_Rcode_to_sjcho/inhouse_PND1_7_from_cellranger/cellranger"

filtered_sce_list = readRDS('/home/sjcho/projects/AT2_MHC2/20241113_after_cellbender/PND_Rcode_to_sjcho/inhouse_PND1_7_from_cellranger/preprocessing/outs/preprocessing/filtered_sce_list.rds')

preprocess <- function(sce, save_path = save_path, sample_name){
  rownames(sce) = uniquifyFeatureNames(rowData(sce)$ID, rowData(sce)$Symbol)
  
  mtgenes = rowData(sce)[grep("^mt-", rowData(sce)$Symbol),]$Symbol
  is.mito = rownames(sce) %in% mtgenes
  
  sce <- addPerCellQC(
    sce,
    subsets = list(MT=mtgenes),
    percent_top = c(50, 100, 200, 500), 
    detection_limit = 0
  )
  
  sce$log10_sum = log10(sce$sum + 1)
  sce$log10_detected = log10(sce$detected + 1)
  
  df <- data.frame(
    log10_sum = sce$log10_sum,
    detected = sce$log10_detected,
    subsets_MT_percent = sce$subsets_MT_percent,
    percent.top_500 = sce$percent.top_500
  )
  p <- ggplot(df, aes(x = log10_sum)) + geom_histogram(bins = 300) + ggtitle(paste0(sample_name, ' log10(nCount)'))
  p <- p + theme_classic()
  p <- p + geom_vline(xintercept = log10(1000), color = "red")
  ggsave(p, file = paste0(save_path, sample_name, '_log10_sum.png'))

  p <- ggplot(df, aes(x = subsets_MT_percent)) + geom_histogram(bins = 300) + ggtitle(paste0(sample_name, ' MT percent'))
  p <- p + theme_classic()
  p <- p + geom_vline(xintercept = 15, color = "red")
  ggsave(p, file = paste0(save_path, sample_name, '_MT_pct.png'))

  p <- ggplot(df, aes(x = detected)) + geom_histogram(bins = 300) + ggtitle(paste0(sample_name, ' log10(nFeature)'))
  p <- p + theme_classic()
  p <- p + geom_vline(xintercept = log10(750), color = "red")
  ggsave(p, file = paste0(save_path, sample_name, '_log10_detected.png'))

  sce <- sce[,sce$sum!=0]
  return(sce)
}

filtering <- function(sce,umi = 1000, mtpct = 5, detect = 750, save_path, sample_name){
  filter_by_total_counts = sce$sum > umi
  filter_by_mt_percent = sce$subsets_MT_percent < mtpct
  filter_by_nfeature = sce$detected > detect
  
  sce <- runColDataPCA(sce, variables = list("sum", "detected", "subsets_MT_percent", "percent.top_500"))
  
  sce$use <- (
    filter_by_total_counts &
      filter_by_mt_percent &
      filter_by_nfeature
  )
  
  p <- plotReducedDim(sce, dimred="PCA_coldata", colour_by="use")
  p <- p + ggtitle(sample_name) + theme_classic()
  ggsave(p, file = paste0(save_path, sample_name, '_PCA_use.png'))

  sce = sce[,sce$use]
  return(sce)
}

for(sample_name in names(filtered_sce_list)){
	sce = filtered_sce_list[[sample_name]]

	### convert esemble to symbol by Read10X_h5
	### we can leap the headache of converting gene symbols by this function
	symbols <- rownames(Read10X_h5(paste0(cellranger_outputs, '/', sample_name, '/outs/raw_feature_bc_matrix.h5')))
	rownames(sce) <- symbols

	sce = preprocess(sce, save_path = save_path, sample_name = sample_name)
	sce = filtering(sce, save_path = save_path, sample_name = sample_name)
	
	filtered_sce_list[[sample_name]] = sce
}

saveRDS(filtered_sce_list, file = paste0(save_path, 'filtered_sce_list_after_basicQC.rds'))

for(sample_name in names(filtered_sce_list)){
  sce = filtered_sce_list[[sample_name]]
  colnames(sce) <- paste0(sample_name, '_', colnames(sce))
  filtered_sce_list[[sample_name]] <- sce
}

sce_merged <- cbind(filtered_sce_list[[1]], filtered_sce_list[[2]])

saveRDS(sce_merged, file = paste0(save_path, 'sce_merged.rds'))
