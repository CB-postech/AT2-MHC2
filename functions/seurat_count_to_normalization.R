library(Seurat)
# library(scater)
library(DropletUtils)
library(scran)
library(magrittr)
# library(vegan)
library(harmony)
library(data.table)
library(ggplot2)
library(BiocParallel)
library(future)

seurat_count_to_normalization <- function(so, metas, fig_save, sample_name, fdr = 0.05, pcs = 20, 
                                            batch = "", lambda = 1, workers = 16, use_vst = FALSE, nfeatures = 2000) {
    so_meta <- so@meta.data[, metas]
    if ('RNA' %in% names(so@assays)) {
        so <- CreateSeuratObject(so@assays$RNA$counts, min.cells = 0, min.features = 0)
    }

    if ('originalexp' %in% names(so@assays)) {
        so <- CreateSeuratObject(so@assays$originalexp$counts, min.cells = 0, min.features = 0)
    }

    sce <- as.SingleCellExperiment(so)
    # sce <- sce[rowSums(counts(sce)) > 0, ]
    # sce <- sce[, colSums(counts(sce)) > 0]
    clusters <- quickCluster(sce)
    # min.size : 100 (default) -> minimum size of each cluster
    # method : hclust, igraph.
    # hclust : hierarchical clustering
    # igraph : graph-based clustering
    # d : dimension
    # 
    param <- MulticoreParam(workers = workers)  # use cores as input

    sce <- computeSumFactors(sce, clusters = clusters, BPPARAM = param)
    sce.norm <- logNormCounts(sce,  pseudo_count = 1)
    print('Scran normalization is Done!')

    gene_variance <- modelGeneVar(sce.norm)
    # , block = so[['Sample']][[1]]
    # Raw and adjusted pvalues for the test against the null hypothesis that bio <= 0
    # fitting 된 value : technical variance로 해석
    # fitting된 curve에서 올라간 부분이 biological variance로 해석

    # block 별로 따로따로 statistics 계산
    # mean, variance -> done by averaging values across levels
    # equiweight=FALSE -> equaly weight가 false -> block마다 cell number를 계산해서 고려해줌
    # Use of block is the recommended approach for accounting for any uninteresting categorical factor of variation

    hvg.norm <- getTopHVGs(gene_variance, fdr.threshold = fdr)
    length(hvg.norm) 

    png(paste0(fig_save, sample_name, '_scran_variance_mean.png'))
    plot(gene_variance$mean, gene_variance$total,
        xlab = "mean_log-expression",
        ylab = "variance",
        main = paste0('genes under FDR 0.05 : ', length(hvg.norm)))
    curve(metadata(gene_variance)$trend(x),
        col = "blue",
        add = T)
    dev.off()

    so <- CreateSeuratObject(assay(sce.norm, "counts"), min.cells = 0, min.features = 0)
    so@assays$RNA$data <- assay(sce.norm, "logcounts")

    all.genes <- rownames(so)
    so <- ScaleData(so, features = all.genes)
    so <- FindVariableFeatures(so, nfeatures = nfeatures)

    top10 <- head(VariableFeatures(so), 10)
    print(top10)
    # plot variable features with and without labels
    plot1 <- VariableFeaturePlot(so)
    plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
    ggsave(plot2+ theme_classic(), file = paste0(fig_save, sample_name, '_sct_VariableFeaturePlot.png'))

    if (!use_vst)
    {
        # default is using HVG detected by scran
        VariableFeatures(so) = hvg.norm
    }

    so <- RunPCA(so, features = VariableFeatures(so), npcs = 50, verbose = FALSE)
    p <- ElbowPlot(so, ndims = 50)
    ggsave(p + theme_classic(), file = paste0(fig_save, sample_name, '_elbow.png'))

    for (meta_features in metas) {
        so[[meta_features]] <- so_meta[[meta_features]]
    }

    if (batch != "") {
        so <- RunHarmony(so, batch, lambda = lambda)
        so <- so %>%
        FindNeighbors(reduction = "harmony", dims = 1:pcs) %>%
        FindClusters(resolution = 1) 

        so <- RunUMAP(so, dims = 1:pcs, reduction = 'harmony')
    }
    else {
        so <- FindNeighbors(so,
                            dims = 1:pcs)
        so <- FindClusters(so,
                        resolution = 1)
        so <- RunUMAP(so, dims = 1:pcs)
    }
    print('Its DONE!')

    return(so)
}

seurat_count_to_normalization_before_harmony <- function(so, metas, fig_save, sample_name, fdr = 0.05, pcs = 20, batch = "", lambda = 1, workers = 16) {
    # so_meta <- so@meta.data[, metas]
    #    so <- CreateSeuratObject(so@assays$RNA@counts, min.cells = 0, min.features = 0)

    # for (meta_features in metas) {
    #     so[[meta_features]] <- so_meta[[meta_features]]
    # }

    sce <- as.SingleCellExperiment(so)
    clusters <- quickCluster(sce)
    # min.size : 100 (default) -> minimum size of each cluster
    # method : hclust, igraph.
    # hclust : hierarchical clustering
    # igraph : graph-based clustering
    # d : dimension
    # 
    library(BiocParallel)

    param <- MulticoreParam(workers = workers)  # use cores as input

    sce <- computeSumFactors(sce, clusters = clusters, BPPARAM = param)
    sce.norm <- logNormCounts(sce,  pseudo_count = 1)

    gene_variance <- modelGeneVar(sce.norm)
    # , block = so[['Sample']][[1]]
    # Raw and adjusted pvalues for the test against the null hypothesis that bio <= 0
    # fitting 된 value : technical variance로 해석
    # fitting된 curve에서 올라간 부분이 biological variance로 해석

    # block 별로 따로따로 statistics 계산
    # mean, variance -> done by averaging values across levels
    # equiweight=FALSE -> equaly weight가 false -> block마다 cell number를 계산해서 고려해줌
    # Use of block is the recommended approach for accounting for any uninteresting categorical factor of variation

    hvg.norm <- getTopHVGs(gene_variance, fdr.threshold = fdr)
    length(hvg.norm) 

    png(paste0(fig_save, sample_name, '_scran_variance_mean.png'))
    plot(gene_variance$mean, gene_variance$total,
        xlab = "mean_log-expression",
        ylab = "variance",
        main = paste0('genes under FDR 0.05 : ', length(hvg.norm)))
    curve(metadata(gene_variance)$trend(x),
        col = "blue",
        add = T)
    dev.off()

    so <- CreateSeuratObject(assay(sce.norm, "counts"), min.cells = 0, min.features = 0)
    so@assays$RNA$data <- assay(sce.norm, "logcounts")

    all.genes <- rownames(so)
    so <- ScaleData(so, features = all.genes)
    so <- FindVariableFeatures(so)

    top10 <- head(VariableFeatures(so), 10)
    # plot variable features with and without labels
    plot1 <- VariableFeaturePlot(so)
    plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
    print(top10)
    ggsave(plot2+ theme_classic(), file = paste0(fig_save, sample_name, '_sct_VariableFeaturePlot.png'))

    VariableFeatures(so) = hvg.norm

    so <- RunPCA(so, features = VariableFeatures(so), npcs = 50, verbose = FALSE)
    p <- ElbowPlot(so, ndims = 50)
    ggsave(p + theme_classic(), file = paste0(fig_save, sample_name, '_elbow.png'))
    print('Its DONE!')
    return(so)
}

seurat_count_to_normalization_before_PCA <- function(so, metas, fig_save, sample_name, fdr = 0.05, pcs = 20, batch = "", lambda = 1, workers = 16) {
    so_meta <- so@meta.data[, metas]

    tryCatch({so <- CreateSeuratObject(so@assays$RNA@counts, min.cells = 0, min.features = 0)},
    error = function(e) {
        so <- CreateSeuratObject(so@assays$RNA$counts, min.cells = 0, min.features = 0)
    }    
    )

    for (meta_features in metas) {
        so[[meta_features]] <- so_meta[[meta_features]]
    }

    sce <- as.SingleCellExperiment(so)
    clusters <- quickCluster(sce)
    # min.size : 100 (default) -> minimum size of each cluster
    # method : hclust, igraph.
    # hclust : hierarchical clustering
    # igraph : graph-based clustering
    # d : dimension
    # 
    library(BiocParallel)

    param <- MulticoreParam(workers = workers)  # use cores as input

    sce <- computeSumFactors(sce, clusters = clusters, BPPARAM = param)
    sce.norm <- logNormCounts(sce,  pseudo_count = 1)

    gene_variance <- modelGeneVar(sce.norm)
    # , block = so[['Sample']][[1]]
    # Raw and adjusted pvalues for the test against the null hypothesis that bio <= 0
    # fitting 된 value : technical variance로 해석
    # fitting된 curve에서 올라간 부분이 biological variance로 해석

    # block 별로 따로따로 statistics 계산
    # mean, variance -> done by averaging values across levels
    # equiweight=FALSE -> equaly weight가 false -> block마다 cell number를 계산해서 고려해줌
    # Use of block is the recommended approach for accounting for any uninteresting categorical factor of variation

    hvg.norm <- getTopHVGs(gene_variance, fdr.threshold = fdr)
    length(hvg.norm) 

    png(paste0(fig_save, sample_name, '_scran_variance_mean.png'))
    plot(gene_variance$mean, gene_variance$total,
        xlab = "mean_log-expression",
        ylab = "variance",
        main = paste0('genes under FDR 0.05 : ', length(hvg.norm)))
    curve(metadata(gene_variance)$trend(x),
        col = "blue",
        add = T)
    dev.off()

    so <- CreateSeuratObject(assay(sce.norm, "counts"), min.cells = 0, min.features = 0)
    so@assays$RNA$data <- assay(sce.norm, "logcounts")

    all.genes <- rownames(so)
    so <- ScaleData(so, features = all.genes)
    so <- FindVariableFeatures(so)

    top10 <- head(VariableFeatures(so), 10)
    # plot variable features with and without labels
    plot1 <- VariableFeaturePlot(so)
    plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
    ggsave(plot2+ theme_classic(), file = paste0(fig_save, sample_name, '_sct_VariableFeaturePlot.png'))

    VariableFeatures(so) = hvg.norm
    print('Its DONE!')
    return(so)
}

seurat_count_to_normalization_stereo <- function(so, metas, fig_save, sample_name, fdr = 0.05, pcs = 20, batch = "", lambda = 1, workers = 16, use_vst = FALSE, nfeatures = 2000) {
    so_meta <- so@meta.data[, metas]
    so <- CreateSeuratObject(so@assays$Spatial$counts, min.cells = 0, min.features = 0)

    for (meta_features in metas) {
        so[[meta_features]] <- so_meta[[meta_features]]
    }
#     so@meta.data[, meta_features] <- so_meta[[meta_features]]

    sce <- as.SingleCellExperiment(so)
    clusters <- quickCluster(sce)
    # min.size : 100 (default) -> minimum size of each cluster
    # method : hclust, igraph.
    # hclust : hierarchical clustering
    # igraph : graph-based clustering
    # d : dimension
    
    param <- MulticoreParam(workers = workers)  # use cores as input
    print('Scran normalization start...')
    sce <- computeSumFactors(sce, clusters = clusters, BPPARAM = param)
    sce.norm <- logNormCounts(sce,  pseudo_count = 1)

    gene_variance <- modelGeneVar(sce.norm)
    # , block = so[['Sample']][[1]]
    # Raw and adjusted pvalues for the test against the null hypothesis that bio <= 0
    # fitting 된 value : technical variance로 해석
    # fitting된 curve에서 올라간 부분이 biological variance로 해석

    # block 별로 따로따로 statistics 계산
    # mean, variance -> done by averaging values across levels
    # equiweight=FALSE -> equaly weight가 false -> block마다 cell number를 계산해서 고려해줌
    # Use of block is the recommended approach for accounting for any uninteresting categorical factor of variation

    hvg.norm <- getTopHVGs(gene_variance, fdr.threshold = fdr)
    print(length(hvg.norm)) 

    png(paste0(fig_save, sample_name, '_scran_variance_mean.png'))
    plot(gene_variance$mean, gene_variance$total,
        xlab = "mean_log-expression",
        ylab = "variance",
        main = paste0('genes under FDR 0.05 : ', length(hvg.norm)))
    curve(metadata(gene_variance)$trend(x),
        col = "blue",
        add = T)
    dev.off()

    so <- CreateSeuratObject(assay(sce.norm, "counts"), min.cells = 0, min.features = 0)
    so@assays$RNA$data <- assay(sce.norm, "logcounts")

    # so <- as.Seurat(sce.norm,
    #             counts = "counts",
    #             data = "logcounts")

    all.genes <- rownames(so)
    so <- ScaleData(so, features = all.genes)
    so <- FindVariableFeatures(so, nfeatures = nfeatures)

    top10 <- head(VariableFeatures(so), 10)
    print(top10)
    # plot variable features with and without labels
    plot1 <- VariableFeaturePlot(so)
    plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
    ggsave(plot2+ theme_classic(), file = paste0(fig_save, sample_name, '_sct_VariableFeaturePlot.png'))

    if (!use_vst)
    {
        # default is using HVG detected by scran
        VariableFeatures(so) = hvg.norm
    }

    so <- RunPCA(so, features = VariableFeatures(so), npcs = 50, verbose = FALSE)
    p <- ElbowPlot(so, ndims = 50)
    ggsave(p + theme_classic(), file = paste0(fig_save, sample_name, '_elbow.png'))

    if (batch != "") {
        so <- RunHarmony(so, batch, lambda = lambda)
        so <- so %>%
        FindNeighbors(reduction = "harmony") %>%
        FindClusters(resolution = 1) 

        so <- RunUMAP(so, dims = 1:pcs, reduction = 'harmony')
    }
    else {
        so <- FindNeighbors(so,
                            dims = 1:pcs)
        so <- FindClusters(so,
                        resolution = 1)
        so <- RunUMAP(so, dims = 1:pcs)
    }
    print('Its DONE!')
    return(so)
}
