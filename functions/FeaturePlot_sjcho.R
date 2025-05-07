fp_sjcho <- function(so, features, order = FALSE, pt.size = 0.001, reduction = 'umap', ncol = 2, min.cutoff = NA, max.cutoff = NA, split.by = NULL) {
    FeaturePlot(so, features = features, order = order, raster = F, reduction = reduction, ncol = ncol, min.cutoff = min.cutoff, max.cutoff = max.cutoff, split.by = split.by,
    cols = c("#e0e0e0", "#b2182b"), pt.size = pt.size) & scale_colour_gradientn(colours = rev(c("#300000", "red","#eeeeee")))
}

fp_sjcho_green <- function(so, features, order = FALSE, pt.size = 0.001, reduction = 'umap', ncol = 2, min.cutoff = NA, max.cutoff = NA, split.by = NULL) {
    FeaturePlot(so, features = features, order = order, raster = F, reduction = reduction, ncol = ncol, min.cutoff = min.cutoff, max.cutoff = max.cutoff, split.by = split.by,
    cols = c("#e0e0e0", "#059c1e"), pt.size = pt.size) & scale_colour_gradientn(colours = rev(c("#001f01", "#04ac09","#eeeeee")))
}

fp_sjcho_graident <- function(so, features, order = FALSE, pt.size = 0.001) {
    FeaturePlot(so, features = features, order = order, raster = F, 
    cols = c("#e0e0e0", "#b2182b"), pt.size = pt.size) & scale_colour_gradientn(colours = rev(c("black", "red","#eeeeee")))
}

fp_sjcho_stereo <- function(so, features, order = FALSE, reduction = 'umap', max.cutoff = NA) {
    FeaturePlot(so, features = features, order = order, raster = F, reduction = reduction, max.cutoff = max.cutoff,
    cols = c("gray", "darkred"), pt.size = 0.3) & scale_colour_gradientn(colours = rev(c("darkred", "lightyellow","gray")))
}
