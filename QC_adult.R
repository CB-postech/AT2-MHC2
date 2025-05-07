# first, cellbender remove-background
# in shell script

# #!/bin/bash

# set -e
# set -u
# set -o pipefail

# samples=("14dpi_floxed", "14dpi_dAT2", "30dpi_floxed", "30dpi_dAT2")

# for sample_ in "${samples[@]}" 
# do
#     mkdir /home/sjcho/projects/AT2_MHC2/20241113_after_cellbender/0.cellbender/outs/"$sample_"
#     cd /home/sjcho/projects/AT2_MHC2/20241113_after_cellbender/0.cellbender/outs/"$sample_"
#     /home/sjcho/.conda/envs/cellbender/bin/cellbender remove-background \
#     --input /home/sjcho/projects/AT2_MHC2/20240929_full_dpi/0.cellragner/"$sample_"/outs/raw_feature_bc_matrix.h5 \
#     --output "$sample_"_output.h5 \
#     --cpu-threads 64
# done

# samples=("naive_dAT2" "naive_floxed" "7dpi_dAT2" "7dpi_floxed")

# for sample_ in "${samples[@]}" 
# do
#     mkdir /home/sjcho/projects/AT2_MHC2/20241113_after_cellbender/0.cellbender/outs/"$sample_"
#     cd /home/sjcho/projects/AT2_MHC2/20241113_after_cellbender/0.cellbender/outs/"$sample_"
#     /home/sjcho/.conda/envs/cellbender/bin/cellbender remove-background \
#     --input /home/sjcho/projects/AT2_MHC2/20240929_full_dpi/0.cellragner/"$sample_"/outs/raw_feature_bc_matrix.h5 \
#     --output "$sample_"_output.h5 \
#     --cpu-threads 64
# done

library(Seurat)
library(magrittr)
library(data.table)
library(ggplot2)
library(stringr)
library(scater)
library(DropletUtils)
library(harmony)
library(pheatmap)

# load custom functions
source('/home/sjcho/yard/functions/R/seurat_count_to_normalization.R')
source('/home/sjcho/yard/functions/R/FeaturePlot_sjcho.R')
source('/home/sjcho/yard/functions/R/save_ggplot2_to_ppt.R')
source('/home/sjcho/projects/AT2_MHC2/20241113_after_cellbender/utils.R')

save_path = '/home/sjcho/projects/AT2_MHC2/20241113_after_cellbender/1.log_normalization/outs/1.1.merge_and_basicQC/'
cellbender_outs = '/home/sjcho/projects/AT2_MHC2/20241113_after_cellbender/0.cellbender/outs'
files <- list.files(cellbender_outs)

counts_list = list()
dpi = c()
condition = c()
for (file in files) {
    fullpath = paste0(cellbender_outs, '/', file)
    counts <- Read10X_h5(paste0(fullpath, '/', file, '_output_filtered_seurat.h5'))
    colnames(counts) <- paste0(file, '_', colnames(counts))
    counts_list[[file]] <- counts
    dpi = c(dpi, rep(str_extract(file, "^[^_]+"), counts %>% ncol))
    condition = c(condition, rep(str_extract(file, "[^_]+$"), counts %>% ncol))
    table(condition)
}

merged <- Reduce(cbind, counts_list)
so <- CreateSeuratObject(counts = merged, min.cells = 0, min.features = 0)
so[['dpi']] = dpi; so[['condition']] = condition; 
so[['log10_UMI']] = log10(so[['nCount_RNA']][[1]] + 1)
so[['mt.pct']] <- PercentageFeatureSet(so, pattern = "^mt-")

### draw QC plot
qc_df = data.frame(log10UMI = so[['log10_UMI']][[1]], 
                   nFeature = so[['nFeature_RNA']][[1]], 
                   mt_pct = so[['mt.pct']][[1]])

p <- ggplot(qc_df, aes(x = log10UMI)) + geom_histogram(fill = 'forestgreen', bins = 200) + NoLegend() +
    geom_vline(xintercept = log10(1500), color = 'red', linetype = 'dashed', size = 0.5) +
    labs(title = 'log10UMI') +  theme_classic()
ggsave(p, file = paste0(save_path, 'QC_log10UMI_hist.png'), width = 5, height = 5)
p <- ggplot(qc_df, aes(x = mt_pct)) + geom_histogram(fill = 'skyblue', bins = 200) + NoLegend() +
    geom_vline(xintercept = 10, color = 'red', linetype = 'dashed', size = 0.5) +
    labs(title = 'mt_pct') +  theme_classic() +
    theme(axis.title.x = element_text(size = 15), 
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 15),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size = 20))
ggsave(p, file = paste0(save_path, 'QC_mt_pct_hist.png'), width = 10, height = 5)

so <- subset(so, subset = log10_UMI > log10(1500) & mt.pct < 10)
saveRDS(so, file = paste0(save_path, 'QC_passed.rds'))
