
log_normalize <- function(so, save_path, data_name, pcs, nfeatures = 2000) {
    so <- NormalizeData(object = so)
    so <- FindVariableFeatures(object = so, nfeatures = nfeatures)
    so <- ScaleData(object = so)
    so <- RunPCA(object = so, npcs = 50)
    p <- ElbowPlot(so, ndims = 50, reduction = "pca")
    ggsave(p, file = paste0(save_path, data_name, '_elbowplot.png'), width = 5, height = 5)
    so <- so %>%
        FindNeighbors(reduction = 'pca') %>%
        FindClusters(resolution = 1) 
    so <- RunUMAP(so, dims = 1:pcs, reduction = 'pca')
    return(so)
}
plot_alveolar_epithelium_features <- function(so, save_path, data_name) {
    p <- fp_sjcho(so, features = c('Sftpc', 'Etv5', 'Abca3', 'Pdpn',
                                    'Aqp5', 'Foxj1', 'Scgb1a1',
                                    'Scgb3a2', 'Krt5', 'Trp63', 'Cldn4', 
                                    'Sprr1a', 'Tnip3', 'Krt8', 'Mki67', 'Top2a'), ncol = 4)
    ggsave(p, file = paste0(save_path, data_name, '_epi_features.png'), width = 13, height = 15)
    p <- fp_sjcho(so, features = c('Sftpc', 'Etv5', 'Foxj1', 
                                    'Pdpn', 'Aqp5', 'Hopx',
                                    'Tnip3', 'Cldn4', 'Krt8',
                                    'Ptprc', 'Pecam1', 'Pdgfra'), ncol = 3)
    ggsave(p, file = paste0(save_path, data_name, '_epi_features2.png'), width = 11, height = 12)
}
plot_lineage_features <- function(so, save_path, data_name) {
    p <- fp_sjcho(so, features = c('Sftpc', 'Etv5', 'Foxj1', 
                                    'Pdpn', 'Aqp5', 'Hopx',
                                    'Tnip3', 'Cldn4', 'Krt8',
                                    'Ptprc', 'Pecam1', 'Pdgfra'), ncol = 3)
    ggsave(p, file = paste0(save_path, data_name, '_lineage_features.png'), width = 11, height = 12)
}
plot_dpi_condition <- function(so, save_path, data_name, group.by = 'seurat_clusters') {
    p <- DimPlot(so, group.by = group.by, split.by = 'condition', label = T, repel = T, ncol = 2)
    ggsave(p, file = paste0(save_path, data_name, '_condition.png'), width = 13.5, height = 6)
    p <- DimPlot(so, group.by = group.by, split.by = 'dpi', label = T, repel = T, ncol = 2)
    ggsave(p, file = paste0(save_path, data_name, '_dpi.png'), width = 13.5, height = 12)
    p <- DimPlot(so, group.by = group.by, label = T, repel = T)
    ggsave(p, file = paste0(save_path, data_name, '_cluster.png'), width = 8, height = 6)
}
plot_tech <- function(so, save_path, data_name) {
  p <- fp_sjcho(so, features = c('log10_UMI', 'nFeature_RNA', 'mt.pct'), ncol = 2)
  # RdBu
  library(RColorBrewer)
  p <- p & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
  ggsave(p, file = paste0(save_path, data_name, '_tech.png'), width = 9, height = 5)
}

colors_dpi = c('#e6b0aa', '#c0392b', '#d7bde2', '#9b59b6', '#a9cce3', '#2980b9', '#abebc6', '#2ecc71')

lineage_cols = c('#bd5e7e', '#292f56', '#faec70', '#9f149d')
names(lineage_cols) = c('Epithelial', 'Immune', 'Endothelial', 'Stromal')

epithelial_cell_order = c('AT1', 'AT2', 'Alveolar_Proliferating', 'Alveolar_Transitory', 'Club', 'Basal', 'Ciliated', 'Airway_Proliferating')
epithelial_cols = c('#950404', '#eb0432', '#240000', '#ff9ed8', '#c29800', '#c14c10', '#ffe27a', '#b9377e'); names(epithelial_cols) = epithelial_cell_order

# alv_cell_order = c('AT2', 'pAT2', 'late_Transitory', 'early_Transitory', 'mature_AT1', 'premature_AT1')
# alv_cols = c('#c93a65', '#e29aad','#a50c16', '#6f1566', 'red', '#babd00')
# names(alv_cols) = alv_cell_order

# alv_cell_full_order = c(alv_cell_order, 'proliferating', 'Interfereon_response')
# alv_cols_full = c(as.vector(alv_cols), 'hotpink', 'black')
# names(alv_cols_full) = alv_cell_full_order

alv_cell_order = c('AT2', 'pAT2', 'late_Transitory', 'early_Transitory', 'mature_AT1', 'premature_AT1')
alv_cols = c('#c93a65', '#e29aad','#a50c16', '#6f1566', 'red', '#babd00')
names(alv_cols) = alv_cell_order

alv_cell_full_order = c(alv_cell_order, 'proliferating', 'Interfereon_response')
alv_cols_full = c(as.vector(alv_cols), 'hotpink', 'black')
names(alv_cols_full) = alv_cell_full_order

Tcell_cols_order = c('naive-like Cd4 T cell', 'naive-like Treg', 'Cxcr3+ Treg', 'Proliferating Treg', 'cycling Cd4 T cell', 'Th1 Cd4 T cell', 'Th17 Cd4 T cell', 'Cd4 Tcm')
Tcell_cols = c('#67cc7c', '#9d2c5c', '#f3b66e', '#016a01', '#caa1dd', '#53b1db', '#f90068', '#7f7f7f')
names(Tcell_cols) = Tcell_cols_order

corrected_deg <- function(DEG_df, celltype_color = NULL) {
  # modeling linear regression
  model <- lm(num_deg ~ num_cell, data = DEG_df)
  
  # calculating residuals for each cell type
  DEG_df$residuals <- residuals(model)
  
  # ranking cell type
  DEG_df_rank <- DEG_df %>%
    arrange(desc(residuals)) %>%
    mutate(rank = row_number())
  
  # plotting rankplot
  ggplot(DEG_df_rank, aes(x = rank, y = residuals, color = celltype)) +
    geom_point(size = 3) +
    geom_line(aes(group = celltype)) +
    geom_text_repel(aes(label = celltype), vjust = -1, hjust = 0.2, size = 3, color = "black") +
    # scale_color_manual(values = celltype_color) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    labs(x = "Rank",
         y = "Corrected Number of DEGs") +
    theme_classic()
}

convert_to_celltypist_input <- function(so, save_path, set_name) {
    cellnames = as.data.frame(Cells(so))
    colnames(cellnames) = 'cellname'
    fwrite(cellnames, paste0(save_path, set_name, 'cellnames.csv'))

    clusters = as.data.frame(so[['seurat_clusters']][[1]])
    colnames(clusters) = 'seurat_clusters'
    rownames(clusters) =  cellnames[['Cells']]
    fwrite(clusters, paste0(save_path, set_name, 'clusters.csv'))

    genenames = as.data.frame(rownames(so))
    colnames(genenames) = 'genename'
    fwrite(genenames, paste0(save_path, set_name, 'genenames.csv'))

    Matrix::writeMM(t(so@assays$RNA$counts), file = paste0(save_path, set_name, 'counts.mtx'))
}



library(ggrepel)
draw_volcano <- function(FM, p_or_q, FC_cutoff, p_val_adj_cutoff, max_p_val_adj = 100, manual_genes = NULL, manual_color = "black") {
  log2FC <- FM[['avg_log2FC']]
  if (p_or_q == 'p_val') {
    p_val_adj <- FM[['p_val']]
  } else {
    p_val_adj <- FM[['p_val_adj']]
  }
  df <- data.frame(FM[['avg_log2FC']], -log10(p_val_adj))
  rownames(df) <- rownames(FM)
  colnames(df) <- c('log2FC', 'log10p_val_adj')
  
  df[['log10p_val_adj']][df[['log10p_val_adj']] > max_p_val_adj] = max_p_val_adj
  
  df_highlight = subset(df, (log2FC > FC_cutoff | log2FC < -(FC_cutoff)) & (-log10(p_val_adj) > -log10(p_val_adj_cutoff)))
  df_highlight[['name']] = rownames(df_highlight)
  df_highlight[['group']][df_highlight[['log2FC']] > FC_cutoff] = 'exp_enrich'
  df_highlight[['group']][df_highlight[['log2FC']] < -(FC_cutoff)] = 'control_enrich'
  
  if (!is.null(manual_genes)) {
    manual_data <- df[rownames(df) %in% manual_genes,]
    manual_data[['name']] <- rownames(manual_data)
  }
  
  p <- ggplot(df, aes(x=log2FC, y=log10p_val_adj)) + 
    geom_point(size = 2) +
    geom_point(data = df_highlight, aes(color = group, size = 2)) +
    scale_color_manual(values = c("control_enrich" = "blue", "exp_enrich" = "red")) +
    geom_hline(yintercept=0, linetype="solid", color="black", size = 0.5) +
    geom_vline(xintercept=0, linetype="solid", color="black", size = 0.5) +
    geom_hline(yintercept=-log10(p_val_adj_cutoff), linetype="dashed", color="blue", size = 0.2) +
    geom_vline(xintercept=c(-FC_cutoff, FC_cutoff), linetype="dashed", color="blue", size = 0.2) +
    theme_classic() +
    ylab('-log10(adjusted p-value)') +
    xlab('log2(average FoldChange)') +
    theme(plot.title = element_text(size = 15, face = 'bold'),
          axis.title = element_text(size = 15),
          axis.text = element_text(size = 10),
          legend.position = 'none')
  
  if (!is.null(manual_genes)) {
    p <- p + geom_text_repel(data = manual_data, aes(label = manual_data$name),
                             box.padding = 2, point.padding = 2, segment.color = 'gray50', 
                             force = 50, size = 4, color = manual_color)
  }
  else {
    p <- p +  geom_text_repel(data = df_highlight, aes(label = df_highlight$name, color = group), 
                    box.padding = 2, point.padding = 2, segment.color = 'gray50', force = 50, size = 4)
  }

  return(p)
}

building_df_geneset_wise <- function(so, features_names, control_cells, experimental_cells) {
    fc_vector = c(); p_value_vector = c(); subtract_vector = c()
    control_mean_vector = c(); experimental_mean_vector = c()
    for (feature in features_names) {
        fc_vector = c(fc_vector, mean(so@meta.data[experimental_cells, paste0(feature, '1')]) / mean(so@meta.data[control_cells, paste0(feature, '1')]))
        p_value_vector = c(p_value_vector, wilcox.test(so@meta.data[experimental_cells, paste0(feature, '1')], so@meta.data[control_cells, paste0(feature, '1')])$p.value)
        subtract_vector = c(subtract_vector, mean(so@meta.data[experimental_cells, paste0(feature, '1')]) - mean(so@meta.data[control_cells, paste0(feature, '1')]))

        control_mean_vector = c(control_mean_vector, mean(so@meta.data[control_cells, paste0(feature, '1')]))
        experimental_mean_vector = c(experimental_mean_vector, mean(so@meta.data[experimental_cells, paste0(feature, '1')]))
    }

    df = data.frame(features = features_names, 
                    subtract_mean = subtract_vector,
                    control_mean = control_mean_vector,
                    exercise_mean = experimental_mean_vector, fc = fc_vector,
                    p_value = p_value_vector)
    return(df)
}

draw_geneset_wise_dotplot <- function(so, features_list, fig_save, set_name, control_cells, experimental_cells, legend = 'Experiment') {
    for (feature_name in names(features_list)) {
        features = features_list[[feature_name]]
        so <- AddModuleScore(so, features = list(features), name = feature_name)
    }

    df <- building_df_geneset_wise(so, names(features_list), control_cells, experimental_cells)

    df_na = df[!(is.na(df[['p_value']])), ]
    df_na[['avg_log2FC']] = log2(df_na[['fc']])
    df_na[['p_val_adj']] = p.adjust(df_na[['p_value']], method = 'BH')
    rownames(df_na) = df_na[['features']]

    df_na[['change']] = ifelse(df_na[['subtract_mean']] > 0, 'Increase with Experiment', 'Decrease with Experiment')
    df_na[['change']][(df_na[['p_val_adj']] > 0.05)] = 'Not significant'

    df_na[['features']] = factor(df_na[['features']], levels = names(features_list))
    p <- ggplot(df_na, aes(x = subtract_mean, y = features)) + 
        geom_point(aes(size = -log10(p_val_adj), color = change)) +
        scale_color_manual(values = c('Increase with Experiment' = 'darkred', 'Decrease with Experiment' = 'darkblue', 'Not significant' = 'darkgrey')) +
        theme_classic() + 
        geom_vline(xintercept = 0, color = 'black', linetype = 'dashed') +
        theme(axis.text.y = element_text(face="bold"),
            axis.title.y = element_text(size = 0),
        ) +
        labs(x = 'Differential Score') + 
        guides(colour = guide_legend(override.aes = list(size=5)))
    ggsave(p, file = paste0(fig_save, 'dot_', set_name, '.png'), width = 6, height = 4)
    ggplot2pptx(p, 6, 4, paste0(fig_save, 'dot_', set_name, '.pptx'))
    return(list(so, p, df_na))
}


with_vlnplot <- function(so, features, fig_save, set_name, width = 8, height = 8) {
    outs <- draw_celltype_wise_dotplot(so, features, fig_save, set_name)
    so <- outs[[1]]
    head(so)
    vln <- VlnPlot(so, features = paste0(set_name, 1), group.by = 'marker_annotation', split.by = 'Condition', cols = c('darkblue', 'darkred'), pt.size = 0) + NoLegend()
    vln <- vln + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank())
    dot <- outs[[2]]
    dot <- dot + coord_flip()
    dot <- dot + theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(vln + dot, file = paste0(fig_save, 'vln_dot_', set_name, '.png'), width = 8, height = 8)
}

# from https://www.sinobiological.com/resource/cytokines/all-proinflammatory-cytokines
cytokine_receptor_lists <- list(IL1_family = c("Il18", "Il18bp", "Il1a", "Il1b", "Il1f10", "Il1rn", "Il1f5", "Il1f6", "Il1f7", "Il1f8", "Il1rl2", "Il1f9", "Il33"), 
                                IL1_receptors = c("Il18r1", "Il18rap", "Il1r1", "Il1r2", "Il1r3", "Il1r8", "Il1r9", "Il1rl1", "Sigirr"), 
                                TNF_family = c("Tnfsf13b", "Tnfsf9", "Tnfsf8", "Cd40lg", "Cd70", "Fasl", "Eda", "Tnfsf14", "Lta", "Ltb", "Tnf", "Tnfsf10", 
                                                "Tnfsf11", "Tnfsf12", "Tnfsf13", "Tnfsf15", "Tnfsf4"), 
                                TNF_receptor = c("Tnfrsf9", "Tnfrsf13c", "Cd27", "Cd40", "Fas", "Tnfrsf6b", "Tnfrsf21", "Eda2r", "Edar", 
                                                "Pglyrp1", "Tnfrsf19l", "Tnfrsf1a", "Tnfrsf1b", "Tnfrsf11a", "Tnfrsf11b", "Tnfrsf12a", "Tnfrsf13b", 
                                                "Tnfrsf14", "Tnfrsf17", "Tnfrsf18", "Tnfrsf19", "Tnfrsf25", "Ltbr", "Tnfrsf4", "Tnfrsf8", "Tnfrsf10b", "Tnfrsf10a", "Tnfrsf10c", "Tnfrsf10d"), 
                                IFN = c("Ifna1", "Ifna10", "Ifna13", "Ifna14", "Ifna2", "Ifna4", "Ifna7", "Ifnb1", "Ifne", "Ifng", "Ifnz", "Ifna8", "Ifna5", "Ifnw1"), 
                                IFN_receptor = c("Ifnar1", "Ifnar2", "Ifngr1", "Ifngr2"), 
                                IL6_family = c("Clcf1", "Cntf", "Il11", "Il31", "Il6", "Lep", "Lif", "Osm"), 
                                IL6_receptor = c("Cntfr", "Il11ra", "Il6r", "Lepr", "Lifr", "Osmr", "Il31ra"),
                                    IL10_family = c("Il10", "Il19", "Il20", "Il22", "Il24", "Il28b", "Il28a", "Il29"),
                                    IL10_receptor = c("Il10ra", "Il10rb", "Il20ra", "Il20rb", "Il22ra2", "Il22ra1"), 
                                    TGFbeta_family = c("Tgfb1", "Tgfb2", "Tgfb3"), 
                                    TGFbeta_receptor = c("Acvr1c", "Atf2", "Eng", "Tgfbr1", "Tgfbr2", "Tgfbr3"), 
                                    chemokine = c("Ccl1", "Ccl11", "Ccl12", "Ccl13", "Ccl14", "Ccl15", "Ccl16", "Ccl17", "Ccl18", "Ccl19", "Ccl2", "Ccl20", 
                                                "Ccl21", "Ccl22", "Ccl23", "Ccl24", "Ccl25", "Ccl26", "Ccl27", "Ccl28", "Ccl3", "Ccl3l3", "Ccl4", "Ccl4l1", 
                                                "Ccl5", "Ccl6", "Ccl7", "Ccl8", "Ccl9", "Cx3cl1", "Cxcl1", "Cxcl10", "Cxcl11", "Cxcl12", "Cxcl13", "Cxcl14", 
                                                "Cxcl15", "Cxcl16", "Cxcl17", "Cxcl2", "Cxcl3", "Pf4", "Cxcl5", "Cxcl6", "Ppbp", "Cxcl9", "Cxcl8", "Xcl1", 
                                                "Xcl2", "Fam19a1", "Fam19a2", "Fam19a3", "Fam19a4", "Fam19a5"), 
                                chemokine_receptor = c("Ccr1", "Ccr2", "Ccr3", "Ccr4", "Ccr5", "Ccr6", "Ccr7", "Ccr8", "Ccrl1", "Cxcr3", "Cxcr4", "Cxcr5", "Cxcr6", "Ackr3", "Cxcr1", "Cxcr2"))

draw_ggplot_from_gprofiler <- function(gp_result, save_path, cutoff, data_name, width = 8, height = 6) {
    gp_result[['log_p']] <- -log10(gp_result$p_value)
    gp_result$term_id <- factor(gp_result$term_id, levels = unique(gp_result$term_id))
    gp_result$term_name <- factor(gp_result$term_name, levels = unique(gp_result$term_name))

    # query is categorical value
    # fill the bars by query
    p <- ggplot(gp_result, aes(x = term_id, y = log_p, fill = query)) +
        geom_bar(stat = "identity", width = 0.7) +
        scale_fill_manual(values = c('PND7' = "darkred", 'PND1' = "darkblue")) +
        labs(title = "", x = "pathway", y = "-log10(p value)") +
        theme(
            plot.title = element_text(hjust = 0.5),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank()
        )
    p <- p + theme_classic() + coord_flip()
    ggsave(p, file= paste0(save_path, paste0('gprofiler2_result_', cutoff, '_', data_name, '.png')), width = width, height = height)
    return(p)
}

### chi-square test

# 각 유전자별로 카이제곱 검정 수행하는 함수
perform_chisq_test <- function(observed_matrix, cell_numbers) {
  # 결과를 저장할 리스트 초기화
  results <- list()
  
  # 각 유전자에 대해 반복
  for(gene in colnames(observed_matrix)) {
    # 관찰값
    observed <- observed_matrix[, gene]
    names(observed) <- rownames(observed_matrix)
    if (sum(observed) != 0) {
      
      # 기대값 계산
      # cell number 비율에 따른 기대값
      expected <- sum(observed) * (cell_numbers / sum(cell_numbers))
      
      # 카이제곱 검정 수행
      chisq_test <- chisq.test(x = observed, p = cell_numbers/sum(cell_numbers))
      
      # 결과 저장
      results[[gene]] <- list(
        statistic = chisq_test$statistic,
        p_value = chisq_test$p.value,
        observed = observed,
        expected = expected
      )
    }
    else {
      results[[gene]] <- list(
        statistic = NA,
        p_value = 1,
        observed = observed,
        expected = observed
      )
    }

  }
  return(results)
}


# 다중검정 보정을 위한 함수
adjust_p_values <- function(results) {
  p_values <- sapply(results, function(x) x$p_value)
  
  adjusted_p <- p.adjust(p_values, method = "BH")
  
  for(gene in names(results)) {
    results[[gene]]$adjusted_p_value <- adjusted_p[gene]
  }
  
  return(results)
}

sort_results_by_pvalue <- function(results) {
  # 각 유전자의 정보를 데이터프레임으로 변환
  summary_df <- data.frame(
    gene = names(results),
    pvalue = sapply(results, function(x) x$p_value),
    adj_pvalue = sapply(results, function(x) x$adjusted_p_value),
    chi_stat = sapply(results, function(x) x$statistic),
    max_contribution = sapply(results, function(x) max(x$contribution)),
    max_contrib_celltype = sapply(results, function(x) {
      names(x$contribution)[which.max(x$contribution)]
    }),
    direction = sapply(results, function(x) {
      # 가장 큰 contribution을 보인 celltype의 fold change 방향
      max_contrib_idx <- which.max(x$contribution)
      if(x$fold_change[max_contrib_idx] > 0) "Up" else "Down"
    }),
    fold_change = sapply(results, function(x) {
      # 가장 큰 contribution을 보인 celltype의 fold change 값
      max_contrib_idx <- which.max(x$contribution)
      x$fold_change[max_contrib_idx]
    }), 
    log2_odds_ratio = sapply(results, function(x) {
     # celltype이 2개인지 확인
     if(length(x$observed) != 2) {
       NA
     }
     # log2 Odds ratio 계산 
     # OR = (obs1/exp1)/(obs2/exp2)
     else {
        log2((x$observed[2]/x$expected[2]) / (x$observed[1]/x$expected[1]))
     }}
  ))
  
  # p-value로 정렬
  summary_df <- summary_df[order(summary_df$adj_pvalue), ]
  
  return(summary_df)
}

# Function to normalize UMI counts across cells in single-cell RNA sequencing data
#
# This function performs random sampling of UMIs to normalize the total UMI count
# per cell to a target value. For cells with fewer UMIs than the target,
# it performs oversampling with replacement. For cells with more UMIs than
# the target, it performs downsampling without replacement.
#
# @param matrix A genes x cells matrix where each element represents UMI counts
# @param target_umi Integer specifying the desired number of UMIs per cell (default: 10000)
# @return A normalized matrix with exactly target_umi counts per cell
# @note Empty cells (zero UMIs) are preserved as is
# @examples
# normalized_matrix <- normalize_umi(raw_matrix, target_umi = 10000)

normalize_umi <- function(matrix, target_umi = 10000) {
  # Calculate total UMIs per cell
  col_sums <- colSums(matrix)
  
  # Initialize result matrix
  result <- matrix(0, nrow = nrow(matrix), ncol = ncol(matrix))
  rownames(result) <- rownames(matrix)
  colnames(result) <- colnames(matrix)
  
  for(i in 1:ncol(matrix)) {
    current_umi <- col_sums[i]
    
    if(current_umi == 0) {
      next  # Skip empty cells
    }
    
    # Calculate sampling probabilities for each gene    
    probs <- matrix[,i] / current_umi
    probs <- probs[probs != 0]

    sampled_counts <- table(sample(
      1:length(probs),
      size = target_umi,
      prob = probs,
      replace = TRUE))
    # if(length(probs) < target_umi) {
    #   sampled_counts <- table(sample(
    #     1:length(probs),
    #     size = target_umi,
    #     prob = probs,
    #     replace = TRUE
    #   ))
    # } else {
    #   sampled_counts <- table(sample(
    #     1:length(probs),
    #     size = target_umi,
    #     prob = probs,
    #     replace = FALSE
    #   ))
    # }
    
    # Store results
    result[names(probs)[as.numeric(names(sampled_counts))], i] <- sampled_counts
  }
  return(result)
}

create_summary_matrix <- function(adjusted_results) {
  # Get all genes (keys from the list)
  genes <- names(adjusted_results)
  
  # Get condition names from the first gene's expected values
  conditions <- names(adjusted_results[[1]]$expected)
  
  # Create column names
  col_names <- c(
    "statistic",
    "p_value",
    "adjusted_p_value",
    paste0(conditions, "_observed"),
    paste0(conditions, "_expected")
  )
  
  # Initialize matrix
  result_matrix <- matrix(
    nrow = length(genes),
    ncol = length(col_names),
    dimnames = list(genes, col_names)
  )
  
  # Fill matrix
  for (i in seq_along(genes)) {
    gene <- genes[i]
    result_matrix[i, "adjusted_p_value"] <- adjusted_results[[gene]]$adjusted_p_value
    result_matrix[i, "p_value"] <- adjusted_results[[gene]]$p_value
    result_matrix[i, "statistic"] <- adjusted_results[[gene]]$statistic
    result_matrix[i, paste0(conditions, "_observed")] <- adjusted_results[[gene]]$observed
    result_matrix[i, paste0(conditions, "_expected")] <- unname(adjusted_results[[gene]]$expected)
  }
  
  return(result_matrix)
}

combine_results <- function(result_list) {
  library(metap)
  genes <- rownames(result_list[[1]])
  columns <- colnames(result_list[[1]])
  p_value_col <- "adjusted_p_value"
  other_cols <- setdiff(columns, p_value_col)
  
  # Combine p-values using selected method
  p_values <- sapply(genes, function(gene) {
    p <- sapply(result_list, function(x) x[gene, p_value_col])
    # if (sum(p == 1) >= 98) {
    #   p_meta = list(fisher = 1, stuffer = 1)
    # } else if (sum(p == 0) >= 98) {
    #   p_meta = list(fisher = 0, stuffer = 0)
    # } else {
    #   # Extract actual p-values from metap objects
    #   p_fisher <- sumlog(p)$p
    #   p_stuffer <- sumz(p)$p
    #   p_meta = list(fisher = p_fisher, stuffer = p_stuffer)
    # }
    # return(p_meta)
    return(median(p))
  })
  
  # Create separate vectors for fisher and stouffer p-values
  # p_fisher <- p_values['fisher', ] %>% unlist
  # p_stuffer <- p_values['stuffer', ] %>% unlist
  
  # Calculate means for other columns
  mean_values <- sapply(other_cols, function(col) {
    sapply(genes, function(gene) {
      mean(sapply(result_list, function(x) x[gene, col]))
    })
  })
  
  # Combine into data frame
  result <- data.frame(
    # p_fisher = p_fisher,
    # p_stuffer = p_stuffer,
    p_median = p_values,
    mean_values,
    row.names = genes
  )
  
  return(result)
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


condition_dpi_levels = c('naive_floxed', '7dpi_floxed', '14dpi_floxed', '30dpi_floxed', 
                                                                'naive_dAT2', '7dpi_dAT2', '14dpi_dAT2', '30dpi_dAT2')

########## contour plot from somi Kim from LCB, postech ############

## so : seurat object
## name_of_dimension : name of dimension that contains coordinates
## condition_name : name of metadata that contains condition information
## ref_condition : name of reference condition (ex : wt)
## exp_condition : name of experimental condition (ex : ko)
## celltype_name : name of metadata that contains celltype information
## nbin : number of bins for density plot
## linewidth : width of the contour line
## type_color : color for each condition
## celltype_color : color for each celltype
## save_path : path to save the plot
## data_name : name of the data
## width : width of the plot
## height : height of the plot
contour_plot <- function(so, name_of_dimension = 'umap', condition_name, ref_condition, exp_condition, celltype_name, nbin = 10, linewidth = 0.25,
                        type_color, celltype_color, save_path, data_name, width, height) {

  # get umap coordinates, condition, and celltype to build densdf
  x <- so@reductions[name_of_dimension][[1]]@cell.embeddings[, 1]
  y <- so@reductions[name_of_dimension][[1]]@cell.embeddings[, 2]
  condition <- so@meta.data[, condition_name]
  celltype <- so@meta.data[, celltype_name]
  densdf <- data.frame(x = x, y = y, type = condition, celltype = celltype)

  p_ref <- ggplot(densdf, aes(x=x, y=y)) +  
      # First add the point layers
      geom_point(data=densdf[densdf$type != ref_condition,], color = 'lightgray', alpha=1) +
      geom_point(data=densdf[densdf$type == ref_condition,], aes(color = celltype), alpha=0.75) +
      scale_color_manual(values = celltype_color) +
      
      # Then add the contour lines (without fill)
      geom_density_2d(data=densdf[densdf$type == ref_condition,], color = type_color[ref_condition], bins=nbin, linewidth  = linewidth) + 
      
      # Set the plot limits
      xlim(min(densdf$x)-1, max(densdf$x)+1) + 
      ylim(min(densdf$y)-1, max(densdf$y)+1) + 
      theme(legend.position = "none",
          text = element_text(size=15),
          panel.background = element_blank(),
          plot.background = element_blank(),
          panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          strip.background = element_rect(fill="lightgray", colour=NA))
  ggsave(p_ref, file = paste0(save_path, data_name, '_', ref_condition, '_contour.png'), width = width, height = height)
  ggsave(p_ref, file = paste0(save_path, data_name, '_', ref_condition, '_contour.pdf'), width = width, height = height)
  ggplot2pptx(p_ref, width, height, paste0(save_path, data_name, '_', ref_condition, '_contour.pptx'))

  p_exp <- ggplot(densdf, aes(x=x, y=y)) +  
      # First add the point layers
      geom_point(data=densdf[densdf$type != exp_condition,], color = 'lightgray', alpha=1) +
      geom_point(data=densdf[densdf$type == exp_condition,], aes(color = celltype), alpha=0.75) +
      scale_color_manual(values = celltype_color) + 
      
      # Then add the contour lines (without fill)
      geom_density_2d(data=densdf[densdf$type == exp_condition,], color = type_color[exp_condition], bins=nbin, linewidth  = linewidth) + 
      
      # Set the plot limits
      xlim(min(densdf$x)-1, max(densdf$x)+1) + 
      ylim(min(densdf$y)-1, max(densdf$y)+1) + 
      theme(legend.position = "none",
          text = element_text(size=15),
          panel.background = element_blank(),
          plot.background = element_blank(),
          panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          strip.background = element_rect(fill="lightgray", colour=NA))
  ggsave(p_exp, file = paste0(save_path, data_name, '_', exp_condition, '_contour.png'), width = width, height = height)
  ggsave(p_exp, file = paste0(save_path, data_name, '_', exp_condition, '_contour.pdf'), width = width, height = height)
  ggplot2pptx(p_exp, width, height, paste0(save_path, data_name, '_', exp_condition, '_contour.pptx'))
}

PND.cols = c('#b00202', '#d04f43', '#e9837e', '#fab7b7', 
        '#6402b1', '#7c33c1', '#9353d1', '#a971e1', '#c08ff0', '#d6adff', 
        '#013707', '#165f1d', '#008a05', '#29b32e', '#61ea66',
        '#6f4a00', '#9c6721', '#cc8640', '#ffa561')
names(PND.cols) = c(c('AT1', 'AT2', 'Ciliated cell', 'Club cell', 
                  'EC.Aplnr', 'EC.artery', 'EC.Car4', 'EC.lymphatic', 'EC.vein', 'EC.prolif.',
                  'Alv.fib', 'Myo.fib', 'Mural cell', 'Mesen.progenitor', 'Mesen.prolif.',
                  'B cell', 'T cell', 'NK cell', 'Myeloid cell'))
