library(officer)
library(rvg)
ggplot2pptx <- function(p, w, h, filename){
  
  editable_graph <- dml(ggobj = p)
  doc <- read_pptx()
  doc <- add_slide(doc)
  doc <- ph_with(x = doc, editable_graph,
                 location = ph_location(width = w, height = h))
  print(doc, target = filename)
}

building_df <- function(so, set_name, control_cells, exp_cells, annotation_name, threshold = 0, is_gene) {
    # initialize vectors
    fc_vector = c(); p_value_vector = c(); subtract_vector = c()
    control_mean_vector = c(); experiment_mean_vector = c()
    control_cell_number = c(); experiment_cell_number = c()
    threshold_status = c()  # New vector to track threshold status

    for (celltype in so[[annotation_name]][[1]] %>% levels) {
        celltype_cells = Cells(so)[so[[annotation_name]][[1]] == celltype]
        control_and_celltype_cells = intersect(control_cells, celltype_cells)
        exp_and_celltype_cells = intersect(exp_cells, celltype_cells)
        control_cell_number = c(control_cell_number, length(control_and_celltype_cells))
        experiment_cell_number = c(experiment_cell_number, length(exp_and_celltype_cells))
        
        if (length(control_and_celltype_cells) > threshold & length(exp_and_celltype_cells) > threshold) {
            threshold_status = c(threshold_status, "above_threshold")
            if (is_gene) {
                con_original = expm1(so@meta.data[control_and_celltype_cells, paste0(set_name, '1')])
                exp_original = expm1(so@meta.data[exp_and_celltype_cells, paste0(set_name, '1')])
                con_mean = mean(con_original)
                exp_mean = mean(exp_original)
                log_con_mean = log(con_mean + 1)
                log_exp_mean = log(exp_mean + 1)
                fc_vector = c(fc_vector, log_exp_mean / log_con_mean)
                p_value_vector = c(p_value_vector, wilcox.test(con_original, exp_original)$p.value)
                subtract_vector = c(subtract_vector, log_exp_mean - log_con_mean)
                control_mean_vector = c(control_mean_vector, log_con_mean)
                experiment_mean_vector = c(experiment_mean_vector, log_exp_mean)
            }
            else {
                fc_vector = c(fc_vector, mean(so@meta.data[exp_and_celltype_cells, paste0(set_name, '1')]) / mean(so@meta.data[control_and_celltype_cells, paste0(set_name, '1')]))
                p_value_vector = c(p_value_vector, wilcox.test(so@meta.data[exp_and_celltype_cells, paste0(set_name, '1')], so@meta.data[control_and_celltype_cells, paste0(set_name, '1')])$p.value)
                subtract_vector = c(subtract_vector, mean(so@meta.data[exp_and_celltype_cells, paste0(set_name, '1')]) - mean(so@meta.data[control_and_celltype_cells, paste0(set_name, '1')]))
                control_mean_vector = c(control_mean_vector, mean(so@meta.data[control_and_celltype_cells, paste0(set_name, '1')]))
                experiment_mean_vector = c(experiment_mean_vector, mean(so@meta.data[exp_and_celltype_cells, paste0(set_name, '1')]))
            }
        }
        else {
            threshold_status = c(threshold_status, "below_threshold")
            fc_vector = c(fc_vector, 1)  # FC of 1 means no change
            p_value_vector = c(p_value_vector, 1)  # p-value of 1 means not significant
            control_mean_vector = c(control_mean_vector, 0)
            experiment_mean_vector = c(experiment_mean_vector, 0)
            subtract_vector = c(subtract_vector, 0)
        }
    }

    df = data.frame(celltype = so[[annotation_name]][[1]] %>% levels,
                   control_cell_number = control_cell_number,
                   subtract_mean = subtract_vector,
                   experiment_cell_number = experiment_cell_number,
                   control_mean = control_mean_vector,
                   experiment_mean = experiment_mean_vector,
                   fc = fc_vector,
                   p_value = p_value_vector,
                   threshold_status = threshold_status)
    return(df)
}

draw_celltype_wise_dotplot <- function(so, set_name, features, condition_name, control, experiment, annotation_name, save_path, cols = c('darkred', 'darkblue', 'darkgray', 'black'), facet_group = NULL, facet_levels = NULL, threshold = 0, width = 6, height = 4) {
    epsilon = 1e-6
    is_gene = FALSE
    if (length(features) == 1) {
        if(features %in% rownames(so)) {
            so[[paste0(set_name, '1')]] = as.vector(so@assays$RNA$data[intersect(rownames(so), features), ])
            is_gene = TRUE
        }
        else {
            so[[paste0(set_name, '1')]] = so[[features]]
        }
    }
    if (length(features) > 1) {
        so <- AddModuleScore(so, features = list(features), name = set_name)
    } 

    control_cells = Cells(so)[so[[condition_name]][[1]] == control]
    exp_cells = Cells(so)[so[[condition_name]][[1]] == experiment]

    df <- building_df(so, set_name, control_cells, exp_cells, annotation_name, threshold, is_gene)
    
    df[['avg_log2FC']] = log2(df[['fc']] + epsilon)
    df[['p_val_adj']] = p.adjust(df[['p_value']], method = 'BH')
    rownames(df) = df[['celltype']]

    # Modified change classification
    df[['change']] = ifelse(df[['threshold_status']] == "below_threshold", 
                           'Not detected',
                           ifelse(df[['subtract_mean']] > 0, 
                                  paste0('Increase with ', experiment),
                                  ifelse(df[['subtract_mean']] < 0,
                                        paste0('Decrease with ', experiment),
                                        'Not significant')))
    
    df[['change']][df[['threshold_status']] == "above_threshold" & df[['p_val_adj']] > 0.05] = 'Not significant'
    df[['celltype']] = factor(df[['celltype']], levels = so[[annotation_name]][[1]] %>% levels)

    # Add 'Lower than threshold' to color scheme
    names(cols) = c(paste0('Increase with ', experiment), 
                   paste0('Decrease with ', experiment), 
                   'Not significant',
                   'Not detected')
    # Add facet group
    if (!(is.null(facet_group))) {
        df$group = facet_group
        df$group = factor(df$group, levels = facet_levels)
    }
    

    p <- ggplot(df, aes(x = subtract_mean, y = celltype)) + 
        geom_point(aes(size = -log10(p_val_adj), color = change)) +
        scale_color_manual(values = cols) +
        theme_classic() + 
        geom_vline(xintercept = 0, color = 'black', linetype = 'dashed') +
        theme(axis.text.y = element_text(face="bold"),
              axis.title.y = element_text(size = 0)) +
        labs(x = 'Differential Score') + 
        guides(colour = guide_legend(override.aes = list(size=5)))

    if (is_gene) {
        p <- p + labs(x = 'Average log2FC')
    }
    
    # # Add facet group
    # if (!(is.null(facet_group))) {
    #     p <- p + facet_wrap(~group, ncol = 1)
    # }

    ggsave(p, file = paste0(save_path, 'dot_', set_name, '.png'), width = width, height = height)
    ggplot2pptx(p, 6, 4, paste0(save_path, 'dot_', set_name, '.pptx'))
    return(list(so, p, df))
}

#### facet example
# dpis = c('naive', ..., 'naive', '7dpi', ..., '7dpi', '14dpi', ..., '14dpi', '30dpi', ..., '30dpi')
# facet_levels = c('naive', '7dpi', '14dpi', '30dpi')
# returns = draw_celltype_wise_dotplot(so.alv, 'HALLMARK_INFLAMMATORY_RESPONSE', subset(HM_gene_sets, gs_name == 'HALLMARK_INFLAMMATORY_RESPONSE')$gene_symbol, 
#                             'condition', 'floxed', 'dAT2', 'annotation_dpi', 
#                             save_path, cols = c('darkred', 'darkblue', 'darkgray', 'black'), facet_group = dpis, facet_levels = levels(so.alv$dpi), threshold = 0, width = 6, height = 12)
# ggsave(returns[[2]] + coord_flip() + theme(axis.text.x = element_text(angle = 90, hjust = 1))+ scale_y_discrete(labels = c(plot_annotations)) + facet_wrap(~group, scales = "free_x", ncol = 4, nrow = 1),
#     filename = paste0(save_path, 'alv_inflammatory_score.png'), width = 6, height = 5)
