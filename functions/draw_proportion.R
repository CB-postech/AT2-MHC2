 make.proportion.to.plot <- function(so, row_meta, col_meta, colors = "", add_blackline = FALSE, flip = FALSE) {
  library(magrittr)
  library(Seurat)
  Idents(so) <- so[[row_meta]][[1]]
  
  print("This function makes a proportion dataframe.")
  print("This will calculate the proportion of coloumn contents in each row contents")
  print("For example, let be row_contents is cluster of Seurat object and col_contents is Normal or Tumor")
  print("Then return df has proportion of normal cell proportion and tumor cell proportion in each cluster")
  print("PLEASE SET ORDER OF CONTENTS BY LEVELS!")
  
  row_contents <- levels(so@meta.data[, row_meta])
  col_contents <- levels(so@meta.data[, col_meta])
  print(row_contents)
  print(col_contents)

  proportion.sub.NorT.df <- as.data.frame(matrix(data = rep(0, row_contents %>% length * col_contents %>% length), nrow = row_contents %>% length, ncol = col_contents %>% length))
  rownames(proportion.sub.NorT.df) = as.character(row_contents)
  colnames(proportion.sub.NorT.df) = c(paste0(col_contents, rep("_proportion", col_contents %>% length)))
  # paste0(col_contents, rep("_proportion", col_contents %>% length)), 

  print('1')

  col_cells <- list()
  for (celltype in col_contents) {
    col_cells[[celltype]] <- Cells(so)[so[[col_meta]] == celltype]
  }
  print(col_cells)

  for (cluster in row_contents) {
    cells <- WhichCells(so, ident = cluster)
    intersect_cells <- list()
    for (celltype in col_contents) {
      intersect_cells[[celltype]] <- intersect(cells, col_cells[[celltype]])
      # proportion.sub.NorT.df[cluster, paste0(celltype, '_proportion')] = intersect_cells[[celltype]] %>% length
      proportion.sub.NorT.df[as.character(cluster), paste0(celltype, '_proportion')] = intersect_cells[[celltype]] %>% length / cells %>% length
    }
  }  
  print('2')
  print(proportion.sub.NorT.df)

  cluster.for.plot <- c()
  proportion.for.plot <- c()
  fill.for.plot <- c()
  for (cluster in row_contents) {
    cluster.for.plot <- c(cluster.for.plot, rep(cluster, col_contents %>% length))
    proportion.for.plot <- c(proportion.for.plot, proportion.sub.NorT.df[cluster, ] %>% unlist %>% as.vector)
    fill.for.plot <- c(fill.for.plot, col_contents %>% as.character)
  }
  print(cluster.for.plot)
  print('3')

  data.frame(cluster.for.plot, proportion.for.plot, fill.for.plot) -> sub.plot.df
  print(sub.plot.df)
  colnames(sub.plot.df) <- c('cluster', 'proportion', 'subtype')

  print('4')

  p <- ggplot(sub.plot.df, aes(fill=factor(subtype, levels = col_contents),x=factor(cluster, levels = row_contents), y=proportion)) + geom_bar(position ='fill', stat='identity')
  p <- p + xlab(row_meta) + ylab(col_meta) + theme_classic() + labs(fill = col_meta)
  # p <- p + theme(aspect.ratio = 2)
  if (add_blackline) {
    p <- p + geom_col(color="black", size = 0.25)
  }
  p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  if (flip) {
    p <- p + coord_flip()
    p <- p + theme(axis.text.x = element_text(angle = 0))
  }
  p <- p + labs(y = 'Proportion')
  if (length(colors) != 1 || colors != "") {
    p <- p + scale_fill_manual(values = colors)
  }
  # plot.margin = unit(c(1, 1, 1, 1), "cm")
  # put margin 1cm to top, bottom, left, right
  # aspect.ratio -> controls the aspect ratio of the plot.
  # ggsave(p, file = paste0(row_meta, '_', col_meta, '.pdf'))
  # x=factor(cluster, level = order.used %>% as.character)
  return(list(sub.plot.df, p, proportion.sub.NorT.df)) 
}
# Msln, Upk3b가 DP이여야 mesothelial cell -> AT1과 같이 clustering이 안됨.

### example code
# source('/home/sjcho/yard/functions/R/make_proportion_to_plot.R')
# Idents(so.epi) <- so.epi[['Condition']][[1]]
# levels(so.epi) <- c('Control', 'Exercise')
# so.epi[['Condition']] <- Idents(so.epi)

# Idents(so.epi) <- so.epi[['marker_annotation2']][[1]]
# levels(so.epi) <- c('AT1', 'DP_AT2_AT1', 'AT2', 'DP_AT2_club', 'club', 'mucous', 'ciliated')
# so.epi[['marker_annotation2']] <- Idents(so.epi)

# Idents(so.epi) <- so.epi[['Sample']][[1]]
# levels(so.epi) <- c('con1', 'con2','con3', 'ex1', 'ex2', 'ex5')
# so.epi[['Sample']] <- Idents(so.epi)

# cols = cols = c("#07105D", "#2D87B1", "#0CA55D", "#42E512", "#A20707", "#DD8215", "#8C30F3")
# p <- make.proportion.to.plot(so.epi, row_meta = 'marker_annotation2', col_meta = 'Condition', add_blackline = TRUE, colors = cols)
# ggsave(p[[2]], file = paste0(fig_save, 'compositional_change.png'), width = 7, height = 7)
