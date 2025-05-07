
library(officer)
library(rvg)
library(cowplot)
# save ggplot2 object as ppt
save_ppt <- function(plotObj, fig_save, fig_name) {
    read_pptx() |> # build ppt
    add_slide() |> # add slide
    ph_with( # add image
        dml(ggobj = plotObj), 
        location = ph_location_fullsize() 
    ) |>
    print(paste0(fig_save, fig_name, '.pptx')) # save ppt
}


library(officer)
library(rvg)
ggplot2pptx <- function(p, width, height, filename){
  
  editable_graph <- dml(ggobj = p)
  doc <- read_pptx()
  doc <- add_slide(doc)
  doc <- ph_with(x = doc, editable_graph,
                 location = ph_location(width = width, height = height))
  print(doc, target = filename)
}

library(grid)
library(gtable)



extract_legend <- function(ggplot_obj) {
  tmp <- ggplot_build(ggplot_obj)
  
  if (!is.null(tmp$plot$scales$get_scales("color"))) {
    color_scale <- tmp$plot$scales$get_scales("color")
    scale_colors <- unique(tmp$data[[1]]$colour)
    scale_labels <- if (!is.null(color_scale$labels)) color_scale$labels else levels(tmp$data[[1]]$group)
    scale_name <- if (!is.null(color_scale$name)) color_scale$name else "color"
  } else {
    scale_colors <- NULL
    scale_labels <- NULL
    scale_name <- NULL
  }
  
  if (!is.null(scale_colors)) {
    dummy_data <- data.frame(
      x = rep(1, length(scale_colors)),
      y = rep(1, length(scale_colors)),
      group = factor(seq_along(scale_colors))
    )
  } else {
    dummy_data <- data.frame(x = 1, y = 1, group = factor(1))
  }
  
  legend_plot <- ggplot() +
    theme_void()
  
  if (!is.null(scale_colors)) {
    legend_plot <- legend_plot +
      geom_point(
        data = dummy_data,
        aes(x = x, y = y, color = group),
        show.legend = TRUE
      ) +
      scale_color_manual(
        values = scale_colors,
        labels = scale_labels,
        name = scale_name
      )
  }
  
  legend_plot <- legend_plot +
    theme(
      legend.position = "center",
      plot.background = element_blank(),
      panel.background = element_blank(),
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank()
    )
    
  return(legend_plot)
}
