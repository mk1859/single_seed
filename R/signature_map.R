# function for signature plots with violin insets
# signatures should be added earlier to Seurat object with AddModuleScore
# it allows remove some of the treatments and set thier order

signature_map <- function (seurat_obj, signature, excluded = NULL, order, column) {
  require (Seurat)
  require (ggplot2)
  
  # export data from Seurat object
  plot <- cbind (as.data.frame (Embeddings(object = seurat_obj, reduction = "pca"))  [,1:2], seurat_obj@meta.data) 
  
  # set order of treatments
  plot$timepoint <- factor(plot$timepoint, levels = order)
  
  # exclude some time points if necessary
  if (!is.null(excluded )){
    plot <-  plot [-grep( excluded, rownames(plot)),]
  }
  
  # PCA plot
  g1 <- ggplot(plot, aes(x=PC_1, y= PC_2, color = .data[[signature]])) +
            geom_point (size = 2) + 
            scale_colour_viridis() +
            theme_classic () +
            theme(axis.line=element_blank(),
                  axis.text.x=element_blank(),
                  axis.text.y=element_blank(),
                  axis.ticks=element_blank(),
                  axis.title.x=element_blank(),
                  axis.title.y=element_blank()) 
                  
  # violin plot for treatments
  g2 <- ggplot(plot, aes(x=.data[[column]], y= .data[[signature]], fill = .data[[column]], color = .data[[column]])) +
          geom_violin () + 
          scale_fill_tableau() +
          scale_color_tableau() +
          theme_classic() +
          theme(legend.position = "none",
              strip.background = element_blank(),
              axis.title.x = element_blank())
  
  # combine plots
  g <- ggdraw() +
         draw_plot(g1) +
         draw_plot(g2, x = 0, y = .75, width = .25, height = .25)
    
  return (g)
}
