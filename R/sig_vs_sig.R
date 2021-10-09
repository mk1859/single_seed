# signature vs signature plot
# it is possible to exclude some treatments and requires adding order of treatments on the legend
# coloumn of treatment to color the point need to be specify

sig_vs_sig <- function (seurat_obj, signature1, signature2, excluded = NULL, order, column) {
  require (Seurat)
  require (ggplot2)
  
  # export data from Seurat object
  plot <- seurat_obj@meta.data
  
  # set order of treatments
  plot$timepoint <- factor(plot$timepoint, levels = order)
  
  # exclude some time points if necessary
  if (!is.null(excluded)){
    plot <-  plot [-grep( excluded, rownames(plot)),]
  }
  
  g <- ggplot(plot, aes(x=.data[[signature1]], y= .data[[signature2]], color = .data[[column]])) +
          geom_point(size = 2) +
          theme_classic() +
          scale_color_tableau()
  
  return (g)

}
