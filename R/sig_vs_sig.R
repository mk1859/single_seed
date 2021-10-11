# signature vs signature plot
# it requires specifying column with treatment names, their order, it is possible to exclude one of the treatments

sig_vs_sig <- function (seurat_obj, signature1, signature2, excluded = NULL, order, column) {
  require (Seurat)
  require (ggplot2)
  
  # export data from the Seurat object
  plot <- seurat_obj@meta.data
  
  # set order of treatments
  plot$timepoint <- factor(plot$timepoint, levels = order)
  
  # exclude treatment if necessary
  if (!is.null(excluded)){
    plot <-  plot [-grep( excluded, rownames(plot)),]
  }
  
  g <- ggplot(plot, aes(x=.data[[signature1]], y= .data[[signature2]], color = .data[[column]])) +
          geom_point(size = 2) +
          theme_classic() +
          scale_color_tableau()
  
  return (g)

}
