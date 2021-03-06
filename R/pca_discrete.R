# PCA plot with discrete values coloring

pca_discrete <- function (seurat_obj, type = "timepoint", tableu ="Tableau 10", excluded = NULL) {
  require (Seurat)
  require (ggplot2)
  require (ggthemes)
  
  # export data from the Seurat object
  plot <- cbind (as.data.frame (Embeddings(object = seurat_obj, reduction = "pca"))  [,1:2], seurat_obj@meta.data)
  
  # exclude treatment if necessary
  if (!is.null(excluded)){
    plot <-  plot [-grep( excluded, rownames(plot)),]
  }
  
  g <- ggplot(plot, aes(x=PC_1, y= PC_2, color = .data[[type]])) +
            geom_point (size = 2) + 
            scale_color_tableau(tableu) +
            theme_classic() +
            theme(axis.line=element_blank(),
                  axis.text.x=element_blank(),
                  axis.text.y=element_blank(),
                  axis.ticks=element_blank(),
                  axis.title.x=element_blank(),
                  axis.title.y=element_blank()) +
            theme(legend.title = element_blank())
  
  return (g)
}
