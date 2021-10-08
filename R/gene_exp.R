# function to create PCA plot with normalized gene expression and violin plot inset showing expression in indicated treatments

gene_exp <- function(seurat_obj, gene, order, column) {
  require (cowplot)
  require (ggplot2)
  require (ggthemes)
  require (viridis)
  
  # export normalized data and filter for gene
  gene_exp <- as.data.frame(GetAssayData (seurat_obj, slot = "data")) 
  gene_exp <- t(gene_exp [which(rownames(gene_exp) == gene),])
  
  # add PCA coordinates and metadata
  gene_exp <- cbind(gene_exp, as.data.frame (Embeddings(object = seurat_obj, reduction = "pca")) [,1:2], seurat_obj@meta.data)       
  
  # set order of treatments
  gene_exp [,eval(column)] <- factor(gene_exp [,eval(column)], levels = order)
  
  # PCA map
  g1 <- ggplot(gene_exp, aes(x=PC_1, y= PC_2, color =	.data[[gene]])) +
          geom_point (size = 2) + 
          scale_colour_viridis(option="C") +
          theme_classic() +
          theme(strip.background = element_blank(),
              axis.title = element_blank(),
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              axis.line = element_blank())
  
  # violin plot for time points
  g2 <- ggplot(gene_exp, aes(x=.data[[column]], y= .data[[gene]], fill = .data[[column]], color = .data[[column]])) +
          geom_violin () + 
          scale_fill_tableau() +
          scale_color_tableau() +
          theme_classic() +
          theme(legend.position = "none",
              strip.background = element_blank(),
              axis.title = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.line.x = element_blank())
  
  # combine plots
  g <- ggdraw() +
          draw_plot(g1) +
          draw_plot(g2, x = 0, y = .7, width = .25, height = .3)
  
  return (g)
}
