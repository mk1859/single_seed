# function to create a multipanel diagram for the fraction of background reads, it takes the order of treatments as input
background_plot <- function(matrix, order, background) {
  require (cowplot)
  require (ggplot2)
  require (ggthemes)
  
  # metadata for seeds
  seed_attr <- data.frame(n_reads = colSums(matrix),                                    # number of reads
                          background = background,                                      # fraction of background reads
                          timepoint = as.factor(gsub('.{0,6}$', '', colnames(matrix)))) # extract information about treatment from seed name
                        
  
  # set order of treatments
  seed_attr$timepoint <- factor(seed_attr$timepoint, levels = order)
  
  # plot for number of reads and fraction of background reads
  g1 <- ggplot(seed_attr, aes(x=log10(n_reads), y=background, color = timepoint)) +
   geom_point(size = 1) + 
   scale_color_tableau() +
   theme_classic() +
   geom_smooth(method='lm', se=FALSE) +
   theme(legend.position = "none")
  
  # boxplot for fraction of background reads
  g2 <- ggplot(seed_attr, aes(x=timepoint, y=background, color = timepoint)) +
   geom_boxplot() + 
   scale_color_tableau() +
   theme_classic() +
   theme(axis.line=element_blank(),
         axis.text.x=element_blank(),
         axis.text.y=element_blank(),
         axis.ticks=element_blank(),
         axis.title.x=element_blank(),
         axis.title.y=element_blank()) +
   theme(legend.position = "none")

  # create a multipanel plot
  g <- plot_grid(g1, g2, ncol = 2, align = "h", 
               rel_widths = c(2,0.8))
               
  return (g)
}
