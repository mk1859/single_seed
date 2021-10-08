# function to create histograms, takes correlation table, order of time points and treshold of correlation as input
corr_hist <- function(corr_table, order, treshold = 0.3) {
  require (reshape2)
  require (ggthemes)
  require (ggplot2)
  
  table_long <- melt (corr_table)                            # long form of table
  table_long$Var2 <- factor(table_long$Var2, levels = order) # correct order

g <- ggplot(table_long, aes(x = value, color = Var2)) +
      geom_histogram (fill = "white", breaks = seq (-.8, .8, .1)) +
      xlab("Correlation with background %") +
      scale_color_tableau() +
      theme_classic() +
      xlim (-1,1)+
      facet_wrap(vars(Var2), nrow = 2) +
      theme(strip.background = element_blank(),
            legend.position = "none") +
      geom_vline(aes(xintercept= treshold ), color="red", size=1)

  return (g)
}
