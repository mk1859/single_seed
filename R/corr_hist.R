# function to create histograms, it takes correlation table, order of the time points and threshold of correlation as input

corr_hist <- function(corr_table, order, threshold = 0.3) {
  require (reshape2)
  require (ggthemes)
  require (ggplot2)
  
  table_long <- melt (corr_table)                            # the long form of the table
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
          geom_vline(aes(xintercept= threshold ), color="red", size=1)

  return (g)
}
