# function to create area plot showing number of genes showing correlation to the background reads
corr_number <- function(corr_table, treshold = 0.3) {
  require (ggthemes)
  require (ggplot2)

  # count genes with treatments above threshold of correlation
  corr_sum <- apply (corr_table,2, function (x) ifelse (x > treshold , 1, 0) )
  
  # number of treatments
  cor_stat <- rowSums(corr_sum [,1:(ncol(corr_sum)-1)])
  
  # number of genes for all seeds
  all_stat <- table (corr_sum [,ncol (corr_sum)] > 0 & cor_stat == 0) [2]
  names(all_stat) <- "all"
  
  # create table for plotting
  cor_stat <- as.data.frame(c (table (as.factor(cor_stat)), all_stat)) %>%
    mutate(., treatments = rownames(.), pos1 = 0, pos2 =0)
    
  cor_stat [1,1] <- cor_stat [1,1] - all_stat # remove genes only in all seeds
  colnames (cor_stat) [1] <- "Freq"
  
  # loop for settings borders of the areas on the plot, starting from 0 it sets borders according to the number of affected genes
  a <- cor_stat$Freq [1]
  b <- 0
  for (i in 1:nrow(cor_stat)) {
    cor_stat$pos1 [i] <- b
    cor_stat$pos2 [i] <- a
    b <- a
    if(i == nrow(cor_stat)) {
      break
    }
    a <- a + cor_stat$Freq [i + 1]
  }
  
  g <- ggplot(cor_stat, aes(x=pos1, y = 1, xmin = pos1, xmax = pos2, ymin = 0, ymax =1, fill = treatments)) +
      scale_x_continuous(name="x") +
      scale_y_continuous(name="y") +
      geom_rect() +
      scale_fill_tableau("Miller Stone") +
      theme_classic() +
      theme(axis.line.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x =element_blank(),
            axis.title.x=element_blank()) +
      coord_flip()
  return (g)
}
