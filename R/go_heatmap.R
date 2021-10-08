# function to create heatmap of enriched GO terms for multiple or single element list
# it allow to choose if term name or id will appear in the plot (term_name = TRUE/FALSE)
# and if parent GO terms and term ontology should be plotted

go_heatmap <- function(res_file, term_name =FALSE, term_category = TRUE, parent_term = TRUE, single_condition = FALSE) {
  require (cowplot)
  require (ggplot2)
  require (ggthemes)
  
  if (single_condition == TRUE) {
    plot <- res_file [order(-log10(res_file$value)),]
    plot$variable <- "x"
    
  } else { # multiple lists
  
  plot <- res_file [order(res_file$category, res_file$parent, -log10(res_file [,ncol(res_file)])),]
  order_go <- as.factor(plot$term_id)
  plot <- melt(plot)
  plot$term_id <- factor(plot$term_id,order_go)
  
  }
  
  if (term_name == TRUE) { # term name or term id on the plot
    g1 <- ggplot(plot, aes(variable, term_name)) + 
            geom_tile(aes(fill = -log10(value))) + 
            scale_fill_gradient2(low = "midnightblue", 
                       mid = "white", 
                       high = "darkred") +
            theme_classic ()
  
  } else {
    g1 <- ggplot(plot, aes(variable, term_id)) + 
            geom_tile(aes(fill = -log10(value))) + 
            scale_fill_gradient2(low = "midnightblue", 
                       mid = "white", 
                       high = "darkred") +
            theme_classic ()
  }
  
  g2 <- ggplot(plot, aes(1, term_id)) + 
            geom_tile(aes(fill = as.factor(category))) + 
            theme_classic () +
            theme(axis.line=element_blank(),
                  axis.text.x=element_blank(),
                  axis.text.y=element_blank(),
                  axis.ticks=element_blank(),
                  axis.title.x=element_blank(),
                  axis.title.y=element_blank()) +
            labs (fill = "ontology")
    
  
  g3 <- ggplot(plot, aes(1, term_id, fill = as.factor(parent))) + 
            geom_tile() + 
            scale_fill_tableau("Tableau 20") +
            theme_classic () +
            theme(axis.line=element_blank(),
                  axis.text.x=element_blank(),
                  axis.text.y=element_blank(),
                  axis.ticks=element_blank(),
                  axis.title.x=element_blank(),
                  axis.title.y=element_blank()) +
            labs (fill = "parent")

   # create multipanel plot
  if (term_category == TRUE & parent_term == TRUE) {
    g <- plot_grid(g3, g2, g1, nrow = 1, align = "h", 
               rel_widths = c(0.75,0.5,1.5))
    return (g)
  }
  
  if (term_category == TRUE & parent_term == FALSE) {
    g <- plot_grid(g2, g1, nrow = 1, align = "h", 
               rel_widths = c(0.5,1.5))
    return (g)
  }
  
  if (term_category == FALSE & parent_term == TRUE) {
    g <- plot_grid(g3, g1, nrow = 1, align = "h", 
               rel_widths = c(0.75,1.5))
    return (g)
  }
  
  if (term_category == FALSE & parent_term == FALSE) {
    return (g1)
  }
}
