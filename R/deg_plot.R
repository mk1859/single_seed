# function to create plot showing number of changed genes in list of comparisons
# direction = TRUE divides genes to upregulated and downregulated 
# limits and by is passed to x axis 

deg_plot <- function(deg_list, direction = TRUE, limits = c(-400,800), by = 200) {
  require (rlist)
  require (ggthemes)
  
  # show separately up- and downregulated genes
  if (direction == TRUE) { 
    
    # check how many genes are up- and downregulated
    degs <- lapply (deg_list, function(x) {
                    as.data.frame(table(x$avg_log2FC > 0))})
    
    names (degs) <- names (deg_list) 
    degs <- list.rbind(degs)
    degs$comparison <- factor (gsub('.{0,2}$', '', rownames(degs)), levels = rev(names (deg_list)))
    degs$Freq [c(T,F)] <- degs$Freq [c(T,F)] * -1
    
    g <- ggplot (degs, aes (x = comparison, y = Freq, fill = Var1)) +
                geom_bar(stat='identity') +
                theme_classic() +
                coord_flip() +
                scale_fill_manual(name="change", 
                    labels = c("downregulated", "upregulated"), 
                    values = c("blue", "red")) +
                scale_y_continuous("genes", breaks= seq (limits [1], limits [2], by), limits = limits)
     
    # all affected genes
    } else { 
   
      # check how many genes are affected
      degs <- lapply (deg_list, function(x) {
                  nrow(x) })
    
      names (degs) <- names (deg_list)
      degs <- as.data.frame(degs)
      degs$comparison = factor (rownames(degs), levels = rev(rownames(degs)))
      
      g <- ggplot (degs, aes (x = comparison, y = degs, fill= comparison)) +
                geom_bar(stat='identity') +
                theme_classic() +
                coord_flip() +
                scale_fill_tableau() +
                scale_y_continuous("genes", breaks= seq (limits [1], limits [2], by), limits = limits)
    }
  
 return (g)
}
