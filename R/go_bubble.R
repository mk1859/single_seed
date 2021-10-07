# function for GO term bubble plots, takes single list with single category results as input

go_bubble <- function(res_file_single) {
   require (ggplot2)
       
  # add fold change enrichemnet of the GO terms
  res_file_single$FC <- (res_file_single$intersection_size / res_file_single$query_size) /
                        (res_file_single$term_size / res_file_single$effective_domain_size)
         
  # add fraction of genes of the term
  res_file_single$frac <- res_file_single$intersection_size / res_file_single$term_size

  g <- ggplot(res_file_single, aes(x=FC, y= -log10(p_value), color = parentTerm, size = frac)) +
          geom_point (alpha =.5) + 
          scale_size(range = c(3, 12)) +
          scale_color_tableau()  +
          theme_classic ()
 
  return (g)
}
