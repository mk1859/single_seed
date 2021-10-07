# analysis of multiple GO categories by warpping go_res.R function on BP, CC and MF ontology categories

multiple_category <- function(genes, background, p_value = 0.05, rrvgo_threshold=0.8) {
    require (gprofiler2)
    require (rrvgo)
    require (org.At.tair.db)
    require (rlist)
    
    # create list with ontologies
    sum_res <- list (BP = NULL, CC = NULL, MF = NULL)
    
    # use go_res for each ontology
    sum_res$BP <- go_res (genes = genes, background = background, p_value = p_value, category = "BP", rrvgo_threshold=rrvgo_threshold)
    sum_res$CC <- go_res (genes = genes, background = background, p_value = p_value, category = "CC", rrvgo_threshold=rrvgo_threshold)
    sum_res$MF <- go_res (genes = genes, background = background, p_value = p_value, category = "MF", rrvgo_threshold=rrvgo_threshold)
    
    # combine to single data frame
    sum_res <- list.rbind(sum_res)
    
    # identify and remove ontologies without terms
    empty_category <- grep("no GO term", sum_res [,1])
    
    if (length(empty_category) > 0 ) {
      sum_res <- sum_res [-grep("no GO term", sum_res [,1]),] 
    }
    
    return (sum_res)
}
