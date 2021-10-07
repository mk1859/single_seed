# analysis of multiple gene lists for GO term enrichment
# when selected category = "all" it wrapps multiple_category.R
# when one of BP, CC or MF selected it wraps go_res.R

multiple_genelists <- function(gene_lists, background, p_value = 0.05, category = "all", rrvgo_threshold=0.9) {
    require (gprofiler2)
    require (rrvgo)
    require (org.At.tair.db)
    require (rlist)

    # multiple categories
    if (category == "all") { 
 
    all_list <- lapply (gene_lists, function (x){ 
          multiple_category (x, background = background, p_value = p_value, rrvgo_threshold = rrvgo_threshold)
         })
    } 
    # selected category
    else { 
    
    all_list <- lapply (gene_lists, function (x){ 
          go_res (x, background = background, p_value = p_value, rrvgo_threshold = rrvgo_threshold)
         })
    }
 
  # combine to data frame
  all_list <- Reduce(function(x, y, ...) merge(x, y, all = TRUE, by = "term_id",suffixes = names (gene_lists),...), all_list)
    
  # data frame require cleaning if list was longer than 1
    
  if (length(gene_lists) > 1) {
  
    # function for cleaning columns of data frame
    clean <- function (filt_frame,term) {
                 sel_column <- filt_frame [,grep(term, colnames(filt_frame))]
                 column_list <- apply (sel_column, 2, function (x) x )
                 coalesce(!!!column_list)
                }
    
    # create new clean data frame
    all_res <- data.frame (term_id = all_list$term_id,
                           term_name = clean(all_list, "term_name"),
                           category = clean(all_list, "source"),
                           parent = clean(all_list, "parentTerm"),
                           p_value = all_list [,"p_value"])
                           
    
    colnames(all_res) [-c(1:4)] <- names (gene_lists)
  
    return (all_res)
  
  } else { #single list
    all_res <- data.frame (term_id = all_list$term_id,
                         term_name = all_list$term_name,
                         category = all_list$source,
                         parent = all_list$parentTerm,
                         value = as.numeric(all_list$p_value))
  return (all_res)
    
  }
  
}
