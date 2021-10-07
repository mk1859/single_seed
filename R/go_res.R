# enrichment of GO terms from single category with calculation of parent GO terms
# requires vector of gene names, vector of background genes, category (BP, CC or MF) and thershold for finding of parent GO terms

go_res <- function(genes, background, p_value = 0.05, category = "BP", rrvgo_threshold=0.8) {
    require (gprofiler2)
    require (rrvgo)
    
    # enrichment calculation
    res <- gost(query = genes, organism = "athaliana", custom_bg = background, user_threshold = p_value, sources = "GO")$result
    
    # selection of GO term category
    category_selected <- res [grep(category, res$source),]
    
    # stop if no go terms at selected category
    if (nrow(category_selected) == 0) {return (print("no GO term at category"))}
    
    else {
    
        # identification of parent GO terms
        scores <- setNames(-log10(category_selected$p_value), category_selected$term_id)
        simMatrix <- calculateSimMatrix(category_selected$term_id, orgdb="org.At.tair.db", ont=category, method="Rel")
        reducedTerms <- reduceSimMatrix(simMatrix, scores, threshold=rrvgo_threshold, orgdb="org.At.tair.db")
        reducedTerms <- reducedTerms [,c(1,2,3,8)]
        colnames (reducedTerms) [1] <- "term_id"
        
        # combine to single data frame
        res <- merge (category_selected [,c(9,10,3,4,5,6,12,11)], reducedTerms, by = "term_id")
    
        return (res)
    }
}
