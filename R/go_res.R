# enrichment of GO terms from a single ontology with identification of parent GO terms
# it requires a vector of gene names, vector of background genes names, ontology (BP, CC or MF) and threshold for finding parent GO terms

go_res <- function(genes, background, p_value = 0.05, category = "BP", rrvgo_threshold=0.8) {
    require (gprofiler2)
    require (rrvgo)
    require (org.At.tair.db)
    
    # enrichment calculation
    res <- gost(query = genes, organism = "athaliana", custom_bg = background, user_threshold = p_value, sources = "GO")$result
    
    # selection of GO term ontology
    category_selected <- res [grep(category, res$source),]
    
    # stop if no go terms at selected ontology
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
