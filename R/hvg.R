# function for identification of most varaible genes in sourat object, with number of top genes to set

hvg <- function(seurat_obj, top =500) {
    require (Seurat)
    require (matrixStats)
    
    # export gnormalized gene ezpression form Seurat object
    norm <- GetAssayData(object = seurat_obj, assay = "SCT", slot = "data")
    
    # callculate variance, crate data frame and order it
    hvg <- rowVars (as.matrix(norm))
    hvg <- data.frame (gene = rownames(norm), var = hvg)
    hvg <- hvg [order (hvg$var, decreasing = TRUE),] [1:top,]
              
  return (hvg)
}
