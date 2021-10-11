# function for identification of the most variable genes in Seurat object, with the top genes threshold to set

hvg <- function(seurat_obj, top =500) {
    require (Seurat)
    require (matrixStats)
    
    # export normalized gene expression from Seurat object
    norm <- GetAssayData(object = seurat_obj, assay = "SCT", slot = "data")
    
    # calculate variance, create data frame and order it
    hvg <- rowVars (as.matrix(norm))
    hvg <- data.frame (gene = rownames(norm), var = hvg)
    hvg <- hvg [order (hvg$var, decreasing = TRUE),] [1:top,]
              
  return (hvg)
}
