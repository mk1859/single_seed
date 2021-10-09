# function for creating Seurat objects for list of treatments 
# with batch removal and background regression duing sctransform normalization

list_seurat <- function (selected_genes, background) {
  require (Seurat)
  require (sctransform)

  # list of seed attributes, take treatment and batch information from seed name
  timepoint_attr <- lapply (selected_genes, function (x) { 
                            data.frame(batch = gsub('.{0,3}$', '', colnames(x)),
                            background = background [names(background) %in% colnames(x)] ) })
                            
  # creating Seurat objects
  timepoint_seurats <- mapply (function (x,y) {
          CreateSeuratObject(counts = x, meta.data = y)},
          matrix_list,timepoint_attr)
  # normalizing Seurat objects
  timepoint_seurats <- lapply (timepoint_seurats, function (x) {
                          Seurat::SCTransform(object = x, verbose =FALSE, 
                                              n_genes = NULL, 
                                              batch_var = "batch",
                                                method = "nb", 
                                              vars.to.regress ='background',
                                              return.only.var.genes = TRUE)})
  # running PCA
  timepoint_seurats <- lapply (timepoint_seurats, function (x) {
  RunPCA(object = x, verbose = FALSE)})
  
  return (timepoint_seurats)
  
}
