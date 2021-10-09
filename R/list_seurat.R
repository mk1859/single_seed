# function for creating Seurat objects for list of treatments 
# with batch removal and background regression duing sctransform normalization

list_seurat <- function (matrix_list, background) {
  require (Seurat)
  require (sctransform)

  # list of seed attributes, take batch information from seed name
  timepoint_attr <- lapply (matrix_list, function (x) { 
                            data.frame(batch = gsub('.{0,3}$', '', colnames(x)),
                            background = background [names(background) %in% colnames(x)] ) })
                            
  # create list of Seurat objects
  timepoint_seurats <- mapply (function (x,y) {
                          CreateSeuratObject(counts = x, meta.data = y)},
                       matrix_list,timepoint_attr)
  
  # normalize Seurat objects
  timepoint_seurats <- lapply (timepoint_seurats, function (x) {
                          Seurat::SCTransform(object = x, verbose =FALSE, 
                                              n_genes = NULL, 
                                              batch_var = "batch",
                                              method = "nb", 
                                              vars.to.regress ='background',
                                              return.only.var.genes = TRUE)})
  # run PCA
  timepoint_seurats <- lapply (timepoint_seurats, function (x) {
                          RunPCA(object = x, verbose = FALSE)})
  
  return (timepoint_seurats) 
}
