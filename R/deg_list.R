# it wrapper function to Seurat FindMarkers that allows to compare gene expression pairwise between conditions named in two vectors
# function generates list as out put which elements contain data frames with significant DEGs

deg_list <- function(seurat_obj, vector1, vector2, padj = 0.05, log2FC_threshold = log2(2)) {
  require (Seurat)
  
  # set
  seurat_obj@active.ident <- as.factor(seurat_obj$timepoint)
  deg_list <- list ()
  for (i in 1:length (vector1)){
  deg_list [[i]] <- FindMarkers(object = seurat_obj, ident.1 = vector2 [i],
                      ident.2 = vector1 [i],logfc.threshold = log2FC_threshold, 
                      test.use = "wilcox",only.pos =FALSE,
                      assay = "SCT", slot ="data") %>%
                    filter(.,p_val_adj < padj)
  }
  names(deg_list) <- paste0(vector1, "_", vector2)
  
  return (deg_list)
}
