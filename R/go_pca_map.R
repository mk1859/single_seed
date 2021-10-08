# function to plot expression of genes belonging to a given GO term on the PCA map
# it uses Seurat object, GO term id, order of treatments to plot on violin plot and name of column with treatments, it also allows to exclude some treatments
# it uses other custom function called signature_map

go_pca_map <-  function (seurat_obj, go_id, excluded = NULL, order, column) {
  require (Seurat)
  require (ggplot2)
  require (org.At.tair.db)
  require (biomaRt)
  require (GO.db)
  require (ggthemes)
  require (viridis)
  
  # export genes included in single seed matrix
  genes <- rownames(GetAssayData(seurat_obj, assay = "SCT", slot = "data"))
  
  # download GO terms for Arabidopsis
  ensembl <- useMart(biomart="plants_mart", dataset = "athaliana_eg_gene", host="plants.ensembl.org")
  
  # assign genes to GO terms
  go <- getBM(attributes = c("ensembl_gene_id", "go_id"), filters = "ensembl_gene_id", values = genes, mart = ensembl)
  
  # pick GO term you are interested in
  go <- go$ensembl_gene_id [go$go_id == go_id]
  
  # create signature using Seurat function
  seurat_obj <- AddModuleScore(seurat_obj, features = list(go), name = "go_term")
  
  # plot signature
  g <- signature_map (seurat_obj, signature = "go_term1", excluded = NULL, order = order , column = column)

  return (g)
}
