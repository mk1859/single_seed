# function to find co-expressed genes in normalized data from Seurat object
# it uses correlatePairs from scran to calculate correlations in expresssion between gene pairs
# which next filter them by set treshold
# next step is cration of graph obcject in which genes are nodes and presence of correltion above treshold is an edge between them
# this graph is used to identify highly connected gene groups with parameters passed to RBGL::highlyConnSG
# finally gene groups are filttred to consist of some minimal number of genes

coexpressed <- function (seurat_obj, 
                          iters=1e8,             # number of iterations to calculate correlation
                          threshold,             # filtering significant correlation
                          sat=1,                 # paramter for highlyConnSG
                          ldv=c(80,40,20,10),    # paramter for highlyConnSG
                          n_genes=10) {             # minimal number of genes in identified groups
                          
  require (scran)
  require (scater)
  require (graph)
  require (RBGL)
  require (dplyr)
  
  # export genes from Seurat object and put counts to SingleCellExperiment
  genes <- as.matrix(GetAssayData(object = seurat_obj, assay = "SCT",slot = "data"))
  genes <- SingleCellExperiment(list(sctransform = genes))
  
  # gene expression correlation and filtering correlated pairs
  genes <- as.data.frame(correlatePairs(genes, iters= iters, assay.type="sctransform")) %>%
             filter(.,limited == TRUE & rho > threshold)                              
  
  # creating a graph from correlating genes
  g <- ftM2graphNEL(cbind(genes$gene1, genes$gene2), W=NULL, V=NULL, edgemode="undirected")
  
  # identification of highly connected gene groups
  clusters <- highlyConnSG(g, sat= sat, ldv=ldv)$clusters
  
  # ordering and filtering gene groups
  clusters <- clusters[order(lengths(clusters), decreasing=TRUE)]
  names(clusters) <- paste0("cluster_", seq(1:length(clusters)))
  clusters <- clusters [lengths(clusters) >= n_genes]
  
  return (clusters)
}
