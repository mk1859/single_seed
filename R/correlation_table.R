# function to calculate correlation of gene count to fraction of background reads
# it takes matrix of counts, normalizes sequencing depth per seed using Seurat 
# then it exptracts information about treatment from seed name 
# finally calculates Pearson correlation to the fraction o background reads for seeds in each traetmant and for all seeds combine

correlation_table <- function(matrix, background) {
  require (Seurat)
  
  #input for normalization
  seed_attr <- data.frame(background = background,                                       # fraction of background
                          timepoint = as.factor(gsub('.{0,6}$', '', colnames(matrix))))  # extract information about treatment from seed name
                          
  # normalize using Seurat
  norm_seeds <- NormalizeData(CreateSeuratObject(counts = matrix, meta.data = seed_attr), verbose = FALSE)
  
  #export normalized data
  norm_seeds <- as.matrix(GetAssayData(object = norm_seeds, slot = "data"))
  
  # correlation of gene's reads with background fraction for each treatment
  timepoints <- levels(seed_attr$timepoint)
  
  background_corr <- lapply (timepoints, function (x) { 
                            apply (norm_seeds [,which(seed_attr$timepoint == x)],1, 
                            function(z) cor.test(z,seed_attr$background
                            [which(seed_attr$timepoint == x)])$estimate)
  })
  
  # global correlation of gene's UMIs with background percent
  background_corr$all <- apply (norm_seeds,1, function(x) {
                                cor.test(x,seed_attr$background)$estimate})
  
  # data frame from the list
  background_corr <- t(matrix(unlist(background_corr), nrow=length(background_corr), byrow=TRUE))
  colnames (background_corr) <- c (timepoints, "all")
  rownames(background_corr) <- rownames (matrix)
                           
  return (background_corr)
}
