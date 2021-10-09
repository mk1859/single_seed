# function to filter genes for each treatment to be above some expression treshold
# it outputs list which elements are gene name vectors

select_genes <- function (matrix, treatments, avg_reads = 1) {
  
  # choose seeds
  selected_seeds <- lapply (treatments, function (x) {
    which(gsub('.{0,6}$', '', colnames(matrix)) == x) })
    
  # select genes
  selected_genes <- lapply (selected_seeds, function (x) {
    matrix [rowMeans(matrix[ ,x]) > avg_reads, x]})
    
  names(selected_genes) <- treatments
  
  return (selected_genes)
}
