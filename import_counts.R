# R function to upload set of matrices contatining read counts (with the same set of genes).
# It assumes each matrix have gene names as row names.


import_counts <- function(
                          dir, #
header = TRUE) {
  require (rlist)
  files <- list.files (dir)
  matrices <- lapply (files, function (x){
    read.csv (paste0(dir, x), sep = "\t", header = header, row.names =1)
  })
  matrices <- list.cbind(matrices)
  return (matrices)
}
