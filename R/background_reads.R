# function that requires both raw and pre-filtered matrices as input and returns fraction of intergenic reads per seed from pre-filtered matrix

background_reads <- function(matrix, prefiltered) {

  #count background reads 
  background_reads <- matrix  %>% 
  filter(rownames(matrix) %in% 
           c(Araport$gene 
              [c(!(Araport$type == "protein_coding"),                      # not encoding proteins
                    grep ("ATC", rownames (matrix)),                       # chloroplast
                    grep ("ATM", rownames (matrix)))],                     # mitochondrial
                    "__no_feature"))  %>%                                  # intergenic
    colSums()
  
  # calculate fraction of intergenic reads
  background_reads <- background_reads/ colSums(matrix)
  
  # filter for seeds remainig after prefiltering
  background_reads <- background_reads [names (background_reads)%in% colnames(prefiltered)]
                           
  return (background_reads)
}
