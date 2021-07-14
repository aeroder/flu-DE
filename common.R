#
# Common
# Function: determine the number of common elements between two vectors in a list, usually called by VennsD.R
# 
# Input: degenes: A list of the genes to compare between each of the sets
#        genes: A full list of possible genes
#        return_list: A boolean value that when TRUE, will return the genes in common between all of the sets.
#                     Defaults to FALSE indicating that the function should create and save a plot
#        samp: A list of the comparisons to make that corresponds to locations within de_list
#
# Author: Allison Roder
# Last updated: 6/15/19
#

common <- function(degenes, genes, return_list = FALSE, samp, ...) {
  # input: 
  # 1. a list of DE genes from the conditions to be compared 
  # 2. a full list of the genes
  # 3. The locations of the conditions to be compared
  
  # set list of DE genes to de
  # set full list of genes to gene_names
  de <- as.list(degenes)
  gene_names <- genes
  
  # get just the list of genes by extracting the row names
  allgenes <- row.names(gene_names)
  
  # Set the total "in common" genes to the full list of genes
  comp <- allgenes
  
  # If only one DE list location is provided, output the number of genes in that list location
  if (length(samp) == 1) {
    comp <- de[[samp[[1]]]]
    
    # If two are given, find the intersect and set length to that number
  } else if (length(samp) == 2) {
    comp <- intersect(de[[samp[1]]], de[[samp[2]]])
    
    # if three are given, first find the intersect of two of them and then the intersect of those with the third
  } else if (length(samp) == 3) {
    comp <- intersect(de[[samp[1]]], 
                      intersect(de[[samp[2]]], 
                                de[[samp[3]]]))
    
    # same as for three but for four
  } else if (length(samp) == 4) {
    comp <- intersect(de[[samp[1]]], 
                      intersect(de[[samp[2]]], 
                                intersect(de[[samp[3]]], 
                                          de[[samp[4]]])))
  }
  # return the length of comp
  if (return_list == TRUE) {
    return(as.list(comp))
  } else {
    length(comp)
  }
  
  

}

