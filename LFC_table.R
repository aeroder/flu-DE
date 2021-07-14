#
# LFC Table
# Function: Input a list of results tables and a set of genes and output a dataframe containing 
#           just the log fold changes (or any column) from those dfs
#           output: a df containg the selected column from each of the input dataframes
#
# Author: AR
# Last updated: 6/28/19
#

LFCtable <- function(gene_subset, col_select = "log2FC", res_tables) {
  
  library(plyr)
  library(dplyr)

  # Initiate a data frame to be filled for the output
  out <- data.frame()
  
  for(i in 1:length(res_tables)) {
    
    # Make a temporary data frame called temp and set the current res_tables value to that
    temp <- res_tables[[i]]
    
    # Check to see if the desired column exists in each data frame, if not - stop and print an error
    if (!col_select %in% colnames(temp)) {
      stop(glue("Column name {col_select} does not exist in one of the results tables"))
    }
    
    # Set the row names (presumably the geneIDs) to be a column called 'gene'
    if (!("gene") %in% colnames(temp)) {
      temp <- rownames_to_column(temp, "gene")
    }
    
    # Filter the columns in the data frame to just select the column with the name given 
    # by the variable 'col_select' - default == "log2FC"
    temp <- dplyr::filter(temp, temp$gene %in% gene_subset)
    rownames(temp) <- temp$gene
    temp <- dplyr::select(temp, col_select)

    # if this is the first run through, set the output df to the values in temp
    if (i == 1) {
      out <- temp
    }
    # on any other run through, bind the column in temp to the existing data frame
    if (i > 1) {
      out <- cbind(out, temp)
    }
    
  }
  
  # Set the names of the DF columns to the names of the input list
  names(out) <- names(res_tables)
  # return the created data frame
  return(out)
  
}













