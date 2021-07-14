# load files in 
# 
# need folder name
# prefixes for each file name 
# load into a list

nextSteps <- function(r_dir, virus=NULL) {
 
    # For DESeq2 analysis
    library("DESeq2")
  
    # For data maniputation and saving
    library("RColorBrewer")
    library("glue")
    library("readr")
    library("readxl")
    library("writexl")
    library("tidyr")
    library("dplyr")
    library("tibble")
  
    # For graphing and dataviz
    library("VennDiagram")
    library("biobroom")
    library("ggplot2")
    library("ggfortify")
    library("pheatmap") 
  
    # set the working director to the input above where the results tables are
    setwd(glue("~/Dropbox/Ghedin_Lab/DARPA_seq/DESeq_runs/{r_dir}"))
        
    # file_prefixes stores the virus names which are the names that the res tables are saved under
    file_prefixes <- c("FluB","H1N1pdm09", "H3N2", 
                       "PR8","FSS13025","MR766","P6740","PRVABC59")
    # if an argument is supplied, divide this into just the flu strains or just the zika strain
    if (virus=="flu") {
        file_prefixes <- c("FluB","H1N1", "H3N2", "PR8")
    } else if (virus=='flu2') {
        file_prefixes <- c("FluB","H1N1pdm09", "H3N2", "PR8")
    } else if (virus=="zika") {
        file_prefixes <- c("FSS13025","MR766","P6740","PRVABC59")
    }
        
    # Load in the full res_tables for each of the viruses indicated
    res_tables <- list()
    for (strain in file_prefixes) {
        res_tables[[length(res_tables) + 1]] <- as.data.frame(read_excel(glue("./{strain}_res.xlsx")))
    }
    names(res_tables) <- file_prefixes
        
    # Load in the q > 005 tables for the viruses indicated
    q_tables <- list()
    for (strain in file_prefixes) {
        q_tables[[length(q_tables) + 1]] <- as.data.frame(read_excel(glue("./{strain}_res_q005.xlsx")))
    }
    names(q_tables) <- file_prefixes
  
    #q_tables <- list()
    #for (strain in file_prefixes) {
    #    temp <- as.data.frame(read_excel(glue("./{strain}_res_q005.xlsx")))
    #    q_tables[[length(q_tables) + 1]] <- as.list(temp$gene)
    #}
        
    de <- lapply(q_tables, '[[', 'gene')
    de <- do.call(c, de)
    names(de) <- NULL  
    de <- unique(de) 

    my_list <- list(res_tables, q_tables, de)
    return(my_list)
    
}
                 