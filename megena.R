# Need to input count tables (normalized? or not) and filter them into separate viruses. All have mock and then separate the other viruses on their own
# Start by making col tables containing just each virus, similar to how we did for each time point

run_megena <- function(count_table_path, colname, prefix, ...) {

    library("glue")
    library("readr")
    library("writexl")
    library("tidyr")
    library("dplyr")
    library("tibble")
    library(MEGENA)
    
    counts_table = read.csv(file = glue("./DESeq_runs/{count_table_path}"), row.names = 1, 
                            stringsAsFactors = FALSE, check.names = FALSE)
    col_data = read.csv(file = glue("./count_Tables/megena_coldata/{colname}"), 
                        header = TRUE, row.names = 1, colClasses = "factor")
    
    # check that all samples from the column data are found in the counts table
    diff_samples = setdiff(rownames(col_data), colnames(counts_table))
    if (length(diff_samples)) {
        stop("some samples not in counts table: ", toString(diff_samples))
    }

    # subset to samples in groups table (also sets samples to be in the same order)
    # allows for use of the same count table with different, smaller/subsetted column data files
    counts_table = counts_table[, rownames(col_data)]
    
    names(counts_table) = glue("{names(counts_table)}_{prefix}")
    
    if (!dir.exists("~/Dropbox/Ghedin_Lab/DARPA_Seq/count_tables/megena_tables")) {
      dir.create("~/Dropbox/Ghedin_Lab/DARPA_Seq/count_tables/megena_tables")
    }
    setwd("~/Dropbox/Ghedin_Lab/DARPA_Seq/count_tables/megena_tables")
    
    file_prefix = gsub(pattern = "_cd.csv", replacement = "", x = colname)
    
    counts_table <- counts_table %>% rownames_to_column("rownames")

    write_excel_csv(counts_table, glue("{file_prefix}_megena.csv"))
       
}    

