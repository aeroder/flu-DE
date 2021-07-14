# Rscript for filtering for low coverage areas of the genome and calculating % coverage:
percentcov <- function(input_file, file_prefix, cutoff = 10, ...) {

    # load libraries in for data maniputation and graphing
    library("RColorBrewer")
    library("glue")
    library("readr")
    library("readxl")
    library("writexl")
    library("tidyr")
    library("dplyr")
    library("tibble")
    library("stringr")
    library('pipeR')
    library('rlist')

    # read in the count table 
    cov_counts <- read.csv(input_file)
    cov_counts <- as.data.frame(cov_counts)

    # find the unique samples in the list
    samples <- as.character(unique(cov_counts$name))

    # initialize a list to put each seperate samples coverage table into
    sample_tables <- list()

    # filter the full count table and make seperate tables for each sample
    for (i in samples) {
      temp_table <- cov_counts %>%
        filter(name == i)
      sample_tables[[length(sample_tables) + 1]] <- temp_table
    }
    names(sample_tables) <- samples

    # clean-up
    rm(i, temp_table)

    # initialize a data frame with column names sample and coverage to put the % coverage into
    percent_coverage <- data.frame("sample"=character(0), "coverage"=numeric(0), stringsAsFactors = FALSE)

    # calculate the percentage of rows < a certain cutoff 
    for (t in 1:length(sample_tables)) {
      # total number of rows
      total <- nrow(sample_tables[[t]])
      # total number of rows with count less than cutoff
      # change the last number in this row if you want a different cutoff
      zero <- nrow(filter(sample_tables[[t]], totalcount <= cutoff))
      # calculate the percent
      percent = 100-(zero/total*100)
      tempname = names(sample_tables[t])
      # populate the % coverage table
      percent_coverage[nrow(percent_coverage) + 1,] = c(tempname,percent)
    }

    # clean-up
    rm(t, total, zero, tempname, percent)

    res_csv <- glue("{file_prefix}_{cutoff}_percentcov.csv")
    write_excel_csv(percent_coverage, path = res_csv)
}
