# ----------------------- general ---------------------------------

# # increase output width
# options(width = 120)
# # print warnings as they occur
# options(warn = 1)
# # java heap size
# options(java.parameters = "-Xmx8G")
# 
# # get scripts directory (directory of this file) and load relevant functions
# args_all = commandArgs(trailingOnly = FALSE)
# scripts_dir = normalizePath(dirname(sub("^--file=", "", args_all[grep("^--file=", args_all)])))
# 
# # enter source statments for functions here once ready
# source(paste0(scripts_dir, "/load-install-packages.R"))
# source(paste0(scripts_dir, "/deseq2-pca.R"))
# source(paste0(scripts_dir, "/deseq2-compare.R"))
# source(paste0(scripts_dir, "/plot-volcano.R"))
# source(paste0(scripts_dir, "/plot-heatmap.R"))

source('~/Dropbox/Ghedin_Lab/DARPA_Seq/R_files/functions/heatmap-AR.R')
source('~/Dropbox/Ghedin_Lab/DARPA_Seq/R_files/functions/deseq2_analysis.R')

# # relevent arguments, only take arguments after the --args flag
# args = commandArgs(trailingOnly = TRUE)
# counts_table_file = args[1]
# col_data_file = args[2]
# 
# # check to make sure there are 2 arguments provided
# if (length(args) < 2) {
#   stop("Must provide a count table and a metadata file")
# }
# 
# # check that input files exist
# if (!file.exists(counts_table_file)) stop("file does not exist: ", counts_table_file)
# if (!file.exists(col_data_file)) stop("file does not exist: ", col_data_file)

# ------------------ Creating directory and loading libraries -----------------------

deseq_start <- function(count_table, column_data, r_dir, group_1 = 1, group_2 = 2, compare = "Mock", returnDDS = FALSE, ...) {

    # Load necessary packages:

    # For DESeq2 analysis
    library("DESeq2")

    # For data maniputation and saving
    library("RColorBrewer")
    library("glue")
    library("readr")
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

    # ----------------- Input importing --------------------

    # import counts table
    setwd("~/Dropbox/Ghedin_Lab/DARPA_Seq")
    counts_table = read.csv(file = glue("./count_Tables/processed_countdata/{count_table}"), row.names = 1, 
                            stringsAsFactors = FALSE, check.names = FALSE)
    message("input counts table gene num:      ", nrow(counts_table))
    message("input counts table sample num:    ", ncol(counts_table))
    message("input counts table sample names:  ", toString(colnames(counts_table)))
    message("")

    # import groups table
    col_data = read.csv(file = glue("./count_Tables/processed_coldata/{column_data}"), 
                        header = TRUE, row.names = 1, colClasses = "factor")
    message("metadata row num:   ", nrow(col_data))
    message("metadata condition names: ", toString(rownames(col_data)))
    message("metadata sample names:  ", toString(colnames(col_data)))
    message("")

    #readline("Enter file prefix for directory name and saving (note: no spaces): ")
    if (!dir.exists(glue("~/Dropbox/Ghedin_Lab/DARPA_Seq/DESeq_runs/{r_dir}"))) {
      dir.create(glue("~/Dropbox/Ghedin_Lab/DARPA_Seq/DESeq_runs/{r_dir}"))
    }
    setwd(glue("~/Dropbox/Ghedin_Lab/DARPA_Seq/DESeq_runs/{r_dir}"))

    # check that all samples from the column data are found in the counts table
    diff_samples = setdiff(rownames(col_data), colnames(counts_table))
    if (length(diff_samples)) {
      stop("some samples not in counts table: ", toString(diff_samples))
    }

    # subset to samples in groups table (also sets samples to be in the same order)
    # allows for use of the same count table with different, smaller/subsetted column data files
    counts_table = counts_table[, rownames(col_data)]

    # group info (use the first column for grouped comparisons)
    loc1 <- grep(group_1, colnames(col_data))
    loc2 <- grep(group_2, colnames(col_data))
    col_data[,loc1] <- as.factor(col_data[,loc1])
    group_a = colnames(col_data)[loc1]
    # also grab the second column just in case
    group_b = colnames(col_data)[loc2]
    message("group name: ", group_a)
    message("secondary group name: ", group_b)
    # find the unique levels in the group data (ie. the different viruses to be compared)
    group_levels = levels(col_data[, group_a])
    gb_levels = levels(col_data[, group_b])

    message("group levels: ", toString(group_levels))
    message("secondary group levels: ", toString(gb_levels))
    message("")
  
    contrast <- c(group_a, toString(group_levels[2]), toString(group_levels[1]))

    # design formula - what should DESeq ultimately compare between?
    # ie. compare between viruses
    design_formula = formula(glue("~ {group_a}"))
    message("design formula: ", design_formula)
    # need to also provide the virus that each should be compared to

    # ------------------ DESeq ----------------------

    # import raw counts and create DESeq object
    # since v1.16 (11/2016), betaPrior is set to FALSE and shrunken LFCs are obtained afterwards using lfcShrink
    dds = DESeqDataSetFromMatrix(countData = counts_table, colData = col_data, design = design_formula)
    dds$Virus <- relevel(dds$Virus, ref = "Mock") # need to figure out how to generalize this

    dds = DESeq(dds, parallel = TRUE, BPPARAM = BiocParallel::MulticoreParam(workers = 4))

    # VST
    vsd = varianceStabilizingTransformation(dds, blind = TRUE)

    # export normalized counts
    norm_counts = counts(dds, normalized = TRUE) %>% 
      round(2) %>% 
      as.data.frame() %>% 
      rownames_to_column("gene")
    write_excel_csv(norm_counts, glue("{r_dir}_normcounts.csv"))
    write_xlsx(list(normalized_counts = norm_counts), path = glue("{r_dir}_normcounts.xlsx"))

    # export variance stabilized counts
    vsd_table = assay(vsd) %>% 
      round(2) %>% 
      as.data.frame() %>% 
      rownames_to_column("gene")
    write_excel_csv(vsd_table, path = glue("{r_dir}_vstcounts.csv"))

    # ----------------------- PCA ------------------------

    # PCA plot
    source('~/Dropbox/Ghedin_Lab/DARPA_Seq/R_files/functions/deseq2-pca.R')
    pca_plot = deseq2_pca(vsd, intgroup = c(group_a, group_b), ntop = 1000)
    save_plot(glue("{r_dir}_pca.pdf"), pca_plot, base_width = 8, base_height = 5, units = "in", dpi = 300)
    save_plot(glue("{r_dir}_pca.png"), pca_plot, base_width = 8, base_height = 5, units = "in", dpi = 300)
    message("Saving pca plot")

    # ---------------differential expression------------------

    # perform comparisons for all results in the results table
    source('~/Dropbox/Ghedin_Lab/DARPA_Seq/R_files/functions/deseq2_analysis.R')
    res_names <- as.list(resultsNames(dds))
    res_names <- res_names[-c(1)]
    for (result in seq_along(res_names)) {
      res <- gsub(pattern = "Virus_", replacement = "", x = res_names[[result]])
      message(glue("Comparison: ", res))
      deseq2_analysis(deseq_data = dds, name = res_names[[result]], contrast = contrast, save_plot = TRUE)
    }
    
    if (returnDDS == TRUE) {
        return(dds)
    }

}