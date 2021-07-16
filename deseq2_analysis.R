
deseq2_analysis <- function(deseq_data, contrast = NULL, name = NULL, r_dir, ...) {

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

  # create a directory to store the output files
  # r_dir <- getwd()
  if(!is.null(r_dir)) {
    if (!dir.exists(r_dir)) dir.create(r_dir)
  }

  # Make a directory to store heatmaps in
  heatmaps_dir = "heatmaps"
  if (!dir.exists(glue("{r_dir}/{heatmaps_dir}"))) {
    dir.create(glue("{r_dir}/{heatmaps_dir}"))
  }

  # Make a directory to store volcano plots in
  volcano_dir = "volcano-plots"
  if (!dir.exists(glue("{r_dir}/{volcano_dir}"))) {
    dir.create(glue("{r_dir}/{volcano_dir}"))
  }

  # Make a directory to store venns diagrams in
  vennD_dir = "venndiagrams"
  if (!dir.exists(glue("{r_dir}/{vennD_dir}"))) {
    dir.create(glue("{r_dir}/{vennD_dir}"))
  }

  # generate a results table based on the supplied name.
  # If no name is given, stop the analysis and return an error
  if(!is.null(name)) {
    res <- results(deseq_data, name = name, cooksCutoff = FALSE)
    res_name <- name
  } else {
    stop("Please supply a name for the results table to be generated.")
  }

  # DESeq results name format : "Virus_FluB_vs_Mock
  # get the file suffix based on comparison name
  file_suffix = gsub(pattern = "Virus_", replacement = "", x = res_name)
  file_suffix = gsub(pattern = "_vs_Mock", replacement = "", x = file_suffix)

  # save unmodified results object
  res_rds = glue("deseq2_res_{file_suffix}.rds")
  saveRDS(res, file = glue("{r_dir}/{res_rds}"))
  message("saving results object: ", res_rds)

  # save unmodified results as csv
  res_df <- as.data.frame(res) %>% rownames_to_column("gene")
  res_csv <- glue("{file_suffix}_res.csv")
  write_excel_csv(res_df, path = glue("{r_dir}/{res_csv}"))
  message("saving results csv: ", res_csv)

  # format results and save as excel spreadsheet
  res_df <- res_df %>%
    dplyr::mutate(
      baseMean       = round(baseMean, 1),
      log2FoldChange = round(log2FoldChange, 2),
      pvalue         = round(pvalue, 15),
      padj           = round(padj, 15)
    ) %>%
    dplyr::rename(
      log2FC = log2FoldChange
    )

    res_df <- res_df %>%
      dplyr::select(gene, baseMean, log2FC, pvalue, padj)

    res_xlsx <- glue("{file_suffix}_res.xlsx")
    write_xlsx(setNames(list(res_df), strtrim(res_name, 31)), path = glue("{r_dir}/{res_xlsx}"))
    message("saving results as excel file: ", res_xlsx)

    # generate volcano plot
    source('~/National Institutes of Health/Ghedin, Elodie (NIH NIAID) [E] - LAB_STUFF/allison/projects/DARPA_seq/R_files/flu-DE/plot-volcano_YH.R')
    plot_volcano(stats_df = res_df,
                 positive_group=contrast[2],
                 negative_group=contrast[3],
                 gene_col = "gene",
                 fc_col = "log2FC",
                 p_col = "padj",
                 p_cutoff = 0.05,
                 fc_label = "Fold Change (log2)", p_label = "Adjusted P Value (-log10)",
                 title = res_name,
                 file_prefix = glue("{r_dir}/{volcano_dir}/plot.volcano2.{file_suffix}"))

    message("# genes padj < 0.1: ", nrow(subset(res_df, padj < 0.1)))
    message("# genes padj < 0.05: ", nrow(subset(res_df, padj < 0.05)))
    message("# genes padj < 0.01: ", nrow(subset(res_df, padj < 0.01)))
    message("# genes padj < 0.001: ", nrow(subset(res_df, padj < 0.001)))

    # save differential expression results in Excel format
    res_q005_xlsx <- gsub(pattern = ".xlsx", replacement = "_q005.xlsx", x = res_xlsx)
    res_q005_df <- subset(res_df, padj < 0.05)
    write_xlsx(setNames(list(res_q005_df), strtrim(res_name, 31)), path = glue("{r_dir}/{res_q005_xlsx}"))
    message("saving DE genes results as excel file: ", res_q005_xlsx)

    # perform vst on data for use in heatmap
    vst <- assay(varianceStabilizingTransformation(deseq_data, blind = TRUE))
    samples_all = colnames(deseq_data)
    betas <- coef(deseq_data)

    # heatmap gene subsets (list with genes, plot title, and file suffix)
    # Heatmap code written by Yuhan!
    hmg = list()
    hmg[[length(hmg) + 1]] = list(genes = res_df %>% dplyr::filter(padj < 0.05) %>% .$gene,
                                  title = "q < 0.05",
                                  file_suffix = "q005")
    message("list generated for q < 0.05: ", length(hmg[[length(hmg)]]))
    hmg[[length(hmg) + 1]] = list(genes = res_df %>% dplyr::filter(padj < 0.01) %>% .$gene,
                                  title = "q < 0.01",
                                  file_suffix = "q001")
    hmg[[length(hmg) + 1]] = list(genes = res_df %>% dplyr::filter(padj < 0.001) %>% .$gene,
                                  title = "q < 0.001",
                                  file_suffix = "q0001")
    hmg[[length(hmg) + 1]] = list(genes = res_df %>% dplyr::filter(pvalue < 0.05) %>% .$gene,
                                  title = "p < 0.05",
                                  file_suffix = "p005")
    hmg[[length(hmg) + 1]] = list(genes = res_df %>% dplyr::filter(pvalue < 0.01) %>% .$gene,
                                  title = "p < 0.01",
                                  file_suffix = "p001")

    # process every gene subset
    for (i in 1:length(hmg)) {

      # generate title and file suffix
      hm_title = glue("{res_name}_{hmg[[i]]$title}")
      hm_file_prefix = glue("{r_dir}/{heatmaps_dir}/hm_{file_suffix}_{hmg[[i]]$file_suffix}")

      # generate heatmaps if gene list is not too small or big
      if (length(hmg[[i]]$genes) > 10 && length(hmg[[i]]$genes) < 10000) {

        mat <- betas[hmg[[i]]$genes,-c(1)]
        # generate heatmap using all samples

        heatmapAR(mat = mat,
                  row_subset = hmg[[i]]$genes,
                  title = hm_title,
                  file_prefix = hm_file_prefix,
                  save_plot = TRUE)

      }
    }


}
