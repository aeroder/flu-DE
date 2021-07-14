miRNA_targets <- function(miR_target_table, miR, cutoff = 0.90) {

  # For data maniputation and saving
  library("RColorBrewer")
  library("glue")
  library("readr")
  library("readxl")
  library("writexl")
  library("tidyr")
  library("plyr")
  library("dplyr")
  library("tibble")
  library("stringr")

  # For graphing and dataviz
  library("VennDiagram")
  library("biobroom")
  library("ggplot2")
  library("ggfortify")
  library("pheatmap")

  targets <- read.csv(glue("./{miR_target_table}"), stringsAsFactors = FALSE)

  target_table <- subset(file, miR > cutoff)

  target_table$tx.ens.id <- str_replace_all(target_table$tx.ens.id, "\\..", "")

  results <- getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id"),
                   filters = "ensembl_transcript_id",
                   values = target_table$tx.ens.id, mart = mart)

  ordered <- match(results$ensembl_transcript_id, target_table$tx.ens.id)

  target_table <- target_table[match(results$ensembl_transcript_id, target_table$tx.ens.id),]
  target_table$ensg.id <- results$ensembl_gene_id

  return(target_table)

}
