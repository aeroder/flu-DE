#
# heatmapAR
# Function: Creates both a png and pdf heatmap using pheatmap of a supplied matrix
# 
# Input: mat: A matrix of values to be plotted in the heatmap
#        row_subset: A subset of the rows contained in mat that should be included (defaults to all rows)
#        col_subset: A subset of the columns contained in mat that should be included (defaults to all cols)
#        title: A title for the output graph
#        file_prefix: A prefix for the output files. Defaults to hmp_{title} unless changed
#        cluster_rows: Whether or not to cluster rows. Defulats to TRU
#        col_ann: A df containing the sample information, if desired, for annotating the columns (optional)
#        ann_colors: A list of lists to assign values in col_ann to specific colors in the key and annotations (optional)
#
# Author: Allison Roder
# updated: 6/15/19
#          7/19/19
#          8/2/19
#

heatmapAR <- function(mat, 
                      row_subset = NULL, 
                      col_subset = NULL,
                      title, 
                      file_prefix = "hmp_", 
                      cluster_rows = TRUE,
                      col_ann = NULL, 
                      ann_colors = NULL,
                      save_plot = FALSE,
                      ...) {
  
    library(glue)
    library(pheatmap)
  
    # confirm that the subsets exist in the full matrix
    if (!is.null(row_subset)) {
        row_subset = intersect(row_subset, rownames(mat))
        if (length(row_subset) < 2) {
            message("too few genes/rows")
        }
    }
    else if (!is.null(col_subset)) {
        col_subset = intersect(col_subset, colnames(mat))
        if (length(col_subset) < 2) {
            message("too few samples/columns")
        }
    } else {
    
    # edit the row names to more simplified versions of miRNA names
    # ie. remove everything after the first _
    rownames(mat) <- sub("\\_.*", "", rownames(mat))
  
    annotationcol <- colnames(mat)
  
    if (!is.null(col_subset)) {
    # edit the column names to just say the virus
    for (i in seq_along(annotationcol)) {
        temp <- gsub(pattern = "Virus_", replacement = "", x = annotationcol[[i]])
        temp <- gsub(pattern = "_vs_Mock", replacement = "", x = temp)
        annotationcol[[i]] <- temp
    }
    colnames(mat) <- annotationcol
    }

    title = glue("{title} - Fold Change over Mock")
    filename = NULL
    
    if (save_plot == TRUE) {
        if (!is.null(file_prefix)) {
        filename = glue("{file_prefix}.pdf")
        }
    }
    
    # Blue - Yellow - Red (might change this)
    cell_colors = colorRampPalette(c("#91bfdb", "#ffffbf", "#fc8d59"))(50)

    # adjust row font size based on number of rows
    fontsize_row = 4
    show_rownames = TRUE
    if ((nrow(mat)) < 80) {
      fontsize_row = 6
    }
    if ((nrow(mat)) > 150) {
      show_rownames = FALSE
    }
  
    annotation_col = NULL
    annotation_colors = NULL
    show_colnames <- TRUE
    if(!is.null(col_ann)) {
      annotation_col <- col_ann
      show_colnames <- FALSE
    }
    if(!is.null(ann_colors)) {
      annotation_colors <- ann_colors
    }
  
    if(cluster_rows == FALSE) {
      cluster_rows = FALSE
    }
    # set the scale of the heatmap to 10

    myBreaks <- c(seq(min(mat), 0, length.out=ceiling(50/2) + 1), 
                  seq(max(mat)/50, max(mat), length.out=floor(50/2)))
  
  
    if (save_plot == TRUE) {
        pheatmap(mat,
                 breaks = myBreaks,
                 color = cell_colors,
                 border_color = NA,
                 cluster_rows = cluster_rows,
                 cluster_cols = FALSE,
                 annotation_col = annotation_col,
                 annotation_colors = annotation_colors,
                 main = title,
                 fontsize_row = fontsize_row,
                 fontsize_col = 12,
                 show_rownames = show_rownames,
                 show_colnames = show_colnames,
                 filename = filename
                 )
    } else {
        pheatmap(mat,
                 breaks = myBreaks,
                 color = cell_colors,
                 cluster_rows = cluster_rows,
                 border_color = NA,
                 cluster_cols = FALSE,
                 annotation_col = annotation_col,
                 annotation_colors = annotation_colors,
                 main = title,
                 fontsize_row = fontsize_row,
                 fontsize_col = 12,
                 show_rownames = show_rownames,
                 show_colnames = show_colnames)
    }  
    # Create pdf heatmap
  
  }
}