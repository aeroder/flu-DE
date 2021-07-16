#
# VennsD
# Function: Use these numbers generate by the common function to make a Venn Diagram
#
# Input: Namevec (optional): A vector containing the names to use for labeling the circles in the VD
#        degenes: A list of the genes to compare between each of the sets
#        genes: A full list of possible genes
#        file_prefix: A prefix to save the file to, shouldn't contain spaces
#        return_list: A boolean value that when TRUE, will return the genes in common between all of the sets.
#                     Defaults to FALSE indicating that the function should create and save a plot
#        a: A list of the comparisons to make that corresponds to locations within de_list
#
# Author: Allison Roder
# Last updated: 6/15/19

vennsD <- function(namevec = NULL, degenes, genes, file_prefix = NULL, return_list = FALSE, save_plots = FALSE, a, ...) {
  namevec <- namevec
  de <- as.list(degenes)
  genes <- genes
  colorvec <- colors


  if (return_list == TRUE) {
    out <- common(de, genes, return_list = TRUE, samp = a)

  } else {
  if (length(a) == 1) {
    out <- draw.single.venn(common(de, genes, samp = c(a)), ...,
                            category = namevec,
                            fill = c("#843b62"))
  }
  if (length(a) == 2) {
    out <- draw.pairwise.venn(common(de, genes, samp = c(a[1])), 
                              common(de, genes, samp = c(a[2])),
                              common(de, genes, samp = a), ...,
                              category = namevec,
                              fill = c("#21A064","#505160"),
                              scaled = FALSE,
                              ext.text = FALSE,
                              cat.pos = 0,
                              cex = rep(1.5,3),
                              cat.cex = rep(1.5,2)
                              )
  }
  if (length(a) == 3) {
    out <- draw.triple.venn(common(de, genes, samp = c(a[1])),
                            common(de, genes, samp = c(a[2])),
                            common(de, genes, samp = c(a[3])),
                            common(de, genes, samp = c(a[1:2])),
                            common(de, genes, samp = c(a[2:3])),
                            common(de, genes, samp = c(a[c(1,3)])),
                            common(de, genes, samp = a), ...,
                            category = namevec,
                            fill = c("#21A064","#32A4C5","#505160"),
                            cex = rep(1.5,7),
                            scaled = TRUE,
                            cat.cex = rep(1.5,3))
  }
  if (length(a) == 4) {
    out <- draw.quad.venn(common(de, genes, samp = c(a[1])),
                          common(de, genes, samp = c(a[2])),
                          common(de, genes, samp = c(a[3])),
                          common(de, genes, samp = c(a[4])),
                          common(de, genes, samp = c(a[1:2])),
                          common(de, genes, samp = c(a[c(1,3)])),
                          common(de, genes, samp = c(a[c(1,4)])),
                          common(de, genes, samp = c(a[2:3])),
                          common(de, genes, samp = c(a[c(2,4)])),
                          common(de, genes, samp = c(a[3:4])),
                          common(de, genes, samp = c(a[1:3])),
                          common(de, genes, samp = c(a[c(1,2,4)])),
                          common(de, genes, samp = c(a[c(1,3,4)])),
                          common(de, genes, samp = c(a[2:4])),
                          common(de, genes, samp = a), ...,
                          category = namevec,
                          fill = c("#6e7f80","#536872","#708090","#536878"),
                          cex = rep(1.5,15),
                          cat.cex = rep(1.5,4))
  }
  if (!exists("out"))
    out <- "Oops"

  #tiff(filename = glue("{file_prefix}.tiff"), units = "in", res = 300)
  #grid.draw(out)

  filename_pdf = glue("{file_prefix}.pdf")
  filename_png = glue("{file_prefix}.png")

    if (save_plots == TRUE) {
      save_plot(filename_png, out, base_width = 8, base_height = 5, units = "in", dpi = 300)
      Sys.sleep(1)
      save_plot(filename_pdf, out, base_width = 8, base_height = 5, units = "in", dpi = 300)
      Sys.sleep(1)
    }
  }
  return(out)
}
