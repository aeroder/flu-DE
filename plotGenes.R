#
# Plot Genes
# Function: Input a DESeq 'dds' object (or several for multiple time points, 
#           and a gene and graph the normalized counts of that gene over time
#
# Author: AR
# Last updated: 8/20/19
#

plotGenes <- function(gene, ddsList, ...) {
  
    library(plyr)
    library(dplyr)
    library(DESeq2)

    if (length(ddsList) == 1) {
        geneCounts <- plotCounts(ddsList[1], gene = gene, intgroup = c("Virus","Time"), returnData = TRUE)
    } else if (length(ddsList) > 1) {
        for (i in seq_along(ddsList)) {
            if (i == 1) {
                geneCounts <- plotCounts(ddsList[[1]], gene = gene, intgroup = c("Virus","Time"), returnData = TRUE)
            } else {
                geneCounts <- rbind(geneCounts, plotCounts(ddsList[[i]], gene = gene, intgroup = c("Virus","Time"), returnData = TRUE))
            }
        }
    }
    
    geneCountsum <- plyr::ddply(geneCounts, c("Virus", "Time"), 
                                summarise,
                                N = length(count),
                                mean = mean(count),
                                sd   = sd(count),
                                se   = sd / sqrt(N)
    )
    
    g <- ggplot(geneCountsum, aes(x = Time, y = mean, group = Virus)) +  
                geom_point(aes(color = Virus)) +
                geom_line(aes(color = Virus)) +
                geom_errorbar(aes(ymin=mean-se, ymax=mean+se, color = Virus), width = 0.1)
    
    return(g)


}













