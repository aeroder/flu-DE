# Try out goseq on the top DE genes

# make a vector of 0s and 1s where 1 is a DE gene and 0 is not DE
# set the names of the vector to the geneIDs


# Function to analyze GO categories of DE genes
#
# Takes as input:
# 1. wald: A list of DESeq Results tables using the wald test
# 2. location: The location within the list to run the goseq on (ie. 3 for the first virus at 24 hours)
# 3. cutoff: the adjusted p value to use for what is considered to be differentially expressed
# 4. numCats: the number of categories to use in the graph. Recommended < 40 or the words will be too small

GOplot <- function(wald, location, cutoff, numCats) {
  wald <- wald
  loc <- location
  pv <- cutoff
  num <- numCats
  
  # create a binary integer with all of the genes in wald
  # if the fold change is less than the cutoff, use a 1 value. If it isn't, use a 0 value
  resSig <- as.integer(wald[[loc]]$padj < pv)
  
  # use the gene ids as the row names for the binary list
  names(resSig) <- rownames(wald[[loc]])
  # get rid of anything that has an NA value
  resSig <- na.omit(resSig)
  
  # create the goseq nullp element
  nulp <- nullp(resSig, 'hg19', 'ensGene', bias.data = NULL, plot.fit = FALSE)
  # perform the goseq to pull out the enriched categories and their associated p values
  gseq <- goseq(nulp, 'hg19', 'ensGene', gene2cat=NULL)
  
  # graph 
  out <- gseq %>% 
    top_n(numCats, wt=-over_represented_pvalue) %>% 
    mutate(hitsPerc=numDEInCat*100/numInCat) %>% 
    ggplot(aes(x=hitsPerc, 
               y=term, 
               colour=over_represented_pvalue, 
               size=numDEInCat)) +
    geom_point() +
    scale_colour_distiller(palette = "PRGn") +
    expand_limits(x=0) +
    labs(x="Hits (%)", y="GO term", colour="P value", size="Count")
  
  return(out)
  
}
