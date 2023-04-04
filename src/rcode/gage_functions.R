# gage_functions - this file contains functions related to GAGE enrichment
require(gage)
require(gageData)
require(pathview)
require(futile.logger)

data(kegg.gs) # the genesets 
data(egSymb) # ids to symbols

getLfcFromDeseq2 <- function(deseq2.res) {
  deseq2.fc=deseq2.res$log2FoldChange
  names(deseq2.fc)=rownames(deseq2.res)
  deseq2.fc
}




#'
#' This function generates gage enrichment result on kegg genesets
#' 
getKeggGsResultsFromDeseqResult <- function(deseq2.res, padjThreshold=0.05, gsets=kegg.gs.sym) { 
  z1 <- deseq2.res[ !is.na(deseq2.res$padj), ]
  z2 <- z1[z1$padj < padjThreshold, ]
  flog.info("padjThreshold: %f dim: %s", padjThreshold, dim(z2))
  exp.fc = getLfcFromDeseq2(z2)
  kegg.gs.sym<-lapply(kegg.gs, eg2sym)
  fc.kegg.p <- gage(exp.fc, gsets = kegg.gs.sym, ref = NULL, samp = NULL)
  fc.kegg.p
}



selectPathways <-function(fc.kegg.p, qValueThreshold) {
  sel <- fc.kegg.p$greater[, "q.val"] < qValueThreshold & + !is.na(fc.kegg.p$greater[, "q.val"])
  path.ids <- rownames(fc.kegg.p$greater)[sel]
  sel.l <- fc.kegg.p$less[, "q.val"] < qValueThreshold &
    + !is.na(fc.kegg.p$less[,"q.val"])
  path.ids.l <- rownames(fc.kegg.p$less)[sel.l]
  path.ids2 <- substr(c(path.ids, path.ids.l), 1, 8)
  path.ids2
} 


