library(Seurat)
library(ggplot2)

quick_cell_count <- function(cCds, cGene) {
  sum(as.matrix(cCds@raw.data[which(rownames(cCds@data)==cGene), ]) > 0)
}

make_feature_plots <- function(cElem, Markers, outDir, cName=NULL, ncol=5, pt.size=1, width=15, height=25) { 
  cList <- list()  
  cCds <- cElem$scds
  cGenes <- Reduce(union, Markers)
  cIdx <- which(cGenes %in% rownames(cCds@data))
  cGenes.keep <- cGenes[cIdx]
  for(cGene in cGenes.keep) {
    flog.info("make_feature_plots() %s", cGene)
    currCds  <- cElem$scds
    p <- FeaturePlot(currCds, cGene, do.return=T, pt.size=pt.size)
    p[[cGene]] <- p[[cGene]] + labs(title=paste0(cGene, " (", quick_cell_count(cCds, cGene), "/", dim(cCds@data)[2],")"), x="", y="")
    cList[[cGene]] <- p[[cGene]]
  }
  plot_grid(plotlist=cList, ncol=ncol)
  if(is.null(cName)) {
    gfile <- file.path(outDir, paste0("markers.pdf"))
  } else {
    gfile <- file.path(outDir, paste0(cName, "-markers.pdf"))
  }
  ggsave(filename = gfile, width=width, height=height, units="in", dpi=300)
}
