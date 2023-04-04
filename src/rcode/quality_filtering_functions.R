## This file includes functions related to quality-filtering
 
library(monocle)
library(Seurat)
library(futile.logger)
library(Matrix)
library(simpleSingleCell)


cleanGeneNamingForCellRanger <- function(monocleCellDataset, columnName="EnsemblID", pattern="premRNA", replacement="") {
    flog.info("changing %s: \"%s\" to \"%s\"", columnName, pattern, replacement)
    a.f <- fData(monocleCellDataset)
    origNames <- a.f[[columnName]]
    changedNames <- gsub(pattern, replacement, a.f[[columnName]])
    flog.info("changed: %d entries", sum(origNames != changedNames))
    a.f[[columnName]] <- changedNames
    fData(monocleCellDataset) <- a.f
    monocleCellDataset
}


#'
#' Function to retrieve id from pData(cds) or stop if it is not present and not given by user
#' 
getIdFromPdata <- function(cds, currId=NULL) {
  a.p <- pData(cds)
  if(is.null(currId) && ("id" %in% colnames(a.p))) {
    currId <-a.p$id[1]
  } 
  if(is.null(currId)) {
    stop("please provide currId")
  }
  flog.info("currId: %s", currId)
  currId
} 

#'
#' This code finds the list of mitochondrial genes (using grepl over `^MT-`) and cells from the data
#' using `nmads=3` mean absolute deviation higher than the median fraction of mitochondrial reads across cells
#' 
#' @return list(mito_genes=, mito_cells=)
#' 
filter_mitochondrial_genes_and_cells <- function(monocleCellDataset, 
                                                 mt_pattern='^MT-', 
                                                 ignore.case=TRUE,
                                                 outDir=getwd(),
                                                 nmads=3, 
                                                 currId=NULL) { 
  flog.info("mt_pattern: %s ignore.case: %s nmads: %d outDir: %s", mt_pattern, ignore.case, nmads, outDir)
  a.f <- fData(monocleCellDataset)
  a.p <- pData(monocleCellDataset)

  currId <- getIdFromPdata(monocleCellDataset, currId)    

  b <- a.f$gene_short_name
  
  idx <- grepl(mt_pattern, b, ignore.case = ignore.case)

  if(sum(idx) > 0) { 
    totalCells <- Matrix::colSums(as.matrix(monocleCellDataset))
    
    mitoCells <- Matrix::colSums(as.matrix(monocleCellDataset[idx, ]))
    mitoFraction <- mitoCells/totalCells
    is.mito <- isOutlier(mitoFraction, nmads=nmads, type="higher", log=FALSE)
    cutoff <- min(mitoFraction[is.mito])
    mitoFile=getFilePath(fileName = paste0(currId, "-mito.eps"), subDir = 'mito', outDir = outDir)
    p<-ggplot(data.frame(mitoFraction=mitoFraction), aes(x=mitoFraction))+geom_histogram() + scale_x_log10()
    p <- p + geom_vline(xintercept = cutoff, color="red", linetype="dashed")
    subtitle1 <- paste0("cutoff: ", floor(1000*cutoff)/1000, " mito_cells: ", sum(is.mito), "/", length(is.mito))  
    p <- p + labs(title=currId, subtitle = subtitle1 )
    ggsave(mitoFile, plot=p, width=7, height=7, dpi=600, units = "in")
    saveRDS(p, file = paste0(mitoFile, ".rds"))
    
    
    a.p["isMito"] <- is.mito
    a.f["isMito"] <- idx
    
    baseColumns <- c("barcode",	"num_genes_expressed",	"sample",	"assembly",	"method",	"id",	"isMito")
    baseColumns.reduced <- Filter(function(x) {  x %in% colnames(a.p) }, baseColumns)
    
    write.table(a.p[a.p$isMito, baseColumns.reduced], 
                getFilePath(paste0(currId, "-mito-cells.tsv"), outDir = outDir, subDir = "mito"), 
                sep="\t", row.names=FALSE, col.names=TRUE, quote = FALSE)
    
    
    geneColumns <-c("EnsemblID",	"gene_short_name", "isMito")
    geneColumns.reduced <- Filter( function(x) { x %in% colnames(a.f) }, geneColumns)
    write.table(a.f[a.f$isMito, ], 
                getFilePath(paste0(currId, "-mito-genes.tsv"), outDir = outDir, subDir = "mito"), 
                sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
    
    
    flog.info("get_mitochondrial_genes_and_cells() %s removed genes: %d cells: %d/%d", currId, sum(idx), sum(is.mito), length(is.mito))
    monocleCellDataset[!idx, !is.mito]
  } else {
      flog.warn("did not find any gene names matching %s pattern. check your gene names for dataset: %s", mt_pattern, currId)
      monocleCellDataset
  }
}


filter_cells_with_low_number_of_genes <- function(cds, 
                                                  nmads=3, 
                                                  type="both", 
                                                  minGenes=0, 
                                                  outDir=getwd(),
                                                  currId=NULL) {
  a.p <- pData(cds)
  id   <- getIdFromPdata(cds, currId)

  flog.info("filter_cells_with_low_number_of_genes(): %s nmads: %d endToCut: %s", id, nmads, type)
  
  lessThanMedian.idx <- a.p$num_genes_expressed < median(a.p$num_genes_expressed)
  
  
  flog.info("\tremoving cells using outlier detection")
  o.idx <- isOutlier(a.p$num_genes_expressed, nmads=nmads, type=type, log=TRUE)
  cutoffThreshold <- max(a.p$num_genes_expressed[o.idx & lessThanMedian.idx])
  
  cutoffHigher <- min(a.p$num_genes_expressed[o.idx & (!lessThanMedian.idx)])
  
  flog.info("\toutlier limit (within %d mean absolute deviation of median): %s", nmads, cutoffThreshold)
  flog.info("\tcells with low number of expressed genes: %d", sum(o.idx))
  
  m.idx <- (a.p$num_genes_expressed <= minGenes) 
  flog.info("\tremoving cells having less than %d expressed genes: %d", minGenes, sum(m.idx))
  
  
  fileName <- file.path(outDir, "cells", paste0(id, '-cells-histogram.eps'))
  dir.create(dirname(fileName), recursive = TRUE, showWarnings = FALSE)
  
  p <- ggplot(a.p, aes(x=num_genes_expressed))+geom_histogram() + scale_x_log10()
  cellsRemoved <- 0
  if(is.numeric(cutoffThreshold) & is.finite(cutoffThreshold)) {
    cellsRemoved = sum(o.idx)
    p <- p + geom_vline(xintercept = cutoffThreshold, color="red", linetype="dashed")
  }
  
  if(is.numeric(cutoffHigher) & is.finite(cutoffHigher)) {
    p <- p + geom_vline(xintercept = cutoffHigher, color="red", linetype="dashed") 
  }
  
  if(minGenes > 0) {
    cellsRemoved = cellsRemoved + sum(m.idx)
    p <- p + geom_vline(xintercept = minGenes, color="magenta", linetype="dashed")
  
  }

  p <- p + labs(title=id, x="genes expressed per cell", y="number of cells", 
                subtitle=paste0("outliers: ", cellsRemoved, " remaining: ",  (dim(a.p)[1] - cellsRemoved) ))
  ggsave(filename = fileName, p, dev="eps", width=7, height=7, units = "in", dpi=600)
  saveRDS(p, file=paste0(fileName, ".rds"))  
  
  flog.info("\tplot: %s", fileName)
  
  fileName2 = file.path(outDir, "cells", paste0(id, '-cells-filtered.tsv'))
  a.p["filter"] <- !(o.idx | m.idx)
  
  write.table(a.p, file=fileName2, sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  flog.info("\tcsvOut: %s", fileName2)
  cds[, !(o.idx | m.idx)]
}

#'
#' Function to merge two data-frames obtained from `get_genes_expressed_in_low_number_of_cells`
#' 
#' Returns merged data.frame(EnsemblID, short_gene_name, filter)
#' 
get_merged_genes_to_filter <- function(x.f, y.f) { 
  
  x = x.f[, c("EnsemblID", "gene_short_name", "filter")]
  y = y.f[, c("EnsemblID", "gene_short_name", "filter")]
  
  flog.info("x: %d/%d y: %d/%d", sum(x$filter), dim(x)[1], sum(y$filter), dim(y)[1])
  z <- merge(x, y, by=c("EnsemblID", "gene_short_name"), all=TRUE)
  
  ## if adding new genes - then we don't know anything about them so far
  z$filter.x[is.na(z$filter.x)] <- FALSE
  
  ## but the old genes should not be perturbed if not present in the newly added dataset
  z$filter.y[is.na(z$filter.y)] <- TRUE
  
  z["filter"] <- z["filter.x"] & z["filter.y"]
  z$filter.x <- NULL
  z$filter.y <- NULL
  
  z
}


#' Function to filter genes expressed in fewer than `cutoff` cells
#'
get_genes_expressed_in_low_number_of_cells <- function(cds, remove=FALSE, minFrac=0.001, outDir=getwd(), currId=NULL) {
  currId <- getIdFromPdata(cds, currId)
  flog.info("Filtering genes expressed in few cells: %s minFrac: %f", currId, minFrac)
  file_name <- file.path(outDir, "genes", paste0(currId, "-genes-histogram.eps"))
  dir.create(dirname(file_name), recursive = TRUE, showWarnings = FALSE)
  
  cds <- detectGenes(cds) ## detect-genes with min_expr=0.1
  a.f <- fData(cds)
  totalCells <- dim(cds)[2]
  cutoff <- floor(minFrac * totalCells)
  flog.info("\ttotalCells: %d cutoff: %d", totalCells, cutoff)
  idx <- a.f$num_cells_expressed < cutoff
  
  p <- ggplot(a.f, aes(x = num_cells_expressed)) + geom_histogram() + scale_x_log10()
  p <- p + labs(title=currId, x='number of cells where gene is expressed', y='number of genes', subtitle=paste0('total: ',  dim(cds)[1] ,' below cutoff: ', sum(idx))) 
  p <- p + geom_vline(xintercept = cutoff, color='red', linetype='dashed')
  p 
  ggsave(filename = file_name, plot = p, width = 7, height= 7, units="in", device="eps", dpi=600)
  saveRDS(p, file = paste0(file_name, ".rds"))
  
  flog.info("\t%s genes expressed in less than %d cells: %d", currId, cutoff, sum(idx))
  flog.info("\tplot: %s", file_name)
  
  a.f[["filter"]] <- idx
  csv_file <- file.path(outDir, "genes", paste0(currId, "-gene-histogram.tsv"))
  write.table(a.f, file=csv_file, quote = FALSE, sep="\t", col.names = TRUE, row.names = FALSE)
  flog.info("\tcsv_file: %s", csv_file)
  a.f
}

find_cells_with_human_to_mouse_ratio <- function(cds1, humanGeneIDs, ratioLow=5.0/6.0, ratioHigh=1.0) {
  a.f <- fData(cds1)
  humanIdx <- match(humanGeneIDs, a.f$EnsemblID)
  humanIdx <- humanIdx[!is.na(humanIdx)]
  ratIdx <- setdiff(1:dim(cds1)[1], humanIdx)
  humanReads <- Matrix::colSums( as.matrix(cds1[humanIdx, ]) )
  ratReads <-Matrix::colSums( as.matrix(cds1[ratIdx, ]) )
  cDf <- data.frame(barcode=names(humanReads), GRCh38=humanReads, Rnor6=ratReads, ratio=humanReads/(humanReads+ratReads))
  cDf$Hs <- (cDf$ratio >= ratioLow) & (cDf$ratio <= ratioHigh)
  flog.debug("find_cells_with_human_to_mouse_ratio() humanIdx: %d ratIdx: %d humanReads: %d ratReads: %d accepted: %d",  
            length(humanIdx), length(ratIdx), sum(humanReads), sum(ratReads), sum(cDf$Hs))
  cDf
}