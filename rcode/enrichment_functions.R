silent_load <- function(PKG) { suppressWarnings(suppressMessages(library(PKG, quietly = T, character.only = T))) }
sapply(c("data.table", 
         "clusterProfiler", 
         "dplyr", 
         "stringr"), silent_load) 


overlap_phyper2<-function(L1,L2,bg=length(unique(c(unlist(L1),unlist(L2))))){
  # takes two list with indices and creates a matrix with pvalue for overlap between elements
  # phyper test uses all entries in L as background if bg is not specified.
  nL1<-length(L1)
  nL2<-length(L2)
  M<-mat.or.vec(nL1,nL2)
  P<-mat.or.vec(nL1,nL2)
  
  nA <- rep(0, nL1)
  nB <- rep(0, nL2)
  
  P[,]<-1
  for (i in 1:nL1){
    nA[i] <- length(L1[[i]])
    for (j in  1:nL2){
      nB[j] <- length(L2[[j]])
      
      M[i,j]<-length(intersect(L1[[i]],L2[[j]]))
      if (length(L1[[i]])==0 | length(L2[[j]])==0 ) {
        P[i,j] <- 1
      }else{
        P[i,j]<-1-phyper(M[i,j],length(L1[[i]]),bg-length(L1[[i]]),length(L2[[j]]))
      }
    }
  }
  colnames(P)<-names(L2)
  rownames(P)<-names(L1)
  colnames(M)<-names(L2)
  rownames(M)<-names(L1)   
  
  return(list(P=P,M=M, A=nA, B=nB))
}



## main code begins here 
getCellTypeFromFileName <- function(currFile) {
  returnValue <- sub(".tsv", "", basename(currFile))
  returnValue
}

readDtFromFile <- function(currFile) {  
  currDt <- fread(currFile)
  currCellType <- getCellTypeFromFileName(currFile)
  flog.info("read %d entries for file: %s", nrow(currDt), currCellType )
  currDt <- currDt[grep("Hs", currDt$species), ]
  flog.info("entries after removing non-Hs genes: %d", nrow(currDt))
  currDt 
}

#'
#' Function that reads a list of tsv files to populate enrichmentCellType for `clusterProfiler::enricher`
#' 
#' @example
#'  baseDir <-'/Users/yo4230sh/November16/'
#'  markerDir <- paste0(baseDir, '/MarkersList_Download/')
#'  flog.info("directory: %s exists: %s", markerDir, dir.exists(markerDir))
#'  fileNames <- Sys.glob(file.path(markerDir, "*.tsv"))
#'  combinedCellTypeDt <- readListOfFilesToGetCombinedDt(fileNames)
#'
readListOfFilesToGetCombinedDt <- function(fileNames) { 
  dtList <- list()
  for(currFile in fileNames) { 
    currName <- getCellTypeFromFileName(currFile)
    dtList[[currName]] <- readDtFromFile(currFile)
  }
  combinedDt <- do.call("rbind", dtList)
  combinedDt
} 

#' Simple function to read the panglodb downloaded files present in `filesDir`
readPangloDbFilesDir <-function(filesDir) {
  flog.info("readPangloDbFilesDir(): %s", filesDir)
  readListOfFilesToGetCombinedDt(Sys.glob(file.path(filesDir, "*.tsv")))
}


####################################################################
## Handling Cloupe Files
####################################################################
getNumClustersInCloupeGeneDt <- function(x) (ncol(x) - 2)/3

extractedSortedGeneOrder <- function(geneDt, clusterId) {
  currAvg <-geneDt[,(2+(clusterId-1)*3 + 1), with=FALSE]
  currLFC <- geneDt[,(2+(clusterId-1)*3 + 2), with=FALSE]
  currPval <- geneDt[,(2+(clusterId-1)*3 + 3), with=FALSE]
  pvalueColumnName <- paste("Cluster", clusterId, "P-Value")
  currOrderedIndex <- order(currPval, decreasing = FALSE)
  sortedGeneOrder <- geneDt[currOrderedIndex, c("GeneName", pvalueColumnName), with=FALSE]
  sortedGeneOrder
}

filterCloupeGeneDt <- function(geneDt) {
  mtGenes <- geneDt$GeneName[grep("^MT-", geneDt$GeneName)]
  geneDt <- geneDt[!grepl("^MT-", geneDt$GeneName), ]
  flog.info("removing %d mitochondrial genes: %s remaining genes: %d", length(mtGenes), paste(mtGenes, collapse=", "), nrow(geneDt))
  geneDt
}

getOrderedGenesPerClusterFromCloupeGeneDt <- function(geneDt) { 
  numGenes <- nrow(geneDt)
  numClusters <- getNumClustersInCloupeGeneDt(geneDt)
  flog.info("processing file: numGenes: %d numClusters: %d", numGenes, numClusters)
  flog.info("read file: %s numGenes: %d numClusters: %d", basename(geneListFile), numGenes, numClusters)
  orderedGenesPerCluster <- list()
  for(clusterId in 1:numClusters) {
    orderedGenesPerCluster[[clusterId]] <- extractedSortedGeneOrder(geneDt, clusterId)
  }
  orderedGenesPerCluster
} 

####################################################################
## Main enrichment functions
####################################################################

#'
#' Function to get enrichment result 
#' 
#' @param currOrderedGenesDt => data.table(geneName, pValue)
#' @param combinedCellTypeDt => standard output read from `readPangloDbFilesDir`
#' @param pValueThreshold => which genes to call as significant
#' @param enrichmentCutoffPValueThreshold cutoff p-value for enricher
#' @param defaultQvalueCutoff (default: 0.2) cutoff q-value for enricher
#' 
getEnrichmentResult <- function(currOrderedGenesDt, combinedCellTypeDt, pValueThreshold, 
                                enrichmentCutoffPValueThreshold, defaultQvalueCutoff=0.2) { 
  currSignificantGenes <- currOrderedGenesDt$GeneName[currOrderedGenesDt[, 2] < pValueThreshold]
  flog.info("geneList pValThreshold: %f numSignificantGenes: %d", pValueThreshold, length(currSignificantGenes))
  term2geneDt <- combinedCellTypeDt[, c("cell type", "official gene symbol"), with=FALSE]
  colnames(term2geneDt) <- c("cellType", "gene")
  myQValueCutoff <-  max(defaultQvalueCutoff, enrichmentCutoffPValueThreshold) 
  flog.info("enricher pvalueCutoff: %f qvalueCutoff: %f", enrichmentCutoffPValueThreshold, myQValueCutoff)
  currResult <- enricher(currSignificantGenes, TERM2GENE = term2geneDt, pvalueCutoff = enrichmentCutoffPValueThreshold, qvalueCutoff = myQValueCutoff)
  currResult
} 

#' Function to return enriched genes per cluster
#' @param orderedGenesPerCluster list of data.table(gene, p-value)
#' @param combinedCellTypeDt => standard output read from `readPangloDbFilesDir`
#' @param pValueThreshold => which genes to call as significant
#' @param enrichmentCutoffPValueThreshold cutoff p-value for enricher
#' 
getEnrichmentPerCluster <- function(orderedGenesPerCluster, 
                                    combinedCellTypeDt, pValueThreshold, enrichmentCutoffPValueThreshold) { 
  enrichmentPerCluster <- list()
  for(clusterId in 1:length(orderedGenesPerCluster)) {
    currInput <- orderedGenesPerCluster[[clusterId]]
    currResult <- getEnrichmentResult(currInput, combinedCellTypeDt, pValueThreshold, enrichmentCutoffPValueThreshold)
    enrichmentPerCluster[[clusterId]] <- currResult
  }
  enrichmentPerCluster
}



getCombinedResult <- function(enrichmentPerCluster) { 
  combinedResult <- list()
  for(clusterId in 1:length(enrichmentPerCluster)) { 
    cResult <- enrichmentPerCluster[[clusterId]]
    innerDt <- as.data.frame(cResult)
    if(nrow(innerDt) > 0) { 
      innerDt$clusterId <- clusterId
      combinedResult[[clusterId]] <- innerDt
    }
  }
  combinedResult <- rbindlist(combinedResult)
  combinedResult
} 

getGeneRatio <- function(geneRatioStr) { 
  countA <- as.numeric( sub("/[0-9]+$", "", geneRatioStr) )
  countB <- as.numeric( sub("^[0-9]+/", "", geneRatioStr) )
  countA/countB
}

getDotPlot <- function(combinedResult, title=NULL) { 
  p <- ggplot(combinedResult, aes(y=ID, x=GeneRatioFloat)) 
  p <- p + geom_point(aes(color=as.factor(clusterId), size=-log10(pvalue) )) 
  if(!is.null(title)) {
    p <-p + ggtitle(title)
  }
  p 
} 

getMatrixFromSigDt <- function(sigDt) { 
  x <- sigDt
  x$GeneName <- NULL
  m <- as.matrix(x)
  rownames(m) <- sigDt$GeneName
  m 
}

getClusteredMatrix <- function(m, clustering='none') {
  if(is.function(clustering))
  {
    m=clustering(m)
  }else
  {
    if(clustering=='row')
      m=m[hclust(dist(m))$order, ]
    if(clustering=='column')
      m=m[,hclust(dist(t(m)))$order]
    if(clustering=='both')
      m=m[hclust(dist(m))$order ,hclust(dist(t(m)))$order]
  }
  m
}

#'
#' Function that takes a processed Seurat Cell Dataset and list of markers to perform panglodb enrichment
#' 
#' @param pbmc - Seurat Cell Dataset
#' @param pbmc.markers - List of marker (genes) identified
#' @param combinedCellTypeDt - Output of \code{\link{readPangloDbFilesDir}}
#' @param pValueThreshold
#' @param enrichmentCutoffPValueThreshold
#' @param cSample - optional sample name for current cell dataset 
#' @param outDir - where results are saved
#'
do_seurat_enrichment <- function(pbmc, pbmc.markers, combinedCellTypeDt, 
                                 pValueThreshold=0.05, 
                                 enrichmentCutoffPValueThreshold = 0.05, 
                                 cSample=NULL, outDir=getwd()) {
  if(is.null(combinedCellTypeDt)) {
    flog.warn("PangloDB combinedCellTypeDt not provided - skipping enrichment step")
    list(scds=pbmc, markers=pbmc.markers)
    
  } else {
    flog.info("do_seurat_enrichment(): sample: %s, outDir: %s markers: %d pValueThreshold: %f enrichmentCutoffPValueThreshold: %f", 
              cSample, outDir, length(pbmc.markers$gene), pValueThreshold, enrichmentCutoffPValueThreshold)
    
    clustIds <- unique(pbmc.markers$cluster)
    mList <- lapply(clustIds, function(x) {
      z <- pbmc.markers[pbmc.markers$cluster == x, c("gene", "p_val")]
      colnames(z) <- c("GeneName", "pValue")
      z
    })
    cpResult <- getEnrichmentPerCluster(orderedGenesPerCluster = mList, 
                                        combinedCellTypeDt = combinedCellTypeDt, 
                                        pValueThreshold = pValueThreshold, 
                                        enrichmentCutoffPValueThreshold = enrichmentCutoffPValueThreshold)
    combinedResult <- getCombinedResult(cpResult)
    combinedResult$clusterId <- clustIds[combinedResult$clusterId]
    combinedResult$GeneRatioFloat <- getGeneRatio(combinedResult$GeneRatio)
    
    if(is.null(cSample)) {
      cSample <- get_sample_name_from_seurat_cds(pbmc)
    } 
    p <- getDotPlot(combinedResult,title=paste0(cSample, " pVal: ", enrichmentCutoffPValueThreshold))
    ggsave(filename = paste0(outDir, "/", cSample,"_dotplot_", enrichmentCutoffPValueThreshold,  ".pdf"))
    combinedResult
  }
} 
