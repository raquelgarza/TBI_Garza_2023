library(stringr)

#'
#' Function that reads cellranger output present in `dataDir` by searching for `searchPattern`
#' 
#' Assumes the structure: `dataDir`/`sequencingMethod`/`sampleId`/`searchPattern`
#'
#' @return data.frame containing `sample`, `method` and `dataDir` 
#' 
read_cellranger_data_files <- function(dataDir, searchPattern="outs/filtered_feature_bc_matrix/") {
  cellrangerOuts <- Sys.glob(file.path(dataDir, "*", "*", searchPattern), dirmark = F)
  a <- cellrangerOuts[1]
  sequencingMethod <- sapply(cellrangerOuts, function (a) gsub("(.*)/([^/]+)/([^/]+)/outs/filtered_feature_bc_matrix/$", "\\2", a))
  xx <- str_split(sequencingMethod, '_')
  assembly <- sapply(xx, function(x) {x[[1]] })
  FACS <- sapply(xx, function(x) { x[2] })
  sampleId <- sapply(cellrangerOuts, function (a) gsub("(.*)/([^/]+)/([^/]+)/outs/filtered_feature_bc_matrix/$", "\\3", a))
  retDf <- data.frame(sample=sampleId, assembly=assembly, method=FACS, dataDir=cellrangerOuts, stringsAsFactors = FALSE)
  retDf$id <- paste(retDf$assembly, retDf$method,retDf$sample, sep='-')
  rownames(retDf) <- retDf$id
  retDf 
} 


#'
#' Function that reads cellranger output present in `dataDir` by searching for `searchPattern`
#' 
#' Assumes the structure: `dataDir`/`sequencingMethod`/`sampleId`/`searchPattern`
#'
#' @return data.frame containing `sample`, `method` and `dataDir` 
#' 
find_cellranger_data_files <- function(dataDir, searchPattern="/outs/filtered_feature_bc_matrix") {

  dataDir = normalizePath(dataDir)
  matchingDirectories = Filter(function (x) grepl(searchPattern, x), list.dirs(dataDir))  
  basePaths = sapply(matchingDirectories, function (a) sub(searchPattern, "", a))

  sampleId <- basename(basePaths)  
  retDf <- data.frame(sample=sampleId, dataDir=matchingDirectories, stringsAsFactors = FALSE)
  retDf 
} 
