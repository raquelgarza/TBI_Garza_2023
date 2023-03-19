source("rcode/pipeline_helper_functions.R")
source("rcode/read_cellranger_data_files.R")

silently_load_package_list(c("monocle", "Seurat", "Matrix", "stringr"))




#'
#' Function to convert duplicated gene names by appending _1, _1_2 to each repeated instance
#' of the name. Used in `convert_monocle_cds_to_seurat_cds`
#' 
fix_duplicated_gene_short_name <- function(mCds) {
  a.f <- fData(mCds)
  while(sum(duplicated(a.f$gene_short_name) > 0)) { 
    dupIdx <- duplicated(a.f$gene_short_name)
    dupGene <- a.f$gene_short_name[dupIdx]
    uniqueNamesOfDupGenes <- unique(dupGene)
    flog.info("found %d duplicated genes", sum(dupIdx))
    for(currName in uniqueNamesOfDupGenes) {
        count <- 1
        currIdx <- which(a.f$gene_short_name == currName)
        flog.info("fixing gene: %s having %d duplicates", currName, length(currIdx))
        for(j in currIdx) {
          a.f$gene_short_name[j] <- paste0(currName, "_", count)
          count <- count+1
        }
    }
  }
  fData(mCds) <- a.f
  mCds
}


#'
#' Get path for where to put the file
#' 
getFilePath <- function(fileName, outDir=NULL, subDir=NULL) {
  if(is.null(outDir)) { 
    outDir <- Sys.Date()
  }
  
  if(is.null(subDir)) { 
    retVal <- file.path(outDir, fileName)  
  } else { 
    retVal <- file.path(outDir, subDir, fileName)
  }
  dir.create(dirname(retVal), showWarnings = FALSE, recursive = TRUE)
  retVal
}


#'
#' Read a file transparently if it is gzipped
#' 
gz_fix <- function(file_name) { 
  if(grepl("^.*.gz$", file_name)) { 
    file_name <- gzfile(file_name, open="r") 
  }
  file_name
} 


#'
#' Convert cellranger v3 filtered outputs to monocle CellDataset
#' 
get_monocle_cell_dataset_from_cellranger_files <- function(matFile, geneFile, barcodeFile) {
  m <- readMM(matFile)
  flog.info("total number of reads: %d", sum(m))
  gene_info <- read.delim(geneFile, stringsAsFactors = FALSE,  header = FALSE)
  #gene_info <- read.delim(geneFile, stringsAsFactors = FALSE,  header = FALSE)
  colnames(gene_info) <- c('EnsemblID', 'gene_short_name', 'biotype')
  
  ## Fix naming for mspecies. 
  for(c_pattern in c("GRCH38_", "Rnor6__")) {
    gene_info$EnsemblID <- gsub(c_pattern, "", gene_info$EnsemblID)
    gene_info$gene_short_name <- gsub(c_pattern, "", gene_info$gene_short_name)
  }
  rownames(gene_info) <- gene_info$EnsemblID
  #rownames(gene_info) <- gene_info$gene_short_name
  barcode_info <- barcodes <- read.delim(barcodeFile, stringsAsFactors = FALSE, header = FALSE)
  barcode_info$barcode <- barcode_info$V1
  barcode_info$V1 <- NULL
  rownames(barcode_info) <- barcode_info$barcode
  
  pd <- new("AnnotatedDataFrame", data=barcode_info)
  fd <- new("AnnotatedDataFrame", data=gene_info)
  
  rownames(m) <- gene_info$EnsemblID
  colnames(m) <- barcode_info$barcode
  
  cds <- newCellDataSet(m, phenoData = pd, featureData = fd)
  cds <- detectGenes(cds)
  cds
}

#'
#' Reads from `outs/filtered_feature_bc_matrix` directory 
#' 
get_monocle_cell_dataset_from_cellranger_filtered_outputs <- function(cellranger_output_dir, 
                                                                      additionalPhenoTypeList = NULL) {
  inner_dir <- paste0(cellranger_output_dir, '/outs/filtered_feature_bc_matrix/')
  if(file.exists(inner_dir)) {
    filtered_data_dir <- paste0(cellranger_output_dir, '/outs/filtered_feature_bc_matrix/')
  } else {
    filtered_data_dir <- cellranger_output_dir
  }
  
  matFile <- file.path(filtered_data_dir, 'matrix.mtx.gz')
  geneFile <- file.path(filtered_data_dir, 'features.tsv.gz')
  barcodeFile <- file.path(filtered_data_dir, 'barcodes.tsv.gz')
  cds <- get_monocle_cell_dataset_from_cellranger_files(matFile, geneFile, barcodeFile)
  
  ## additional phenotypic data to be added to monocle CellDataset
  pd <- phenoData(cds)
  if(!is.null(additionalPhenoTypeList)) {
    for(name in names(additionalPhenoTypeList)) {
      pd[[name]] <- additionalPhenoTypeList[[name]] 
    }
    phenoData(cds) <- pd
  }
  
  cds  
} 


check_directory_is_cellranger_output_directory <- function(filtered_data_dir) {
  matFile <- file.path(filtered_data_dir, 'matrix.mtx.gz')
  geneFile <- file.path(filtered_data_dir, 'features.tsv.gz')
  barcodeFile <- file.path(filtered_data_dir, 'barcodes.tsv.gz')
  file.exists(matFile) && file.exists(geneFile) && file.exists(barcodeFile)  
}

find_cellranger_outputs <-function(dataDir) {
  outDirs  <- Filter(check_directory_is_cellranger_output_directory, list.dirs(dataDir))
  names(outDirs) <- gsub(dataDir, "", outDirs)
  outRes <- lapply(outDirs, get_monocle_cell_dataset_from_cellranger_filtered_outputs)
}

#'
#' Function to get monocle cds to Seurat cds
#' 
convert_monocle_cds_to_seurat_cds <- function(monocle_cds) {
  monocle_cds <- fix_duplicated_gene_short_name(monocle_cds)
  m <- as.matrix(monocle_cds)
  names1 <- fData(monocle_cds)$gene_short_name 
  rownames(m) <- names1
  CreateSeuratObject(m, meta.data = pData(monocle_cds))
}

#'
#' Read Monocle cell datasets based on output of `read_cellranger_data_files`
#' 
#' @param sampleDf data.frame containing `sample`, `method` and `dataDir` 
#' @return list of monocle CellDataset
#' 
get_monocle_cds_from_sample_df <- function(sampleDf) {
  inner_function <- function(currDataDir, currPhenoData) {
    get_monocle_cell_dataset_from_cellranger_filtered_outputs(
      currDataDir, 
      additionalPhenoTypeList = currPhenoData
    )
  }
  result <- list()
  ## This is really bad code (should be done with apply-like functions but works for now)
  for(i in 1:nrow(sampleDf)) {
    currRow <- sampleDf[i,]
    currDataDir<- currRow[['dataDir']]
    sampleName <- currRow[['id']]
    flog.info("reading dataset: %s", currDataDir)
    result[[sampleName]] <- inner_function(currDataDir, currRow)
  }
  result
}

#' 
#' Function for retrieving a CellDataset from output of `get_monocle_cds_from_sample_df`
#' 
peek_monocle_cds_list <- function(monocleCdsList, i) {
  curr_cds <- monocleCdsList[[i]]
  curr_name <- names(monocleCdsList)[[i]]
  flog.info("sample: %s genes: %d barcodes: %d", curr_name, dim(curr_cds)[1], dim(curr_cds)[2])
  curr_cds
}

#' 
#' Function to find number of reads for each barcode in each sample 
#' 
get_sample_counts_from_seurat_cds_list <- function(seuratCdsList) { 
  sampleCountsList <-  lapply(seuratCdsList, function(x) { 
    a <- data.frame(barcode=x$barcode, reads=x$nCount_RNA) # sample=x$sample, assembly=x$assembly, method=x$method) 
    colnames(a) <- c('barcode', x$id[[1]])
    a
  })
  sampleCountsDf <- Reduce(function(x, y) merge(x, y, by="barcode", all=TRUE), sampleCountsList)
  rownames(sampleCountsDf) <- sampleCountsDf$barcode
  sampleCountsDf$barcode  <- NULL
  sampleCountsDf[is.na(sampleCountsDf)] <- 0
  sampleCountsDf
} 

get_cell_reads_from_monocle_cds <- function(x) { 
  z <- Matrix::colSums(as.matrix(x))
  a <- data.frame(barcode=names(z), reads=z) # sample=x$sample, assembly=x$assembly, method=x$method) 
  colnames(a) <- c('barcode', x$id[[1]])
  a
}

#' 
#' Function to find number of reads for each barcode in each sample 
#' 
get_cell_reads_from_monocle_cds_list <- function(monocleCdsList) { 
  sampleCountsList <-  lapply(monocleCdsList, get_cell_reads_from_monocle_cds)
  sampleCountsDf <- Reduce(function(x, y) merge(x, y, by="barcode", all=TRUE), sampleCountsList)
  rownames(sampleCountsDf) <- sampleCountsDf$barcode
  sampleCountsDf$barcode  <- NULL
  sampleCountsDf[is.na(sampleCountsDf)] <- 0
  sampleCountsDf
} 


#'
#' Function that returns samples in `sampleCountsDf` matching `c_pattern` restriction 
#' 
#' @return list(counts, correlation, df)
#' 
getCountsForSampleRestriction <- function(sampleCountsDf, c_pattern) {
  cDf <- sampleCountsDf[,str_detect(colnames(sampleCountsDf), c_pattern)]
  cDf <- cDf[rowSums(cDf) > 0,]
  # Do the same, but with colors corresponding to value
  a <- as.matrix(cDf)
  plot_colored_correlation(cor(a))
  c <- t(a) %*% a
  list(counts=c, correlation=cor(a), df=cDf)
}


