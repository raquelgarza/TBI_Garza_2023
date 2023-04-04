#!/usr/bin/env Rscript

# This file runs panglodb enrichment by running do_seurat_enrichment function and providing
# panglodb enrichment

# FIXME This needs to be separated from Seurat run

# Usage: sc_pipeline/do_panglodb_enrichment.R [-[-PROCESSED_CDS_LIST|i] <character>] [-[-OUT_DIR|o] <character>] [-[-PANGLODB_DIR|p] <character>] [-[-LOG_FILE|l] [<character>]] [-[-help|h]]
#
# -i|--PROCESSED_CDS_LIST    RDS_File containing processed CDS list of which each item has the following attributes 
#                                 * scds: Seurat Cell Dataset 
#                                 * markers: List of significant markers
# -o|--OUT_DIR               where results should be saved
# -p|--PANGLODB_DIR          panglodb_directory
# -l|--LOG_FILE              where output should be logged
# -h|--help                  show this help message

ENRICHMENT_OUTPUT_FILE_NAME='enrichment.rds' 

thisFile <- function() {
        cmdArgs <- commandArgs(trailingOnly = FALSE)
        needle <- "--file="
        match <- grep(needle, cmdArgs)
        if (length(match) > 0) {
                # Rscript
                return(normalizePath(sub(needle, "", cmdArgs[match])))
        } else {
                # 'source'd via R console
                return(normalizePath(sys.frames()[[1]]$ofile))
        }
}
baseDirectory <- dirname(dirname(thisFile()))
cat("Library source directory:", baseDirectory, "\n")

source("rcode/pipeline_helper_functions.R")
spec = matrix(c(
  'PROCESSED_CDS_LIST' , 'i', 1, "character", "RDS_File containing processed CDS list of which each item has the following attributes - scds: Seurat Cell Dataset - markers: List of significant markers",
  'OUT_DIR', 'o', 1, "character", "where results should be saved", 
  'PANGLODB_DIR', 'p', 1, "character", "panglodb_directory containing PangloDB *.tsv files ", 
  'LOG_FILE', 'l', 2, "character", "where output should be logged",
  'help' , 'h', 0, "logical", "show this help message"
), byrow=TRUE, ncol=5)
opt = getopt(spec)

DEFAULT_OPTIONS = list()
DEFAULT_OPTIONS[["LOG_FILE"]]="do_panglodb_enrichment.log"

opt <- set_default_options(opt, DEFAULT_OPTIONS, get_Rscript_filename())

if(is.null(opt$PROCESSED_CDS_LIST) || is.null(opt$OUT_DIR) || is.null(opt$PANGLODB_DIR)) {
  cat(getopt(spec, usage=TRUE))
  q(status=1)
}
source("rcode/plotting_functions.R")
source("rcode/enrichment_functions.R")

flog.info("opt: %s", toString(opt))
processedCdsList <- readRDS(opt$PROCESSED_CDS_LIST)

# Read PangloDB directory containing *.tsv files 
combinedCellTypeDt <- readPangloDbFilesDir(opt$PANGLODB_DIR)

doEnrichment <- function( processedCdsList, opt, combinedCellTypeDt) {
  enrichmentResult <- lapply(processedCdsList, function(result) { 
    cSample <- result$scds@project.name
    currOutDir <- file.path(opt$OUT_DIR, cSample)
    dir.create(currOutDir, showWarnings = F, recursive = T)
    flog.info("cSample: %s outDir: %s", cSample, currOutDir)
    do_seurat_enrichment(result$scds, result$markers, combinedCellTypeDt, outDir=currOutDir, cSample=cSample)
  })
  saveRDS(enrichmentResult, file.path(opt$OUT_DIR, "enrichmentResult.rds"))
  enrichmentResult
} 

doEnrichment(processedCdsList, opt, combinedCellTypeDt)
