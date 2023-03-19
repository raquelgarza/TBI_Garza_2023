# Usage: scripts/do_quality_filtering.R [-[-verbose|v] [<integer>]] [-[-help|h]] [-[-DATA_DIR|i] <character>] [-[-OUT_DIR|o] <character>] [-[-LOG_FILE|l] [<character>]] [-[-FILTERED_SAMPLEDF_FILE|f] [<character>]] [-[-RESULTS_FILE|r] [<character>]]
# -v|--verbose                   more logging (DEBUG) [default: FALSE]
# -h|--help                      show this help message
# -i|--DATA_DIR                  input directory containing ./{method}/{sample}/outs/filtered_feature_bc_matrix/
# -o|--OUT_DIR                   output directory where {sampleDf.tsv, qa-results.rda}
# -l|--LOG_FILE                  log-file [default: OUT_DIR/cell_quality_analysis.log]
# -f|--FILTERED_SAMPLEDF_FILE    filtered sample file [default: OUT_DIR/sampleDf-quality-filtered.tsv]
# -r|--RESULTS_FILE              results file [default: OUT_DIR/qa-results.rda]

# @see https://f1000research.com/articles/5-2122/v2
# @see https://stackoverflow.com/questions/51295402/r-on-macos-error-vector-memory-exhausted-limit-reached 

source("rcode/pipeline_helper_functions.R")

DEFAULT_OPTIONS=list()
DEFAULT_OPTIONS["FILTER_MITO_NMADS"]=3
DEFAULT_OPTIONS["FILTER_LOW_NUM_GENES_NMADS"]=3
DEFAULT_OPTIONS["FILTER_LOW_NUM_GENES_TYPE"]="both"
DEFAULT_OPTIONS["FILTER_MIN_NUM_GENES_EXPRESSED"]=0
DEFAULT_OPTIONS["RESULTS_FILE"]="qa-results.rds"
DEFAULT_OPTIONS["FILTERED_SAMPLEDF_FILE"]="qa-SampleDf.rds"
DEFAULT_OPTIONS['LOG_FILE'] <- "do_quality_filtering.log"

#get options, using the spec as defined by the enclosed list.
#we read the options from the default: commandArgs(TRUE).
spec = matrix(c(
  'SAMPLEDF_FILE' , 'idf', 1, "character", "input sampleDf file containing full-paths",
  'RDS_FILE', 'irds', 1, 'character', 'input rds file', 
  'OUT_DIR' , 'o', 1, "character", "output directory where {sampleDf.tsv, qa-results.rda}", 
  'LOG_FILE' , 'l' , 2, "character", paste0("log-file default: OUT_DIR/", DEFAULT_OPTIONS$LOG_FILE), 
  'FILTERED_SAMPLEDF_FILE', 'f', 2, "character", paste0("filtered sample file default: OUT_DIR/",DEFAULT_OPTIONS$FILTERED_SAMPLEDF_FILE), 
  'RESULTS_FILE', 'r', 2, "character", paste0("results file default: OUT_DIR/", DEFAULT_OPTIONS$RESULTS_FILE),
  'FILTER_MITO_NMADS', 'm', 2, "integer", paste0("Filter cells having number mean absolute deviation (default: ", 
                                                         DEFAULT_OPTIONS$FILTER_MITO_NMADS,
                                                         ") higher than the median fraction of mitochondrial reads across cells using pattern ^MT in gene names. Set this to 0 or below to skip this step"),
  'FILTER_LOW_NUM_GENES_NMADS', 'n', 2, "integer", paste0("Filter cells having nmads (default: ", 
                                                                    DEFAULT_OPTIONS$FILTER_LOW_NUM_GENES_NMADS,
                                                                    ") higher than the median number of expressed genes on",
                                                                    DEFAULT_OPTIONS$FILTER_LOW_NUM_GENES_TYPE, " side"),
  'FILTER_MIN_NUM_GENES_EXPRESSED', 'min_genes', 2, "integer", "Filter to remove cells having less than this many genes (default: 0)",
  'help' , 'h', 0, "logical", "show this help message"
), byrow=TRUE, ncol=5)
opt = getopt(spec)


opt <- set_default_options(opt, DEFAULT_OPTIONS, get_Rscript_filename())

PKGS = c("monocle", "stringr", "reshape2", "data.table", "ggplot2", "gplots", "simpleSingleCell", "scater")
silently_load_package_list(PKGS)


outDir <- opt$OUT_DIR
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)


## Load helper  functions 
source("rcode/cellranger_helper_functions.R")
source("rcode/quality_filtering_functions.R")
source("rcode/plotting_functions.R")

## Main script starts here
sampleDf.minimal <- read_sampledf(opt$SAMPLEDF_FILE) 
monocleCdsList <- readRDS(opt$RDS_FILE)

flog.info("monocleCdsList: %s", paste(names(monocleCdsList), collapse=","))
saveRDS(sampleDf.minimal, file=opt$FILTERED_SAMPLEDF_FILE)

## do mito-filtering
if(opt$FILTER_MITO_NMADS <= 0) { 
    flog.warn("\nskipping filtering of mito genes FILTER_MITO_NMADS: %d\n", opt$FILTER_MITO_NMADS)
} else {
    flog.info("Filtering mito genes with nmads=%s",opt$FILTER_MITO_NMADS)

    monocleCdsList <- lapply(names(monocleCdsList), function(n) {
     x = monocleCdsList[[n]]
     filter_mitochondrial_genes_and_cells(x, outDir=outDir, nmads=as.numeric(opt$FILTER_MITO_NMADS), currId=n)
    })
    names(monocleCdsList) <- sampleDf.minimal$id
    sampleDf.minimal$isNonMitoNumCells <- sapply(monocleCdsList, function(cds1) { dim(cds1)[2]} )
    write_sampledf(sampleDf.minimal, file=opt$FILTERED_SAMPLEDF_FILE)
    saveRDS(monocleCdsList, file=opt$RESULTS_FILE)
    
    flog.info("saving minimal plot")
    make_merged_plot(outDir, "mito")
    
}
 
## do low-num-genes filtering
flog.info("filter_cells_with_low_number_of_genes()")
monocleCdsList <- sapply(names(monocleCdsList), function(n) {
  x = monocleCdsList[[n]]
  filter_cells_with_low_number_of_genes(x, outDir=outDir,
                                        nmads = as.numeric(opt$FILTER_LOW_NUM_GENES_NMADS),
                                        type = opt$FILTER_LOW_NUM_GENES_TYPE,
                                        minGenes = as.numeric(opt$FILTER_MIN_NUM_GENES_EXPRESSED),
                                        currId=n
                                        ) })
names(monocleCdsList) <- sampleDf.minimal$id
sampleDf.minimal$lowFilteredNumCells <- sapply(monocleCdsList, function(cds1) { dim(cds1)[2]} )
saveRDS(sampleDf.minimal, file=opt$FILTERED_SAMPLEDF_FILE)
saveRDS(monocleCdsList, file=opt$RESULTS_FILE)


flog.info("remove genes that are low-expressed across all samples")
genesToFilter <- lapply(names(monocleCdsList), function(n) {
    x = monocleCdsList[[n]]
    get_genes_expressed_in_low_number_of_cells(x, outDir=outDir, currId=n) 
})
reducedSetOfGenes <- Reduce(get_merged_genes_to_filter, genesToFilter)

#' Internal function that reduces specific cds in monocleCdsList
#'
#' @param n index of cds to reduce
#' @param reduced_gene_list limit cds to genes present in this list
#' @return cds with genes present in \code{reduced_gene_list}
#'
get_reduced_cds <- function(n, reduced_gene_list) {
 cds1 = monocleCdsList[[n]]
 flog.info("id: %s", n)
 c_genes <- fData(cds1)
 c_genes_1 <- merge(c_genes, reduced_gene_list, by=c("EnsemblID", "gene_short_name"))
 zz <- as.logical(c_genes_1$filter)
 flog.info("filtering genes: %d remaining: %d", sum(zz), sum(!zz))
 cds1[!zz, ]
}

monocleCdsList <- sapply(names(monocleCdsList), function(x) get_reduced_cds(x, reducedSetOfGenes) )
sampleDf.minimal$reducedNumGenes <- sapply(monocleCdsList, function(cds1) { dim(cds1)[1]} )
names(monocleCdsList) <- sampleDf.minimal$id
write_sampledf(sampleDf.minimal, file=opt$FILTERED_SAMPLEDF_FILE)
saveRDS(monocleCdsList, file=opt$RESULTS_FILE)

flog.info("saving gene plot")
make_merged_plot(outDir, "genes")

#flog.info("saving cell plot")
#make_merged_plot(outDir, "cells")


