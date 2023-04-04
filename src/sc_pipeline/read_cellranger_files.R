source("rcode/pipeline_helper_functions.R")
source("rcode/cellranger_helper_functions.R")

#get options, using the spec as defined by the enclosed list.
#we read the options from the default: commandArgs(TRUE).
spec = matrix(c(
    'SAMPLE_DF_FILE' , 'i', 1, "character", "input sampleDf file containing full-paths",
    'LOG_FILE', 'l', 2, "character", "where output should be logged (default: OUT_DIR/read_cellranger_files.log)",
    'OUT_RDS_FILE', 'f', 1, "character", "output rds file containing read cellranger datasets (default: OUT_DIR/unprocessed.rds)", 
    'OUT_SAMPLEDF_FILE', 'o', 2, "character", "where to save outputSampleDf file with cells/reads info",
    'help' , 'h', 0, "logical", "show this help message"
), byrow=TRUE, ncol=5)
opt <- getopt(spec)

opt$OUT_DIR <- dirname(opt$OUT_RDS_FILE)
dir.create(opt$OUT_DIR, recursive = T, showWarnings = F)

DEFAULT_OPTIONS <- list()
DEFAULT_OPTIONS["LOG_FILE"] <- file.path(opt$OUT_DIR, "read_cellranger_files.log")
DEFAULT_OPTIONS["OUT_SAMPLEDF_FILE"] <- opt$SAMPLE_DF_FILE

opt <- set_default_options(opt, DEFAULT_OPTIONS, get_Rscript_filename())

## Main script starts here
flog.info("trying to read: %s", opt$SAMPLE_DF_FILE)
sampleDf <- read.csv(opt$SAMPLE_DF_FILE, sep="\t", header=T)
flog.info("sampleDf has %d samples", dim(sampleDf)[1])

monocleCdsList.orig <- list() 
for(i in 1:nrow(sampleDf)) {
    cId <- sampleDf$id[i]
    flog.info("reading sample - id: %s sample: %s", cId, sampleDf$sample[i])
    cds <- get_monocle_cell_dataset_from_cellranger_filtered_outputs(sampleDf$dataDir[i])
    pData(cds)$sample <- cId
    flog.info("cds sample: %s", cId)
    monocleCdsList.orig[[cId]] <- cds 
}
names(monocleCdsList.orig) <- sampleDf$id

sampleDf$origNumCells <- sapply(monocleCdsList.orig, function(x) { dim(x)[2] })
sampleDf$origNumReads <- sapply(monocleCdsList.orig, function(x) { sum(as.matrix(x)) })

write_sampledf(sampleDf, opt$OUT_SAMPLEDF_FILE)
saveRDS(monocleCdsList.orig, file=opt$OUT_RDS_FILE)
