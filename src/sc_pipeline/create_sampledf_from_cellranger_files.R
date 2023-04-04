#!/usr/bin/env Rscript
# Usage: sc_pipeline/create_sampledf_from_cellranger_files.R [-[-help|h]] [-[-DATA_DIR|i] <character>] [-[-SAMPLEDF_FILE|o] [<character>]] [-[-PATTERN|p] [<character>]]
#     -h|--help             show this help message
#     -i|--DATA_DIR         input directory containing ./{method}/{sample}/outs/filtered_feature_bc_matrix/
#     -o|--SAMPLEDF_FILE    output file [default: sampleDf.tsv]
#     -p|--PATTERN          search pattern [default: outs/filtered_feature_bc_matrix/]

source("rcode/pipeline_helper_functions.R")

#get options, using the spec as defined by the enclosed list.
#we read the options from the default: commandArgs(TRUE).
spec = matrix(c(
  'help' , 'h', 0, "logical", "show this help message",
  'DATA_DIR' , 'i', 1, "character", "input directory containing ./{method}/{sample}/outs/filtered_feature_bc_matrix/",
  'SAMPLEDF_FILE' , 'o', 2, "character", "output file [default: sampleDf.tsv]",
  'PATTERN', 'p', 2, "character", "search pattern [default: outs/filtered_feature_bc_matrix/]"
), byrow=TRUE, ncol=5)
opt = getopt(spec)

DEFAULT_OPTS=list()
DEFAULT_OPTS[['SAMPLEDF_FILE']] = "sampleDf.tsv"
DEFAULT_OPTS[['PATTERN']]="outs/filtered_feature_bc_matrix"

opt <- set_default_options(opt, DEFAULT_OPTS)

# # if help was asked for print a friendly message
# # and exit with a non-zero error code
# if ( !is.null(opt$help) ) {
#   cat(getopt(spec, usage=TRUE))
#   q(status=1)
# }

# #set some reasonable defaults for the options that are needed,
# #but were not specified.
# if ( is.null(opt$verbose ) ) { flog.threshold(INFO) } else { flog.threshold(DEBUG) }
# if ( is.null(opt$SAMPLEDF_FILE ) ) { opt$SAMPLEDF_FILE="sampleDf.tsv" }
# if ( is.null(opt$PATTERN) ) { opt$PATTERN = "outs/filtered_feature_bc_matrix/" }

silently_load_package_list(c("getopt", "stringr", "data.table", "futile.logger"))

source("rcode/read_cellranger_data_files.R")

flog.info("looking for pattern: %s in dir: %s", opt$PATTERN, opt$DATA_DIR)
sampleDf <- find_cellranger_data_files(opt$DATA_DIR, opt$PATTERN)
sampleDf$id <- sampleDf$sample
write_sampledf(sampleDf, opt$SAMPLEDF_FILE)
