#!/usr/bin/env Rscript

source("rcode/pipeline_helper_functions.R")

silently_load_package_list(c(
  "argparse", "futile.logger", "Seurat", "tools"
))


get_parser <- function() {
  # initial reading of arguments
  parser <- ArgumentParser()
  parser$add_argument("sample_df", type="character", help="Path to sample_df file")
  parser$add_argument("out_rds", type="character", help="Where to save the output")
  parser$add_argument("--out_sample_df", type="character", help="Where to save the output sample_df file", default=NULL)
  parser
}

write_sampledf <- function(sampleDf, file) {
    flog.info("written %d samples columns: (%s) to file: %s", dim(sampleDf)[1], paste(colnames(sampleDf), collapse=','), file)
    switch(file_ext(file),
           "rds"=saveRDS(sampleDf, file=file),
           "tsv"=write.table(sampleDf, file=file, sep="\t", row.names=F, col.names=T, quote=F)
    )
}

#' Read list of seurat objects from cellranger files
#'
#' @param sample_df_file Input file to read containing sample information (should have columns `dataDir` and `id`)
#' @param out_rds_file Where to write the RDS file containing list of Seurat objects
#' @param out_sample_df Where to write the output sampleDf file
#'
#' @returns List of Seurat Objects
#'
read_10x_using_seurat <- function(sample_df_file, out_rds_file, out_sample_df=NULL) {
  out_dir <- dirname(out_rds_file)
  dir.create(out_dir, recursive = T, showWarnings = F)
  flog.info("created %s: %d", out_dir, dir.exists(out_dir))
  # send everything to log file as well.
  log_file <- file.path(out_dir, "read_10x_using_seurat.log")
  flog.appender(name='ROOT', appender.tee(log_file))


  flog.info("trying to read: %s", sample_df_file)
  sampleDf <- read.csv(sample_df_file, sep="\t", header=T)
  flog.info("sampleDf has %d samples", dim(sampleDf)[1])


  monocleCdsList.orig <- list()
  for(i in 1:nrow(sampleDf)) {
    cId <- sampleDf$id[i]
    flog.info("reading sample: %s ", cId, sampleDf$sample[i])

    filtered_data_dir <- sampleDf$dataDir[i]
    inner_dir <- file.path(filtered_data_dir, '/outs/filtered_feature_bc_matrix/')
    if(file.exists(inner_dir)) {
      filtered_data_dir <- inner_dir
    }
    counts <- Read10X(data.dir = filtered_data_dir)
    cds <- CreateSeuratObject(counts = counts)
    cds@project.name <- cId
    cds@meta.data$sample <- cId

    flog.info("cds sample: %s", cId)
    monocleCdsList.orig[[cId]] <- cds
  }
  names(monocleCdsList.orig) <- sampleDf$id


  sampleDf$origNumCells <- sapply(monocleCdsList.orig, function(x) { dim(GetAssayData(x))[2] })
  sampleDf$origNumReads <- sapply(monocleCdsList.orig, function(x) { sum(x@meta.data$nCount_RNA) })

  if(is.null(out_sample_df)) {
    out_sample_df <- file.path(out_dir, "sampleDf.tsv")
  }
  write_sampledf(sampleDf, file = out_sample_df)
  saveRDS(monocleCdsList.orig, file=out_rds_file)

  monocleCdsList.orig
}


## Here is the main flow - this part is only run when called from terminal (using Rscript) - not when sourced
if(!interactive()) {

  layout <- layout.format(paste0('~l [~t] :~f() ~m'))
  flog.layout(layout)

  parser <- get_parser()
  parsed_args <- parser$parse_args()
  outCdsList <- read_10x_using_seurat(parsed_args$sample_df, parsed_args$out_rds, parsed_args$out_sample_df)
}
