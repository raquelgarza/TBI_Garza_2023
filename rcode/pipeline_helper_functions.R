silent_load <- function(PKG) {
  suppressWarnings(suppressMessages(suppressPackageStartupMessages(library(PKG, quietly = T, character.only = T))))
}

# pipeline helper functions
suppressMessages(library("getopt"))
suppressMessages(library("futile.logger"))
suppressMessages(library("stringr"))
suppressMessages(library("tools"))

update_logger <- function(calling_script, logfile=NULL) {
  calling_script = basename(calling_script)
  layout <- layout.format(paste0('~l [~t] ', calling_script, ':~f() ~m'))
  flog.layout(layout)
  if(!is.null(logfile)) {
    dir.create(dirname(logfile), showWarnings = F, recursive = T)
    flog.appender(name='ROOT', appender.tee(logfile))
  }
} 

set_default_options <- function(opt, DEFAULT_OPTIONS, calling_script=get_Rscript_filename()) {
  
  # if help was asked for print a friendly message
  # and exit with a non-zero error code
  if ( !is.null(opt$help) ) {
    cat(getopt(spec, usage=TRUE))
    q(status=1)
  }
  
  for( key in names(DEFAULT_OPTIONS)) {
    if(is.null(opt[[key]])) {
      opt[[key]] = DEFAULT_OPTIONS[[key]]
      flog.info("setting %s to: %s", key, opt[[key]])
    }
  }
  
  ## better logging
  if ( is.null(opt$verbose ) ) { flog.threshold(INFO) } else { flog.threshold(DEBUG) }
  if ( is.null(opt$LOG_FILE) ) { opt$LOG_FILE = paste0( calling_script, ".log") }
  update_logger(calling_script, opt$LOG_FILE)
  opt
}

read_sampledf <- function(fileName) {
    flog.info("reading file: %s", fileName)
    sampleDf <- switch(file_ext(fileName), 
                     "rds" = readRDS(fileName),
                     "tsv" = read.csv(fileName, sep="\t", header=T))
    flog.info("read %d samples columns: (%s)", dim(sampleDf)[1], paste(colnames(sampleDf), collapse = ','))
    sampleDf
}

write_sampledf <- function(sampleDf, file) {
    flog.info("written %d samples columns: (%s) to file: %s", dim(sampleDf)[1], paste(colnames(sampleDf), collapse=','), file)  
    switch(file_ext(file),
           "rds"=saveRDS(sampleDf, file=file),
           "tsv"=write.table(sampleDf, file=file, sep="\t", row.names=F, col.names=T, quote=F)
    )
}

silently_load_package_list <- function(PKGS) {
  for(PKG in PKGS) { 
    silent_load(PKG)
    flog.info("loading library: %s", PKG)
  }
}

