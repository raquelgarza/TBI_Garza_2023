# This script removes cells from GRCh38 samples

source("rcode/pipeline_helper_functions.R")

silently_load_package_list(c("getopt", "tools"))

ASSEMBLY=c("GRCh38", "Hg38")

spec = matrix(c(
  'IN_SAMPLEDF_FILE' , 'idf', 1, "character", "input sampleDf file",
  'IN_RDS_FILE', 'irds', 1, "character", "input RDS file containing list of monocleDatasets", 
  'HUMAN_GENE_IDS_FILE', 'g', 1, "character", "file naming human Gene IDS", 
  'HUMAN_READS_MIN_RATIO', 'minR', 2, 'numeric', "min-threshold for fraction of human-reads/total (default: 0.9)", 
  'HUMAN_READS_MAX_RATIO', 'maxR', 2, 'numeric', "max-threshold for fraction of human-reads/total (default: 1.0)", 
  'OUT_DIR' , 'odir', 1, "character", "output directory", 
  'OUT_SAMPLEDF_FILE', 'odf', 1, "character", "output sampledf rds file", 
  'OUT_RDS_FILE', 'ords', 1, "character", "output cellDataset list rds file",
  'FILTER_MODE', 'fm', 2, "logical", "whether to only return Hg38/GRCh38 assembly as output.\n\tIf T, filter to Hg38 using sample as the merge column,\n\tElse, return full list with entries filtered for that dataset.",
  'LOG_FILE', 'l', 2, "character", "where output should be logged",
  'help' , 'h', 0, "logical", "show this help message"
), byrow=TRUE, ncol=5)
opt = getopt(spec)

DEFAULT_OPTIONS=list()
DEFAULT_OPTIONS["LOG_FILE"]="remove_contaminated_cells.log"
DEFAULT_OPTIONS["HUMAN_READS_MIN_RATIO"]=0.9
DEFAULT_OPTIONS["HUMAN_READS_MAX_RATIO"]=1.0 
DEFAULT_OPTIONS["FILTER_MODE"]=FALSE
DEFAULT_OPTIONS["OUT_RDS_FILE"] ="remove_contaminated_cells.rds"
DEFAULT_OPTIONS["OUT_SAMPLEDF_FILE"] = "remove_contaminated_cells.tsv"

opt <- set_default_options(opt, DEFAULT_OPTIONS, get_Rscript_filename())



PKGS = c("stringr", "monocle", "cowplot", "reshape2", "data.table", "ggplot2", "gplots", "simpleSingleCell")
silently_load_package_list(PKGS)

flog.info("Executing: %s", get_Rscript_filename())
source("rcode/cellranger_helper_functions.R")
source("rcode/plotting_functions.R")
source("rcode/quality_filtering_functions.R")


# ### Current Testing Paths
# opt$OUT_DIR <- "/Users/yo4230sh/AllAnalysis_Files/asyn_project_out/remove_contaminated"
# opt$IN_RDS_FILE <- "/Users/yo4230sh/AllAnalysis_Files/asyn_project_out/unprocessed.rds"
# opt$IN_SAMPLEDF_FILE <- "/Users/yo4230sh/AllAnalysis_Files/asyn_project_out/sampleDf.tsv"
# opt$HUMAN_GENE_IDS_FILE <-  "/Users/yo4230sh/DropboxDNLP/GoogleDrive_Backup/Yogita/graftsCellAnalysis/DataFiles/humanGeneIDs.rds"

dir.create(opt$OUT_DIR, showWarnings = F, recursive = T)




## Inputs
humanGeneIDs <- readRDS(opt$HUMAN_GENE_IDS_FILE)

sampleDf <- read_sampledf(opt$IN_SAMPLEDF_FILE)

# ## check that sampleDf has the necessary columns
# for(necessaryColNames in c("sample", "id", "assembly")) {
#     if (opt$FILTER_MODE && (!(necessaryColNames %in% colnames(sampleDf)))) {
#         flog.error("In FILTER_MODE, sampleDf should have column: %s", necessaryColNames)
#         quit(status=1)
#     }
# }



monocleCdsList <- readRDS(opt$IN_RDS_FILE)


samples <- unique(sampleDf$sample)

flog.info("colNames: %s", paste(colnames(fData(monocleCdsList[[1]])), collapse=","))

monocleCdsList <- lapply(monocleCdsList, function(x) {
    x <- cleanGeneNamingForCellRanger(x, columnName="EnsemblID", pattern="premRNA", replacement="")
    x <- cleanGeneNamingForCellRanger(x, columnName="gene_short_name", pattern="premRNA", replacement="")
    x
    }) 

## using multi-assembly reads to identify human only cells

if(opt$FILTER_MODE) {
    flog.info("finding Hs_vs_nonHs only on \"Mspecies\" samples")
    mDf <- sampleDf[sampleDf$assembly=="Mspecies", ]
} else {
    flog.info("sampleDf$assembly missing or not containing `Mspecies` - proceeding with all samples")
    mDf <- sampleDf
}

mrList <- list()
mId <- mDf$id[[1]]




for(mId in mDf$id) {
  cSample <- mId
  cds1 <- monocleCdsList[[mId]]
  a.f <- fData(cds1)
  flog.info("mId: %s cSample: %s cds: %dx%d", mId, cSample, dim(cds1)[1], dim(cds1)[2])

  # try prefix information
  humanIdx <- which(stringr::str_starts(a.f$EnsemblID, "grch38"))
  flog.info("humanGenes found based on prefix search: %d", length(humanIdx))

  if(length(humanIdx) == 0) {
    flog.info("writing to match to humanGeneIDs")
    humanIdx <- match(humanGeneIDs, a.f$EnsemblID)
    humanIdx <- humanIdx[!is.na(humanIdx)]
  }

  ratIdx <- setdiff(1:dim(cds1)[1], humanIdx)

  flog.info("humanIdx: %d ratIdx: %d",  length(humanIdx), length(ratIdx))
  humanReads <- Matrix::colSums( as.matrix(cds1[humanIdx, ]) )
  ratReads <-Matrix::colSums( as.matrix(cds1[ratIdx, ]) )
  flog.info("%s total reads: human: %d rat: %d ", cSample, sum(humanReads), sum(ratReads))
  cDf <- data.frame(barcode=names(humanReads), GRCh38=humanReads, Rnor6=ratReads, ratio=humanReads/(humanReads+ratReads))
  cDf$Hs <- (cDf$ratio >= as.numeric(opt$HUMAN_READS_MIN_RATIO)) & (cDf$ratio <= as.numeric(opt$HUMAN_READS_MAX_RATIO))
  flog.info("%s num_human_cells: %d/%d", cSample, sum(cDf$Hs), length(cDf$Hs))
  mrList[[cSample]] <- cDf
}

mrPlots <- mapply(function(x, n) {
  cTitle = paste0(n, " (",  sum(x$Hs), "/",length(x$Hs), ")")
  p <- ggplot(x, aes(x=GRCh38, y=Rnor6, color=ratio, shape=Hs))+geom_point()
  p <- p + scale_x_continuous(trans = "log10") + scale_y_continuous(trans="log10")
  p <- p + labs(title=cTitle)
  ggsave(plot=p, filename = file.path(opt$OUT_DIR, paste0(n, ".eps")))
  p
}, mrList, names(mrList), SIMPLIFY = FALSE)

plot_grid(plotlist=mrPlots, ncol=1, labels = "auto", align='v')
ggsave(file.path(opt$OUT_DIR, "mspecies-rat-vs-human.eps"), units="in", dpi=600)
saveRDS(list("mrList"=mrList, "mrPlots"=mrPlots), file=file.path(opt$OUT_DIR, "mspecies-analysis.rds"))

filter_to_accepted_cells <- function(cds1, hsDf) {
  b_cds1 <- pData(cds1)$barcode
  b_accepted <- rownames(hsDf)[hsDf$Hs]
  idx<- b_cds1 %in% b_accepted
  flog.info("cds1: %d mspecies-accepted: %d final: %d", length(b_cds1), length(b_accepted), sum(idx))
  cds1[, idx]
}

if(opt$FILTER_MODE) { ## filter to Hg38 using sample as the merge column
    names(mrList) <- mDf$sample
    filteredSampleDf <- sampleDf[sampleDf$assembly %in% ASSEMBLY, ]
    get_mr_list_name <- function(i) { filteredSampleDf$sample[i] }
} else { ## return full list with entries filtered for that dataset
    filteredSampleDf <- sampleDf
    get_mr_list_name <-  function(i) { filteredSampleDf$id[i] }
}

filteredCdsList <- list()
for(i in 1:nrow(filteredSampleDf)) {
    # flog.warn("Filtering to assembly: %s", filteredSampleDf$assembly[i])
    cId <- filteredSampleDf$id[i]
    cSample <- get_mr_list_name(i)
    hsDf <- mrList[[as.character(cSample)]]
    cds1 <- monocleCdsList[[cId]]
    flog.info("cId: %s cSample: %s accepted: %d", cId, cSample, sum(hsDf$Hs))
    filteredCdsList[[cId]] <- filter_to_accepted_cells(cds1, hsDf)
} 

filteredSampleDf$decontaminatedNumCells <- sapply(filteredCdsList, function(x) { dim(x)[2] })
filteredSampleDf$decontaminatedNumReads <- sapply(filteredCdsList, function(x) { sum(as.matrix(x)) })

dfm <- filteredSampleDf[,c("sample",  "origNumCells", "decontaminatedNumCells")]
dfm$origNumCells <- dfm$origNumCells - dfm$decontaminatedNumCells
colnames(dfm) <- c("Sample",  "Rat", "Human")
dfm <- melt(dfm, id.vars = 1)
p <- ggplot(data=dfm, aes(y=Sample, x=value)) + geom_bar(aes(fill = variable),stat = "identity",position = "dodge")
p <- p +  theme(axis.text.y = element_text(size = 6)) + labs(x="Number of cells")
p
ggsave(file.path(opt$OUT_DIR, "num_decontaminated_cells.png"), p, limitsize = FALSE, width=15, height=8, units="cm")
names(filteredCdsList) <- filteredSampleDf$id
write_sampledf(filteredSampleDf, opt$OUT_SAMPLEDF_FILE)
flog.info("Finishing %s", get_Rscript_filename())



## commenting this out here to handle rnor6 removal for the asyn project and saving downstream
## for new workflows this would happen here instead of glue code at end of project
# saveRDS(filteredCdsList, file=opt$OUT_RDS_FILE)


## glue code to remove rnor6 genes from downstream processing
# filteredCdsList <- readRDS("~/AllAnalysis_Files/asyn_project_out/decontaminated/decontaminated.rds")

remove_rnor6_genes <- function(test_cds) {
  humanIdx <- which(stringr::str_starts(fData(test_cds)$EnsemblID, "grch38"))
  cellData <- exprs(test_cds)
  cellData <- cellData[humanIdx, ]
  new_rownames <- stringr::str_replace(rownames(cellData), "grch38_", "")
  rownames(cellData) <- new_rownames
  newFeatureData <- featureData(test_cds)
  newFeatureData <- newFeatureData[humanIdx, ]
  rownames(newFeatureData) <- new_rownames
  newFeatureData$EnsemblID <- new_rownames
  newFeatureData$gene_short_name <- stringr::str_replace(newFeatureData$gene_short_name, "grch38_", "")
  new_cds <- newCellDataSet(
    cellData,
    phenoData = phenoData(test_cds),
    featureData = newFeatureData,
  )
  flog.info("removed rnor6 here")
  new_cds
}

updated_cds_list <- list()
for(i in 1:length(filteredCdsList)) {
  flog.info("%d/%d %s", i, length(filteredCdsList), names(filteredCdsList)[[i]])
  updated_cds_list[[i]] <- remove_rnor6_genes(filteredCdsList[[i]])
}
names(updated_cds_list) <- names(filteredCdsList)

# saveRDS(updated_cds_list, file=file.path("~/AllAnalysis_Files/asyn_project_out/decontaminated", "decontaminated_rnor6_genes_removed.rds"))
saveRDS(updated_cds_list, file=opt$OUT_RDS_FILE)
