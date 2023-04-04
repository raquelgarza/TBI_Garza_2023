start_time <- Sys.time()

source("rcode/pipeline_helper_functions.R")

silently_load_package_list(c(
  "monocle", "Seurat", "ggplot2", "cowplot", "data.table", "reshape2", "stringr",
  "dplyr", "jsonlite", "harmony" # , "clusterProfiler"
))


spec <- matrix(c(
  'IN_SAMPLEDF_FILE' , 'idf', 1, "character", "input sampleDf file",
  'IN_RDS_FILE', 'irds', 1, "character", "input RDS file containing list of monocleDatasets", 
  'OUT_RDS_FILE', 'ords', 1, "character", "output cellDataset list rds file",
  'MERGE_JSON_FILE', 'm', 1, "character", "JSON input",
  'LOG_FILE', 'l', 2, "character", "where output should be logged",
  'RESOLUTION', 'r', 1, "double", "resolution for Seurat Clustering",
  'help' , 'h', 0, "logical", "show this help message"
), byrow=TRUE, ncol=5)
opt <- getopt(spec)


update_logger(get_Rscript_filename(), logfile = opt$LOG_FILE)


source("rcode/cellranger_helper_functions.R")
source("rcode/plotting_functions.R")
source("rcode/seurat_functions.R")
source("rcode/enrichment_functions.R")

inputCdsList <- readRDS(opt$IN_RDS_FILE)
sampleDf <- read_sampledf(opt$IN_SAMPLEDF_FILE)
mergeJson <- read_json(opt$MERGE_JSON_FILE)

OUT_DIR <- dirname(opt$OUT_RDS_FILE)

if( is.null(opt$RESOLUTION)) { opt$RESOLUTION = 0.8 }
flog.info("mergeJson: %s", toString(mergeJson))


default_pipeline_options <- get_seurat_pipeline_options()

updateProjectNamesForSeuratCdsList <- function(seuratCdsList) {

  cNames <- names(seuratCdsList)
  seuratCdsList <- lapply(cNames, function (x) {
  seuratCdsList[[x]]@project.name <- x
  seuratCdsList[[x]]
  })
  names(seuratCdsList) <- cNames
  seuratCdsList
}



selectValidSamples <- function(x) {
  validIds <- mergeSamples[[x]]
  p1 <- sampleDf %>% filter(id %in% validIds) %>% select(sample, id)
  p1
  cdsList <- seuratCdsList[ p1$id ]
  names(cdsList) <- p1$id
  cdsList
}

#'
#' Plot and save cluster+sample information side-by-side
#'
my_dim_plot <- function(pbmc, currSample='', seuratOutDir='.', reduction='umap',  l_size=7, width=16, height=8) {
  p <- DimPlot(pbmc, reduction = reduction, group.by = 'sample') + theme(legend.position="bottom" , legend.text=element_text(size=l_size)) + ggtitle("sample")
  q <- DimPlot(pbmc, reduction = reduction, group.by = 'ident', label = T) + theme(legend.position="bottom", legend.text=element_text(size=l_size))+ggtitle("cluster")
  r <- CombinePlots(list(q, p), ncol = 2)
  ggsave(file.path(seuratOutDir, paste0(reduction, "-", currSample ,".pdf")), plot=r, width=width, height=height)
  r
}


get_cds_and_markers <- function(pbmc, currSample, resolution=opt$RESOLUTION, seuratOutDir=OUT_DIR, pipeline_options=default_pipeline_options) {

  currentClusterName = paste0("umap_", resolution)
  pbmc@meta.data[[currentClusterName]] <- pbmc@meta.data[["seurat_clusters"]]
  my_dim_plot(pbmc, reduction="umap", seuratOutDir=seuratOutDir, currSample = currSample)

  flog.info("do differential expression for each cluster")
  pbmc.markers <- FindAllMarkers(object = pbmc,
                                 only.pos = pipeline_options[["FindAllMarkers"]][["only.pos"]],
                                 min.pct = pipeline_options[["FindAllMarkers"]][["min.pct"]],
                                 thresh.use = pipeline_options[["FindAllMarkers"]][["thresh.use"]],
                                 verbose = T)

  reduced.markers <- pbmc.markers %>% group_by(cluster)

  markersFile <- file.path(seuratOutDir, paste(currSample,  "markers.rds", sep="-"))
  saveRDS(reduced.markers, markersFile)

  write.table(reduced.markers,
              file=file.path(seuratOutDir, paste(currSample,  "markers.tsv", sep="-")),
              sep="\t", quote=F, row.names=F, col.names = T)


  result.orig <- list(scds=pbmc, markers=reduced.markers)

  flog.info("finished processing sample: %s", currSample)
  result.orig
}

# For merging strategies
# - https://htmlpreview.github.io/?https://github.com/immunogenomics/harmony/blob/master/docs/SeuratV3.html
# - https://htmlpreview.github.io/?https://github.com/satijalab/seurat.wrappers/blob/master/docs/harmony.html

seurat_v2_merge_fn <- function(cdsList, mergedSampleName) {
  result <- merge(x=cdsList[[1]], y=cdsList[2:length(cdsList)], add.cell.ids=names(cdsList))
  result <- result %>%
    Seurat::NormalizeData(verbose = T) %>%
    do_seurat_dispersion(pipeline_options = default_pipeline_options, seuratOutDir = OUT_DIR, currSample = mergedSampleName) %>%
    ScaleData(verbose = T) %>%
    do_seurat_pca(genes=VariableFeatures(result), pipeline_options = default_pipeline_options, seuratOutDir = OUT_DIR, currSample=mergedSampleName) %>%
    RunUMAP(reduction="pca", dims=1:30) %>%
    FindNeighbors(dims=1:30) %>%
    FindClusters(resolution=opt$RESOLUTION) %>%
    get_cds_and_markers(currSample = mergedSampleName)
}


#' Integration (Seurat v4) with cca reduction - rpca crashes out
#'
#' https://satijalab.org/seurat/articles/integration_rpca.html
#'
seurat_v3_merge_fn <- function(cdsList, mergedSampleName) {
  cdsListNames <- names(cdsList)
  cdsList <- mapply(function(x, n) {
    flog.info("performing normalization and variable feature selection for sample: %s", n)
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method="vst", nfeatures=2000)
  }, cdsList, cdsListNames)

  flog.info("finding integration features for %d", length(cdsList))
  features <- SelectIntegrationFeatures(object.list = cdsList)

  # cdsList <- lapply(X = cdsList, FUN = function(x) {
  #   x <- ScaleData(x, features = features, verbose = T)
  #   x <- RunPCA(x, features = features, verbose = T)
  # })

  flog.info("found integration features: %d", length(features))
  anchors <- FindIntegrationAnchors(object.list = cdsList, anchor.features = features, reduction="rpca")

  result <- IntegrateData(anchorset = anchors)
  flog.info("Finished integration")
  DefaultAssay(result, "integrated")

  result <- result %>%
    ScaleData(verbose = T) %>%
    RunPCA(npcs = 20, verbose = T) %>%
    RunUMAP(reduction = "pca", dims = 1:20)

  p1 <- DimPlot(result, reduction="umap", group.by = "sample")
  ggsave(p1, filename =file.path(dirname(opt$OUT_RDS_FILE), paste0(mergedSampleName,"_seurat-v3-dimplot.pdf")))

  result
}

harmony_merge_fn <- function(cdsList, mergedSampleName, seed=2021) {

  # set a seed (default: 2021) so that plots are reproducible
  set.seed(seed)

  # start by merging using SeuratV2 and then do the fancy stuff
  result <- merge(x=cdsList[[1]], y=cdsList[2:length(cdsList)], add.cell.ids=names(cdsList)) %>%
    Seurat::NormalizeData(verbose = T) %>%
    do_seurat_dispersion(pipeline_options = default_pipeline_options, seuratOutDir = OUT_DIR, currSample = mergedSampleName) %>%
    ScaleData(verbose = T) %>%
    do_seurat_pca(genes=VariableFeatures(result), pipeline_options = default_pipeline_options, seuratOutDir = OUT_DIR, currSample=mergedSampleName) %>%
    RunHarmony("sample") %>%
    RunUMAP(reduction="harmony", dims=1:30) %>%
    FindNeighbors(reduction="harmony", dims=1:30) %>%
    FindClusters(resolution=opt$RESOLUTION) %>%
    get_cds_and_markers(currSample = mergedSampleName)
}


##########################################################################
## Run basic pipeline
##########################################################################
## Merge samples
seuratCdsList <-lapply(inputCdsList, convert_monocle_cds_to_seurat_cds)
names(seuratCdsList) <- names(inputCdsList)


seuratCdsList <- updateProjectNamesForSeuratCdsList(seuratCdsList)

flog.info("seuratCdsList: %s", toString(names(seuratCdsList)))

mergeSamples <- mergeJson[["sample"]]

flog.info("mergeSamples: %s", toString(mergeSamples))

mergeOptions <- mergeJson[["strategy"]]

listOfCdsListsToBeMerged <- lapply(names(mergeSamples), selectValidSamples)
# saveRDS(listOfCdsListsToBeMerged, file="/tmp/testing.RDS")

MERGE_ALGORITHM <- mergeJson[["strategy"]]
flog.info("MERGE_ALGORITHM IS: %s", MERGE_ALGORITHM)

if(MERGE_ALGORITHM == "SeuratV2") {
  flog.info("merging using Seurat V2")
  outputCdsList <- mapply( seurat_v2_merge_fn,listOfCdsListsToBeMerged, names(mergeSamples))
} else if(MERGE_ALGORITHM == "SeuratV3") {
  flog.info("merging using SeuratV3 with rpca reduction")
  outputCdsList <- mapply(seurat_v3_merge_fn,listOfCdsListsToBeMerged, names(mergeSamples))

} else if(MERGE_ALGORITHM == "Harmony") {
  flog.info("merging using Harmony")
  outputCdsList <- mapply(harmony_merge_fn,listOfCdsListsToBeMerged, names(mergeSamples))
}

# names(outputCdsList) <- names(mergeSamples)
flog.info("saving merged_rds to %s", opt$OUT_RDS_FILE)
saveRDS(outputCdsList, file=opt$OUT_RDS_FILE)

elapsed_time <- Sys.time() - start_time
flog.info("Finished successfully in %f seconds", elapsed_time)