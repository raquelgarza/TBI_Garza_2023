source("rcode/pipeline_helper_functions.R")
PKGS <- c("monocle", "Seurat", "cowplot", "data.table", "reshape2", "stringr", "dplyr", "clusterProfiler", "jsonlite", "getopt")
silently_load_package_list(PKGS)

spec <- matrix(c(
  'IN_SAMPLEDF_FILE' , 'idf', 1, "character", "input sampleDf file",
  'IN_RDS_FILE', 'irds', 1, "character", "input RDS file containing list of monocleDatasets", 
  'OUT_RDS_FILE', 'ords', 1, "character", "Final seurat processed output", 
  'OUT_SAMPLEDF_FILE', 'ordf', 2, "character", "Copy out sampleDf file to seuratOutDir", 
  'OUT_DIR', 'odir', 1, "character", "output directory for saving stuff",
  'OPTIONS_JSON_FILE', 'j', 2, "character", "json file containing seurat options", 
  'RESOLUTION', 'r', 1, "double", "resolution for clustering",
  'PANGLODB_DIR', 'p', 1, "character", "panglodb_directory", 
  'MARKERS_FILE', 'mfile', 1, 'character', 'markers R file', 
  'LOG_FILE', 'l', 2, "character", "where output should be logged",
  'help' , 'h', 0, "logical", "show this help message"
), byrow=TRUE, ncol=5)

opt <- getopt(spec)

DEFAULT_OPTIONS <- list()
DEFAULT_OPTIONS[["TOP_K"]] <- c(10, 20, 50, 100)
DEFAULT_OPTIONS[["OUT_DIR"]] <- "/tmp"
DEFAULT_OPTIONS[["OUT_RDS_FILE"]] <- file.path(DEFAULT_OPTIONS$OUT_DIR, "results-processed.rds")
DEFAULT_OPTIONS[["LOG_FILE"]] <- file.path(DEFAULT_OPTIONS$OUT_DIR, "run_seurat.log")
DEFAULT_OPTIONS[['OPTIONS_JSON_FILE']] <- 'seurat_defaults.json'
DEFAULT_OPTIONS[["MARKERS_FILE"]] <- NULL 
DEFAULT_OPTIONS[["PANGLODB_DIR"]] <- NULL
DEFAULT_OPTIONS[["OUT_SAMPLEDF_FILE"]] <- file.path(opt$OUT_DIR, "sampleDf_processed.tsv")

opt <- set_default_options(opt, DEFAULT_OPTIONS, get_Rscript_filename())

source("rcode/cellranger_helper_functions.R")
source("rcode/plotting_functions.R")
source("rcode/seurat_functions.R")
source("rcode/enrichment_functions.R")


if(!dir.exists(opt$OUT_DIR)) { 
  dir.create(opt$OUT_DIR, recursive = TRUE, showWarnings = FALSE)
} 

sampleDf <- read_sampledf(opt$IN_SAMPLEDF_FILE)
seuratCdsList <- readRDS(opt$IN_RDS_FILE)
seuratOptions <- get_seurat_pipeline_options() # seuratOptions <- read_json(opt$OPTIONS_JSON_FILE)

if(!is.null(opt$MARKERS_FILE)) {
  flog.info("reading markers file: %s", opt$MARKERS_FILE)
  Markers <- readRDS(opt$MARKERS_FILE)
} else {
  flog.info("markers not provided - not doing feature plots")
  Markers <- NULL
}

if(!is.null(opt$PANGLODB_DIR)) { 
  flog.info("reading panglodb: %s", format(opt$PANGLODB_DIR))
  combinedCellTypeDt <- readPangloDbFilesDir(opt$PANGLODB_DIR)
} else {
  flog.info("panglodb files not provided - no enrichment")
  combinedCellTypeDt <- NULL
} 


flog.info("samples: %d file: %s", length(seuratCdsList), opt$IN_RDS_FILE)

###### BEGINNING OF FUNCTIONS 

quick_cell_count <- function(cCds, cGene) {
  cData <- GetAssayData(cCda, slot="data")
  cCountData <- GetAssayData(cData, slot="counts")
  sum(as.matrix(cCountData[which(rownames(cData)==cGene), ]) > 0)
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


plot_top_k_markers <- function(x, currSample, k=5, localOutDir="/tmp") {
  pbmc <- x$scds
  flog.info("processing reduced gene list for k: %d", k)
  reduced.markers <- x$markers %>% group_by(cluster) %>% top_n(k, p_val_adj)
  genes.use <- unique(reduced.markers$gene)
  z <- toString(k)
  sampleName1 <- paste0(currSample, "_top_", z)
  
  result <- list(scds=pbmc, markers=reduced.markers)
  
  ## make heatmap
  make_topk_heatmap(result$scds, reduced.markers, localOutDir, saveHeatmap=T, cex.row=3, currSample=sampleName1)
  
  ## make violin plots
  violinPlots <- list()
  for(cGene in genes.use) {
    p <- VlnPlot(pbmc, features=cGene)
    violinPlots[[cGene]] <- p 
    ggsave(file.path(localOutDir, paste0(currSample, "_violin_", cGene,".pdf")), plot=p)
  } 
  saveRDS(violinPlots, file=file.path(localOutDir, paste0(sampleName1, "_violinPlots.rds")))
  
  ## make feature plots
  featurePlots <- list()
  for(cGene in genes.use) {
    p <- FeaturePlot(pbmc, features=cGene, shape.by="sample")
    violinPlots[[cGene]] <- p 
    ggsave(file.path(localOutDir, paste0(currSample, "_feature_", cGene,".pdf")), plot=p)
  } 
  saveRDS(violinPlots, file=file.path(localOutDir, paste0(sampleName1, "_featurePlots.rds")))
  
  result  
}


run_seurat_pipeline_local <- function(pbmc, 
                                seuratOutDir, 
                                currSample=NULL, 
                                resolution=1, 
                                pipeline_options=get_seurat_pipeline_options(NULL),
                                combinedCellTypeDt=NULL) {
  
  if(!dir.exists(seuratOutDir)) { 
    dir.create(seuratOutDir, recursive = T, showWarnings = F)
  }
  
  
  flog.info("processing: %s", currSample)
  
  pbmc <- NormalizeData(object = pbmc, 
                        normalization.method = pipeline_options[["NormalizeData"]][["normalization.method"]],
                        scale.factor = pipeline_options[["NormalizeData"]][["scale.factor"]])
  
  
  pbmc <- do_seurat_dispersion(pbmc, pipeline_options, seuratOutDir, currSample)
  
  pbmc <- ScaleData(object = pbmc)
  
  pc.genes <- VariableFeatures(pbmc)
  flog.info("pc.genes: %d", length(pc.genes))
  
  pbmc <- do_seurat_pca(pbmc, pc.genes, pipeline_options, seuratOutDir, currSample)
  
  flog.info("running UMAP")
  pbmc <- RunUMAP(pbmc, dims = 1:30, verbose = FALSE)
  
  pbmc <- FindNeighbors(pbmc, dims = 1:30, verbose = T)
  pbmc <- FindClusters(pbmc, verbose = T, resolution = resolution)
  
  currentClusterName = paste0("umap_", resolution)
  pbmc@meta.data[[currentClusterName]] <- pbmc@meta.data[["seurat_clusters"]]
  
  my_dim_plot(pbmc, reduction="umap", seuratOutDir=seuratOutDir, currSample = currSample)
  
  flog.info("Finished running UMAP")
  
  # flog.info("run tsne and clustering")
  # pbmc <- do_seurat_clustering_and_tsne(pbmc, resolution, pipeline_options, seuratOutDir, currSample, genes.use = NULL)
 
  # currentClusterName = paste0("tsne_", resolution)
  # pbmc@meta.data[[currentClusterName]] <- pbmc@meta.data[["seurat_clusters"]]
  # my_dim_plot(pbmc, reduction="tsne", seuratOutDir=seuratOutDir, currSample = currSample)


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



##### END OF FUNCTIONS


processedCdsList <- lapply(names(seuratCdsList), function(n) {
  x <- seuratCdsList[[n]]
  cSample <- n
  currSampleOutDir <- file.path(opt$OUT_DIR, x@project.name)
  
  flog.info("running %s => %s", x@project.name, currSampleOutDir)
  
  result <- run_seurat_pipeline_local(x, 
                      seuratOutDir = currSampleOutDir, 
                      currSample=n, 
                      resolution = opt$RESOLUTION,
                      pipeline_options = seuratOptions, 
                      combinedCellTypeDt=NULL)

  flog.info("making top feature plots")
  # make_feature_plots_local(result, Markers, currSampleOutDir, cSample)
  
  k=5
  # plot_top_k_markers(result, cSample, 5, opt$OUT_DIR)
  
  result 
})

names(processedCdsList) <- sapply(seuratCdsList, function(x) { x@project.name })
write_json(seuratOptions, file.path(opt$OUT_DIR, "seurat_options.json")) 
saveRDS(processedCdsList, opt$OUT_RDS_FILE)
write_sampledf(sampleDf, opt$OUT_SAMPLEDF_FILE)
flog.info("========================= DONE ====================")

