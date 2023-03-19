library(Seurat, quietly = T)
# source("enrichment_functions.R")

get_sample_name_from_seurat_cds <- function(pbmc) { str_split(pbmc@meta.data$id[1], "-")[[1]][3] }

make_topk_heatmap <-function(pbmc, pbmc.markers, outDir, k=10, saveHeatmap=T, cex.row=7, currSample=NULL) {
  if(is.null(currSample)) {
    currSample <- get_sample_name_from_seurat_cds(pbmc)
  } 
  flog.info("make_topk_heatmap() sample: %s k: %d outDir: %s", currSample, k, outDir)

  # setting slim.col.label to TRUE will print just the cluster IDS instead of
  # every cell name
  p  <- DoHeatmap(object = pbmc, features = pbmc.markers$gene)
  p <- p + theme(axis.text.y = element_text(size=cex.row), axis.text.x= element_text(size=cex.row))
  topk_heatmapFile <- file.path(outDir, paste0(currSample, "-heatmap", ".eps"))
  if(saveHeatmap) {
    p
    ggsave(p, filename = topk_heatmapFile, width=7, height = 11, units="in")
  } 
}

#'
#' Function to quickly filter a Seurat CellDataset for genes appearing in less than `min.cells`
#' 
my_filter <- function(object, min.cells) {
  genes.use <- rownames(object@data)
  if (min.cells > 0) {
    num.cells <- Matrix::rowSums(object@data > 0)
    genes.use <- names(num.cells[which(num.cells >= min.cells)])
    object@raw.data <- object@raw.data[genes.use, ]
    object@data <- object@data[genes.use, ]
  }
  object 
} 

get_seurat_pipeline_options <- function(curr_options=NULL) { 
  pipeline_options=list()

  pipeline_options[["NormalizeData"]] = list("normalization.method"="LogNormalize", "scale.factor"=1000)
  
  pipeline_options[["dispersionFile.height"]]=7
  pipeline_options[["dispersionFile.width"]]=7

  pipeline_options[["FindVariableGenes"]]=list("x.low.cutoff"=0.0125, "x.high.cutoff"=3, "y.cutoff"=0.5)
  
  pipeline_options[["RunPCA"]]=list("pcs.print"=1:9, "genes.print"=5)
  
  pipeline_options[["VizPCA"]]=list("height"=11, "width"=7, "pcs.use"=1:9, "filename"="vizPCA.pdf")
  
  pipeline_options[["PCHeatmap"]]=list("height"=11, "width"=7, "pc.use"=1:10, "do.balanced" = T, "label.columns" = F, "use.full"   = F) 

  pipeline_options[["FindClusters"]]=list("dims.use"=1:20, "reduction.type"="pca")
  
  pipeline_options[["FindAllMarkers"]]=list("only.pos"=F, "min.pct"=0.25, "thresh.use"=0.25)
  
  pipeline_options[["MergeSeurat"]]= list("do.normalize"=F, "do.scale"=F, "do.center"=F, "min.genes"=10, "min.cells"=10)
  
  if(!is.null(curr_options)) {
    for(key in names(curr_options)) {
      flog.info("overriding %s: %s (original: %s)", key, 
                toString(curr_options[[key]]), 
                toString(pipeline_options[[key]]))
      pipeline_options[[key]] = curr_options[[key]]
    }
  }
  
  for(key in names(pipeline_options)) {
    flog.info("pipeline_option %s -> %s", toString(key), toString(pipeline_options[[key]]))
  }
  
  pipeline_options
} 

do_seurat_dispersion <- function(pbmc, pipeline_options, seuratOutDir, currSample) {
  ## Find the variable genes and then? scale the data!  
  
  flog.info("do_seurat_dispersion() ")
  dispersionFile <- file.path(seuratOutDir, paste0(currSample, "-dispersion.pdf"))

  # pbmc <- FindVariableFeatures(object = pbmc,
  #                           mean.function = ExpMean,
  #                           dispersion.function = LogVMR,
  #                           x.low.cutoff = pipeline_options[["FindVariableGenes"]][["x.low.cutoff"]],
  #                           x.high.cutoff = pipeline_options[["FindVariableGenes"]][["x.high.cutoff"]],
  #                           y.cutoff = pipeline_options[["FindVariableGenes"]][["y.cutoff"]],
  #                           display.progress = T,
  #                           sort.results = T)
  pbmc <- pbmc %>%
    FindVariableFeatures(selection.method = "vst", nfeatures=2000)

  flog.info("number of variable genes: %d", length(VariableFeatures(pbmc)))
  
  
  
  # Identify the 10 most highly variable genes
  top10 <- head(VariableFeatures(pbmc), 20)
  
  # plot variable features with and without labels
  plot1 <- VariableFeaturePlot(pbmc)
  plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

  ggsave(dispersionFile, plot=plot2,
      height=pipeline_options[["dispersionFile.height"]], 
      width=pipeline_options[["dispersionFile.width"]])
  
  pbmc 
} 

do_seurat_pca<- function(pbmc, genes, pipeline_options, seuratOutDir, currSample) { 
  
  vizFile <- file.path(seuratOutDir, paste(currSample,   "vizPCA.pdf", sep="-"))
  pcaPlotFile <- file.path(seuratOutDir,paste(currSample,  "pca.pdf", sep="-" ))
  elbowFile <- file.path(seuratOutDir, paste(currSample, "pcElbowPlot.pdf", sep="-"))

  pbmc <- RunPCA(object = pbmc, 
                 do.print = F, 
                 pc.genes = genes,
                 pcs.print = pipeline_options[["RunPCA"]][["pcs.print"]], 
                 genes.print = pipeline_options[["RunPCA"]][["genes.print"]])
  

  p <- VizDimLoadings(object = pbmc, reduction = "pca")
  ggsave(vizFile, plot=p,
      width=pipeline_options[["VizPCA"]][["width"]], 
      height=pipeline_options[["VizPCA"]][["height"]])
  
  if(length(unique(pbmc@meta.data$sample)) < 6) { 
  p <- DimPlot(object=pbmc, group.by = 'sample', shape.by = 'sample') 
  } else {
    p <- DimPlot(object=pbmc, group.by = 'sample')
  }
  ggsave(pcaPlotFile, plot=p)
  
  pbmc
}


do_seurat_clustering_and_tsne <- function(pbmc, resolution, pipeline_options, seuratOutDir, currSample, genes.use) {
  
  flog.info("do_seurat_clustering_and_tsne() resolution: %f outDir: %s sample: %s", 
            resolution, seuratOutDir, currSample)
  
  flog.info("%s", toString( pipeline_options[["FindClusters"]]))
  if(is.null(genes.use)) {
    pbmc <- FindNeighbors(object = pbmc, 
                          # reduction.type = pipeline_options[["FindClusters"]][["reduction.type"]], 
                          # dims.use = pipeline_options[["FindClusters"]][["dims.use"]]
                          )
    pbmc <- FindClusters(object = pbmc, 
                         resolution = resolution, 
                         # print.output = 0, 
                         # save.SNN = T
                         )
    pbmc <- RunTSNE(object = pbmc, 
                    # dims.use = pipeline_options[["FindClusters"]][["dims.use"]], 
                    do.fast = TRUE)
    fileName1 <- file.path(seuratOutDir, paste(currSample, "tsne.pdf", sep="-"))
  } else {
    flog.info("do_seurat_clustering_and_tsne() - using genes instead of PC")
    pbmc <- FindClusters(object = pbmc, 
                         force.recalc = T,
                         resolution = resolution, print.output = 0, save.SNN = F,
                         genes.use=genes.use)
    pbmc <- RunTSNE(object = pbmc, 
                    do.fast = TRUE, genes.use=genes.use)
    fileName1 <- file.path(seuratOutDir, paste(currSample, "_top_", length(genes.use), "tsne.pdf", sep="-"))
  }
  
  flog.info("saving tsne_plot: %s", fileName1)
  p <- DimPlot(object = pbmc, reduction="tsne")
  ggsave(fileName1, plot=p)
  
  flog.info("do_seurat_clustering_and_tsne() -- end")
  pbmc
} 

run_seurat_pipeline <- function(pbmc, 
                                seuratOutDir, 
                                currSample=NULL, 
                                topK=NULL, 
                                resolution=0.6, 
                                pipeline_options=get_seurat_pipeline_options(NULL),
                                combinedCellTypeDt=NULL) {
  if(is.null(currSample)) { 
    currSample <- str_split(pbmc@meta.data$id[1], '-')[[1]][3]
  }
  
  dir.create(file.path(seuratOutDir), recursive = T, showWarnings = F)
  
  
  flog.info("processing: %s", currSample)
  
  pbmc <- NormalizeData(object = pbmc, 
                        normalization.method = pipeline_options[["NormalizeData"]][["normalization.method"]],
                        scale.factor = pipeline_options[["NormalizeData"]][["scale.factor"]])
 
  
  pbmc <- do_seurat_dispersion(pbmc, pipeline_options, seuratOutDir, currSample)
  
  pbmc <- ScaleData(object = pbmc, vars.to.regress = c("nUMI"))
  
  pc.genes <- pbmc@var.genes
  flog.info("pc.genes: %d", length(pc.genes))
  
  pbmc <- do_seurat_pca(pbmc, pc.genes, pipeline_options, seuratOutDir, currSample)

  pbmc <- do_seurat_clustering_and_tsne(pbmc, resolution, pipeline_options, seuratOutDir, currSample, genes.use = NULL)
  
  pbmc.markers <- FindAllMarkers(object = pbmc, 
                                 only.pos = pipeline_options[["FindAllMarkers"]][["only.pos"]], 
                                 min.pct = pipeline_options[["FindAllMarkers"]][["min.pct"]], 
                                 thresh.use = pipeline_options[["FindAllMarkers"]][["thresh.use"]])
  
  reduced.markers <- pbmc.markers %>% group_by(cluster)
  genes.use=NULL
  
  markersFile <- file.path(seuratOutDir, paste(currSample,  "markers.rds", sep="-"))
  saveRDS(reduced.markers, markersFile)
  result.orig <- list(scds=pbmc, markers=reduced.markers)
  result.orig <- do_seurat_enrichment(result.orig$scds, result.orig$markers, 
                                        combinedCellTypeDt, outDir=seuratOutDir, cSample=currSample)

  saveRDS(result.orig, file=file.path(seuratOutDir, paste(currSample, "processed.rds", sep="-")))
  
  
  if(!is.null(topK)) {
    for(k in topK) {
      localOutDir <- file.path(seuratOutDir, paste0("top_", k))
      dir.create(localOutDir, showWarnings = F, recursive = T)
      
      flog.info("processing reduced gene list for k: %d", k)
      reduced.markers <- pbmc.markers %>% group_by(cluster) %>% top_n(k, p_val_adj)
      genes.use <- reduced.markers$genes
      z <- toString(k)
      sampleName1 <- paste0(currSample, "_top_", z)
      rdsFile1 <- file.path(localOutDir, paste(sampleName1,  "processed.rds", sep="-"))
      
      pbmc <- do_seurat_pca(pbmc, genes.use, pipeline_options, localOutDir, sampleName1)
      pbmc <- do_seurat_clustering_and_tsne(pbmc, resolution, pipeline_options, localOutDir, sampleName1, genes.use = genes.use)

      result <- list(scds=pbmc, markers=reduced.markers)
      result <- do_seurat_enrichment(result$scds, result$markers, combinedCellTypeDt, outDir=localOutDir, cSample=sampleName1)
      
      make_topk_heatmap(pbmc, reduced.markers, localOutDir, saveHeatmap=T, cex.row=3, currSample=sampleName1) 
      saveRDS(result, file=rdsFile1)
      
    }
  }
  flog.info("finished processing sample: %s", currSample)
  result.orig
} 


make_violin_plots <- function(cElem, outDir, nCol=5, min.cells=10, cSample=NULL, markerPValueAdjThreshold=NULL) { 
  pbmc <- cElem$scds
  pbmc.markers <- cElem$markers
  pbmc <- my_filter(pbmc, min.cells)
  if(is.null(cSample)) { 
    # for backward compatibility
    cSample <- str_split(pbmc@meta.data$id[1], "-")[[1]][3]
  } 
  plotList <- list()
  cGenes <- Reduce(union, Markers)
  cIdx <- which(cGenes %in% rownames(pbmc@data))
  cGenes.keep <- cGenes[cIdx]
  flog.info("sample: %s genes: %s", cSample, paste(cGenes.keep, collapse = ", "))
  nClusters <- length(unique(pbmc.markers$cluster))
  cWidth  <- ceiling(nClusters/10)*2*nCol
  cHeight <- ceiling(length(cGenes.keep)/nCol)*1.5
  cXsize  <- 0.5
  p <- VlnPlot(object = pbmc, features.plot = cGenes.keep, y.log=T, nCol = nCol, size.x.use=cXsize, size.title.use=10, point.size.use=0.25, do.return = T)+ggtitle(cSample)
  ggsave(filename = paste0(outDir, "/", cSample, "-vln.pdf"), width=cWidth, height=cHeight, units="in", dpi=600)
  p
} 


getGeneSpecificCdsList <- function(processedCdsList, geneName="TH") { 
  processedCdsList <- lapply(processedCdsList, function(x) { 
    x$scds <- my_filter(x$scds, min.cells=10) 
    x
  })
  cdsList.geneSpecific <- lapply(names(processedCdsList), function(sampleName) {
    sCds <- processedCdsList[[sampleName]]$scds
    sCds.geneSpecific <- FilterCells(sCds, geneName, 0)
  })
  names(cdsList.geneSpecific) <- names(processedCdsList)
  cdsList.geneSpecific
}

updateProjectNamesForSeuratCdsList <- function(seuratCdsList) {
  cNames <- names(seuratCdsList)
  seuratCdsList <- lapply(cNames, function (x) {
  seuratCdsList[[x]]@project.name = x
  seuratCdsList[[x]]
  })
  names(seuratCdsList) <- cNames
  seuratCdsList
}


#' Function to merge multiple seurat objects
MergeSeuratMultiple <- function(cdsList, merge_seurat_options, project.name=NULL) {
  flog.info("MergeSeuratMultiple: (%s)", toString(names(merge_seurat_options))) 
  flog.info("MergeSeuratMultiple: (%s)", toString(merge_seurat_options)) 
  
  cdsIds = names(cdsList)
  
  for(i in 1:length(cdsList)) {
    cCds = cdsList[[i]]
    cId = cdsIds[[i]]
    flog.info("%d/%d merging %s", i, length(cdsList), cId)
    if(i==1) {
      result = cCds
    } else if(i==2) {
      result = merge(
        result, cCds, 
        add.cell.id1 = cdsIds[[1]],
        add.cell.id2 = cId, 
        do.normalize = merge_seurat_options[["do.normalize"]],
        do.scale = merge_seurat_options[["do.scale"]],
        do.center= merge_seurat_options[["do.center"]],
        min.genes= merge_seurat_options[["min.genes"]],
        min.cells= merge_seurat_options[["min.cells"]]
      )
    } else {
      result <- merge(
        result, cCds, 
        add.cell.id2 = cId,
        do.normalize = merge_seurat_options[["do.normalize"]],
        do.scale = merge_seurat_options[["do.scale"]],
        do.center= merge_seurat_options[["do.center"]],
        min.genes= merge_seurat_options[["min.genes"]],
        min.cells= merge_seurat_options[["min.cells"]]
      )
    }
  }
  if(!is.null(project.name)) {
    flog.info("setting merged project.name: %s", project.name)
    result@project.name = project.name
  }
  result
}



