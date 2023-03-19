# This file has plotting related functions
source("rcode/pipeline_helper_functions.R")

silently_load_package_list(c("ellipse", "ggplot2", "gplots", "cowplot"))


#'
#' Nice colored correlation plot
#' 
plot_colored_correlation <- function(ctab) {
  colorfun <- colorRamp(c("#CC0000","white","#3366CC"), space="Lab")
  plotcorr(ctab, col=rgb(colorfun((ctab+1)/2), maxColorValue=255), mar = c(0.1, 0.1, 0.1, 0.1))
}

#'
#' Function that makes single A4 merged plot of all the `*.eps.rda` plots
#' present in `outDir/subDir` and saves resulting plot in `outDir/merged-subDir.eps`
#' 
make_merged_plot <- function(outDir, subDir, ncol=3) { 
  plotDir <- file.path(outDir, subDir)
  outPlot <- paste0(outDir, "/merged-", subDir, ".eps")
  theme_set(theme_cowplot(font_size = 5))
  plotFiles <- Sys.glob(file.path(plotDir, "/*.eps.rds"))
  pl <- list()
  for(currPlotFile in plotFiles) {
    currId <- gsub(".eps.rds", "", basename(currPlotFile))
    flog.info("currId: %s", currId)
    p <- readRDS(currPlotFile)
    pl[[currId]] <- p
  }
  nr = ceiling( length(plotFiles) / ncol)  
  plot_grid(plotlist = pl, ncol=ncol, labels="AUTO")
  ggsave(file=outPlot, width=5*ncol, height=5*nr, units="in", dpi=600)
} 

# 
# get_human_vs_rat_plot(zz)
#
get_human_vs_rat_plot <- function(z2, ratio10x=0.0) {
  idList <- str_split(colnames(z2)[1], '-')[[1]]
  colnames(z2) <- sapply(colnames(z2), function(x) str_split(x, '-')[[1]][1])
  sampleName <- paste(idList[2], idList[3], sep='-')
  z2$ratio <- z2$GRCh38 / (z2$GRCh38 + z2$Rnor6) 
  flog.info('%s', sampleName)
  p <- ggplot(z2, aes(x=ratio)) + geom_histogram() 
  p <- p + geom_vline(xintercept = 5/6, linetype="dashed", color="red")
  p <- p + labs(x='', y='#cells', title=sampleName) 
  p 
}

#' 
#'  Function to make histogram per column
#'  
#'  @return ggplot histogram
#'   
make_histogram_plot_per_column <- function(x) {
  x %>% gather() %>% head()
  ggplot(gather(x), aes(value)) + 
    geom_histogram(bins = 10) + 
    facet_wrap(~key, scales = 'free_x')
}