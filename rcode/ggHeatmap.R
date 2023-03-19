require(ggplot2)
require(reshape2)

getClusteredMatrix <- function(m, clustering='none') {
  if(is.function(clustering))
  {
    m=clustering(m)
  }else
  {
    if(clustering=='row')
      m=m[hclust(dist(m))$order, ]
    if(clustering=='column')
      m=m[,hclust(dist(t(m)))$order]
    if(clustering=='both')
      m=m[hclust(dist(m))$order ,hclust(dist(t(m)))$order]
  }
  m
}


ggheatmap <- function(m, clustering='both', outputPlotFile='~/test.pdf', meltedColNames=c('GeneName', 'Cluster', 'LFC')) { 
  m <- getClusteredMatrix(m, clustering)
  melted.m <- melt(m)
  colnames(melted.m) <- meltedColNames
  gg <- ggplot(melted.m, aes_string(y=meltedColNames[[1]], x=meltedColNames[[2]], fill=meltedColNames[[3]]))
  gg <- gg + geom_tile( colour = "white")
  gg <- gg + scale_fill_gradient2(low = "blue", mid = "white", high = "red")
  gg <- gg + labs(x="", y="")
  gg <- gg + theme_bw()
  gg <- gg + theme(panel.grid=element_blank(), panel.border=element_blank())
  y.axis.text.style <- element_text(face = "bold", color = "black", size = 6)
  gg <- gg + theme(axis.text.y = y.axis.text.style)
  gg
  nx <- ncol(m)
  ny <- nrow(m)
  if(is.null(outputPlotFile)) { 
    ggsave(outputPlotFile, plot=gg, scale=2.0, width=1*nx, height=ceiling(ny*0.10), units="cm", dpi=600)
}
gg
} 
