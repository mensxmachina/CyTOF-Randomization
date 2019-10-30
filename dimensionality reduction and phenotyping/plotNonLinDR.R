## Functions for plotting the results from the nonlinear dimensionality reduction (tSNE) analysis

plot.tSNE <- function(Y,title = NULL) {
  df <- data.frame(tSNE1 = Y[,1],tSNE2 = Y[,2])
  
  glim <- 30#ceiling(max(abs(df)))+1
  r <- rev(RColorBrewer::brewer.pal(11, "Spectral"))
  g <- ggplot(df, aes(tSNE1, tSNE2)) +  stat_binhex(binwidth=0.3) +
    scale_fill_gradientn(colours=r, limits=c(0,30)) + geom_density2d(colour = "black", size = 0.1) +
    ylim(c(-glim,glim)) + xlim(c(-glim,glim)) + ggtitle(title) 
  
  return(g)
}
plot.Marker.on.tSNE <- function(df.tSNE,marker.name) {
  r <- rev(RColorBrewer::brewer.pal(11, "Spectral"))
  
  df.tSNE[,3] = df.tSNE[,3]/quantile(df.tSNE[,3],.99)
  g <- ggplot(data = df.tSNE, aes(x = tSNE1, y = tSNE2, z = df.tSNE[,3])) + 
    stat_summary_hex(fun = mean, bins=200) + ggtitle(marker.name) +
    labs(fill=gsub('\\s','\n','normalized abundance')) +
    scale_fill_gradientn(colours = r, limits = c(0,1),oob=scales::squish) +
    theme( 
      axis.text = element_blank(),
      axis.ticks=element_blank(),
      axis.title=element_blank(),
      axis.line = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      plot.title = element_text(size=8),
      legend.title = element_text(size = 6),
      legend.text = element_text(size = 6),
      legend.key.width = unit(3,'mm'),
      legend.key.height = unit(1,'mm'),
      legend.margin = margin(0,0,0,0), 
      legend.title.align=0.5,
      legend.direction = "horizontal",
      legend.justification="center" ,
      legend.box.just = "bottom",
      legend.position = "bottom") + 
    guides(fill = guide_legend(label.position = "top"))
  
  return(g)
}