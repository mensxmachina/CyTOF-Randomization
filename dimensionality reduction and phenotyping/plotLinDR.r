## Functions for plotting the results from the linear dimensionality reduction analyses

plot.pca <- function(Y,title = NULL) {
  df <- data.frame(PC1 = Y[,1],PC2 = Y[,2])
  
  r <- rev(RColorBrewer::brewer.pal(11, "Spectral"))
  g <- ggplot(df, aes(PC1, PC2)) +  stat_binhex(binwidth = 0.2) + ggtitle(title) + 
    scale_fill_gradientn(colours=r, limits=c(0,70)) + geom_density2d(colour = "black", size = 0.1) +
    ylim(c(-6,6)) + xlim(c(-6,6)) 
  
  return(g)
}

explained.variance <- function(sdev){
  return(100*sdev^2/sum(sdev^2))
}


plot.jive <- function(f,titles) {
  if(!require(r.jive)){
    install.packages("r.jive")
    require(r.jive)
  }
  
  l <- length(f$data)
  VarJoint = rep(0, l)
  for (i in 1:l) VarJoint[i] = norm(f$joint[[i]], type = "F")^2/norm(f$data[[i]], type = "F")^2
  
  VarIndiv = rep(0, l)
  for (i in 1:l) VarIndiv[i] = norm(f$individual[[i]], type = "F")^2/norm(f$data[[i]], type = "F")^2
  
  VarResid = 1 - VarJoint - VarIndiv
  
  df = data.frame( case = factor(as.vector(matrix(sapply(titles,rep,3),nrow = 1)), levels = titles),
                   variable = factor(rep(c("Joint","Individual","Residual"),l), levels = c("Joint","Individual","Residual")),
                   value = matrix(rbind(VarJoint, VarIndiv, VarResid),6,1) )
  
  # df$case = factor(df$case, levels = title)
  
  r <- rev(RColorBrewer::brewer.pal(11, "Spectral"))
  g <- ggplot(df,aes(fill = variable, y = value, x = case)) + ggtitle("Variation Explained") + 
    geom_bar(stat="identity", position = position_fill(reverse = TRUE)) + scale_fill_manual(values = r[1:3])+
    theme(axis.title.x = element_blank(), 
          axis.title.y = element_blank(),
          legend.title = element_blank(),
          axis.text = element_text(color="gray30", size=6),
          axis.text.x = element_text(angle=0,vjust=0.5),
          plot.title = element_text(size=8),
          text = element_text(size=8),
          legend.position="bottom",
          legend.justification="center",
          legend.key.size = unit(3,"mm"),
          legend.key.width = unit(3,'mm'),
          legend.text = element_text(margin = margin(r = 2, unit = "mm"),size = 6),
          legend.margin = margin(0,0,0,0)
    )
  return(g)
}


plot.MDS <- function(df){
  if(!require(ggrepel)){
    install.packages('ggrepel')
    require(ggrepel)
  }
  
  ggplot(df, aes(x = MDS1, y = MDS2, color = condition)) + 
    geom_point(size = 2, alpha = 0.8) + 
    ggrepel::geom_label_repel(aes(label = sample_id), point.padding = 0.2) + 
    theme_bw() +
    scale_color_manual(values = c("#1F78B4","#B15928"))
}