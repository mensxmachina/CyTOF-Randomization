library(flowCore)
library(RColorBrewer)
library(tidyverse)
library(cowplot)

## Load the data
nonrand = read.FCS('non_randomized.fcs', transform=FALSE)
rand.type1 = read.FCS('randomized_type1.fcs', transform=FALSE)
rand.type2 = read.FCS('randomized_type2.fcs', transform=FALSE)
rand.extreme = read.FCS('randomized_maximal.fcs', transform=FALSE)
rand.extreme.repl1 = read.FCS('randomized_maximal_replicate1.fcs', transform=FALSE)
rand.extreme.repl2 = read.FCS('randomized_maximal_replicate2.fcs', transform=FALSE)

# asinh transforms
nonrand@exprs[,3:58] <- asinh(nonrand@exprs[,3:58]/5)
rand.type1@exprs[,3:58] <- asinh(rand.type1@exprs[,3:58]/5)
rand.type2@exprs[,3:58] <- asinh(rand.type2@exprs[,3:58]/5)

## PANEL ANNOTATIONS
allMarkers <- read.csv('markers.csv',stringsAsFactors = F)

## Exclude non Marker and failed Marker variables
nonMarkers.Cytk <- as.numeric(unlist(sapply(c("Barcode","Time","cell length","unkn","EQBeads","DNA"),grep, allMarkers$CytkPanel)))
failedMarkers <- as.numeric(unlist(sapply(c("TCRab", "\\bCD4\\b", "ICOS", "CCR6", "CD57", "CD28"),grep, allMarkers$CytkPanel)))
excludedMarkers <- c(failedMarkers,nonMarkers.Cytk)

## Identify the cytokine and lineage markers
cytokine.Markers <- as.numeric(unlist(sapply(c("IL","IFNg","Granzyme","TNFa","GM-CSF"),grep, allMarkers$CytkPanel)))
lineage.Markers <- setdiff(c(1:NROW(allMarkers)),c(excludedMarkers,cytokine.Markers))

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
#                                                         #
#  Perform PCA (we subsample without loss of generality)  #
#                                                         #
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

## Functions to use
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
  ggplot(df, aes(x = MDS1, y = MDS2, color = condition)) + 
    geom_point(size = 2, alpha = 0.8) + 
    ggrepel::geom_label_repel(aes(label = sample_id), point.padding = 0.2) + 
    theme_bw() +
    scale_color_manual(values = c("#1F78B4","#B15928"))
}

set.seed(123)
subsample <- sample(seq(1:dim(nonrand)[1]),dim(nonrand)[1]/10)

## PCA
pca_nonrand = prcomp(nonrand@exprs[subsample,-excludedMarkers], center = T, scale. = F) # non-randomized
pca_rand.type1 = prcomp(rand.type1@exprs[subsample,-excludedMarkers], center = T, scale. = F) # randomized type 1
pca_rand.type2 = prcomp(rand.type2@exprs[subsample,-excludedMarkers], center = T, scale. = F) # randomized type 2
pca_rand.extreme = prcomp(rand.extreme@exprs[subsample,-excludedMarkers], center = T, scale. = F) # extreme randomization


## how much variance each PC is explaining ?
expVar.type1 = data.frame(nonrand = cumsum(explained.variance(pca_nonrand$sdev)),
                          rand = cumsum(explained.variance(pca_rand.type1$sdev)) )

expVar.type2 = data.frame(nonrand = cumsum(explained.variance(pca_nonrand$sdev)),
                          rand = cumsum(explained.variance(pca_rand.type2$sdev)) )

expVar.extreme = data.frame(nonrand = cumsum(explained.variance(pca_nonrand$sdev)),
                            rand = cumsum(explained.variance(pca_rand.extreme$sdev)) )


## perform JIVE
J.type1 <- r.jive::jive(list(t(nonrand@exprs[subsample,-excludedMarkers]),
                             t(rand.type1@exprs[subsample,-excludedMarkers])),scale = F, center = T, maxiter = 100)

J.type2 <- r.jive::jive(list(t(nonrand@exprs[subsample,-excludedMarkers]),
                             t(rand.type2@exprs[subsample,-excludedMarkers])),scale = F, center = T, maxiter = 100)

J.extreme <- r.jive::jive(list(t(nonrand@exprs[subsample,-excludedMarkers]),
                               t(rand.extreme@exprs[subsample,-excludedMarkers])),scale = F, center = T, maxiter = 100)

## Perform MDS (median marker values)
data.MDS = data.frame( nonrand = apply(exprs(nonrand)[,-excludedMarkers],2,median),
                       type1 = apply(exprs(rand.type1)[,-excludedMarkers],2,median),
                       type2 = apply(exprs(rand.type2)[,-excludedMarkers],2,median),
                       maximal = apply(exprs(rand.extreme)[,-excludedMarkers],2,median))
rownames(data.MDS) = allMarkers$CytkPanel[-excludedMarkers]

mds = limma::plotMDS(data.MDS,plot = FALSE)

ggMDS <- data.frame(MDS1 = mds$x, MDS2 = mds$y, sample_id = colnames(data.MDS),
                    condition = factor(c("Non-Random.",rep("Randomized",3))))


## plot the 2D PCA embeddings
plot.pca(pca_nonrand$x[,1:2],'Non-randomized')
plot.pca(pca_rand.type1$x[,1:2],'Randomized, type 1')
plot.pca(pca_rand.type2$x[,1:2],'Randomized, type 2')
plot.pca(pca_rand.extreme$x[,1:2],'Randomized, extreme')

## plot the explained variance per PC
ggplot(expVar.type1, aes(nonrand, rand)) + geom_point(color = "tomato", alpha = 0.7) + 
  labs(x="non-randomized (%)",y="randomized, type 1 (%)") + xlim(c(0,100)) + ylim(c(0,100)) + geom_abline(intercept = 0, lty = 2, size = 0.25)

ggplot(expVar.type2, aes(nonrand, rand)) + geom_point(color = "tomato", alpha = 0.7) + 
  labs(x="non-randomized (%)",y="randomized, type 2 (%)") + xlim(c(0,100)) + ylim(c(0,100)) + geom_abline(intercept = 0, lty = 2, size = 0.25)

ggplot(expVar.extreme, aes(nonrand, rand)) + geom_point(color = "tomato", alpha = 0.7) + 
  labs(x="non-randomized (%)",y="randomized, extreme (%)") + xlim(c(0,100)) + ylim(c(0,100)) + geom_abline(intercept = 0, lty = 2, size = 0.25)

## plot the common, individual and residual variances
plot.jive(J.type1,titles = c(gsub('\\-','\n','non-randomized'), "randomized, type 1"))
plot.jive(J.type2,titles = c(gsub('\\-','\n','non-randomized'), "randomized, type 2"))
plot.jive(J.extreme,titles = c(gsub('\\-','\n','non-random.'), gsub(' ','\n','extreme')))


## plot the multidimensional scaling difference
plot.MDS(ggMDS)

## ## ## ## ## ## ## ## ## ## ## ##
#                                 #
#  Perform tSNE (no subsampling)  #
#                                 #
## ## ## ## ## ## ## ## ## ## ## ##

## Functions to use
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
# plot.Marker.on.tSNE <- function(df.tSNE,marker.name) {
#   r <- rev(RColorBrewer::brewer.pal(11, "Spectral"))
#   
#   df.tSNE[,3] = df.tSNE[,3]/quantile(df.tSNE[,3],.99)
#   g <- ggplot(data = df.tSNE, aes(x = tSNE1, y = tSNE2, z = df.tSNE[,3])) + 
#     stat_summary_hex(fun = mean, bins=200) + ggtitle(marker.name) +
#     labs(fill=gsub('\\s','\n','normalized abundance')) +
#     scale_fill_gradientn(colours = r, limits = c(0,1),oob=scales::squish)
#   
#   return(g)
# }


## Perform tSNE
set.seed(123)
tsne_out_nonrand <- Rtsne.multicore::Rtsne.multicore(nonrand@exprs[,-excludedMarkers], dims = 2)

set.seed(123)
tsne_out_rand_type1 <- Rtsne.multicore::Rtsne.multicore(rand.type1@exprs[,-excludedMarkers], dims = 2)

set.seed(123)
tsne_out_rand_type2 <- Rtsne.multicore::Rtsne.multicore(rand.type2@exprs[,-excludedMarkers], dims = 2)

set.seed(123)
tsne_out_rand_extreme <- Rtsne.multicore::Rtsne.multicore(rand.extreme@exprs[,-excludedMarkers], dims = 2)

set.seed(123)
tsne_out_rand_extreme_repl1 <- Rtsne.multicore::Rtsne.multicore(rand.extreme.repl1@exprs[,-excludedMarkers], dims = 2)

set.seed(123)
tsne_out_rand_extreme_repl2 <- Rtsne.multicore::Rtsne.multicore(rand.extreme.repl2@exprs[,-excludedMarkers], dims = 2)

## Plot the tSNE embeddings
plot.tSNE(tsne_out_nonrand$Y,'Non-randomized')
plot.tSNE(tsne_out_rand_type1$Y,'Randomized, type 1')
plot.tSNE(tsne_out_rand_type2$Y,'Randomized, type 2')
plot.tSNE(tsne_out_rand_extreme$Y,'Randomized, extreme')
plot.tSNE(tsne_out_rand_extreme_repl1$Y,'Randomized, extreme, replicate')
plot.tSNE(tsne_out_rand_extreme_repl2$Y,'Randomized, extreme, replicate')

## Plot the distribution of each marker on the tSNE embedding
tSNE.nonrand <- data.frame(tSNE1 = tsne_out_nonrand$Y[,1], 
                           tSNE2 = tsne_out_nonrand$Y[,2],
                           as.matrix(nonrand@exprs[,-excludedMarkers]))

tSNE.rand.type1 <- data.frame(tSNE1 = tsne_out_rand_type1$Y[,1], 
                              tSNE2 = tsne_out_rand_type1$Y[,2],
                              as.matrix(rand.type1@exprs[,-excludedMarkers]))


tSNE.rand.type2 <- data.frame(tSNE1 = tsne_out_rand_type2$Y[,1], 
                              tSNE2 = tsne_out_rand_type2$Y[,2],
                              as.matrix(rand.type2@exprs[,-excludedMarkers]))

tSNE.rand.extreme <- data.frame(tSNE1 = tsne_out_rand_extreme$Y[,1], 
                                tSNE2 = tsne_out_rand_extreme$Y[,2],
                                as.matrix(rand.extreme@exprs[,-excludedMarkers]))

tSNE.rand.extreme.repl1 <- data.frame(tSNE1 = tsne_out_rand_extreme_repl1$Y[,1], 
                                     tSNE2 = tsne_out_rand_extreme_repl1$Y[,2],
                                     as.matrix(rand.extreme.repl1@exprs[,-excludedMarkers]))

tSNE.rand.extreme.repl2 <- data.frame(tSNE1 = tsne_out_rand_extreme_repl2$Y[,1], 
                                      tSNE2 = tsne_out_rand_extreme_repl2$Y[,2],
                                      as.matrix(rand.extreme.repl2@exprs[,-excludedMarkers]))


M <- allMarkers[-excludedMarkers,2]
colnames(tSNE.rand.extreme.repl1)[3:ncol(tSNE.rand.extreme.repl1)] <- colnames(tSNE.rand.extreme)[3:ncol(tSNE.rand.extreme)] <- colnames(tSNE.rand.type2)[3:ncol(tSNE.rand.type2)] <- colnames(tSNE.rand.type1)[3:ncol(tSNE.rand.type1)] <- colnames(tSNE.nonrand)[3:ncol(tSNE.nonrand)] <- M

nrand_label <- ggdraw() + draw_label("non-rand.", size = 6)
rand_label_type1 <- ggdraw() + draw_label("type 1", size = 6)
rand_label_type2 <- ggdraw() + draw_label("type 2", size = 6)
rand_label_extreme <- ggdraw() + draw_label("maximal", size = 6)

p.all <- list()
for (m in 1:length(M)){
  
  g1 <- plot.Marker.on.tSNE(tSNE.nonrand[,c(1,2,(2+m))], marker.name = "")
  g2 <- plot.Marker.on.tSNE(tSNE.rand.type1[,c(1,2,(2+m))], marker.name = "")
  g3 <- plot.Marker.on.tSNE(tSNE.rand.type2[,c(1,2,(2+m))], marker.name = "")  
  g4 <- plot.Marker.on.tSNE(tSNE.rand.extreme[,c(1,2,(2+m))], marker.name = "")
  
  marker_label <- ggdraw() + draw_label(M[m], angle = 90, size = 6)
  
  p1 <- plot_grid(nrand_label,rand_label_type1, rand_label_type2,rand_label_extreme,
                  g1 + theme(legend.position="none",plot.margin=margin(t=-0.5,b=0.1, unit="cm")), 
                  g2 + theme(legend.position="none",plot.margin=margin(t=-0.5,b=0.1,unit="cm")),
                  g3 + theme(legend.position="none",plot.margin=margin(t=-0.5,b=0.1,unit="cm")),
                  g4 + theme(legend.position="none",plot.margin=margin(t=-0.5,b=0.1,unit="cm")),
                  nrow = 2, rel_heights = c(.1,1))
  
  p3 <- plot_grid(marker_label, p1, rel_widths = c(.1,1))
  
  p.all[[m]] <- p3 + annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1, color = "gray95", fill = NA)
}
p.all[[m+1]] <- get_legend(g1)

p <- plot_grid(plotlist = p.all, nrow = 8, ncol = 4)
p <- add_sub(p,expression(paste(bold("Supplementary Figure 4")," tSNE analysis results for all markers. Each subplot displays the marker distribution when the")), 
             x = 0, y = 0, hjust = 0, vjust = 0, size = 8)
p <- add_sub(p,"non-randomized data (left) and the randomized data are employed",
             x = 0, y = 1, hjust = 0, vjust = 1.3, size = 8)

ggsave(paste0("supplementary_figure4_full.tif"),plot = p,device = "tiff",dpi = 300,width = 20,height = 18, units ="cm")
