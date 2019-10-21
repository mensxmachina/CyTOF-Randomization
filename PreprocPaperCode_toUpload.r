##
## Preprocessing paper code 
##
## editor: Giorgos Papoutsoglou
##
## date: 07/05/2019

library(flowCore)
library(tidyverse)
library(cowplot)
library(FlowSOM)
library(RColorBrewer)
# packages also required 'r.jive', 'limma', 'ggrepel','fossil','clue'


## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
#                                                #
#  Read the non-randomized and randomized files  #
#                                                #
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

nonrand = read.FCS('nonrandomized.fcs', transform=FALSE)
rand.type1 = read.FCS('randomized_type1.fcs', transform=FALSE)
rand.type2 = read.FCS('randomized_type2.fcs', transform=FALSE)
rand.extreme = read.FCS('randomized_maximal.fcs', transform=FALSE)

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

J <- r.jive::jive(list(t(nonrand@exprs[subsample,-excludedMarkers]),
                       t(rand.type1@exprs[subsample,-excludedMarkers]),
                       t(rand.type2@exprs[subsample,-excludedMarkers]),
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
plot.jive(J, titles = c(gsub('\\-','\n','non-randomized'), "type 1", "type 2", " extreme"))

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


## Perform tSNE
set.seed(123)
tsne_out_nonrand <- Rtsne.multicore::Rtsne.multicore(nonrand@exprs[,-excludedMarkers], dims = 2)

set.seed(123)
tsne_out_rand_type1 <- Rtsne.multicore::Rtsne.multicore(rand.type1@exprs[,-excludedMarkers], dims = 2)

set.seed(123)
tsne_out_rand_type2 <- Rtsne.multicore::Rtsne.multicore(rand.type2@exprs[,-excludedMarkers], dims = 2)

set.seed(123)
tsne_out_rand_extreme <- Rtsne.multicore::Rtsne.multicore(rand.extreme@exprs[,-excludedMarkers], dims = 2)

## Plot the tSNE embeddings
plot.tSNE(tsne_out_nonrand$Y,'Non-randomized')
plot.tSNE(tsne_out_rand_type1$Y,'Randomized, type 1')
plot.tSNE(tsne_out_rand_type2$Y,'Randomized, type 2')
plot.tSNE(tsne_out_rand_extreme$Y,'Randomized, extreme')

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


M <- allMarkers[-excludedMarkers,2]
colnames(tSNE.rand.extreme)[3:ncol(tSNE.rand.extreme)] <- colnames(tSNE.rand.type2)[3:ncol(tSNE.rand.type2)] <- colnames(tSNE.rand.type1)[3:ncol(tSNE.rand.type1)] <- colnames(tSNE.nonrand)[3:ncol(tSNE.nonrand)] <- M

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
  
  p1 <- plot_grid(nrand_label,rand_label_type1, rand_label_type2,rand_label_extreme, # rand_label_type2,
                  g1 + theme(legend.position="none",plot.margin=margin(t=-0.5,b=0.1, unit="cm")), 
                  g2 + theme(legend.position="none",plot.margin=margin(t=-0.5,b=0.1,unit="cm")),
                  g3 + theme(legend.position="none",plot.margin=margin(t=-0.5,b=0.1,unit="cm")),
                  g4 + theme(legend.position="none",plot.margin=margin(t=-0.5,b=0.1,unit="cm")),
                  nrow = 2, rel_heights = c(.1,1)) #
  
  p3 <- plot_grid(marker_label, p1, rel_widths = c(.1,1))
  
  p.all[[m]] <- p3 + annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1, color = "gray95", fill = NA)
  
  # g1 <- plot.Marker.on.tSNE(tSNE.nonrand[,c(1,2,(2+m))], marker.name = paste0(M[m],"Non-randomized")) 
  # g2 <- plot.Marker.on.tSNE(tSNE.rand[,c(1,2,(2+m))], marker.name = paste0(M[m],", Randomized"))
  # plot_grid(g1 + theme(legend.position="none"), g2 + theme(legend.position="none"), get_legend(g1), nrow = 1, rel_widths = c(1,1,0.3))
  
}
p.all[[m+1]] <- get_legend(g1)

p <- plot_grid(plotlist = p.all, nrow = 8, ncol = 4)
p <- add_sub(p,expression(paste(bold("Supplementary Figure 4")," tSNE analysis results for all markers. Each subplot displays the marker distribution when the")), 
             x = 0, y = 0, hjust = 0, vjust = 0, size = 8)
p <- add_sub(p,"non-randomized data (left) and the randomized data are employed",
             x = 0, y = 1, hjust = 0, vjust = 1.3, size = 8)


## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## #
#                                                              #
# Ensemble clustering using SOM (Nowicka et al., F1000, 2017)  #
#                                                              #
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## #

## Functions to use
ensemble_clustering_SOM <- function(data, col.used, nclust = 20, seed = 0){
  
  if (seed == 1) {
    set.seed(1234)
  }
  
  som <- BuildSOM(data, colsToUse = col.used)
  
  ## Get the cell clustering into 100 SOM codes
  cell_clustering_som <- som$map$mapping[,1]
  
  ##
  ## Metaclustering using ConsensuClusterPlus
  ##
  
  codes <- som$map$codes
  plot_outdir <- "consensus_plots"
  nmc <- nclust
  if (seed == 1) {
    mc <- ConsensusClusterPlus::ConsensusClusterPlus(t(codes), maxK = nmc, reps =100,pItem =0.9, pFeature =1, title = plot_outdir, plot ="png",clusterAlg ="hc", innerLinkage ="average",finalLinkage ="average",distance ="euclidean", seed =1234)
  } else {
    mc <- ConsensusClusterPlus::ConsensusClusterPlus(t(codes), maxK = nmc, reps =100,pItem =0.9, pFeature =1, title = plot_outdir, plot ="png",clusterAlg ="hc", innerLinkage ="average",finalLinkage ="average",distance ="euclidean")
  }
  
  ## Get cluster ids for each cell
  code_clustering <- mc[[nmc]]$consensusClass
  clusterID <- code_clustering[cell_clustering_som]
  
  return(clusterID)
}
plot_clustering_heatmap_wrapper <- function(expr, expr01,  cell_clustering, color_clusters, cluster_merging =NULL, annot = TRUE){
  # Calculate the median expression
  expr_median <- data.frame(expr, cell_clustering = cell_clustering) %>% group_by(cell_clustering) %>% summarize_all(funs(median))
  expr01_median <- data.frame(expr01, cell_clustering = cell_clustering) %>% group_by(cell_clustering) %>% summarize_all(funs(median))
  # Calculate cluster frequencies
  clustering_table <- as.numeric(table(cell_clustering))
  clustering_prop <- round(clustering_table /sum(clustering_table) *100, 2)
  # Sort the cell clusters with hierarchical clustering  
  d <- dist(expr_median[, colnames(expr)], method ="euclidean")
  cluster_rows <- hclust(d, method ="average")
  expr_heat <- as.matrix(expr01_median[, colnames(expr01)])
  rownames(expr_heat) <- expr01_median$cell_clustering
  # Colors for the heatmap
  color_heat <- colorRampPalette(rev(brewer.pal(n =9, name ="RdYlBu")))(50)
  heat_breaks = seq(from =0, to =1, by =0.02)
  legend_breaks = seq(from =0, to =1, by =0.2)
  labels_row <- paste0(expr01_median$cell_clustering, " (", clustering_prop ,"%)")
  # Annotation for the original clusters
  annotation_row <- data.frame(Cluster = factor(expr01_median$cell_clustering))
  rownames(annotation_row) <- rownames(expr_heat)  
  color_clusters1 <- color_clusters[1:nlevels(annotation_row$Cluster)]
  names(color_clusters1) <- levels(annotation_row$Cluster)
  annotation_colors <- list(Cluster = color_clusters1)
  # Annotation for the merged clusters
  if(!is.null(cluster_merging)){
    cluster_merging$new_cluster <- factor(cluster_merging$new_cluster)    
    annotation_row$Cluster_merging <- cluster_merging$new_cluster
    color_clusters2 <- color_clusters[1:nlevels(cluster_merging$new_cluster)]
    names(color_clusters2) <- levels(cluster_merging$new_cluster)
    annotation_colors$Cluster_merging <- color_clusters2
  }
  
  p<- pheatmap::pheatmap(expr_heat, 
                         breaks = heat_breaks,
                         color = color_heat, 
                         cluster_cols =FALSE, 
                         cluster_rows = cluster_rows, 
                         labels_row = labels_row,
                         display_numbers =TRUE, 	
                         number_color ="black",
                         fontsize =8, fontsize_number =6,  
                         legend_breaks = legend_breaks,
                         annotation_legend = annot,
                         annotation_row = annotation_row, 
                         annotation_colors = annotation_colors)
  return(p)
}
scale01 <- function(data){
  rng <- matrixStats::colQuantiles(data, probs = c(0.01, 0.99))
  data01 <- t((t(data) - rng[, 1]) / (rng[, 2] - rng[, 1]))
  data01[data01 <0] <- 0
  data01[data01 >1] <- 1
  return(data01)
}
cluster.matches <- function(ID1,ID2){
  C <- table(ID1,ID2) # Create the contigency table
  R <- C/rowSums(C) # Recall
  P <- t(t(C)/colSums(C)) # Precision
  FM1 <- 2*P*R/(P+R) # F measure matrix
  FM1[is.na(FM1)] <- 0
  cm <- clue::solve_LSAP(FM1,maximum = TRUE)
  return(cm)
}
metric100 <- function(dataSOM1,dataSOM2){
  temp.nr <- temp.r <- matrix(rep(0,nrow(dataSOM1$data)*100),nrow(dataSOM1$data),100)
  FM.100 <- ARI.100 <- rep(0,100)
  for (i in 1:100){
    temp.nr[,i] <- ensemble_clustering_SOM(dataSOM1, lineage.Markers, 20)
    temp.r[,i] <- ensemble_clustering_SOM(dataSOM2, lineage.Markers, 20)
    FM.100[i] <- FMeasure(temp.nr[,i],temp.r[,i])
    ARI.100[i] <- fossil::adj.rand.index(temp.nr[,i],temp.r[,i])
  }
  (mean(FM.100))
  (mean(ARI.100))
  return(list(FM.100, ARI.100))
}

## Transform the data into a FlowSOM object
nonrand.SOM <- ReadInput('nonrandomized.fcs', transform =FALSE, scale =FALSE)
nonrand.SOM$data[,3:58] <- asinh(nonrand.SOM$data[,3:58]/5) # arcsinh transform

rand.type1.SOM <- ReadInput('randomized_type1.fcs', transform =FALSE, scale =FALSE)
rand.type1.SOM$data[,3:58] <- asinh(rand.type1.SOM$data[,3:58]/5) # arcsinh transform

rand.type2.SOM <- ReadInput('randomized_type2.fcs', transform =FALSE, scale =FALSE)
rand.type2.SOM$data[,3:58] <- asinh(rand.type2.SOM$data[,3:58]/5) # arcsinh transform

rand.extreme.SOM <- ReadInput('randomized_maximal.fcs', transform =FALSE, scale =FALSE)

## Scale marker abundances to values between 0 and 1 (only for visualization)
nonrand01 <- scale01(nonrand.SOM$data)
rand01.type1 <- scale01(rand.type1.SOM$data)
rand01.type2 <- scale01(rand.type2.SOM$data)
rand01.extreme <- scale01(rand.extreme.SOM$data)

## Name the columns with the names of the markers
MarkerNames <- allMarkers$CytkPanel
MarkerNames[c(16,19,26,30,38,45,47)]<-c("CCR6","CD195","CD278","CD194","HLA_DR","CD152","CD279")
colnames(nonrand.SOM$data) <- colnames(nonrand01) <- MarkerNames
colnames(rand.type1.SOM$data) <- colnames(rand01.type1) <- MarkerNames
colnames(rand.type2.SOM$data) <- colnames(rand01.type2) <- MarkerNames
colnames(rand.extreme.SOM$data) <- colnames(rand01.extreme) <- MarkerNames

## Perform Ensemble clustering using SOM
somID.nonrand <- ensemble_clustering_SOM(nonrand.SOM, lineage.Markers, 20, seed = 1)
somID.rand.type1 <- ensemble_clustering_SOM(rand.type1.SOM, lineage.Markers, 20, seed = 1)
somID.rand.type2 <- ensemble_clustering_SOM(rand.type2.SOM, lineage.Markers, 20, seed = 1)
somID.rand.extreme <- ensemble_clustering_SOM(rand.extreme.SOM, lineage.Markers, 20, seed = 1)

## Compare the two clusterings using the F1 measure
FM.type1 <- FMeasure(somID.nonrand,somID.rand.type1)
FM.type2 <- FMeasure(somID.nonrand,somID.rand.type2)
FM.extreme <- FMeasure(somID.nonrand,somID.rand.extreme)

## Compare the two clusterings using the Adjusted Rand Index (0..1)
ARI.type1 <- fossil::adj.rand.index(somID.nonrand,somID.rand.type1)
ARI.type2 <- fossil::adj.rand.index(somID.nonrand,somID.rand.type2)
ARI.extreme <- fossil::adj.rand.index(somID.nonrand,somID.rand.extreme)

## Employ the Hungarian algorithm to find the best matching between the two clusterings
cluster.matches.type1 <- cluster.matches(somID.rand.type1,somID.nonrand)
cluster.matches.type2 <- cluster.matches(somID.rand.type2,somID.nonrand)
cluster.matches.extreme <- cluster.matches(somID.rand.extreme,somID.nonrand)

## Repeat 100 times to check the average discrepancy
metric.100.type1 <- metric100(nonrand.SOM,rand.type1.SOM)
metric.100.type2 <- metric100(nonrand.SOM,rand.type2.SOM)
metric.100.extreme <- metric100(nonrand.SOM,rand.extreme.SOM)

## Create clustering heatmaps
color_clusters <- c("#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72","#B17BA6", "#FF7F00", "#FDB462", "#E7298A", "#E78AC3",
                    "#33A02C", "#B2DF8A", "#55A1B1", "#8DD3C7", "#A6761D","#E6AB02", "#7570B3", "#BEAED4", "#666666", "#999999")

plot_clustering_heatmap_wrapper(expr = nonrand.SOM$data[, lineage.Markers],expr01 = nonrand01[, lineage.Markers],
                                cell_clustering = somID.nonrand, color_clusters = color_clusters, annot = FALSE)

plot_clustering_heatmap_wrapper(expr = rand.type1.SOM$data[, lineage.Markers],expr01 = rand01.type1[, lineage.Markers],
                                cell_clustering = cluster.matches.type1[somID.rand.type1], color_clusters = color_clusters, annot = FALSE)

plot_clustering_heatmap_wrapper(expr = rand.type2.SOM$data[, lineage.Markers],expr01 = rand01.type1[, lineage.Markers],
                                cell_clustering = cluster.matches.type2[somID.rand.type2], color_clusters = color_clusters, annot = FALSE)

plot_clustering_heatmap_wrapper(expr = rand.extreme.SOM$data[, lineage.Markers],expr01 = rand01.type1[, lineage.Markers],
                                cell_clustering = cluster.matches.extreme[somID.rand.extreme], color_clusters = color_clusters, annot = FALSE)


## Plot the Clusters on the tSNE map (run tSNE in previous section before this one)
SOM.nonrand <- data.frame(tSNE1 = tsne_out_nonrand$Y[,1],
                          tSNE2 = tsne_out_nonrand$Y[,2],  
                          nonrand.SOM$data[, lineage.Markers],
                          clusterID = factor(somID.nonrand, levels = 1:20))

SOM.rand.type1 <- data.frame(tSNE1 = tsne_out_rand_type1$Y[,1],
                             tSNE2 = tsne_out_rand_type1$Y[,2],  
                             rand.type1.SOM$data[, lineage.Markers],
                             clusterID = factor(somID.rand.type1, levels = 1:20))

SOM.rand.type2 <- data.frame(tSNE1 = tsne_out_rand_type2$Y[,1],
                             tSNE2 = tsne_out_rand_type2$Y[,2],  
                             rand.type2.SOM$data[, lineage.Markers],
                             clusterID = factor(somID.rand.type2, levels = 1:20))

SOM.rand.extreme <- data.frame(tSNE1 = tsne_out_rand_extreme$Y[,1],
                               tSNE2 = tsne_out_rand_extreme$Y[,2],  
                               rand.extreme.SOM$data[, lineage.Markers],
                               clusterID = factor(somID.rand.extreme, levels = 1:20))


ggplot(SOM.nonrand,  aes(x = tSNE1, y = tSNE2, color = clusterID)) +
  geom_point(size = 0.1, shape = 20) + theme_bw() + scale_color_manual(values = color_clusters[1:20]) + 
  guides(color = guide_legend(override.aes = list(size =4), ncol =2))

# don't forget to match the colors before ploting the next figures
ggplot(SOM.rand.type1, aes(x = tSNE1, y = tSNE2, color = clusterID)) +
  geom_point(size = 0.1, shape = 20) + theme_bw() + scale_color_manual(values = color_clusters[cluster.matches.type1]) + 
  guides(color = guide_legend(override.aes = list(size =4), ncol =2))

ggplot(SOM.rand.type2, aes(x = tSNE1, y = tSNE2, color = clusterID)) +
  geom_point(size = 0.1, shape = 20) + theme_bw() + scale_color_manual(values = color_clusters[cluster.matches.type2]) + 
  guides(color = guide_legend(override.aes = list(size =4), ncol =2))

ggplot(SOM.rand.extreme, aes(x = tSNE1, y = tSNE2, color = clusterID)) +
  geom_point(size = 0.1, shape = 20) + theme_bw() + scale_color_manual(values = color_clusters[cluster.matches.extreme]) + 
  guides(color = guide_legend(override.aes = list(size =4), ncol =2))
