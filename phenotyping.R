library(flowCore)
library(FlowSOM)
library(RColorBrewer)
library(tidyverse)
library(cowplot)


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

nonrand.SOM <- FlowSOM::ReadInput('non_randomized.fcs', transform =FALSE, scale =FALSE)
nonrand.SOM$data[,3:58] <- asinh(nonrand.SOM$data[,3:58]/5) # arcsinh transform

rand.type1.SOM <- FlowSOM::ReadInput('type1.fcs', transform =FALSE, scale =FALSE)
rand.type1.SOM$data[,3:58] <- asinh(rand.type1.SOM$data[,3:58]/5) # arcsinh transform

rand.type2.SOM <- FlowSOM::ReadInput('type2.fcs', transform =FALSE, scale =FALSE)
rand.type2.SOM$data[,3:58] <- asinh(rand.type2.SOM$data[,3:58]/5) # arcsinh transform

rand.extreme.SOM <- FlowSOM::ReadInput('randomized_maximal.fcs', transform =FALSE, scale =FALSE)
rand.extreme.repl1.SOM <- FlowSOM::ReadInput('randomized_maximal_replicate1.fcs', transform =FALSE, scale =FALSE)
rand.extreme.repl2.SOM <- FlowSOM::ReadInput('randomized_maximal_replicate2.fcs', transform =FALSE, scale =FALSE)

## Scale marker abundances to values between 0 and 1 (only for visualization)
nonrand01 <- scale01(nonrand.SOM$data)
rand01.type1 <- scale01(rand.type1.SOM$data)
rand01.type2 <- scale01(rand.type2.SOM$data)
rand01.extreme <- scale01(rand.extreme.SOM$data)
rand01.extreme.repl1 <- scale01(rand.extreme.repl1.SOM$data)
rand01.extreme.repl2 <- scale01(rand.extreme.repl2.SOM$data)

## PANEL ANNOTATIONS
allMarkers <- read.csv('markers.csv',stringsAsFactors = F)

## Exclude non Marker and failed Marker variables
nonMarkers.Cytk <- as.numeric(unlist(sapply(c("Barcode","Time","cell length","unkn","EQBeads","DNA"),grep, allMarkers$CytkPanel)))
failedMarkers <- as.numeric(unlist(sapply(c("TCRab", "\\bCD4\\b", "ICOS", "CCR6", "CD57", "CD28"),grep, allMarkers$CytkPanel)))
excludedMarkers <- c(failedMarkers,nonMarkers.Cytk)

## Identify the cytokine and lineage markers
cytokine.Markers <- as.numeric(unlist(sapply(c("IL","IFNg","Granzyme","TNFa","GM-CSF"),grep, allMarkers$CytkPanel)))
lineage.Markers <- setdiff(c(1:NROW(allMarkers)),c(excludedMarkers,cytokine.Markers))


## Name the columns with the names of the markers
MarkerNames <- allMarkers$CytkPanel
MarkerNames[c(16,19,26,30,38,45,47)]<-c("CCR6","CD195","CD278","CD194","HLA_DR","CD152","CD279")
colnames(nonrand.SOM$data) <- colnames(nonrand01) <- MarkerNames
colnames(rand.type1.SOM$data) <- colnames(rand01.type1) <- MarkerNames
colnames(rand.type2.SOM$data) <- colnames(rand01.type2) <- MarkerNames
colnames(rand.extreme.SOM$data) <- colnames(rand01.extreme) <- MarkerNames
colnames(rand.extreme.repl1.SOM$data) <- colnames(rand01.extreme) <- MarkerNames
colnames(rand.extreme.repl2.SOM$data) <- colnames(rand01.extreme) <- MarkerNames

## Perform Ensemble clustering using SOM
somID.nonrand <- ensemble_clustering_SOM(nonrand.SOM, lineage.Markers, 20, seed = 1)
somID.rand.type1 <- ensemble_clustering_SOM(rand.type1.SOM, lineage.Markers, 20, seed = 1)
somID.rand.type2 <- ensemble_clustering_SOM(rand.type2.SOM, lineage.Markers, 20, seed = 1)
somID.rand.extreme <- ensemble_clustering_SOM(rand.extreme.SOM, lineage.Markers, 20, seed = 1)
somID.rand.extreme.repl1 <- ensemble_clustering_SOM(rand.extreme.repl1.SOM, lineage.Markers, 20, seed = 1)
somID.rand.extreme.repl2 <- ensemble_clustering_SOM(rand.extreme.repl2.SOM, lineage.Markers, 20, seed = 1)

## Compare the two clusterings using the F1 measure
FM.type1 <- FMeasure(somID.nonrand,somID.rand.type1)
FM.type2 <- FMeasure(somID.nonrand,somID.rand.type2)
FM.extreme <- FMeasure(somID.nonrand,somID.rand.extreme)
FM.extreme.repl1 <- FMeasure(somID.nonrand,somID.rand.extreme.repl1)
FM.extreme.repl2 <- FMeasure(somID.nonrand,somID.rand.extreme.repl2)

## Compare the two clusterings using the Adjusted Rand Index (0..1)
ARI.type1 <- fossil::adj.rand.index(somID.nonrand,somID.rand.type1)
ARI.type2 <- fossil::adj.rand.index(somID.nonrand,somID.rand.type2)
ARI.extreme <- fossil::adj.rand.index(somID.nonrand,somID.rand.extreme)
ARI.extreme.repl1 <- fossil::adj.rand.index(somID.nonrand,somID.rand.extreme.repl1)
ARI.extreme.repl2 <- fossil::adj.rand.index(somID.nonrand,somID.rand.extreme.repl2)

## Employ the Hungarian algorithm to find the best matching between the two clusterings
cluster.matches.type1 <- cluster.matches(somID.rand.type1,somID.nonrand)
cluster.matches.type2 <- cluster.matches(somID.rand.type2,somID.nonrand)
cluster.matches.extreme <- cluster.matches(somID.rand.extreme,somID.nonrand)
cluster.matches.extreme.repl1 <- cluster.matches(somID.rand.extreme.repl1,somID.nonrand)
cluster.matches.extreme.repl2 <- cluster.matches(somID.rand.extreme.repl2,somID.nonrand)

## Repeat 100 times to check the average discrepancy
metric.100.type1 <- metric100(nonrand.SOM,rand.type1.SOM)
metric.100.type2 <- metric100(nonrand.SOM,rand.type2.SOM)
metric.100.extreme <- metric100(nonrand.SOM,rand.extreme.SOM)
metric.100.extreme.repl1 <- metric100(rand.extreme.SOM,rand.extreme.repl1.SOM)
metric.100.extreme.repl2 <- metric100(rand.extreme.SOM,rand.extreme.repl2.SOM)


## Create clustering heatmaps
color_clusters <- c("#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72","#B17BA6", "#FF7F00", "#FDB462", "#E7298A", "#E78AC3",
                    "#33A02C", "#B2DF8A", "#55A1B1", "#8DD3C7", "#A6761D","#E6AB02", "#7570B3", "#BEAED4", "#666666", "#999999")

plot_clustering_heatmap_wrapper(expr = nonrand.SOM$data[, lineage.Markers],expr01 = nonrand01[, lineage.Markers],
                                cell_clustering = somID.nonrand, color_clusters = color_clusters, annot = FALSE)

plot_clustering_heatmap_wrapper(expr = rand.type1.SOM$data[, lineage.Markers],expr01 = rand01.type1[, lineage.Markers],
                                cell_clustering = cluster.matches.type1[somID.rand.type1], color_clusters = color_clusters, annot = FALSE)

plot_clustering_heatmap_wrapper(expr = rand.type2.SOM$data[, lineage.Markers],expr01 = rand01.type2[, lineage.Markers],
                                cell_clustering = cluster.matches.type2[somID.rand.type2], color_clusters = color_clusters, annot = FALSE)

plot_clustering_heatmap_wrapper(expr = rand.extreme.SOM$data[, lineage.Markers],expr01 = rand01.extreme[, lineage.Markers],
                                cell_clustering = cluster.matches.extreme[somID.rand.extreme], color_clusters = color_clusters, annot = FALSE)

plot_clustering_heatmap_wrapper(expr = rand.extreme.repl1.SOM$data[, lineage.Markers],expr01 = rand01.extreme.repl1[, lineage.Markers],
                                cell_clustering = cluster.matches.extreme.repl1[somID.rand.extreme.repl1], color_clusters = color_clusters, annot = FALSE)

plot_clustering_heatmap_wrapper(expr = rand.extreme.repl2.SOM$data[, lineage.Markers],expr01 = rand01.extreme.repl2[, lineage.Markers],
                                cell_clustering = cluster.matches.extreme.repl2[somID.rand.extreme.repl2], color_clusters = color_clusters, annot = FALSE)


## Plot the Clusters on the tSNE map (from the dimensionality reduction results)
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

SOM.rand.extreme.repl1 <- data.frame(tSNE1 = tsne_out_rand_extreme_repl1$Y[,1],
                                     tSNE2 = tsne_out_rand_extreme_repl1$Y[,2],  
                                     rand.extreme.repl1.SOM$data[, lineage.Markers],
                                     clusterID = factor(somID.rand.extreme.repl1, levels = 1:20))

SOM.rand.extreme.repl2 <- data.frame(tSNE1 = tsne_out_rand_extreme_repl2$Y[,1],
                                     tSNE2 = tsne_out_rand_extreme_repl2$Y[,2],  
                                     rand.extreme.repl2.SOM$data[, lineage.Markers],
                                     clusterID = factor(somID.rand.extreme.repl2, levels = 1:20))


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

ggplot(SOM.rand.extreme.repl1, aes(x = tSNE1, y = tSNE2, color = clusterID)) +
  geom_point(size = 0.1, shape = 20) + theme_bw() + scale_color_manual(values = color_clusters[cluster.matches.extreme.repl1]) + 
  guides(color = guide_legend(override.aes = list(size =4), ncol =2))

ggplot(SOM.rand.extreme.repl2, aes(x = tSNE1, y = tSNE2, color = clusterID)) +
  geom_point(size = 0.1, shape = 20) + theme_bw() + scale_color_manual(values = color_clusters[cluster.matches.extreme.repl2]) + 
  guides(color = guide_legend(override.aes = list(size =4), ncol =2))
