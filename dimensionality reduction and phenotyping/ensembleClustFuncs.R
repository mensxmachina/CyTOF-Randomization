# Functions for performing and plotting results on cell phenotyping using ensemble clustering


# Ensemble clustering
# -------------------
# step 1: perform clusteting using self-organizing maps
# step 2: perform mataclusteting using consensus clustering
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


# Plot cluster heatmaps
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


# Scale the cluster data for plotting
scale01 <- function(data){
  rng <- matrixStats::colQuantiles(data, probs = c(0.01, 0.99))
  data01 <- t((t(data) - rng[, 1]) / (rng[, 2] - rng[, 1]))
  data01[data01 <0] <- 0
  data01[data01 >1] <- 1
  return(data01)
}


# Match clusters between randomization analysis results
cluster.matches <- function(ID1,ID2){
  C <- table(ID1,ID2) # Create the contigency table
  R <- C/rowSums(C) # Recall
  P <- t(t(C)/colSums(C)) # Precision
  FM1 <- 2*P*R/(P+R) # F measure matrix
  FM1[is.na(FM1)] <- 0
  cm <- clue::solve_LSAP(FM1,maximum = TRUE)
  return(cm)
}


# Perform ensemble clustering 100 times
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