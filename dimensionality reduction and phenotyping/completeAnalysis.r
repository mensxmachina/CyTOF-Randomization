##
## DIMENSIONALITY REDUCTION AND CELL PHENOTYPING
##
## This script will load the data and apply linear dimensionality reduction 
## analysis (PCA, JIVE, MDS), nonlinear dimensionality reduction analysis (tSNE),
## and cell phenotyping using ensemble clustering (see paper for more details)
##
## date: 07/05/2019

## libraries
source("libraries.R")

##  Data and annotations
source("readData.R")

## ## ## ## ## ## ## ## ## ## ## ## #
#                                   #
#  Linear dimensionality reduction  #
#                                   #
## ## ## ## ## ## ## ## ## ## ## ## #

## Functions to use for plotting linear dimensionality results
source('plotLinDR.r')

## we subsample without loss of generality. set seed for reproducibility
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

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
#                                             #
#  Non linear dimensionality reduction (tSNE) #
#                                             #
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

## Functions to use
source("plotNonLinDR.R")

## Perform tSNE (no subsampling)
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

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## #
#                                                              #
# Ensemble clustering using SOM (Nowicka et al., F1000, 2017)  #
#                                                              #
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## #

## Functions to use
source("ensembleClustFuncs.R")

## Transform the data into a FlowSOM object
nonrand.SOM <- ReadInput('../data/nonrandomized.fcs', transform =FALSE, scale =FALSE)
nonrand.SOM$data[,3:58] <- asinh(nonrand.SOM$data[,3:58]/5) # arcsinh transform

rand.type1.SOM <- ReadInput('../data/randomized_type1.fcs', transform =FALSE, scale =FALSE)
rand.type1.SOM$data[,3:58] <- asinh(rand.type1.SOM$data[,3:58]/5) # arcsinh transform

rand.type2.SOM <- ReadInput('../data/randomized_type2.fcs', transform =FALSE, scale =FALSE)
rand.type2.SOM$data[,3:58] <- asinh(rand.type2.SOM$data[,3:58]/5) # arcsinh transform

rand.extreme.SOM <- ReadInput('../data/randomized_maximal.fcs', transform =FALSE, scale =FALSE)

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
