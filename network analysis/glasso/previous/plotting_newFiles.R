
#### Reconstructing network in randomized vs non-randomized data on several files ####

# This script will draw the results for a specific experiment

#### 

#cleaning memory
rm(list =ls())

#set seed
set.seed(12345)

#experiment to plot
toPlot <- 'PBMC_randomized_extreme.fcs'

#loading the results
load('../../Results/glasso/completeAnalysisResults_newFiles.RData')
load('newMarkerNames.RData')

#selecting the results to plot
idx <- which(names(randomResults) == toPlot)
nonRandomResult <- nonRandomResults
randomResult <- randomResults[[idx]]

# new names
rownames(nonRandomResult) <- colnames(nonRandomResult) <- newMarkerNames
rownames(randomResult) <- colnames(randomResult) <- newMarkerNames

#let's create a cytoscape file
markerNames <- rownames(nonRandomResult)
nMarkers <- length(markerNames)
sink(paste0('../../Results/glasso/', toPlot, '.sif'))
for(i in 1:(nMarkers - 1)){
  for(j in (i + 1):nMarkers){
    if(nonRandomResult[i, j] == 1 && randomResult[i, j] == 1){
      cat(markerNames[i], ' common ', markerNames[j], '\n')
    }
    if(nonRandomResult[i, j] == 1 && randomResult[i, j] == 0){
      cat(markerNames[i], ' NRnotR ', markerNames[j], '\n')
    }
    if(nonRandomResult[i, j] == 0 && randomResult[i, j] == 1){
      cat(markerNames[i], ' RnotNR ', markerNames[j], '\n')
    }
  }
}
sink()
