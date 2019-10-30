
#### Reconstructing network in randomized vs non-randomized data on several files ####

# This script will apply the glasso method on all non-randomized and randomized files, 
# comparing the results between the two different preprocessing

#### Let's set up the analysis ####

#cleaning memory
rm(list =ls())

#set seed
set.seed(12345)

#libraries

#libraries
if(!require(flowCore)){
  source('http://www.bioconductor.org/biocLite.R')
  biocLite('flowCore')
  require(flowCore)
}
if(!require(huge)){
  install.packages('huge')
  require(huge)
}

#small / large correlation threshold
smallCorr <- 0.01;
largeCorr <- 0.1;

#characteristics of the files to read
annotations <- read.csv('../../Data/annotation.csv', header = TRUE, sep = ';', stringsAsFactors = FALSE)
nFiles <- dim(annotations)[1];
nonRandomFiles <- dir('../../Data/gated data/', pattern = '.fcs')
randomFiles <- dir('../../Data/gated_randomized/', pattern = '.fcs')

#markers names
markerNames <- read.csv('../../Data/markers.csv', header = TRUE, stringsAsFactors = FALSE)

#columns to eliminate
toEliminate <- c("Barcode", "Time", "cell length", "unkn", "EQBeads", "DNA");
toEliminatePhos <- c("TCRab", "\\bCD4\\b", "ICOS", "CCR6", "CD57", "CD28")

#lists to save the graphs
nonRandomResults <- vector('list', nFiles);
names(nonRandomResults) <- annotations$exp;
randomResults <- nonRandomResults

#matrix with summary statistics. A = nonRandom, B = random
toCompute <- c('exp', 'panel', 'num_samples', 'num_markers', 'num_possible_edges', 
               'average_abs_corr_A', 'average_abs_corr_B', 
               'average_corr_diff', 'per_small_corr_diff', 'per_large_corr_diff',
               'num_edges_A', 'num_edges_B', 
               'num_diff_edges', 'edges_in_A_not_in_B', 'edges_in_B_not_in_A', 
               'average_abs_corr_diff_edges_A', 'average_abs_corr_diff_edges_B')
nToCompute <- length(toCompute)
summaryStats <- data.matrix(matrix(NA, nFiles, nToCompute));
colnames(summaryStats) <- toCompute;

#### Performing the analyses ####

#looping over the files
for(i in 1:nFiles){
  
  #current experiment 
  currentExp <- annotations$exp[i];
  currentPanel <- annotations$panel[i];
  summaryStats[i, 'exp'] <- currentExp
  summaryStats[i, 'panel'] <- currentPanel
  message(paste0('Current experiment: ', currentExp, ', number ', i , ' of ', nFiles))
  
  #finding out the markers to be kept
  if(currentPanel == 1){ # 1: Phos; 2: Cytk
    
    #looking in the phos panel
    markersToEliminate <- as.numeric(unlist(sapply(c(toEliminate, toEliminatePhos), grep, markerNames$PhosPanel)))
    markersToKeep <- markerNames$channel[setdiff(1:(dim(markerNames)[1]), markersToEliminate)]
    newMarkerNames <- markerNames$PhosPanel[setdiff(1:(dim(markerNames)[1]), markersToEliminate)];
    
  }else{
    
    #looking in the cytk and phos markers
    markersToEliminate1 <- as.numeric(unlist(sapply(toEliminate, grep, markerNames$CytkPanel)))
    markersToEliminate2 <- as.numeric(unlist(sapply(toEliminatePhos, grep, markerNames$PhosPanel)))
    markersToEliminate <- c(markersToEliminate1, markersToEliminate2)
    markersToKeep <- markerNames$channel[setdiff(1:(dim(markerNames)[1]), markersToEliminate)]
    newMarkerNames <- markerNames$CytkPanel[setdiff(1:(dim(markerNames)[1]), markersToEliminate)];
    
  }
  
  #storing info
  summaryStats[i, 'num_markers'] <- length(markersToKeep)
  summaryStats[i, 'num_possible_edges'] <- length(markersToKeep) * (length(markersToKeep) - 1) / 2
  
  #loading non-random data
  fileName <- nonRandomFiles[grep(currentExp, nonRandomFiles)]
  nonRandom <- read.FCS(file.path('../../Data/gated data', fileName))
  
  #matrix of measurements
  nonRandom <- exprs(nonRandom)
 
  #info on samples
  summaryStats[i, 'num_samples'] <- dim(nonRandom)[1]
   
  #eliminating columns related to technical measurements
  nonRandom <- nonRandom[ , markersToKeep]
  colnames(nonRandom) <- newMarkerNames;
  nVars <- dim(nonRandom)[2]
  
  # transforming with the arcsinh
  nonRandom <- asinh(nonRandom/5)
  
  #reading random data
  fileName <- randomFiles[grep(currentExp, nonRandomFiles)]
  random <- read.FCS(file.path('../../Data/gated_randomized', fileName))
  
  #measurement matrix
  random <- exprs(random)
  
  #selecting columns 
  random <- random[ , markersToKeep]
  colnames(random) <- newMarkerNames

  #correlation matrices
  nonRandomCor <- cor(nonRandom)
  randomCor <- cor(random)
  
  #difference between unique values (upper matrix excluding the diagonal)
  idx <- upper.tri(nonRandomCor, diag = FALSE)
  diffCor <- as.numeric(abs(nonRandomCor[idx] - randomCor[idx]))
  
  #storing correlations statistics
  summaryStats[i, 'average_abs_corr_A'] <- mean(abs(nonRandomCor[idx]))
  summaryStats[i, 'average_abs_corr_B'] <- mean(abs(randomCor[idx]))
  summaryStats[i, 'average_corr_diff'] <- mean(diffCor)
  summaryStats[i, 'per_small_corr_diff'] <- sum(diffCor <= smallCorr) / length(diffCor) * 100
  summaryStats[i, 'per_large_corr_diff'] <- sum(diffCor >= largeCorr) / length(diffCor) * 100
  
  #applying the glasso method to the non-randomized data
  nonRandomRes <- huge(nonRandom, method = 'glasso')
  nonRandomRes <- huge.select(nonRandomRes)
  nonRandomGraph <- as.matrix(nonRandomRes$refit)
  rownames(nonRandomGraph) <- colnames(nonRandomGraph) <- newMarkerNames
  nonRandomResults[[i]] <- nonRandomGraph;
  
  #applying the glasso method to the randomized data
  randomRes <- huge(random, method = 'glasso')
  randomRes <- huge.select(randomRes)
  randomGraph <- as.matrix(randomRes$refit)
  rownames(randomGraph) <- colnames(randomGraph) <- newMarkerNames
  randomResults[[i]] <- randomGraph;
  
  #storing glasso statistics
  summaryStats[i, 'num_edges_A'] <- sum(nonRandomGraph[idx])
  summaryStats[i, 'num_edges_B'] <- sum(randomGraph[idx])
  summaryStats[i, 'num_diff_edges'] <- sum(abs(nonRandomGraph[idx] - randomGraph[idx]))
  summaryStats[i, 'edges_in_A_not_in_B'] <- sum(nonRandomGraph[idx] == 1 & randomGraph[idx] == 0)
  summaryStats[i, 'edges_in_B_not_in_A'] <- sum(nonRandomGraph[idx] == 0 & randomGraph[idx] == 1)
  summaryStats[i, 'average_abs_corr_diff_edges_A'] <- mean(abs(nonRandomCor[nonRandomGraph != randomGraph]))
  summaryStats[i, 'average_abs_corr_diff_edges_B'] <- mean(abs(randomCor[nonRandomGraph != randomGraph]))
  
}

#### saving ####
write.csv(summaryStats, file = '../../Results/glasso/completeAnalysisResults.csv', row.names = FALSE, quote = FALSE)
save(nonRandomResults, randomResults, summaryStats, file = '../../Results/glasso/completeAnalysisResults.RData')
#save(summaryStats, file = 'completeAnalysisResults.RData')
