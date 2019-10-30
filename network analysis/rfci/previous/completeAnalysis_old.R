
#### Reconstructing network in randomized vs non-randomized data on several files ####

# This script will apply the fci method on all non-randomized and randomized files, 
# comparing the results between the two different preprocessing

#### Let's set up the analysis ####

#cleaning memory
rm(list =ls())

#set seed
set.seed(12345)

#libraries
if(!require(flowCore)){
  install.packages('flowcore')
  require(flowCore)
}
if(!require(pcalg)){
  install.packages('pcalg')
  require(pcalg)
}

#indipendence test
source('independenceTest.R')

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
               'num_edges_A', 'num_edges_B', 'SHD')
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
    
    #Adding the failed markers and looking in the phos panel
    markersToEliminate <- as.numeric(unlist(sapply(c(toEliminate, toEliminatePhos), grep, markerNames$PhosPanel)))
    markersToKeep <- markerNames$channel[setdiff(1:(dim(markerNames)[1]), markersToEliminate)]
    newMarkerNames <- markerNames$PhosPanel[setdiff(1:(dim(markerNames)[1]), markersToEliminate)];
    
  }else{
    
    #only cytk markers
    markersToEliminate <- as.numeric(unlist(sapply(toEliminate, grep, markerNames$PhosPanel)))
    markersToKeep <- markerNames$channel[setdiff(1:(dim(markerNames)[1]), markersToEliminate)]
    newMarkerNames <- markerNames$PhosPanel[setdiff(1:(dim(markerNames)[1]), markersToEliminate)];
    
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
  
  #applying the fci method to the non-randomized data
  suffStat <- list()
  suffStat$varNames <- colnames(nonRandom)
  suffStat$dataset <- nonRandom
  timeNonRandom <- system.time(nonRandomRes <- fci(suffStat = suffStat, indepTest = bicTest, alpha = 0.1,
                      p = length(newMarkerNames), maj.rule = TRUE))
  nonRandomGraph <- nonRandomRes@amat
  rownames(nonRandomGraph) <- colnames(nonRandomGraph) <- newMarkerNames
  nonRandomResults[[i]] <- nonRandomGraph;
  
  #applying the fci method to the randomized data
  suffStat <- list()
  suffStat$varNames <- colnames(random)
  suffStat$dataset <- random
  timeRandom <- system.time(randomRes <- fci(suffStat = suffStat, indepTest = bicTest, alpha = 0.1,
                                             p = length(newMarkerNames), maj.rule = TRUE))
  randomGraph <- randomRes@amat
  rownames(randomGraph) <- colnames(randomGraph) <- newMarkerNames
  randomResults[[i]] <- randomGraph;
  
  #storing fci statistics
  summaryStats[i, 'num_edges_A'] <- sum(nonRandomGraph[idx] != 0)
  summaryStats[i, 'num_edges_B'] <- sum(randomGraph[idx] != 0)
  summaryStats[i, 'SHD'] <- shd(as(getGraph(nonRandomRes), 'graphNEL'), 
                                as(getGraph(randomRes), 'graphNEL'))
  
  #total time
  message(paste0('Total time: ', round((timeNonRandom['elapsed'] + timeRandom['elapsed'])/60, 3), ' minutes'))
  
}

#### saving ####
write.csv(summaryStats, file = '../../Results/fci/completeAnalysisResults.csv', row.names = FALSE, quote = FALSE)
save(nonRandomResults, randomResults, summaryStats, file = '../../Results/fci/completeAnalysisResults.RData')
