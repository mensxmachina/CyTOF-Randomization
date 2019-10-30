
#### Reconstructing network in randomized vs non-randomized data on several files ####

# This script will apply the rfci method on all non-randomized and randomized files, 
# comparing the results between the two different preprocessing

#### Let's set up the analysis ####

#cleaning memory
rm(list =ls())

#set seed
set.seed(12345)

#libraries
if(!require(flowCore)){
  source("https://bioconductor.org/biocLite.R")
  biocLite("flowCore")
  require(flowCore)
}
if(!require(pcalg)){
  install.packages('pcalg')
  require(pcalg)
}
if(!require(parallel)){
  install.packages('parallel')
  require(parallel)
}
if(!require(ParallelPC)){
  install.packages('ParallelPC')
  require(ParallelPC)
}

#indipendence test
source('independenceTest.R')
source('functionRfciParallel.R')
source('functionSkeletonParallel.R')
source('PAG_SHD.R')

#small / large correlation threshold
smallCorr <- 0.01;
largeCorr <- 0.1;

#percentage of samples to keep
percSamples <- 0.05

#characteristics of the files to read
nonRandomFile <- 'nonrandomized.fcs'
randomFiles <- setdiff(dir('../../Data/', pattern = '.fcs'), nonRandomFile)
nFiles <- length(randomFiles);

#markers names
markerNames <- read.csv('../../Data/markers.csv', header = TRUE, stringsAsFactors = FALSE)

#columns to eliminate
toEliminate <- c("Barcode", "Time", "cell length", "unkn", "EQBeads", "DNA");
toEliminatePhos <- c("TCRab", "\\bCD4\\b", "ICOS", "CCR6", "CD57", "CD28")

#lists to save the graphs
randomResults <- vector('list', nFiles);
names(randomResults) <- randomFiles;
randomGraphs <- randomResults

#matrix with summary statistics. A = nonRandom, B = random
toCompute <- c('exp', 'num_samples', 'num_markers', 'num_possible_edges', 
               'average_abs_corr_A', 'average_abs_corr_B', 
               'average_corr_diff', 'per_small_corr_diff', 'per_large_corr_diff',
               'num_edges_A', 'num_edges_B', 'num_diff_edges', 
               'edges_in_A_not_in_B', 'edges_in_B_not_in_A', 'SHD')
nToCompute <- length(toCompute)
summaryStats <- data.matrix(matrix(NA, nFiles, nToCompute));
colnames(summaryStats) <- toCompute;

#### Performing the analyses ####

#looping over the files
for(i in 1:nFiles){
  
  #current experiment 
  currentExp <- randomFiles[i];
  summaryStats[i, 'exp'] <- currentExp
  message(paste0('Current experiment: ', currentExp, ', number ', i , ' of ', nFiles))
  
  #finding out the markers to be kept
  markersToEliminate1 <- as.numeric(unlist(sapply(toEliminate, grep, markerNames$CytkPanel)))
  markersToEliminate2 <- as.numeric(unlist(sapply(toEliminatePhos, grep, markerNames$PhosPanel)))
  markersToEliminate <- c(markersToEliminate1, markersToEliminate2)
  markersToKeep <- markerNames$channel[setdiff(1:(dim(markerNames)[1]), markersToEliminate)]
  newMarkerNames <- markerNames$PhosPanel[setdiff(1:(dim(markerNames)[1]), markersToEliminate)];

  #storing info
  summaryStats[i, 'num_markers'] <- length(markersToKeep)
  summaryStats[i, 'num_possible_edges'] <- length(markersToKeep) * (length(markersToKeep) - 1) / 2
  
  #loading non-random data
  nonRandom <- read.FCS(file.path('../../Data', nonRandomFile))
  
  #matrix of measurements
  nonRandom <- exprs(nonRandom)
  
  #keeping only a portion of the data
  samplesToKeep <- sample(1:dim(nonRandom)[1], dim(nonRandom)[1] * percSamples)
  if(length(samplesToKeep) <= 1){
    next
  }
  nonRandom <- nonRandom[samplesToKeep, ]
 
  #info on samples
  summaryStats[i, 'num_samples'] <- dim(nonRandom)[1]
   
  #eliminating columns related to technical measurements
  nonRandom <- nonRandom[ , markersToKeep]
  nVars <- dim(nonRandom)[2]
  
  # transforming with the arcsinh
  nonRandom <- asinh(nonRandom/5)
  
  #reading random data
  random <- read.FCS(file.path('../../Data', randomFiles[i]))
  
  #measurement matrix
  random <- exprs(random)
  
  #keeping only a portion of the data
  random <- random[samplesToKeep, ]
  
  #selecting columns 
  random <- random[ , markersToKeep]
  
  # transforming with the arcsinh
  if(!grepl('extreme', randomFiles[i])){
    random <- asinh(random/5)
  }  

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
  
  #applying the rfci method to the non-randomized data
  suffStat <- list()
  suffStat$varNames <- colnames(nonRandom)
  suffStat$dataset <- nonRandom
  if(i == 1){
    timeNonRandom <- system.time(nonRandomRes <- rfci_parallel(suffStat = suffStat, indepTest = bicTest, 
                                                               alpha = 0.1, p = length(newMarkerNames), 
                                                               maj.rule = TRUE, num.cores = 4))
  }
  nonRandomGraph <- nonRandomRes@amat
  rownames(nonRandomGraph) <- colnames(nonRandomGraph) <- newMarkerNames
  nonRandomResult <- nonRandomRes;
  
  #applying the fci method to the randomized data
  suffStat <- list()
  suffStat$varNames <- colnames(random)
  suffStat$dataset <- random
  timeRandom <- system.time(randomRes <- rfci_parallel(suffStat = suffStat, indepTest = bicTest, 
                                                       alpha = 0.1, p = length(newMarkerNames), 
                                                       maj.rule = TRUE, num.cores = 4))
  randomGraph <- randomRes@amat
  rownames(randomGraph) <- colnames(randomGraph) <- newMarkerNames
  randomGraphs[[i]] <- randomGraph
  randomResults[[i]] <- randomRes;
  
  #storing fci statistics
  summaryStats[i, 'num_edges_A'] <- sum(nonRandomGraph != 0)
  summaryStats[i, 'num_edges_B'] <- sum(randomGraph != 0)
  summaryStats[i, 'SHD'] <- pagSHD(nonRandomGraph, randomGraph)
  summaryStats[i, 'num_diff_edges'] <- sum((nonRandomGraph != 0) !=
                                             (randomGraph != 0))
  summaryStats[i, 'edges_in_A_not_in_B'] <- sum((nonRandomGraph != 0) &
                                             (randomGraph == 0))
  summaryStats[i, 'edges_in_B_not_in_A'] <- sum((nonRandomGraph == 0) &
                                                  (randomGraph != 0))
  
  #total time
  message(paste0('Total time: ', round((timeNonRandom['elapsed'] + timeRandom['elapsed'])/60, 3), ' minutes'))
  
}

#### saving ####
write.csv(summaryStats, file = '../Results/rfci/completeAnalysisResults.csv', row.names = FALSE, quote = FALSE)
save(nonRandomResult, nonRandomGraph, randomResults, randomGraphs, summaryStats, file = '../Results/rfci/completeAnalysisResults.RData')

