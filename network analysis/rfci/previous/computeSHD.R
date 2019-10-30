
#### Reconstructing network in randomized vs non-randomized data on several files ####

# This script will apply the SHD metric to the results of the rfci method

#### Let's set up the analysis ####

#cleaning memory
rm(list =ls())

# loading the SHD function
source('PAG_SHD.R')

# loading the previous results
# save(nonRandomResults, nonRandomGraphs, randomResults, randomGraphs, summaryStats, 
# file = '../../Results/rfci/completeAnalysisResults.RData')
load('../../Results/rfci/completeAnalysisResults.RData')
summaryStats <- data.frame(summaryStats, stringsAsFactors = FALSE)
nFiles <- dim(summaryStats)[1]

# adding the SHD to the results matsrix
summaryStats$SHD <- NA;

#### Performing the analyses ####

#looping over the files
for(i in 1:nFiles){
  
  # computing the SHD
  G1 <- nonRandomGraphs[[i]]
  G2 <- randomGraphs[[i]]
  summaryStats$SHD[i] <- pagSHD(nonRandomGraphs[[i]], randomGraphs[[i]])

}

#### saving ####
write.csv(summaryStats, file = '../../Results/rfci/completeAnalysisResults.csv', row.names = FALSE, quote = FALSE)
save(nonRandomResults, nonRandomGraphs, randomResults, randomGraphs, summaryStats, file = '../../Results/rfci/completeAnalysisResults.RData')

