##  Read the non-randomized and randomized files and the marker annotations

nonrand = read.FCS('../data/nonrandomized.fcs', transform=FALSE)
rand.type1 = read.FCS('../data/randomized_type1.fcs', transform=FALSE)
rand.type2 = read.FCS('../data/randomized_type2.fcs', transform=FALSE)
rand.extreme = read.FCS('../data/randomized_maximal.fcs', transform=FALSE)

# asinh transforms
nonrand@exprs[,3:58] <- asinh(nonrand@exprs[,3:58]/5)
rand.type1@exprs[,3:58] <- asinh(rand.type1@exprs[,3:58]/5)
rand.type2@exprs[,3:58] <- asinh(rand.type2@exprs[,3:58]/5)

## PANEL ANNOTATIONS
allMarkers <- read.csv('../data/markers.csv',stringsAsFactors = F)

## Exclude non Marker and failed Marker variables
nonMarkers.Cytk <- as.numeric(unlist(sapply(c("Barcode","Time","cell length","unkn","EQBeads","DNA"),grep, allMarkers$CytkPanel)))
failedMarkers <- as.numeric(unlist(sapply(c("TCRab", "\\bCD4\\b", "ICOS", "CCR6", "CD57", "CD28"),grep, allMarkers$CytkPanel)))
excludedMarkers <- c(failedMarkers,nonMarkers.Cytk)

## Identify the cytokine and lineage markers
cytokine.Markers <- as.numeric(unlist(sapply(c("IL","IFNg","Granzyme","TNFa","GM-CSF"),grep, allMarkers$CytkPanel)))
lineage.Markers <- setdiff(c(1:NROW(allMarkers)),c(excludedMarkers,cytokine.Markers))
