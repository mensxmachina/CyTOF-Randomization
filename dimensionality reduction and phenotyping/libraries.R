# install and load the required libraries
# libraries fossil and ConsensusClusterPlus conflict. we will avoid loading them together.

## CRAN libraries
if(!require(tidyverse)){
  install.packages('tidyverse')
  require(tidyverse)
}

if(!require(cowplot)){
  install.packages('cowplot')
  require(cowplot)
}

if(!require(RColorBrewer)){
  install.packages('RColorBrewer')
  require(RColorBrewer)
}

if(!require(r.jive)){
  install.packages('r.jive')
  require(r.jive)
}

if(!require(fossil)){
  install.packages('fossil')
  require('fossil')
}
detach("package:fossil", unload=TRUE)

if(!require(clue)){
  install.packages('clue')
  require(clue)
}

## Bioconductor libraries
if (!requireNamespace('BiocManager', quietly = TRUE)){
  install.packages('BiocManager')
  BiocManager::install()
}

if(!require(flowCore)){
  BiocManager::install('flowCore')
  require(flowCore)
}

if(!require(limma)){
  BiocManager::install('limma')
  require(limma)
}

if(!require(FlowSOM)){
  BiocManager::install('FlowSOM')
  require(FlowSOM)
}

if(!require(ConsensusClusterPlus)){
  BiocManager::install('ConsensusClusterPlus')
  require(ConsensusClusterPlus)
}

if(!require(pheatmap)){
  BiocManager::install('pheatmap')
  require(pheatmap)
}

if(!require(devtools)) install.packages("devtools")
if(!require(Rtsne.multicore)){
  devtools::install_github("RGLab/Rtsne.multicore")
  require(Rtsne.multicore)
}
