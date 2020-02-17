library(flowCore)

nonrand = read.FCS('non_randomized.fcs', transform=FALSE)

## RANDOMIZE THE DATA
rand.type1 <- rand.type2 <- rand.extreme <- nonrand

set.seed(123)
for (i in 3:53 ){
  rand.type1@exprs[,i] <- rand.type1@exprs[,i] - runif(dim(rand.type1@exprs)[1]) # type 1 randomization
  rand.type2@exprs[,i] <- rand.type2@exprs[,i] + rnorm(dim(rand.type2@exprs)[1], mean = 0, sd = 1) # type 2 randomization
}

# maximal randomization
rand.extreme@exprs[,3:58] <- asinh(rand.extreme@exprs[,3:58]/5) # arcsinh transform
set.seed(123)
for (i in 3:53 ){
  rand.extreme@exprs[,i] <- rand.extreme@exprs[,i] + rnorm(dim(rand.extreme@exprs)[1], mean = 0, sd = 0.3) # extreme randomization
}

# asinh transforms
nonrand@exprs[,3:58] <- asinh(nonrand@exprs[,3:58]/5)
rand.type1@exprs[,3:58] <- asinh(rand.type1@exprs[,3:58]/5)
rand.type2@exprs[,3:58] <- asinh(rand.type2@exprs[,3:58]/5)


write.FCS(rand.type1,'type1.fcs',endian = "little")
write.FCS(rand.type2,'type2.fcs',endian = "little")
write.FCS(rand.extreme.repl2,'maximal3.fcs',endian = "little")
