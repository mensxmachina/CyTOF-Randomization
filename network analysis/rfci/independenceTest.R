#bic-based test of conditional independence for continous measurements
#the test first fits two models:
# m0: y~s
# m1: y~x+s
#and check which model is better with the BIC criterion
#Then the following models are fit
# m0: x~s
# m1: x~y+s
#and the best model is again found with the BIC criterion
#If any of the two m1 models comes out to be better than 
#its corresponding null m0 model, then the null assumption 
#of independence is rejected (return 0)
#Otherwise, the independence is accepted (return 1)

bicTest <- function(x,y, S, suffStat){
  
  #single tests
  res1 <- singleBicTest(x, y, S, suffStat);
  res2 <- singleBicTest(y, x, S, suffStat);
  
  #return
  if(any(c(res1, res2) == 0)){
    return(0)
  }else{
    return(1)
  }
  
}

#single test
singleBicTest <- function(x, y, S, suffStat){

  # bic1 < bic0
  # 
  # ln(n)(k + 1) - 2 ln(L1) < ln(n)(k) - 2 ln(L0)
  # ln(n) < 2 (ln(L1) - ln(L0))
    
  #collecting the suffStat
  varNames <- suffStat$varNames
  dataset <- as.data.frame(suffStat$dataset)
  nSamples <- dim(dataset)[1]
  
  #fitting the null model
  if(length(S) > 0){
    f0 <- as.formula(paste0(varNames[y], ' ~ ',
                            paste(varNames[S], 
                                  collapse = '+')))
  }else{
    f0 <- as.formula(paste0(varNames[y], ' ~ 1'))
  }
  m0 <- lm(formula = f0, data = dataset)
  
  #fitting the m1 model
  if(length(S) > 0){
    f1 <- as.formula(paste0(varNames[y], ' ~ ',
                            paste(c(varNames[x],
                                  varNames[S]), 
                                  collapse = '+')))
  }else{
    f1 <- as.formula(paste0(varNames[y], ' ~ ', varNames[x]))
  }
  m1 <- lm(formula = f1, data = dataset)
  
  #checking if m1 is better
  isM1Better <- log(nSamples) < 2 * (logLik(m1) - logLik(m0))

  #if m1 is better, then we return 0 (i.e., dependence between x and y given S)
  if(isM1Better){
    return(0)
  }else{
    return(1)
  }
  
}

