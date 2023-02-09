model_gridsParameters = function(data,
                                 learner, 
                                 metric,
                                 samplingmethod, 
                                 numberFOLD, 
                                 numberFoldRepeats,
                                 searchMethod){
  tuneLength = 2000
  
  if(learner == 'glmnet'){
    
    grids = expand.grid(alpha =  c(seq(0, 1, length.out = 5), 0.05, 0.025), 
                        lambda =  c(exp(seq(log(2e-3), log(1e0), length.out = 5)), 0.35, 0.45, 0.67, 0.78, 0.89)
    )
    
    fittingParams <- list(metric = metric,
                          
                          family = "binomial",
                          
                          importance = TRUE,
                          
                          standardize = FALSE,
                          
                          tuneGrid = grids, 
                          
                          tuneLength = tuneLength)
    
    
  }else if (learner == 'lda'){
    
    grids = NULL
    fittingParams = list(metric = metric,
                         tuneLength = tuneLength)
    
  }  else if (learner == 'bayesglm'){
    
    grids = NULL
    fittingParams = list(metric = metric,
                         tuneLength = tuneLength)
    
  }else if (learner == 'svmRadial'){
    
    grids =  expand.grid(sigma = c(0.001, 
                                   0.1, 0.25, 0.5, 0.75,0.9),
                         C = 10^seq(-3, 3, length =15))
    
    fittingParams = list(metric = metric,
                         tuneGrid = grids,
                         prob.model = TRUE,
                         tuneLength = tuneLength)
    
  }else if(learner ==  'svmLinear'){
    
    grids =  expand.grid(C = 10^seq(-3, 3, length = 15))
    
    fittingParams = list(metric = metric,
                         tuneGrid = grids,
                         prob.model = TRUE,
                         tuneLength = tuneLength)
    
    
  }else if(learner ==  'knn'){
    
    
    grids =  expand.grid(k = seq(3,  round(sqrt(nrow(data))), 2))
    
    fittingParams = list(metric = metric,
                         tuneGrid = grids,
                         tuneLength = tuneLength)
    
    
  }else if (learner == 'nnet'){
    
    grids = expand.grid(decay = seq(0.001, 0.9999, length = 5),
                        size = c(1, 2, 3, 5, 10, 7, 12))
    
    fittingParams = list(metric = metric,
                         tuneGrid = grids,
                         tuneLength = tuneLength,
                         trace = FALSE,
                         verbosity = 0)
    
  } else if (learner == "rf"){
    
    sqtmtry<- round(sqrt(ncol(data)))
    grids <- expand.grid(mtry = seq(3, sqtmtry, by=2))
    
    fittingParams = list(metric = metric,
                         tuneGrid = grids,
                         tuneLength = tuneLength,
                         importance=TRUE,
                         ntree=1500)
    
    
  }else if (learner == 'rda'){
    grids = expand.grid(gamma = seq(0.001, 0.99, length = 5), 
                        lambda =  seq(0.001, 0.99,by = 0.2))
    
    fittingParams = list(metric = metric,
                         tuneGrid = grids,
                         tuneLength = tuneLength)
    
  }
  
  else if (learner == 'xgbLinear'){
    grids = expand.grid(lambda = c(0,0.0001,0.001),
                        nrounds= seq(50, 200, by= 50),
                        alpha =  c(0,0.0001,0.001),
                        eta= c(0.1, 0.3, 0.9))
    
    fittingParams = list(metric = metric,
                         tuneGrid = grids,
                         tuneLength = tuneLength)
    
  } else if (learner == 'qda'){
    
    grids = NULL
    fittingParams = list(metric = metric,
                         tuneLength = tuneLength)
    
  }else if (learner == 'rpart1SE'){
    
    grids = NULL
    fittingParams = list(metric = metric,
                         tuneGrid = grids,
                         tuneLength = tuneLength)
    
    
  }else if (learner == 'rpart'){
    grids = expand.grid(cp = seq(0.1, 0.99, length = 10))
    
    fittingParams = list(metric = metric,
                         tuneGrid = grids,
                         tuneLength = tuneLength)
    
  }
  
  
  
  if(samplingmethod == 'LOOCV'){
    
    numberFOLD= dim(data)[1] ;  numberFoldRepeats = 1
    
  } else if((samplingmethod == "boot632") | (samplingmethod == "upSampleboot632") | (samplingmethod == "downsampledboot632") | (samplingmethod == "roseboot632") ){
    
    numberFOLD= numberFOLD ;  numberFoldRepeats =  1
    
  } else{
    
    numberFOLD= numberFOLD ;  numberFoldRepeats =  numberFoldRepeats
  }
  
  
  
  
  if ((learner ==  'lda') |(learner ==  'qda') | (learner == 'bayesglm') | (learner == 'rpart1SE')) {
    ngridsDim= 1
    
  }else {
    
    ngridsDim= dim(grids)[1]
  }
  
  
  if((samplingmethod == "boot632") | (samplingmethod == "upSampleboot632") | (samplingmethod == "downsampledboot632") | (samplingmethod == "roseboot632") ) {
    
    set.seed(1234)
    seeds <- vector(mode = "list", length = (numberFOLD*numberFoldRepeats)+2)
    for(i in 1:((numberFOLD*numberFoldRepeats)+1) ) seeds[[i]] <- sample.int(1000000, ngridsDim)
    ## For the last model:
    seeds[[((numberFOLD*numberFoldRepeats)+2)]] <- sample.int(1000000, 1)
    
  }else{
    
    set.seed(1234)
    seeds <- vector(mode = "list", length = (numberFOLD*numberFoldRepeats)+1)
    for(i in 1:(numberFOLD*numberFoldRepeats)) seeds[[i]] <- sample.int(1000000, ngridsDim)
    ## For the last model:
    seeds[[((numberFOLD*numberFoldRepeats)+1)]] <- sample.int(1000000, 1)
    
  }
  
  
  
  if((samplingmethod == 'stratifiedRepeatedCV')| (samplingmethod =='upSamplestratifiedRepeatedCV') | (samplingmethod == 'downsampledstratifiedRepeatedCV')| (samplingmethod == 'rosestratifiedRepeatedCV')  ){
    
    samplingmethod2 = 'repeatedcv'
    
    set.seed(1234)
    cvIndex <- createMultiFolds(data$response,
                                k = numberFOLD,
                                times = numberFoldRepeats)
    
    
  } else if((samplingmethod == "boot632") | (samplingmethod == "upSampleboot632") | (samplingmethod == "downsampledboot632") | (samplingmethod == "roseboot632")){
    
    samplingmethod2 = "boot632"
    cvIndex = NULL
    
  }else{
    
    samplingmethod2 =  samplingmethod
    cvIndex = NULL
    
  }
  
  
  
  if( (samplingmethod == "downsampledboot632") |  (samplingmethod == 'downsampledstratifiedRepeatedCV') ){
    
    sampling = 'down'
    
  }else if ((samplingmethod == "upSampleboot632") | (samplingmethod =='upSamplestratifiedRepeatedCV')){
    
    sampling = 'up'
    
  }else if ((samplingmethod == "roseboot632") | (samplingmethod == 'rosestratifiedRepeatedCV') ){
    
    sampling = 'rose'
    
  }else {
    
    sampling = NULL
  }
  
  
  
  if((metric == 'Accuracy') | (metric == 'Kappa') ){
    
    if((learner ==  'svmLinear')| (learner == 'svmRadial')){
      
      resamplingParams <- list(index = cvIndex,
                               
                               method = samplingmethod2, 
                               
                               number = numberFOLD,
                               
                               repeats = numberFoldRepeats,
                               
                               verbose = FALSE, 
                               
                               classProbs = TRUE,
                               
                               search = searchMethod,
                               
                               savePredictions = "final",
                               
                               allowParallel = TRUE,
                               
                               sampling = sampling,
                               
                               seeds = seeds)
      
    }else{
      
      resamplingParams <- list(index = cvIndex,
                               
                               method = samplingmethod2, 
                               
                               number = numberFOLD,
                               
                               repeats = numberFoldRepeats,
                               
                               verbose = FALSE, 
                               
                               classProbs = TRUE,
                               
                               search = searchMethod,
                               
                               savePredictions = "final",
                               
                               allowParallel = TRUE,
                               
                               sampling = sampling,
                               
                               seeds = seeds)
    }
  }else{
    if((learner ==  'svmLinear')| (learner == 'svmRadial')){
      
      resamplingParams <- list(index = cvIndex,
                               
                               method = samplingmethod2, 
                               
                               number = numberFOLD,
                               
                               repeats = numberFoldRepeats,
                               
                               verbose = FALSE, 
                               
                               classProbs = TRUE,
                               
                               summaryFunction = twoClassSummary,
                               
                               search = searchMethod,
                               
                               savePredictions = "final",
                               
                               allowParallel = TRUE,
                               
                               sampling = sampling,
                               
                               seeds = seeds)
      
    }else{
      
      resamplingParams <- list(index = cvIndex,
                               
                               method = samplingmethod2, 
                               
                               number = numberFOLD,
                               
                               repeats = numberFoldRepeats,
                               
                               verbose = FALSE, 
                               
                               classProbs = TRUE,
                               
                               summaryFunction = twoClassSummary,
                               
                               search = searchMethod,
                               
                               savePredictions = "final",
                               
                               allowParallel = TRUE,
                               
                               sampling = sampling,
                               
                               seeds = seeds)
    }
  }
  
  list(learner = learner,
       fittingParams = fittingParams,
       resamplingParams = resamplingParams)
}
