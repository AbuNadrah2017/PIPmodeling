

wrapperEvaluator <- function(learner, 
                             resamplingParams=list(), 
                             fittingParams=list()) {
  
  wrapperEvaluatorFunction <- function(original_data, 
                                       
                                       class,
                                       
                                       features) {
    
    # Check for empty set of features
    if( length(features) == 0 || features[1] == ''){
      stop('An empty set of features has been received for evaluation')
    }
    # Check for missing data
    if( any( apply(original_data, 2, function(x) { any(is.na(x)) } ) ) ){
      stop('Feature selection cannot be performed with missing values. Try to impute them previously with the preProcces function of the caret package')
    }
    # Obtain only the desired columns
    
    if ((learner == 'glmnet') & (length(features) == 1)){
    
      train_data <- subset(original_data, select = c(features,class))
      train_data['ones'] = rep(1, nrow(train_data))

     }else{
      train_data <- subset(original_data, select = c(features,class))
    }
    
       
    # trainControl
    ctrl <- do.call(caret::trainControl,resamplingParams)
    
    set.seed(1234)
    
      # train
      modelFit <- do.call(caret::train, append(list(
        form = as.formula(paste( class, ".", sep = "~")),
        data = train_data,
        method = learner,
        trControl = ctrl),
        fittingParams))
  
    
    # Row number of the best parameters tuned (these parameters have achieved the best result)
    
    best_tune = modelFit$bestTune
    
    rowBestTune <- as.numeric(rownames(modelFit$bestTune))
    
    
    # Best result for the selected metric
    best_metric_result = modelFit$results[rowBestTune,modelFit$metric]
  
  
    result = list(best_tune = best_tune, 
                  
                  best_metric_result = best_metric_result, 
                  
                  model=modelFit, 
                  
                  featuresnames= features)
  
    result
  }
  
  
  # Determine according to the selected metric whether maximization or minimization is to be done
  if (is.null(fittingParams$metric)) { # Metric not specified in the parameters
    
    target <- "unspecified"
    
  } else {# Metric specified in the parameters
    
    if (fittingParams$metric %in% c("Accuracy","Kappa","ROC","Sens","Spec","Rsquared","F","Precision","Recall")) {
      target <- "maximize"
    } else if (fittingParams$metric %in% c("RMSE","MAE","logLoss")) {
      target <- "minimize"
    } else {
      stop("Metric not supported")
    }
    
  }
  
  # Added as an attribute
  attr(wrapperEvaluatorFunction, 'target') <- target
  attr(wrapperEvaluatorFunction, 'shortName') <- paste(learner, "wrapper")
  attr(wrapperEvaluatorFunction, 'name') <- paste(learner, "Wrapper")
  attr(wrapperEvaluatorFunction, 'kind') <- "Set measure"
  attr(wrapperEvaluatorFunction, 'needsDataToBeDiscrete') <- FALSE
  attr(wrapperEvaluatorFunction,'needsDataToBeContinuous') <- FALSE
  
  return( wrapperEvaluatorFunction )
}

#https://stackoverflow.com/questions/32099991/seed-object-for-reproducible-results-in-parallel-operation-in-caret

