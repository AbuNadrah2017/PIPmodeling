#learner='knn'

#data = scaled_traindata
#featureSetEval = evaluator

#library(pROC)

sequentialBackwardSearch <- function(data,
                                     class,
                                     featureSetEval,
                                     samplingmethod,
                                     numberFOLD, 
                                     numberFoldRepeats,
                                     searchMethod,
                                     combined_variables){
  
  stopCriterion=-1
  stop=FALSE
  
  if (attr(featureSetEval, 'kind') == "Individual measure") {
    stop('Only feature set measures can be used');
  }
  # Extract and eliminate the class to have only the features in the variable 'features'
  column.names <- names(data) 
  class.position <- which(column.names == class) 
  features <- column.names[-class.position] 
  
  
 if((combined_variables == 'clinical+TMB') | (combined_variables == 'clinical+rna+TMB')| (combined_variables == 'clinical+ihc+TMB')| (combined_variables == 'clinical+rna+ihc+TMB')) {
    
    fixed_variables = c("Treatment", 'ECOG' , "TMB 3 classes")
    
  }else{
    
    fixed_variables = c("Treatment", 'ECOG')
  }
  
  
  feat.sub <-  as.vector(features)[!(as.vector(features) %in% fixed_variables) ]#as.vector(features)  # here
  
  excluded.features <- NULL
  
  if (stopCriterion != -1) {
    maxIterations <- min(length(feat.sub)-1, stopCriterion)
  } else {
    maxIterations <- length(feat.sub)-1
  }
  
  # Check for maximization-minimization
  metricTarget <- attr(featureSetEval,'target')
  if(metricTarget=="maximize"){
    max <- TRUE
  }else if(metricTarget=="minimize"){
    max <- FALSE
  }else{ # Metric is not specified
    # Wrapper methods use by default RMSE for regression and Acuraccy for classification (in filter methods the metric is always specified)
    max <- ifelse(is.factor(data[,class]), TRUE, FALSE)
  }
  
  for (i in seq(maxIterations)) {
    
    #print("####################")
    #print(i)
    
    best.value <- featureSetEval(data, class, c(feat.sub, fixed_variables))  #### here
    
    best_tuned_result = best.value$best_metric_result
   
  
    caretparametersInternal = Nestedmodel_FinalgridsParam(data=data,
                                                          learner=best.value$model$method,
                                                          metric=best.value$model$metric,
                                                          samplingmethod  = samplingmethod,
                                                          numberFOLD = numberFOLD , 
                                                          numberFoldRepeats = numberFoldRepeats,
                                                          searchMethod=searchMethod, 
                                                          finalTune =  best.value$best_tune)
    
    
      
      fittingParamsInternal = caretparametersInternal$fittingParams
      resamplingParamsInternal = caretparametersInternal$resamplingParams

    
    evaluatorInternal <- wrapperEvaluator(learner=best.value$model$method,
                                  resamplingParams=resamplingParamsInternal,
                                  fittingParams=fittingParamsInternal)
    
    ########################################################
    
    best.feat <- NULL
    best.feat.value <- NULL
    
    # Step 1 (Exclusion): Eliminate a feature in each step, if with it, we can get a better evaluation (if not, end of the algorithm)
    for (i in seq(along = feat.sub)) {
      feat <- feat.sub[i]
      
      start.time <- Sys.time()
      
    #print(feat)
      
      feat.prueba <- feat.sub
      feat.prueba <- feat.prueba[feat.prueba != feat]
      
      
      value <- evaluatorInternal(data, class, c(feat.prueba, fixed_variables))  ### this is used in the corrected
      inner_best_result = value$best_metric_result
      
      # Find the feature that removing it, we can get a better evaluation
      if(max){ # Classification -> maximize
        if (is.null(best.feat.value) || best.feat.value < inner_best_result) {
          best.feat.value <- inner_best_result
          best.feat <- feat
        }
      }else{ # Regression -> minimize
        if (is.null(best.feat.value) || best.feat.value > inner_best_result) {
          best.feat.value <- inner_best_result
          best.feat <- feat
        }
      }
      
      end.time <- Sys.time()
      difference <- difftime( end.time, start.time, units='mins')
      #print(difference)
    }
    
    # If removing it we can get a better of even evaluation, we remove it
    if(max){ # Classification -> maximize
      if (is.null(best_tuned_result) || best_tuned_result <= best.feat.value) {
        best_tuned_result <- best.feat.value
        best.value <- value
        
        # Remove the selected feature
        feat.sub <- feat.sub[feat.sub != best.feat]
      } else if (stop) {
        break
      }
    }else{ # Regression -> minimize
      if (is.null(best_tuned_result) || best_tuned_result >= best.feat.value) {
        best_tuned_result <- best.feat.value
        best.value <- value
        
        # Remove the selected feature
        feat.sub <- feat.sub[feat.sub != best.feat]
      } else if (stop) {
        break
      }
    }    
    
  }
  
  ## confusion matrix
  set.seed(1234)
  filt_wrapper_xtrain = data[, c(class, best.value$featuresnames)]

  predicted = predict(best.value$model, newdata = filt_wrapper_xtrain )
  predictions = factor(predicted, levels = levels(as.factor(filt_wrapper_xtrain [,class])))
  confusionmatrix = confusionMatrix(as.factor(predictions), as.factor(filt_wrapper_xtrain[,class]))
  

  # List with results
  res <- list(NULL)
  best.set.aux <- matrix(rep(0,(ncol(data)-1)), ncol=(ncol(data)-1), byrow=FALSE, dimnames=list(c(),column.names[-class.position]))
  best.set.aux[which(column.names[-class.position]%in%feat.sub)] <- 1
  
  
  res[[1]] <- best.set.aux
  res[[2]] <-  best.value$featuresnames
  res[[3]] <- best_tuned_result
  res[[4]] <-  best.value$best_tune
  res[[5]] <-  best.value$model
  res[[6]] <-  confusionmatrix
  
  ## variables importance
  #res[[7]] <- roc_nb
  
  names(res) <- c("bestFeaturesID", 
                  'bestFeaturesnames',
                  "bestFitness", 
                  'best parameter',
                  'best model',
                  'best confusion matrix'#,
                  #'best auc value'
  ) 
  
  return(list(best_result = res,
              all_models_list = NULL, 
              all_model_table_summary =  NULL))
}


















# SFBS
sequentialFloatingBackwardSelectionSearch <- function(data,
                                                      class,
                                                      featureSetEval,
                                                      samplingmethod,
                                                      numberFOLD, 
                                                      numberFoldRepeats,
                                                      searchMethod,
                                                      combined_variables) {
  
  if (attr(featureSetEval, 'kind') == "Individual measure") {
    stop('Only feature set measures can be used');
  }
  
  # Extract and eliminate the class to have only the features in the variable 'features'
  column.names <- names(data) 
  class.position <- which(column.names == class) 
  features <- column.names[-class.position] 
  
  
  if((combined_variables == 'clinical+TMB') | (combined_variables == 'clinical+rna+TMB')| (combined_variables == 'clinical+ihc+TMB')| (combined_variables == 'clinical+rna+ihc+TMB')) {
    
    fixed_variables = c("Treatment", "TMB 3 classes", "TCGA mutation subtype")
    
  }else{
    
    fixed_variables = c("Treatment")#, 'ECOG'
  }
  
  
  feat.sub <-  as.vector(features)[!(as.vector(features) %in% fixed_variables) ]#as.vector(features)  # here
  
  excluded.features <- c()
  
  # Check for maximization-minimization
  metricTarget <- attr(featureSetEval,'target')
  if(metricTarget=="maximize"){
    max <- TRUE
  }else if(metricTarget=="minimize"){
    max <- FALSE
  }else{ # Metric is not specified
    # Wrapper methods use by default RMSE for regression and Acuraccy for classification (in filter methods the metric is always specified)
    max <- ifelse(is.factor(data[,class]), TRUE, FALSE)
  }
  
  for (i in seq(length(features) - 1 - length(fixed_variables))) {
    
    #print("####################")
    #print(i)
    
    best.value <- featureSetEval(data, class, c(feat.sub, fixed_variables))
    best_tuned_result = best.value$best_metric_result
    
    
    # Step 1 (Exclusion): Eliminate a feature in each step, if with it, we can get a better evaluation (if not, end of the algorithm)
    best.feat <- NULL
    best.feat.value <- NULL
  
     ########################## This is added to the corrected but not in the original 
     caretparametersInternal = Nestedmodel_FinalgridsParam(data=data,
                                                          learner=best.value$model$method,
                                                          metric=best.value$model$metric,
                                                          samplingmethod  = samplingmethod,
                                                          numberFOLD = numberFOLD , 
                                                          numberFoldRepeats = numberFoldRepeats,
                                                          searchMethod=searchMethod, 
                                                          finalTune =  best.value$best_tune)
    
      
    fittingParamsInternal = caretparametersInternal$fittingParams
    resamplingParamsInternal = caretparametersInternal$resamplingParams
    
    
    evaluatorInternal <- wrapperEvaluator(learner=best.value$model$method,
                                          resamplingParams=resamplingParamsInternal,
                                          fittingParams=fittingParamsInternal,
                                          modelCrossValType = modelCrossValType)
    ########################################################
    
    
    
    for (i in seq(along = feat.sub)) {
      feat <- feat.sub[i]
      #print(feat)
      
      feat.prueba <- feat.sub
      feat.prueba <- feat.prueba[feat.prueba != feat]
      
      value <- evaluatorInternal(data, class, c(feat.prueba, fixed_variables))  ###
      inner_best_result = value$best_metric_result
      
      # Find the feature that removing it, we can get a better evaluation
      if(max){ # Classification -> maximize
        if (is.null(best.feat.value) || best.feat.value < inner_best_result) {
          best.feat.value <- inner_best_result
          best.feat <- feat
        }
      }else{ # Regression -> minimize
        if (is.null(best.feat.value) || best.feat.value > inner_best_result) {
          best.feat.value <- inner_best_result
          best.feat <- feat
        }
      }
      
    }
    
    
    # If removing it we can get a better of even evaluation, we remove it
    if(max){ # Classification -> maximize
      best <- best_tuned_result <= best.feat.value
    }else{ # Regression -> minimize
      best <- best_tuned_result >= best.feat.value
    }
    
    if (is.null(best_tuned_result) || best) {
      best_tuned_result <- best.feat.value
      
      # Remove the selected feature
      feat.sub <- feat.sub[feat.sub != best.feat]
      
      # Save the excluded feature for the future
      excluded.features[length(excluded.features) + 1] <- best.feat
      
      
      # Step 2: Conditional inclusion.
      crit.func.max <- evaluatorInternal(data, class, c(feat.sub, fixed_variables))   
      crit.func.max.value = crit.func.max$best_metric_result
      continue <- TRUE
      
      #print("####################")
      #print(crit.func.max.value)
      
      
      # We can include 1 or more features in each step
      while (continue == TRUE) {
        worst.feat.value <- FALSE
        
        #print("####################")
        #print('inner forward')
        
        #print(excluded.features)
        
        # See if including a feature of the removed, we can get a better evaluation
        for (i in seq(along = excluded.features)) {
          
          feat <- excluded.features[i] 
          
          
          feat.prueba <- feat.sub
          feat.prueba[length(feat.prueba) + 1] <- feat
          
          #print(as.vector(feat.prueba))
          
          crit.func.eval <- evaluatorInternal(data, class,c(feat.prueba, fixed_variables))  
          crit.func.eval.value  <- crit.func.eval$best_metric_result
          
          
          #print(crit.func.eval.value)
          #print("####################")
          
          
          if(max){ # Classification -> maximize
            best.critic <- crit.func.eval.value > crit.func.max.value
          }else{ # Regression -> minimize
            best.critic <- crit.func.eval.value < crit.func.max.value
          }
          
          
          if (best.critic) {
            worst.feat <- feat
            crit.func.max.value <- crit.func.eval.value
            worst.feat.value <- TRUE
            
            # Do not include the feature that was just removed
            if(worst.feat == best.feat) 
              worst.feat.value <- FALSE
          }
        }
        # Include the feature in the result set of features
        if (worst.feat.value == TRUE) {
          feat.sub[length(feat.sub) + 1] <- worst.feat
          excluded.features <- excluded.features[excluded.features != worst.feat]  
        } else {
          continue <- FALSE
        }
      }
    }
  }
  
  
  # List with results
  feat.sub = c(feat.sub, fixed_variables)
  
  res <- list(NULL)
  best.set.aux <- matrix(rep(0,(ncol(data)-1)), ncol=(ncol(data)-1), byrow=FALSE, dimnames=list(c(),column.names[-class.position]))
  best.set.aux[which(column.names[-class.position]%in%feat.sub)] <- 1
  res[[1]] <- best.set.aux
  
  ## best model
  best.model <- evaluatorInternal(data, class, feat.sub) 
  
  ## confusion matrix
  set.seed(1234)
  filt_wrapper_xtrain = data[, c(class, best.value$featuresnames)]
  
  predicted = predict(best.value$model, newdata = filt_wrapper_xtrain )
  predictions = factor(predicted, levels = levels(as.factor(filt_wrapper_xtrain [,class])))
  confusionmatrix = confusionMatrix(as.factor(predictions), as.factor(filt_wrapper_xtrain[,class]))
  
  
  
  res[[2]] <-  best.model$featuresnames
  res[[3]] <- best_tuned_result
  res[[4]] <-  best.model$best_tune
  res[[5]] <-  best.model$model
  res[[6]] <-  confusionmatrix
  
  
  names(res) <- c("bestFeaturesID", 
                  'bestFeaturesnames',
                  "bestFitness", 
                  'best parameter',
                  'best model',
                  'best confusion matrix'#,
                  #'best auc value'
  ) 
  
  return(list(best_result = res,
              all_models_list = NULL, 
              all_model_table_summary =  NULL))
  
}









# exhaustive search function
exhaustiveSearch <- function(data, class, featureSetEval) {
  
  
  if (attr(featureSetEval, 'kind') == "Individual measure") {
    stop('Only feature set measures can be used');
  }
  # Extract and eliminate the class to have only the features in the variable 'features'
  column.names <- names(data)
  class.position <- which(column.names == class)
  features <- column.names[-class.position]
  
  # Check for maximization-minimization
  metricTarget <- attr(featureSetEval,'target')
  if(metricTarget=="maximize"){
    max <- TRUE
  }else if(metricTarget=="minimize"){
    max <- FALSE
  }else{ # Metric is not specified
    # Wrapper methods use by default RMSE for regression and Acuraccy for classification (in filter methods the metric is always specified)
    max <- ifelse(is.factor(data[,class]), TRUE, FALSE)
  }
  
  
  # In best.set, we store the features that are part of the solution
  best.set <- NULL
  best.value <- NULL
  best.model <- NULL
  best.parameter <- NULL
  best.features_names <- NULL
  confusionmatrix <- NULL
  roc_nb <- NULL
  
  # Queue the root of the tree
  queue <- NULL
  queue[[length(queue)+1]] <- list(list(), as.list(features))
  
  
  
  valuess2 = NULL
  table_summary2 = NULL
  
  # Visite each feature
  while (length(queue) > 0) {
    
    # Pop (introduce) an unvisited node 
    node <- queue[[1]] 
    trunk <- node[[1]] 
    branches <- node[[2]]
    queue[[1]] <- NULL
    
    # visit each branch of the current node
    values  = NULL
    table_summary = NULL
    
    for (i in seq(along=branches)) {
      set <- c(trunk, branches[[i]])
      
      #print(set)
      
      # Evaluate and check if better
      value <- featureSetEval(data, class, unlist(set))
      
      tune_parameter = value$best_tune
      tune_model = value$model
      tuned_best_result = value$best_metric_result
      features_name = value$featuresnames
      
      
      
      #print(c(unlist(set), unlist(tune_parameter), unlist(tuned_best_result)))
      
      
      ## predictions
      if ((tune_model$method == 'glmnet') & (length(unlist(set)) == 1)){
        train_data <- subset(data, select = c(unlist(set),class))
        train_data['ones'] = rep(1, nrow(data))
        
      }else{
        train_data <- subset(data, select = c(unlist(set),class))
      }
      
      ## confusion matrix
      set.seed(1234)
      predicted = predict(tune_model, newdata = train_data)
      predictions = factor(predicted, levels = levels(as.factor(train_data[,class])))
      confusionmatrix = confusionMatrix(as.factor(predictions), as.factor(train_data[,class]))
      
      ## predictions and roc curve
      predictions_nb <- predict(tune_model, newdata = train_data, type="prob")[,confusionmatrix$positive]
      roc_nb <- pROC::auc(pROC::roc(train_data[,class], predictions_nb, quiet = TRUE))
      
      
      value = list(value$featuresnames, value$model, value$best_tune, value$best_metric_result,   confusionmatrix, roc_nb)
      valuess = list(values, value)
      values = valuess
      
      
      table_summ = c(paste(unlist(features_name), collapse = ",") , 
                     
                     paste(unlist(tune_parameter), collapse = ","),
                     
                     round(tuned_best_result,3) ,
                     
                     confusionmatrix$overall[1], 
                     
                     round(confusionmatrix$overall[2], 3),  
                     
                     round(as.numeric(roc_nb), 3))
      
      
      table_summary_in  = rbind(table_summary, table_summ )
      table_summary = table_summary_in
      
      
      # Store the new feature if is better
      if(max){ # Classification -> Maximize
        if (is.null(best.value) || tuned_best_result > best.value) {
          best.value <- tuned_best_result
          best.set <- set
          best.model <-  tune_model
          best.parameter <-  tune_parameter
          best.features_names <- features_name
          
          best.confusion_matrix <- confusionmatrix
          best.auc_value <- roc_nb
        }
      }else{ # Regression -> Minimize
        if (is.null(best.value) || tuned_best_result < best.value) {
          best.value <- tuned_best_result
          best.set <- set
          best.model <-  tune_model
          best.parameter <-  tune_parameter
          best.features_names <- features_name
          
          best.confusion_matrix <- confusionmatrix
          best.auc_value <- roc_nb
        }
      }
      
      # Generate branch nodes if there are remaining features to combine. In breadth
      n <- length(branches)
      if (i < n) {
        queue[[length(queue)+1]] <- list(set, branches[(i+1):n])
      }
    }
    
    valuess22 = list(valuess2, valuess) 
    valuess2 = valuess22
    
    table_summary2_out  = rbind(table_summary2, table_summary_in)
    table_summary2 = table_summary2_out
    
    
    
  }
  
  # List with results
  res <- list(NULL)
  
  best.set.aux <- matrix(rep(0,length(features)), ncol=length(features), byrow=FALSE, dimnames=list(c(),features))
  
  pos <- match(unlist(best.set),features)
  
  best.set.aux[pos] <- 1
  
  res[[1]] <- best.set.aux
  res[[2]]  <- best.features_names
  res[[3]] <- best.value
  
  res[[4]] <-  best.parameter
  res[[5]] <-  best.model
  
  res[[6]] <-  best.confusion_matrix
  
  ## variables importance
  res[[7]] <- best.auc_value
  
  names(res) <- c("bestFeaturesID", 
                  'bestFeaturesnames',
                  "bestFitness", 
                  'best parameter',
                  'best model',
                  'best confusion matrix',
                  'best auc value') 
  
  table_summary2 = as.data.frame( table_summary2)
  row.names( table_summary2) = NULL
  colnames(table_summary2)  = c('features', 
                                paste( names(tune_parameter), collapse = ","),
                                'metric result',
                                'accuracy',
                                'kappa',
                                'auc value')
  
  
  return(list(best_result = res,
              all_models_list = valuess22, 
              all_model_table_summary =  table_summary2))
  
}


