

# final training model wrapper
trainwrapperFunc = function(data,
                            class,
                            metric,
                            learner,
                            VariablesearchAlgorithm,
                            samplingmethod, 
                            numberFOLD, 
                            numberFoldRepeats, 
                            searchMethod,
                            scale_method,
                            combined_variables,
                            path,
                            clincal_model_type,
                            modelCrossValType = 'NSRKF',
                            outer_method = 'SKFold', 
                            n_outer_folds = 10,
                            outer_repeatize = 10,
                            cv.cores,
                            finalCV = TRUE,
                            nested_features_selection_search = TRUE,
                            clinical_variables_names,
                            TMB_variables_names,
                            path1){ 
  

  
  clinicalModel = NULL
  
  # train dataset processed
  proprecedTrainTestData = data_preprocesed_DataExtract(combined_variables= combined_variables)
  
  selected_clinical_variables = proprecedTrainTestData$selected_clinical_variables
  
  
  ################################## train, test and all datasets
  # train data
  train_data = proprecedTrainTestData$trainData
  
  ### test data
  test_data = proprecedTrainTestData$testData
  
  # all data
  all_dataDF = proprecedTrainTestData$processed_data
  
  
  
  ################################## scalings 
  # train data scaling
  train_results = traindata_scaling_func(train_data, 
                                         scale_method= scale_method, 
                                         criterion = "AICc")
  scaled_traindata = train_results$scaled_data
  
  
  # test data scaling
  scaled_testdata1 = testdata_scaling_func(train_results$train_scaling_objects,
                                           test_data, 
                                           scale_method= train_results$scale_method)
  scaled_testdata = na.omit(scaled_testdata1)
  

  #All data scaling
  scaled_all_dataDF = testdata_scaling_func(train_results$train_scaling_objects,
                                            all_dataDF,
                                           scale_method= train_results$scale_method)
  
  ##################################
  caretparameters = model_gridsParameters(data=scaled_all_dataDF,
                                          learner=learner,
                                          metric=metric,
                                          samplingmethod = samplingmethod, 
                                          numberFOLD = numberFOLD , 
                                          numberFoldRepeats = numberFoldRepeats,
                                          searchMethod = searchMethod)
  
  
  fittingParams = caretparameters$fittingParams
  resamplingParams = caretparameters$resamplingParams
  
  evaluator <- wrapperEvaluator(learner=learner,
                                resamplingParams=resamplingParams,
                                fittingParams=fittingParams)
  

  if((modelCrossValType == 'SRKF')){
    
    if((combined_variables != 'clinical') | (combined_variables != 'clinical+TMB')){
      
      filt_xtrain = scaled_all_dataDF
      filt_xtest  = scaled_testdata
      
    }else{
      
      df_x = scaled_all_dataDF[, -which(names(scaled_all_dataDF) %in% c('response', 
                                                                        clinical_variables_names,
                                                                        TMB_variables_names))]
      
      args <- list(y = scaled_all_dataDF[,class], 
                   x = df_x,
                   p_cutoff = 0.01, 
                   rsq_cutoff = 0.85^2, 
                   type = "names")
      
      fset <- do.call(ttest_filter, args)
      
      filt_xtrain <- scaled_all_dataDF[, c(class, clinical_variables_names, fset, TMB_variables_names)]
      filt_xtest <- scaled_testdata[, c(class, clinical_variables_names, fset, TMB_variables_names), drop = FALSE]
      
    }
    
    
    ### balance out the training part
    resmaple_out <- randomsample(y = filt_xtrain[,class], 
                                 x = filt_xtrain[, -which(names(filt_xtrain) %in% c(class))]  )
    
    resmaple_y2 <- as.data.frame(resmaple_out$y)
    colnames(resmaple_y2) = class
    resmaple_x2 <- resmaple_out$x
    filt_resmaple_train_data = cbind(resmaple_y2, resmaple_x2)
    
    
    trainedModelR =  trainclassificationFunc(scaled_all_dataDF=filt_resmaple_train_data, 
                                             scaled_testdata=filt_xtest,
                                             class=class, 
                                             evaluator=evaluator, 
                                             VariablesearchAlgorithm=VariablesearchAlgorithm,
                                             learner=learner, 
                                             metric = metric,
                                             samplingmethod  = samplingmethod,
                                             numberFOLD = numberFOLD , 
                                             numberFoldRepeats = numberFoldRepeats,
                                             searchMethod=searchMethod,
                                             combined_variables = combined_variables)
    
  }else if( (modelCrossValType == 'NSRKF') ){             # & (leaner != 'lda') & (learner != 'bayesglm')
    
    trainedModelR =  nestcv_train(class = class,
                                  train_data =  scaled_all_dataDF,
                                  outer_method = outer_method, 
                                  n_outer_folds = n_outer_folds,
                                  outer_repeatize = outer_repeatize,
                                  cv.cores = cv.cores,
                                  learner = learner,
                                  fittingParams = fittingParams,
                                  resamplingParams = resamplingParams,
                                  VariablesearchAlgorithm = VariablesearchAlgorithm,
                                  combined_variables = combined_variables,
                                  clinical_variables_names = clinical_variables_names,
                                  TMB_variables_names = TMB_variables_names,
                                  finalCV = finalCV,
                                  nested_features_selection_search = nested_features_selection_search)
    
    
  } 


  if(modelCrossValType == 'NSRKF'){
    
    filen=paste0(modelCrossValType,  "_", 
                 outer_method,  "_", 
                 n_outer_folds,  "_",  
                 outer_repeatize, "_",
                 VariablesearchAlgorithm, "_",
                 clincal_model_type, "_", 
                 combined_variables, "_",  
                 learner, "_",  
                 metric, "_", 
                 samplingmethod ,"_", 
                 numberFOLD ,"_",  
                 numberFoldRepeats, "_", 
                 searchMethod, "_", 
                 scale_method,  "_", 
                 finalCV, "_",
                 nested_features_selection_search) 
  }else{
    
    filen=paste0(modelCrossValType,  "_", 
                 VariablesearchAlgorithm, "_",
                 clincal_model_type, "_", 
                 combined_variables, "_",  
                 learner, "_",  
                 metric, "_", 
                 samplingmethod ,"_", 
                 numberFOLD ,"_",  
                 numberFoldRepeats, "_", 
                 searchMethod, "_", 
                 scale_method) 
  }
 
  
  dir.create(file.path(path, 'results', combined_variables, filen), showWarnings = FALSE)
  
  setwd(file.path(path, 'results', combined_variables, filen))
  
  file1= paste0("model_result.rds")
  
  saveRDS(list(train_data = train_data, 
               test_data =  test_data,
               scaled_traindata = scaled_traindata,
               scaled_testdata = scaled_testdata,
               DataPreproceed = proprecedTrainTestData,
               selected_clinical_variables = selected_clinical_variables,
               train_dataScaledObject =  train_results,  
               modelResults = trainedModelR,
               clinicalModel = clinicalModel), file1)
  
}








# final training model
trainclassificationFunc <- function(scaled_all_dataDF,
                                    scaled_testdata,
                                    class, 
                                    evaluator, 
                                    VariablesearchAlgorithm = 'exhaustive', 
                                    learner, 
                                    metric,
                                    samplingmethod,
                                    numberFOLD, 
                                    numberFoldRepeats,
                                    searchMethod,
                                    combined_variables){
  
    # Check that the name of the dependent variable exists
  if(!class%in%colnames(scaled_all_dataDF)){
    stop('The dependent variable does not exist in the dataset.')
  }
  # Check the number of columns
  if((ncol(scaled_all_dataDF)-1)<2){
    stop('The feature selection process requires more than 1 feature')
  }
  # Check for missing data
  if( any( apply(scaled_all_dataDF, 2, function(x) { any(is.na(x)) } ) ) ){
    stop('Feature selection cannot be performed with missing values. Try to impute them previously with the preProcces function of the caret package')
  }
  
  
  # The search algorithm is called with the evaluation method
  if(VariablesearchAlgorithm == 'exhaustive'){
    
    t <- proc.time() 
    exhaustiveSearchResults <- exhaustiveSearch(data=scaled_all_dataDF,
                                                class=class,
                                                featureSetEval = evaluator) 
    time <- proc.time()-t  
  }else if((VariablesearchAlgorithm == 'backward')| (VariablesearchAlgorithm == 'GAthenBackward'))   {
    
    t <- proc.time() 
    exhaustiveSearchResults <- sequentialBackwardSearch(data = scaled_all_dataDF, 
                                                        class=class,
                                                        featureSetEval=evaluator,
                                                        samplingmethod = samplingmethod, 
                                                        numberFOLD = numberFOLD , 
                                                        numberFoldRepeats = numberFoldRepeats,
                                                        searchMethod = searchMethod,
                                                        combined_variables = combined_variables) 
    time <- proc.time()-t  
  } else if (VariablesearchAlgorithm ==  'SFBS'){
    
    t <- proc.time() 
    
    exhaustiveSearchResults <- sequentialFloatingBackwardSelectionSearch(data=scaled_all_dataDF, 
                                              class=class,
                                              featureSetEval=evaluator,
                                              samplingmethod = samplingmethod, 
                                              numberFOLD = numberFOLD , 
                                              numberFoldRepeats = numberFoldRepeats,
                                              searchMethod = searchMethod,
                                              combined_variables = combined_variables) 
      
    time <- proc.time()-t  
  }
  
  else if (VariablesearchAlgorithm == 'GA nonparsimony'){
    
    exhaustiveSearchResults <- geneticAlgorithm(data = scaled_all_dataDF, 
                                                class=class, 
                                                featureSetEval= evaluator,
                                                combined_variables=combined_variables)
    
  }else if (VariablesearchAlgorithm == 'GA parsimony'){
    
    exhaustiveSearchResults = FinalGA_model(trainData = scaled_all_dataDF,
                                            testData = scaled_testdata,
                                            learner=learner, 
                                            metric = metric,
                                            samplingmethod  = samplingmethod,
                                            numberFOLD = numberFOLD , 
                                            numberFoldRepeats = numberFoldRepeats,
                                            combined_variables = combined_variables)
  }
  
  
  ## get mode coeficients
  
  result = list(best_model = exhaustiveSearchResults$best_result,
                
                table_all_models_summaries = exhaustiveSearchResults$all_model_table_summary,
                
                list_all_models = exhaustiveSearchResults$all_models_list )
  
  return(result)
}


