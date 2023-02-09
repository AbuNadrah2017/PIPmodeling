source('MIA model.R')

predict_func = function(trainObject, 
                        traindata,
                        testdata,
                        testdata2,
                        test1 = 'no',
                        test2 = 'no',
                        class,
                        clincal_model_type = 'newClinical'){
  
  
  # train prediction
  if(clincal_model_type == 'jcoClinical'){
    
    train_featureNames =  c("ECOG",
                            "Baseline LDH" ,
                            "Neutro Lympho ratio",
                            "Lung metastasis",
                            "Liver metastasis",
                            "Treatment" ,
                            "Line")
    
    trainModels = NULL
    trainModel_name = 'glmnet'
    
  }else if(clincal_model_type == 'newClinical'){
    
    trainObject2 = trainObject$modelResults$best_model
    trainModels = trainObject2$`best model`
    trainModel_name = trainModels$method
    train_featureNames = trainObject2$bestFeaturesnames   # train selected features
  }
  
  
  if ((trainModel_name == 'glmnet') & (length(train_featureNames) == 1)){
    train_selecteddata =  traindata[train_featureNames]
    train_selecteddata['ones'] = rep(1, nrow(train_selecteddata))
    
  }else{
    train_selecteddata =  traindata[train_featureNames]
  }
  
  
  #################################
  # train prediction
  if(clincal_model_type == 'jcoClinical'){
    miaprecition =  MIAPrediction(dataN=train_selecteddata)
    train_predicted = miaprecition$predicted
    train_predictionsProb <- miaprecition$predictedprob
    confusionmatrix = confusionMatrix(as.factor(train_predicted), as.factor(traindata[,class]))
    train_confusionmatrixSUMMARY = confusionmatrix   # train confusion matrix 
    
  }else{
    
    train_predicted = predict(trainModels, newdata = train_selecteddata)
    train_predictionsProb <- predict(trainModels,
                                     newdata = train_selecteddata, 
                                     type="prob")
    row.names(train_predictionsProb) = row.names(train_selecteddata)
    confusionmatrix = confusionMatrix(as.factor(train_predicted), as.factor(traindata[,class]))
    train_confusionmatrixSUMMARY = confusionmatrix    # train confusion matrix 
  }
  
  
  train_confusionmatrix = train_confusionmatrixSUMMARY$table
  
  train_ROC <- roc(traindata[,class],
                   train_predictionsProb[,train_confusionmatrixSUMMARY$positive], 
                   quiet = TRUE, ci=TRUE,
                   direction=">")
  
  train_auc = auc(train_ROC)# train auc
  
  
  # test prediction 1
  if(test1 == 'no'){
    
    test_selecteddata = NULL
    test_predicted = NULL
    test_predictionsProb = NULL
    test_confusionmatrixSUMMARY= NULL
    test_confusionmatrix = NULL
    test_ROC= NULL
    test_auc =  NULL
    
  }else{
  test_selecteddata = testdata[train_featureNames]
  if ((trainModel_name == 'glmnet') & (length(train_featureNames) == 1)){
    test_selecteddata =  testdata[train_featureNames]
    test_selecteddata['ones'] = rep(1, nrow(test_selecteddata))
    
  }else{
    test_selecteddata =  testdata[train_featureNames]
  }
  
  
  if(('MStage'  %in%  names(test_selecteddata))){
    test_selecteddata$MStage = as.factor(test_selecteddata$MStage)  
  }
  
  # test prediction
  if(clincal_model_type == 'jcoClinical'){
    miaprecition =  MIAPrediction(dataN=test_selecteddata)
    test_predicted = miaprecition$predicted
    test_predictionsProb <- miaprecition$predictedprob
    test_confusionmatrixSUMMARY = confusionMatrix(as.factor(test_predicted), as.factor(testdata[,class]))
    test_confusionmatrix = test_confusionmatrixSUMMARY$table
    
  }else{
    test_predicted = predict(trainModels, newdata = test_selecteddata)
    test_predictionsProb <- predict(trainModels, newdata = test_selecteddata, type="prob")
    row.names(test_predictionsProb) = row.names(test_selecteddata)
    test_confusionmatrixSUMMARY = confusionMatrix(as.factor(test_predicted), as.factor(testdata[,class]))
    test_confusionmatrix = test_confusionmatrixSUMMARY$table
  }
  
  test_ROC <- roc(testdata[,class],
                  test_predictionsProb[,test_confusionmatrixSUMMARY$positive], 
                  quiet = TRUE, ci=TRUE,
                  direction=">")
  test_auc = auc(test_ROC)
  }
  
  if(test2 == 'no'){
    
    test_selecteddata2 = NULL
    test_predicted2 = NULL
    test_predictionsProb2 = NULL
    test_confusionmatrixSUMMARY2= NULL
    test_confusionmatrix2 = NULL
    test_ROC2= NULL
    test_auc2 =  NULL
    
  }else{
    
  # test prediction 1
  test_selecteddata2 = testdata2[train_featureNames]
  if ((trainModel_name == 'glmnet') & (length(train_featureNames) == 1)){
    test_selecteddata2 =  testdata2[train_featureNames]
    test_selecteddata2['ones'] = rep(1, nrow(test_selecteddata2))
    
  }else{
    test_selecteddata2 =  testdata2[train_featureNames]
  }
  
  
  if(('MStage'  %in%  names(test_selecteddata2))){
    test_selecteddata2$MStage = as.factor(test_selecteddata2$MStage)  
  }
  
  # test prediction
  if(clincal_model_type == 'jcoClinical'){
    miaprecition2 =  MIAPrediction(dataN=test_selecteddata2)
    test_predicted2 = miaprecition2$predicted
    test_predictionsProb2 <- miaprecition2$predictedprob
    test_confusionmatrixSUMMARY2 = confusionMatrix(as.factor(test_predicted2),
                                                   as.factor(testdata2[,class]))
    test_confusionmatrix2 = test_confusionmatrixSUMMARY2$table
    
  }else{
    test_predicted2 = predict(trainModels, 
                              newdata = test_selecteddata2)
    
    test_predictionsProb2 <- predict(trainModels, 
                                     newdata = test_selecteddata2, type="prob")
    
    test_confusionmatrixSUMMARY2 = confusionMatrix(as.factor(test_predicted2),
                                                   as.factor(testdata2[,class]))
    
    test_confusionmatrix2 = test_confusionmatrixSUMMARY2$table
  }
  
  test_ROC2 <- roc(testdata2[,class],
                   test_predictionsProb2[,test_confusionmatrixSUMMARY2$positive], 
                   quiet = TRUE, ci=TRUE,
                   direction=">")
  test_auc2 = auc(test_ROC2)
  
  }
  
  
  
  # roc curve
  
  prediction_resultlist = list(
    trainModel_name = trainModel_name,
    train_selecteddata = train_selecteddata ,
    train_predicted = train_predicted,
    train_predictionsProb = train_predictionsProb,
    train_predictedProb_names = cbind(train_predictionsProb, train_predicted),
    
    train_confusionmatrixSUMMARY = train_confusionmatrixSUMMARY,
    train_confusionmatrix = train_confusionmatrix,
    train_ROC = train_ROC,
    train_auc =train_auc,
    
    test_selecteddata = test_selecteddata,
    test_predicted =test_predicted,
    test_predictionsProb =test_predictionsProb,
    test_predictedProb_names = cbind(test_predictionsProb, test_predicted),
    test_confusionmatrixSUMMARY= test_confusionmatrixSUMMARY,
    test_confusionmatrix =test_confusionmatrix,
    test_ROC=test_ROC,
    test_auc =  test_auc,
    
    test_selecteddata2 = test_selecteddata2,
    test_predicted2 =test_predicted2,
    test_predictionsProb2 =test_predictionsProb2,
    test_confusionmatrixSUMMARY2= test_confusionmatrixSUMMARY2,
    test_confusionmatrix2 =test_confusionmatrix2,
    test_ROC2=test_ROC2,
    test_auc2 =  test_auc2)
  
  return(prediction_resultlist)
}






roccurvePlot = function(predictedObjects, combined_variables, test1 = 'no', test2 = 'no', metric = 'Spec'){
  
  train_ROC = predictedObjects$train_ROC
  test_ROC = predictedObjects$test_ROC
  
  model_name = predictedObjects$trainModel_name
  
  if(model_name == 'rpart'){
    model_name = 'decision tree'
  }else if (model_name == 'glmnet'){
    model_name = 'penalised logistic regression'
  }else if (model_name == 'bayesglm'){
    model_name = 'Bayesian generalised model'
  }else if (model_name == 'lda'){
    model_name = 'linear discriminant analysis'
  }else{
    model_name
  }
  
  
  
  if((test2 == 'no') & (test1 == 'yes')){
  
    plot(train_ROC, col="green", 
         main = paste0('Multi omic: ', combined_variables, '; Classifier: ', model_name),
         legacy.axes = TRUE, 
         xlim=c(1, 0),
         #print.auc=TRUE, 
         #ci = TRUE,
         print.auc.pattern = "%.2f (%.2f-%.2f)")
    
    lines(test_ROC, col="blue")
    
    legend(x=0.6, y=022, #"bottomright",
           col=c("green", "blue"),
           legend=c("discovery", 
                    "validation"), 
           lty=1)
    
    
    legend_train <- paste0(
      "discovery: AUC = ", round(auc(train_ROC), 2), " (95% CI = ", round(ci.auc(train_ROC)[1], 2), " - ",  round(ci.auc(train_ROC)[3], 2), ")"
    )

    
    legend_test <- paste0(
      "validation: AUC = ", round(auc(test_ROC), 2), " (95% CI = ",  round(ci.auc(test_ROC)[1], 2), " - ",  round(ci.auc(test_ROC)[3], 2), ")"
    )
    
    
#    legend_train <- sprintf("discovery cohort (AUC: %.2f)", auc(train_ROC))
 #   legend_test <- sprintf("validation cohort I (AUC: %.2f)", auc(test_ROC))
  
 
    legend(x=0.5, y=0.1, #"bottomright",
           col=c("green", 
                 "blue"), lty=1,
           legend=c(legend_train, legend_test))
    
  } else if((test1 == 'no') & (test1 == 'no')){
    
    plot(train_ROC, col="green", 
         main = paste0('Multi omic: ', combined_variables, '; Classifier: ', model_name),
         legacy.axes = TRUE, 
         xlim=c(1, 0),
         print.auc.pattern = "%.2f (%.2f-%.2f)")
    
    legend(x=0.6, y=022, #"bottomright",
           col=c("green"),
           legend=c("discovery cohort"), 
           lty=1)
    
    
    
    legend_train <- paste0(
      "discovery cohort: AUC = ", round(auc(train_ROC), 2), " (95% CI = ", round(ci.auc(train_ROC)[1], 2), " - ",  round(ci.auc(train_ROC)[3], 2), ")"
    )
    
    
    
    #legend_train <- sprintf("discovery cohort (AUC: %.2f)", auc(train_ROC))
    
    legend(x=0.6, y=0.2, #"bottomright",
           col=c("green"), lty=1,
           legend=c(legend_train))
    
    
  } else{
    
    test_ROC2 = predictedObjects$test_ROC2
    plot(train_ROC, col="green", 
         main = paste0('Multi omic: ', combined_variables, '; Classifier: ', model_name),
         legacy.axes = TRUE, 
         xlim=c(1, 0),
         print.auc.pattern = "%.2f (%.2f-%.2f)")
    
    lines(test_ROC, col="blue")
    lines(test_ROC2, col="red")
    
    legend(x=0.6, y=022, #"bottomright",
           col=c("green", "blue", 'red'),
           legend=c("discovery", 
                    "validation", 
                    "validation"), 
           lty=1)
    
    
    
    legend_train <- paste0(
      "discovery cohort: AUC = ", round(auc(train_ROC), 2), " (95% CI = ", round(ci.auc(train_ROC)[1], 2), " - ",  round(ci.auc(train_ROC)[3], 2), ")"
    )
    
    
    legend_test <- paste0(
      "validation cohort I: AUC = ", round(auc(test_ROC), 2), " (95% CI = ",  round(ci.auc(test_ROC)[1], 2), " - ",  round(ci.auc(test_ROC)[3], 2), ")"
    )
    
    
    legend_test2 <- paste0(
      "validation cohort II: AUC = ", round(auc(test_ROC2), 2), " (95% CI = ",  round(ci.auc(test_ROC2)[1], 2), " - ",  round(ci.auc(test_ROC2)[3], 2), ")"
    )
    
    
    #legend_train <- sprintf("discovery cohort (AUC: %.2f)", auc(train_ROC))
    #legend_test <- sprintf("validation cohort I (AUC: %.2f)", auc(test_ROC))
    #legend_test2 <- sprintf("validation cohort II (AUC: %.2f)", auc(test_ROC2))
    
    legend(x=0.6, y=0.2, #"bottomright",
           col=c("green", 
                 "blue",
                 'red'), lty=1,
           legend=c(legend_train, legend_test, legend_test2))
    
  }
  
}


roc_cruve=function(predictedObjects, combined_variables, test1, test2, metric){roccurvePlot(predictedObjects, combined_variables, test1, test2, metric)}



TrainModelPredictionResult = function(trainObject, test_data=NULL, class, clinicalModel = TRUE){
  
  
  traindata = trainObject$train_dataScaledObject$scaled_data 
  
  if(clinicalModel == TRUE){
    ##### clinical Model Result
    
    trainModel = trainObject$clinicalModel$modelResults$best_model
    trainModel_name = trainModel$`best model`$method
    
    train_featureNames = trainObject$selected_clinical_variables
  }else{
    
    trainModel = trainObject$modelResults$best_model
    trainModel_name = trainModel$`best model`$method
    
    train_featureNames = trainModel$bestFeaturesnames   # train selected features
  }
  
  # train prediction
  
  
  if ((trainModel_name == 'glmnet') & (length(train_featureNames) == 1)){
    
    train_selecteddata =  traindata[train_featureNames]
    train_selecteddata['ones'] = rep(1, nrow(train_selecteddata))
    
  }else{
    train_selecteddata =  traindata[train_featureNames]
  }
  
  
  train_predicted = predict(trainModel$`best model`, newdata = train_selecteddata)
  train_predictionsProb <- predict(trainModel$`best model`, newdata = train_selecteddata, type="prob")
  
  predictions = factor(train_predicted, levels = levels(as.factor(traindata[,class])))
  confusionmatrix = confusionMatrix(as.factor(predictions), as.factor(traindata[,class]))
  
  
  train_confusionmatrixSUMMARY =   confusionmatrix   # train confusion matrix 
  train_confusionmatrix = train_confusionmatrixSUMMARY$table
  
  train_ROC <- roc(traindata[,class],
                   train_predictionsProb[,train_confusionmatrixSUMMARY$positive], 
                   quiet = TRUE, ci=TRUE,
                   direction=">")
  
  train_auc = auc(train_ROC)# train auc
  
  
  
  
  traindata['Predicted class'] = train_predicted 
  
  # Comp. btw Predicted and True class
  compResults= train_predicted == traindata[,class]
  traindata['Comp. btw Predicted and True class'] =  compResults
  
  
  prediction_resultlist = list(
    trainModel_name = trainModel_name,
    train_selecteddata = train_selecteddata ,
    train_predicted = train_predicted,
    train_predictionsProb = train_predictionsProb,
    train_confusionmatrixSUMMARY = train_confusionmatrixSUMMARY,
    train_confusionmatrix = train_confusionmatrix,
    train_ROC = train_ROC,
    train_auc =train_auc,
    train_featureNames=train_featureNames,
    traindata = traindata)
  
  return(prediction_resultlist)
}



CombinedModelroccurvePlot = function(clincailtrain_ROC1, 
                                     ClinicPathtrain_ROC1, 
                                     ClinicSpatPathtrain_ROC1,
                                     ClinicPathSpathtrain_ROC1 ){
  
  train_ROC1 = clincailtrain_ROC1
  train_ROC2 = ClinicPathtrain_ROC1
  train_ROC3 = ClinicSpatPathtrain_ROC1
  train_ROC4 = ClinicPathSpathtrain_ROC1 
  
  
  plot(train_ROC1, col="green", 
       #main = predictedObjects1$trainModel_name,
       legacy.axes = TRUE, 
       xlim=c(1, 0),
       print.auc.pattern = "%.2f (%.2f-%.2f)")
  
  lines(train_ROC2, col="blue")
  lines(train_ROC3, col="red")
  lines(train_ROC4, col="black")
  
  
  legend(x=0.6, y=022, #"bottomright",
         col=c("green", "blue", "red", "grey"), 
         legend=c("clinical", 
                  "Clinical & path",
                  "Clinical & Spath", 
                  "Clinical, path & Spath"), lty=1)
  
  legend_clinical <- sprintf("clinical (AUC: %s)", paste0(round(auc(train_ROC1), 2),', ', '95% CI: ', round(ci.auc(train_ROC1)[1],2), '-', round(ci.auc(train_ROC1)[3], 2)))
  legend_ClinicPathtrain <- sprintf("Clinical & path (AUC: %s)", paste0(round(auc(train_ROC2), 2),', ', '95% CI: ', round(ci.auc(train_ROC2)[1],2), '-', round(ci.auc(train_ROC2)[3], 2)))
  legend_ClinicSpatPathtrain <- sprintf("Clinical & Spath (AUC: %s)", paste0(round(auc(train_ROC3), 2),', ', '95% CI: ', round(ci.auc(train_ROC3)[1],2), '-', round(ci.auc(train_ROC3)[3], 2)))
  legend_ClinicPathSpathtrain <- sprintf("Clinical, path & Spath (AUC: %s)", paste0(round(auc(train_ROC4), 2),', ', '95% CI: ', round(ci.auc(train_ROC4)[1],2), '-', round(ci.auc(train_ROC4)[3], 2)))
  
  legend(x=0.6, y=0.3, #"bottomright",
         col=c("green", "blue", "red", "grey"), lty=1,
         legend=c(legend_clinical, legend_ClinicPathtrain,
                  legend_ClinicSpatPathtrain, 
                  legend_ClinicPathSpathtrain ))
  
}





VarImportantPlot = function(varImpDF){
  
  varImpDF = as.data.frame(varImpDF$importance[1])
  df2 <- do.call(cbind.data.frame, varImpDF)
  colnames(df2) = 'R'
  rownames(df2) = rownames(varImpDF)
  
  p2=ggplot2::ggplot(df2) +
    geom_col(aes(x=reorder(rownames(df2),R), y=R), col = "blue", show.legend = F) +
    coord_flip() +
    scale_fill_grey() +
    #scale_fill_manual(values = c("blue"))+
    xlab('Features')+
    ylab('Importance')+
    theme_bw()
  
  return(p2)
}


##############################################

predict_funcPIPTesting = function(trainObject,  traindata,  testdata, class){
  
  trainModel = trainObject$best_model
  trainModel_name = trainModel$`best model`$method
  
  # train prediction
  train_featureNames = trainModel$bestFeaturesnames   # train selected features
  
  
  if ((trainModel_name == 'glmnet') & (length(train_featureNames) == 1)){
    
    train_selecteddata =  traindata[train_featureNames]
    train_selecteddata['ones'] = rep(1, nrow(train_selecteddata))
    
  }else{
    train_selecteddata =  traindata[train_featureNames]
  }
  
  
  train_predicted = predict(trainModel$`best model`, newdata = train_selecteddata)
  train_predictionsProb <- predict(trainModel$`best model`, newdata = train_selecteddata, type="prob")
  train_confusionmatrixSUMMARY = trainModel$`best confusion matrix`   # train confusion matrix 
  train_confusionmatrix = train_confusionmatrixSUMMARY$table
  
  train_ROC <- roc(traindata[,class],
                   train_predictionsProb[,train_confusionmatrixSUMMARY$positive], 
                   quiet = TRUE, ci=TRUE,
                   direction=">")
  
  train_auc = auc(train_ROC)# train auc
  
  
  # test prediction
  test_selecteddata = testdata[train_featureNames]
  
  if ((trainModel_name == 'glmnet') & (length(train_featureNames) == 1)){
    
    test_selecteddata =  testdata[train_featureNames]
    test_selecteddata['ones'] = rep(1, nrow(test_selecteddata))
    
  }else{
    test_selecteddata =  testdata[train_featureNames]
  }
  
  
  if(('MStage'  %in%  names(test_selecteddata))){
    test_selecteddata$MStage = as.factor(test_selecteddata$MStage)  
  }
  
  
  test_predicted = predict(trainModel$`best model`, newdata = test_selecteddata)
  test_predictionsProb <- predict(trainModel$`best model`, newdata = test_selecteddata, type="prob")
  
  
  # roc curve
  
  prediction_resultlist = list(
    trainModel_name = trainModel_name,
    train_selecteddata = train_selecteddata ,
    train_predicted = train_predicted,
    train_predictionsProb = train_predictionsProb,
    train_confusionmatrixSUMMARY = train_confusionmatrixSUMMARY,
    train_confusionmatrix = train_confusionmatrix,
    train_ROC = train_ROC,
    train_auc =train_auc,
    
    test_selecteddata = test_selecteddata,
    test_predicted =test_predicted,
    test_predictionsProb =test_predictionsProb)
  
  return(prediction_resultlist)
}









# CinicalPlusRNApredict_func = function(trainObject, 
#                                       traindata,
#                                       testdata,
#                                       class,
#                                       clincal_model_type = 'newClinical'){
#   
#   
#   # train prediction
#   if(clincal_model_type == 'jcoClinical'){
#     
#     train_featureNames =  c("ECOG",
#                             "Baseline LDH" ,
#                             "Neutro Lympho ratio",
#                             "Lung metastasis",
#                             "Liver metastasis",
#                             "Treatment" ,
#                             "Line")
#     
#     trainModel = NULL
#     trainModel_name = 'glmnet'
#     
#   }else{
#     
#     trainObject = trainObject$modelResults
#     trainModel = trainObject$best_model
#     trainModel_name = trainModel$`best model`$method
#     train_featureNames = trainModel$bestFeaturesnames   # train selected features
#   }
#   
#   
#   if ((trainModel_name == 'glmnet') & (length(train_featureNames) == 1)){
#     train_selecteddata =  traindata[train_featureNames]
#     train_selecteddata['ones'] = rep(1, nrow(train_selecteddata))
#     
#   }else{
#     train_selecteddata =  traindata[train_featureNames]
#   }
#   
#   
#   #################################
#   # train prediction
#   if(clincal_model_type == 'jcoClinical'){
#     miaprecition =  MIAPrediction(dataN=train_selecteddata)
#     train_predicted = miaprecition$predicted
#     train_predictionsProb <- miaprecition$predictedprob
#     confusionmatrix = confusionMatrix(as.factor(train_predicted), as.factor(traindata[,class]))
#     train_confusionmatrixSUMMARY = confusionmatrix   # train confusion matrix 
#     
#   }else{
#     
#     train_predicted = predict(trainModel$`best model`, newdata = train_selecteddata)
#     train_predictionsProb <- predict(trainModel$`best model`,
#                                      newdata = train_selecteddata, 
#                                      type="prob")
#     confusionmatrix = confusionMatrix(as.factor(train_predicted), as.factor(traindata[,class]))
#     train_confusionmatrixSUMMARY = confusionmatrix    # train confusion matrix 
#   }
#   
#   
#   train_confusionmatrix = train_confusionmatrixSUMMARY$table
#   
#   train_ROC <- roc(traindata[,class],
#                    train_predictionsProb[,train_confusionmatrixSUMMARY$positive], 
#                    quiet = TRUE, ci=TRUE,
#                    direction=">")
#   
#   train_auc = auc(train_ROC)# train auc
#   
#   
#   # test prediction 1
#   test_selecteddata = testdata[train_featureNames]
#   if ((trainModel_name == 'glmnet') & (length(train_featureNames) == 1)){
#     test_selecteddata =  testdata[train_featureNames]
#     test_selecteddata['ones'] = rep(1, nrow(test_selecteddata))
#     
#   }else{
#     test_selecteddata =  testdata[train_featureNames]
#   }
#   
#   
#   if(('MStage'  %in%  names(test_selecteddata))){
#     test_selecteddata$MStage = as.factor(test_selecteddata$MStage)  
#   }
#   
#   # test prediction
#   if(clincal_model_type == 'jcoClinical'){
#     miaprecition =  MIAPrediction(dataN=test_selecteddata)
#     test_predicted = miaprecition$predicted
#     test_predictionsProb <- miaprecition$predictedprob
#     test_confusionmatrixSUMMARY = confusionMatrix(as.factor(test_predicted), as.factor(testdata[,class]))
#     test_confusionmatrix = test_confusionmatrixSUMMARY$table
#     
#   }else{
#     test_predicted = predict(trainModel$`best model`, newdata = test_selecteddata)
#     test_predictionsProb <- predict(trainModel$`best model`, newdata = test_selecteddata, type="prob")
#     test_confusionmatrixSUMMARY = confusionMatrix(as.factor(test_predicted), as.factor(testdata[,class]))
#     test_confusionmatrix = test_confusionmatrixSUMMARY$table
#   }
#   
#   test_ROC <- roc(testdata[,class],
#                   test_predictionsProb[,test_confusionmatrixSUMMARY$positive], 
#                   quiet = TRUE, ci=TRUE,
#                   direction=">")
#   test_auc = auc(test_ROC)
#   
#   
#   
#   # roc curve
#   
#   prediction_resultlist = list(
#     trainModel_name = trainModel_name,
#     train_selecteddata = train_selecteddata ,
#     train_predicted = train_predicted,
#     train_predictionsProb = train_predictionsProb,
#     train_confusionmatrixSUMMARY = train_confusionmatrixSUMMARY,
#     train_confusionmatrix = train_confusionmatrix,
#     train_ROC = train_ROC,
#     train_auc =train_auc,
#     
#     test_selecteddata = test_selecteddata,
#     test_predicted =test_predicted,
#     test_predictionsProb =test_predictionsProb,
#     test_confusionmatrixSUMMARY= test_confusionmatrixSUMMARY,
#     test_confusionmatrix =test_confusionmatrix,
#     test_ROC=test_ROC,
#     test_auc =  test_auc)
#   
#   return(prediction_resultlist)
# }




