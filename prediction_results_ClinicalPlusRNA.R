
setwd("C:/Users/nurudeen.adegoke/OneDrive - Melanoma Institute Australia/PIP project/PIP model in R/Multi omics code")

path = getwd()

library(caret)
library(pROC)
library(gamlss)
##library(doMC)
library(xgboost)
library(plyr)
library(e1071)
library(class)
##library(arm)
library(ROSE)


library(dplyr)

source('prediction_funtion.R')
source('bootstrapp_odd_ratios_estimates.R')
source('Calibration_plots.R')
source('Kaplan_Meier_Plot.R')

#source('data merged function.R')
#source('data_precprocesss.R')
#source('cdf_scaling.R')
#source('data_scaling.R')




result = function(){
  
  
  
  if(modelCrossValType == 'NSRKF'){
    
    filen=paste0(modelCrossValType,  "_", 
                 outer_method,  "_", 
                 n_outer_folds,  "_", 
                 outer_repeatize, '_',
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
                 finalCV, '_',
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
  
  
  setwd(file.path(path, 'results', filen))  #
  
  #setwd(file.path(path, 'results', combined_variables, filen))
  
  file1= paste0("model_result.rds")
  trainedModel = readRDS(file1)
  
  
  combined_variables2 = c('clinical+ihc',
                         'clinical+ihc+TMB',
                         'clinical+rna+ihc',
                         'clinical+rna+ihc+TMB')
  
  if (any(combined_variables2 == combined_variables)){
    
    
    test_data = NULL
    
    test_data = NULL   ## remove NA
    
    #train_data = data_preprocesedTrain(trainedModel$train_data, combined_variables)
    
    scaled_traindata = trainedModel$train_dataScaledObject$scaled_data
    
    #scaled_traindata = testdata_scaling_func(trainedModel$train_dataScaledObject$train_scaling_objects,
    #                                         train_data, 
    #                                         scale_method= trainedModel$train_dataScaledObject$scale_method)
    
    # # # test I
    scaled_testdata = NULL
    
    test1 = 'no'
    
  }else{

    # ## test data
    # train_data = trainedModel$DataPreproceed$trainData
    # 
    # scaled_traindata = testdata_scaling_func(trainedModel$train_dataScaledObject$train_scaling_objects,
    #                                          train_data,
    #                                          scale_method= trainedModel$train_dataScaledObject$scale_method)
    # 
    # 
    # 
    # # # # test I
    # test_data = trainedModel$DataPreproceed$testData
    # test_data = na.omit(test_data)   ## remove NA
    # 
    # scaled_testdata = testdata_scaling_func(trainedModel$train_dataScaledObject$train_scaling_objects,
    #                                          test_data,
    #                                          scale_method= trainedModel$train_dataScaledObject$scale_method)
    
  
    scaled_traindata = trainedModel$scaled_traindata
    scaled_testdata = trainedModel$scaled_testdata
    #scaled_testdata$ECOG= rep('>=1',length(scaled_testdata$ECOG))
    
    test1 = 'yes'
    
  }
  
  
  
  prediction =  predict_func(trainObject= trainedModel,
                             traindata= scaled_traindata, 
                             testdata=  scaled_testdata, 
                             testdata2= NULL, 
                             test1 = test1,
                             test2 = 'no',
                             class = class,
                             clincal_model_type = clincal_model_type)
  
  
  print(prediction$train_confusionmatrixSUMMARY)
  print(prediction$test_confusionmatrixSUMMARY)
  
  #plot(trainedModel$modelResults$best_model$`best model`)
  #trainedModel$modelResults$best_model$`best model`
  
  cal_trainData = cbind(scaled_traindata, prediction$train_predictionsProb['R'])
  cal_testData = cbind(scaled_testdata, prediction$test_predictionsProb['R'])
  
  cal_trainData$y = ifelse(cal_trainData$response == 'R', 1, 0)
  cal_testData$y = ifelse(cal_testData$response == 'R', 1, 0)
  
  names(cal_trainData)[dim(cal_trainData)[2]-1] = 'pred'
  names(cal_testData)[dim(cal_testData)[2]-1] = 'pred'
  
  cal_trainData['train_predicted'] =  prediction$train_predictedProb_names$train_predicted
  cal_testData['test_predicted'] =  prediction$test_predictedProb_names$test_predicted  
  
  
  print(c(min(prediction$train_predictionsProb), max(prediction$train_predictionsProb)))
  print(c(min(prediction$test_predictionsProb), max(prediction$test_predictionsProb)))
  
  print(prediction$test_predictionsProb2)
  print(prediction$test_predicted2)
  print(prediction$test_predictedProb_names)
  #par(mfrow=c(1,2))
  
  # roc curve
  roc_cruve(prediction, combined_variables=combined_variables, test1 = test1, test2='no', metric=metric)
  
  p1 <- recordPlot()
  
  if(clincal_model_type == 'newClinical'){
    
    ## variables importance
    V <- varImp(trainedModel$modelResults$best_model$`best model`, scale = FALSE)
    #plot(test_class_imp_rec)
    V2 = as.data.frame(V$importance[1])
    
    df2 <- do.call(cbind.data.frame, V2)
    colnames(df2) = 'R'
    rownames(df2) = rownames(V2)
    
    
    p2=ggplot2::ggplot(df2) +
      geom_col(aes(x=reorder(rownames(df2),R), y=R), col = "blue", show.legend = F) +
      coord_flip() +
      scale_fill_grey() +
      #scale_fill_manual(values = c("blue"))+
      xlab('Features')+
      ylab('Importance')+
      
      theme_bw() +
      
      #theme_light() +
      
      theme(plot.title = element_text(colour = "black", face = "bold", size = 15, hjust = 0.5),
                            legend.text = element_text(colour = "black", face = "bold", size = 12),
                            legend.title = element_text(colour = "black", face = "bold", size = 12),
                            axis.text = element_text(colour = "black", face='bold', size = 12),
                            axis.text.y = element_text(colour = "black", face='bold', size = 12),
                            axis.text.x = element_text(colour = "black", face='bold', size = 12),
                            axis.title.y = element_text(colour = "black", face='bold', size = 12),
                            axis.title.x = element_text(colour = "black", face='bold', size = 12))  # customize plot and risk table with a theme

    
    #dev.off()
    library(gridGraphics)
    library(cowplot)
    a=plot_grid(p1, p2,
                nrow = 2,
                #rel_widths = c(1.2, 1),
                labels = 'AUTO',
                hjust = 0, vjust = 1)
    
    p1=plot_grid(p1,
               # labels = 'AUTO',
                hjust = 0, vjust = 1)
    
    p2=plot_grid(p2,
                #labels = 'AUTO',
                hjust = 0, vjust = 1)
    p11 = p1
    p22 = p2
  }else{
    
    a= p1
    p2 = NULL
  }
  
  return(list(plt =a,
              ROC_curvePlot = p1,
              feaselectPlot = p2,
              p11 = p1,
              p22 = p2,
              prediction_r = prediction,
              trainedModel_r = trainedModel,
              cal_trainData = cal_trainData,
              cal_testData = cal_testData,
              filen = filen)) 
}



learner = c('bayesglm', 
            'glmnet', 
            'knn',
            'rf',
            'lda',
            'nnet',
            'svmRadial',
            'svmLinear')[2]

metric= c("ROC",
          "Spec",
          "Kappa")[2]


scale_method = c('zero', 'none')[1]

class = 'response'

samplingmethod = c('LOOCV', 
                   'stratifiedRepeatedCV',
                   "boot632",
                   'upSamplestratifiedRepeatedCV',
                   "upSampleboot632")[2]


searchMethod = 'grid'


combined_variables = c('clinical', 
                       'clinical+rna',
                       'clinical+TMB',
                       'clinical+ihc',
                       'clinical+rna+TMB',
                       'clinical+ihc+TMB',
                       'clinical+rna+ihc',
                       'clinical+rna+ihc+TMB')[4]


VariablesearchAlgorithm = c('backward', 
                            'SFBS',
                            'GA parsimony')[1]

clincal_model_type = c('newClinical')[1]


if(samplingmethod == 'LOOCV'){
  
  numberFOLD = 0; numberFoldRepeats = 0
  
}else if((samplingmethod == "boot632") | (samplingmethod == "upSampleboot632")| (samplingmethod == "downsampledboot632")){
  
  numberFOLD= 10 ;  numberFoldRepeats =  500
  
}else{
  numberFOLD = 20 ;  numberFoldRepeats = 20
}


modelCrossValType = c('NSRKF', 'SRKF')[1]
outer_method = c( 'SKFold',  'cv', 'LOOCV')[1]  ### outer fould validations


if(outer_method == 'SKFold'){
  n_outer_folds =  20
  outer_repeatize = 5
  
}else if(outer_method == 'cv'){
  n_outer_folds =  20
  outer_repeatize = 0
  
}else {
  
  n_outer_folds =  NULL
  outer_repeatize = NULL
}



finalCV = FALSE  # see nestedCV_copied.R for usage
nested_features_selection_search = TRUE


#########################################
DataImport = read.csv(paste0('Preprocessed_', 
                             combined_variables, '_ModelCall.csv'), 
                      row.names = 1)
names(DataImport) <- gsub("\\.", " ", names(DataImport))

train_data = DataImport[DataImport$Cohort == "Discovery", ]
test_data = DataImport[DataImport$Cohort != "Discovery", ]

table(train_data$response)
table(test_data$response)

sum(table(train_data$response))
sum(table(test_data$response))
###################################


aa=result()
prediction = aa$prediction_r
print(prediction$train_confusionmatrixSUMMARY)
print(prediction$test_confusionmatrixSUMMARY)
b = aa$trainedModel_r$modelResults$best_model$`best model`
coef(b$finalModel, b$bestTune$lambda)

#dev.off()

setwd(file.path(path, 'slide','2023', 'plots'))  

file12= paste0(aa$filen, '_RP_', ".png")
png(file = file12, width = 800, height = 600)
aa$ROC_curvePlot
dev.off()


file1= paste0(aa$filen, '_FP_', ".png")
png(file = file1, width = 800, height = 600)
aa$feaselectPlot
dev.off()


library('lares')
#mplot_cuts(aa$cal_trainData$pred, splits = 4)

plot_quantiles(aa$cal_trainData$pred, fill = "deepskyblue",  subtitle = 'Discovery Cohort', splits = 4)
plot_quantiles(aa$cal_testData$pred, fill = 'red',  subtitle = 'Validation Cohort', splits = 4)

#########
setwd(file.path(path, 'slide','2023', 'OR'))  
### odd ratios, CI and p-values
if(learner == 'glmnet'){
  odd_ratio_results = ODD_ratio_estimates(aa,
                                          learner = 'glmnet',
                                          boot.size = 1000L)

  write.csv(odd_ratio_results, paste0(aa$filen, '_OR.csv') )

}else{

  odd_ratio_results = NULL

}


setwd(file.path(path))  

response_prop_sumr  = PFS_and_response_prop_sumr(aa, 
                                                 percentiles = c(0.25, 0.75), 
                                                 Kaplan_mair_title = '',
                                                 title = 'Calibration plot')

write.csv(response_prop_sumr$res.prop.table, 'aaresponder.csv')

write.csv(response_prop_sumr$nonres.prop.table, 'aanonresponder.csv')







##########################
DataImport = read.csv(paste0('Preprocessed_', combined_variables, '_ModelCall.csv'), row.names = 1)
names(DataImport) <- gsub("\\.", " ", names(DataImport))

train_data = DataImport[DataImport$Cohort == "Discovery", ]
test_data = DataImport[DataImport$Cohort != "Discovery", ]
  
table(train_data$response)
sum(table(train_data$response))
table(test_data$response)
sum(table(test_data$response))

# cal_result_train = calibrationPlot(data = aa$cal_trainData, 
#                                    obs = "y", 
#                                    nTiles = 10,
#                                    pred = "pred",
#                                    xlab = "Predicted probability",
#                                    ylab = "Observed probability", 
#                                    title = "Calibration plot for development data",
#                                    data_summary = TRUE)
# 
# cal_result_test = calibrationPlot(data = aa$cal_testData, 
#                                    obs = "y", 
#                                    nTiles = 10,
#                                    pred = "pred",
#                                    xlab = "Predicted probability",
#                                    ylab = "Observed probability", 
#                                    title = "Calibration plot for development data",
#                                    data_summary = TRUE)
# 
# 
# cal_result_train$data_summary
# cal_result_test$data_summary
# 
# 

#https://datascienceplus.com/machine-learning-results-one-plot-to-rule-them-all/
# https://onlinelibrary.wiley.com/doi/full/10.1002/imt2.43
# https://jokergoo.github.io/ComplexHeatmap-reference/book/more-examples.html


