library(caret)
library(Rfast)
options(warn = -1)


source('data merged function.R')
source('data_precprocesss.R')
source('cdf_scaling.R')
source('data_scaling.R')
source('base_model_notNested.R')
source('base_feature_selectionCorrected_notNested.R')
source('caret_train_grids_parameters.R')
source('inner_filter.R')
source('nestedCV_copied.R')
source('trrain_model_function_Modified.R')


library(doParallel)
library(doMC)
library(parallel)
library(foreach)


start.time <- Sys.time()


#ncores <-  getOption("mc.cores", 2L)
ncores <- getOption("mc.cores", detectCores())

#registerDoMC(cores = ncores)

# combined_variables = c('clinical', 
#                        'clinical+rna',
#                        'clinical+TMB',
#                        'clinical+ihc',
#                        'clinical+rna+TMB',
#                        'clinical+ihc+TMB',
#                        'clinical+rna+ihc',
#                        'clinical+rna+ihc+TMB')[2]

VariablesearchAlgorithm = c('backward', 
                            'SFBS', 
                            'GA parsimony')[1]

clincal_model_type = c('newClinical', 'jcoClinical')[1]

class = 'response'
samplingmethod = c('LOOCV',
                    'stratifiedRepeatedCV',
                   "boot632",
                   'upSamplestratifiedRepeatedCV',
                   "upSampleboot632")[2]

searchMethod = 'grid'
path = getwd()
path1 = file.path(path, 'clinical')


modelCrossValType = c('NSRKF', 'SRKF')[2]  ###NSRKF: nested stratified k-fold; SRKF: stratified k-fold 
outer_method = c( 'SKFold',  'cv', 'LOOCV')[2]  ### outer fould validations

cv.cores = ncores

if(outer_method == 'SKFold'){
  n_outer_folds =  20
  outer_repeatize = 2
  
}else if(outer_method == 'cv'){
  n_outer_folds =  20
  outer_repeatize = 0
  
}else {
  
  n_outer_folds =  NULL
  outer_repeatize = NULL
}




clinical_variables_names = c("Sex", 
                             "Age", 
                             "Treatment", "Line", 
                             "Primary melanoma",
                             "Cutaneous primary",
                             "Previous MAPKi",
                             "Baseline LDH", 
                             "M Stage", 
                             "ECOG" ,
                             "Hemoglobin" , 
                             "Neutro Lympho ratio",
                             "Brain metastasis",
                             "Lung metastasis", 
                             "Liver metastasis")

if((combined_variables == 'clinical+TMB') | (combined_variables == 'clinical+rna+TMB')| (combined_variables == 'clinical+ihc+TMB')| (combined_variables == 'clinical+rna+ihc+TMB')) {
  
  TMB_variables_names = c('TCGA mutation subtype', 
                          'TMB 3 classes',
                          'TMB  mutations per Mb  Non synonomous')
  
}else{
  
  TMB_variables_names = NULL
}



finalCV = FALSE  # see nestedCV_copied.R for usage
nested_features_selection_search = TRUE


for(t in 1:length(combined_variables)){
  
  for (i in 1:length(metric)){
    
      for(j in 1:length(scale_method)){
        
        for(r in 1:length(samplingmethod)){
          
          if(samplingmethod == 'LOOCV'){
            
            numberFOLD = 0; numberFoldRepeats = 0
            
          }else if((samplingmethod == "boot632") | (samplingmethod == "upSampleboot632")| (samplingmethod == "downsampledboot632")){
            
            numberFOLD= 5000 ;  numberFoldRepeats =  1
            
          }else{
            
            numberFOLD = 10 ;  numberFoldRepeats = 10
            
          }
                    
          for (k in 1:length(numberFoldRepeats)){
          
          trainwrapperFunc(data=data,
                           class=class,
                           metric=metric[i],
                           learner=learner,
                           VariablesearchAlgorithm=VariablesearchAlgorithm,
                           samplingmethod=samplingmethod[r], 
                           numberFOLD=numberFOLD, 
                           numberFoldRepeats=numberFoldRepeats[k], 
                           searchMethod = searchMethod,
                           scale_method = scale_method[j],
                           combined_variables = combined_variables[t],
                           path=path,
                           clincal_model_type=clincal_model_type,
                           modelCrossValType = modelCrossValType,
                           outer_method = outer_method, 
                           n_outer_folds = n_outer_folds,
                           outer_repeatize = outer_repeatize,
                           cv.cores = cv.cores,
                           finalCV = finalCV,  # see nestedCV_copied.R for usage
                           nested_features_selection_search = nested_features_selection_search,  # if feature selection should be done in inner nested fold
                           clinical_variables_names = clinical_variables_names,
                           TMB_variables_names = TMB_variables_names,
                           path1 = path1)
        }
      }
    }
  }
}

options(warn=0)

end.time <- Sys.time()
difference <- difftime( end.time, start.time, units='mins')
print(difference)



#coef(mymodel$finalModel, mymodel$bestTune$lambda)


