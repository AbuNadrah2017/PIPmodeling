

traindata_scaling_func = function(data, scale_method,  criterion = "AICc"){
  
  if(scale_method == 'zero'){
    
    pP <- preProcess(data, method = c('center', 'scale'))
    
    scaled_data <- predict(pP , data)
    
  }else if(scale_method == 'range'){
    
    pP <- preProcess(data, method = c('range'))
    
    scaled_data <- predict( pP , data)
    
    
  }else if(scale_method == 'cdf'){
    
    pP <- preprocess_CDF(data,  criterion = criterion)
    
    scaled_data <- pP$train_data
    
  }else if(scale_method == 'none'){
    
    pP <- NULL
    
    scaled_data <- data
  }
  
  return(list(train_scaling_objects=pP, scaled_data=scaled_data, scale_method= scale_method))
}


testdata_scaling_func = function(train_scaling_objects, 
                                 
                                 test_data, 
                                 
                                 scale_method){
  
  
  if(scale_method == 'zero'){
    
    pP <- predict(train_scaling_objects, test_data) 
    
    
  }else if(scale_method == 'range'){
    
    pP <- predict(train_scaling_objects, test_data) 
    
  }else if(scale_method == 'cdf'){
    
    pP <- predict_CDF(train_scaling_objects, test_data)
  
  }else if(scale_method == 'none'){
    
    pP <- test_data
  }
  
  return(pP)
}





# 
# ## test and train dataset
# data <- read.csv("combined_dataframe_with_progression_and_molecular_dataset.csv", row.names=1)
# #test_data = read.csv('Validation Set PIP_Clinical Data_final.csv')[, -1]
# 
# 
# # train and test datasets processed
# train_data = train_clinicaldata_func(data)
# #test_data = test_clinicaldata_func(train_data, test_data)
# 
# # train and test scaling
# # train
# train_results = traindata_scaling_func(train_data, 
#                                        scale_method= 'cdf', 
#                                        criterion = "AICc")
# 
# scaled_traindata = train_results$scaled_data
# 
# # test 
# #scaled_testdata= testdata_scaling_func(train_results$train_scaling_objects,
# #                                       test_data, 
# #                                       scale_method= train_results$scale_method)
# 
# #head(scaled_traindata)
# #head(scaled_testdata)