get_typesMine <- function(x, coarse = TRUE) {
  if(is.null(colnames(x)))
    stop("`x` must have column names")
  if(is.matrix(x)) {
    out <- rep(class(x[1,1]), ncol(x))
  } else {
    if(is.data.frame(x)) {
      out <- unlist(lapply(x, function(x) class(x)[1]))
    }
  }
  
  if(coarse) {
    num_classes <- c("integer", "numeric", "double")
    str_classes <- c("factor", "character")
    out <- ifelse(out %in% num_classes, "numeric", out)
    out <- ifelse(out %in% str_classes, "string", out)
    out <- ifelse(out %in% c("numeric", "string"), out, "other")
  }
  names(out) <- colnames(x)
  out
}

data_preproces = function(data, combined_variables){
  
  clinical_variables1a = c("response",
                           "Sex",
                           "Age",
                           "Treatment" ,
                           "Line" ,
                           "Primary melanoma" ,
                           "Cutaneous primary" ,
                           "Previous MAPKi",
                           "Baseline LDH" ,
                           "M Stage" ,
                           "ECOG",
                           "Hemoglobin",
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
  

  
  data2 = data[,c(clinical_variables1a, TMB_variables_names )]
  
  data2$response = factor(data2$response , levels = c('R', 'NR'))
  data2$Sex = factor(data2$Sex , levels = c('Male', 'Female'))
  data2$Treatment = factor(data2$Treatment , levels = c('PD1', 'IPI+PD1'))
  data2$Line = factor(data2$Line , levels = c('1st Line', '>1st line'))
  data2$`Primary melanoma` = factor(data2$`Primary melanoma` , levels = c('Occult', 
                                                                          'Scalp.Face.Neck', 
                                                                          'other'))
  data2$`Cutaneous primary` = factor(data2$`Cutaneous primary` , levels = c('YES', 'No'))
  data2$`Previous MAPKi` = factor(data2$`Previous MAPKi` , levels = c('Yes', 'No'))
  data2$`Baseline LDH` = factor(data2$`Baseline LDH` , levels = c('Normal', 'Elevated'))
  data2$`M Stage` = factor(data2$`M Stage`, levels = c('0', '1'))
  data2$ECOG = factor(data2$ECOG, levels = c( '0', '>=1'))
  data2$`Brain metastasis` = factor(data2$`Brain metastasis`, levels = c('Yes', 'No'))
  data2$`Lung metastasis` = factor(data2$`Lung metastasis` , levels = c('Yes', 'No'))
  data2$`Liver metastasis` = factor(data2$`Liver metastasis`, levels = c('Yes', 'No'))
  
  ############ check string and numerical features
  vars = as.vector(get_typesMine(data2, coarse = TRUE))
  orig_vars <- vars
  vars <- vars %in% c("string")
  string_features = names(data2)[vars]
  numfeat =  names(data2)[!vars]
  
  data3 = data2
  for(i in string_features){
    data3[,i]  = as.factor(data3[,i])
    
  }
  
  
  
  for(i in  numfeat){
    data3[,i]  = as.numeric(unlist(data3[,i]))
  }
  
  
  return(list(data=data3,
              string_features = string_features,
              numfeat = numfeat))
  
}

################## import dataset
 combined_variables = c('clinical', 
                        'clinical+rna',
                        'clinical+TMB',
                        'clinical+ihc',
                        'clinical+rna+TMB',
                        'clinical+ihc+TMB',
                        'clinical+rna+ihc',
                        'clinical+rna+ihc+TMB')[3]

DataImport = read.csv(paste0('Preprocessed_', 
                             combined_variables, '_ModelCall.csv'),
                      row.names = 1)

names(DataImport) <- gsub("\\.", " ", names(DataImport))


train_data = DataImport[DataImport$Cohort == "Discovery", ]
#test_data = DataImport[DataImport$Cohort != "Discovery", ]



data_preproces_res = data_preproces(train_data)
data = data_preproces_res$data
string_features = data_preproces_res$string_features
numfeat  = data_preproces_res$numfeat

clinical_summary = function(data){
  
    clinclica_details = function(df_use, string_features){
    
    for (i in 1:length(string_features)){
      
      A2 = table(as.vector(unlist(df_use[string_features[i]])), df_use$response)
      
      B2 = round(A2 / colSums(A2)[col(A2)] * 100, 1)
      
      print('####################################')
      print('####################################')
      print(string_features[i])
 
      if(dim(B2)[1] == 2){
        
        print(list(A2, B2,  chisq.test(A2, correct = TRUE)))
      }else{
        
        print(list(A2, B2, chisq.test(A2, correct = FALSE)))
      }
      
      A22 = table(as.vector(unlist(df_use[string_features[i]])))
      

      print(A22)
      
      print(round(A22 /sum(A22)* 100, 1))
      
    }
    
  }
  
    
    
  num_function = function(data2){
    ### continuous
    for(i in numfeat){
      
      print('####################################')
      print('####################################')
      print(i)
      
      x1 = data2[,i][data2$response == levels(data2$response)[1]]
      x2 = data2[,i][data2$response == levels(data2$response)[2]]
      
      
      print(list(T2_range = range(x1, na.rm = TRUE),
                 T1_range = range(x2, na.rm = TRUE), 
                 T2_median = median(x1, na.rm = TRUE),
                 T1_median = median(x2, na.rm = TRUE), 
                 tst_comp = wilcox.test(x1, x2),
                 total_mdeian = median(data2[,i]),
                 toal_range = range(data2[,i], na.rm = TRUE)))
      
    }
  
  }
  
  
  result = list(categ =   clinclica_details(df_use=data, string_features),
                nume = num_function(data))
  
  return(result)
  
}


clinical_summary(data)


# library('pmsampsize')
# 
# pmsampsize(type = "b",
#            rsquared = 0.6, 
#            parameters = 8, 
#            prevalence = 0.18)





