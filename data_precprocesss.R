

data_preprocesed_DataExtract = function(combined_variables= 'clinical+pathology'){
  

  DataImport = read.csv(paste0('Preprocessed_', combined_variables, '_ModelCall.csv'), row.names = 1)
  names(DataImport) <- gsub("\\.", " ", names(DataImport))
  
  
  if((combined_variables == 'clinical+TMB') | (combined_variables == 'clinical+rna+TMB')| (combined_variables == 'clinical+ihc+TMB')| (combined_variables == 'clinical+rna+ihc+TMB')) {
    
    DataImport$`TMB 2 classes` = as.factor(DataImport$`TMB 2 classes`)
    DataImport$`TCGA mutation subtype` = as.factor(DataImport$`TCGA mutation subtype`)
    DataImport$`TMB 3 classes` = as.factor(DataImport$`TMB 3 classes`)
    
    DataImport =  subset(DataImport, select=-c(`TMB 2 classes`,
                                               `Variants inside target regions and after quality filters`,
                                               `Germline variants`,
                                               `Somatic variants`,
                                               `Non coding somatic variants`,
                                               `Synonymous somatic variants`,
                                               `Non synonymous somatic variants`,
                                               `TMB  mutations per Mb  all somatic`))

  }else{
    
    DataImport = DataImport
  }

  
  DataImport$response = factor( DataImport$response,  levels = c("R", 'NR'))
  levels(DataImport$response) <- c("R", "NR") 
  
  #DataImport$response = factor( DataImport$response,  levels = c("NR", 'R'))
  #levels(DataImport$response) <- c("NR", "R") 

    
  DataImport$Sex = as.factor(DataImport$Sex)
  DataImport$Treatment = as.factor(DataImport$Treatment)
  DataImport$Line = as.factor(DataImport$Line)
  DataImport$`Primary melanoma` = as.factor(DataImport$`Primary melanoma`)
  DataImport$`Cutaneous primary` = as.factor(DataImport$`Cutaneous primary`)
  DataImport$`Previous MAPKi` = as.factor(DataImport$`Previous MAPKi`)
  DataImport$`Baseline LDH` = as.factor(DataImport$`Baseline LDH`)
  # m stage at entry
  DataImport$`M Stage` = as.factor(DataImport$`M Stage`)
  DataImport$ECOG = as.factor(DataImport$ECOG) 
  DataImport$`Brain metastasis` = as.factor(DataImport$`Brain metastasis`)
  DataImport$`Lung metastasis` = as.factor(DataImport$`Lung metastasis`)
  DataImport$`Liver metastasis` = as.factor(DataImport$`Liver metastasis`)
  

  train_data = DataImport[DataImport$Cohort == "Discovery", ]
  test_data = DataImport[DataImport$Cohort != "Discovery", ]
  
  
  processed_data =  subset(DataImport, select=-c(Cohort))
  trainData =  subset(train_data, select=-c(Cohort))
  testData =  subset(test_data, select=-c(Cohort))
  
  
  selected_clinical_variables = NULL  
  
  return(list(processed_data=processed_data,
              testData = testData,
              trainData = trainData,
              all_Data = DataImport,  
              selected_clinical_variables=selected_clinical_variables))
}






data_preprocesed22 = function(combined_variables= 'clinical+pathology',
                             VariablesearchAlgorithm = 'backward'){

  
  # excluded rows from clinical data
  if(combined_variables == 'clinical+TMB'){
    data1 = clinicalDataImport %>% filter(Comment != "(Clinical + TMB ok)") 
  } else {   
    data1 = clinicalDataImport %>% filter((Comment != "Exclude: Adjuvant") & (Comment != "(Clinical + TMB ok)") )
  }
  
  
  # hemoglobin
  clinicalDataomit =  data1[!is.na(data1$Hemoglobin), ]
  
  # Neutro.Lympho.ratio 
  clinicalDataomit$`Neutro Lympho ratio`[ clinicalDataomit$`Neutro Lympho ratio`<0.1] = 0.1 # round ratio < 0.1 to 0.1
  clinicalDataomit$`Neutro Lympho ratio`[ clinicalDataomit$`Neutro Lympho ratio`>11] = 11   # round ratio < 0.1 to 0.1
 
 
  # new factor
  clinicalDataomit$Response_PIP = factor( clinicalDataomit$Response_PIP,  levels = c("R", 'NR'))
  levels( clinicalDataomit$Response_PIP) <- c("R", "NR") 
  
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
  
  ##########
  names(clinicalDataomit)[c(7:11, 14:24)]  = clinical_variables1a  
  
  
  clinicalDataomit$Sex = as.factor(clinicalDataomit$Sex)
  clinicalDataomit$Treatment = as.factor(clinicalDataomit$Treatment)
  clinicalDataomit$Line = as.factor(clinicalDataomit$Line)
  clinicalDataomit$`Primary melanoma` = as.factor(clinicalDataomit$`Primary melanoma`)
  clinicalDataomit$`Cutaneous primary` = as.factor(clinicalDataomit$`Cutaneous primary`)
  clinicalDataomit$`Previous MAPKi` = as.factor(clinicalDataomit$`Previous MAPKi`)
  clinicalDataomit$`Baseline LDH` = as.factor(clinicalDataomit$`Baseline LDH`)
  # m stage at entry
  clinicalDataomit$`M Stage` = as.factor(clinicalDataomit$`M Stage`)
  clinicalDataomit$ECOG = as.factor(clinicalDataomit$ECOG) 
  clinicalDataomit$`Brain metastasis` = as.factor(clinicalDataomit$`Brain metastasis`)
  clinicalDataomit$`Lung metastasis` = as.factor(clinicalDataomit$`Lung metastasis`)
  clinicalDataomit$`Liver metastasis` = as.factor(clinicalDataomit$`Liver metastasis`)
  
  
  RNA_variables_names1 = names(SelectedRnaFeatures)[-c((length(names(SelectedRnaFeatures))-2):length(names(SelectedRnaFeatures)) )]
  TMB_variables_names1 = names(TMBData)[c(8:18)]
  IHC_variables_names1 = names(SelectedIHCFeatures)[-c( (length(names(SelectedIHCFeatures)) -1):length(names(SelectedIHCFeatures)))]
  
  
  All_RNA_variables_names = names(RNaFeaturesAll)[-c((length(names(RNaFeaturesAll))-2):length(names(RNaFeaturesAll)) )]  
  All_IHC_variables_names = names(IHCFeaturesAll)[-c(1, 2)]  
  
  clinical_variables1 = clinical_variables1a 
  clinical_variables  =  clinical_variables1[-1]
  
  if (combined_variables == 'clinical'){  # clinical 
  
    all_Data = clinicalDataomit[c("response", clinical_variables, 'Cohort')]
    
    all_Data  = na.omit(all_Data)

  } else if (combined_variables == 'clinical+rna'){ # combined RNA and clinical Features
   
    if  ( VariablesearchAlgorithm == 'GA parsimony'){
      
      TrainData = combinedClinicalandRNA(RNaFeaturesAll,  clinicalDataomit)
      
      RNA_variables_names = All_RNA_variables_names
      
    }else{
      
      TrainData = combinedClinicalandRNA(SelectedRnaFeatures,  clinicalDataomit)
      
      RNA_variables_names = RNA_variables_names1
    }
  
    all_Data = TrainData[c("response", clinical_variables , RNA_variables_names, 'Cohort')]
  
    all_Data  = na.omit(all_Data)
    
  }else if (combined_variables == 'clinical+ihc'){   # combined IHC and clinical Features
    
    if  ( VariablesearchAlgorithm == 'GA parsimony'){
      
      TrainData = combinedClinicalandIHC(IHCFeaturesAll,  clinicalDataomit)
      
      IHC_variables_names = All_IHC_variables_names
      
    }else{
      
      TrainData = combinedClinicalandIHC(SelectedIHCFeatures,  clinicalDataomit)
      
      IHC_variables_names = IHC_variables_names1
    }
    
    all_Data = TrainData[c("response", clinical_variables , IHC_variables_names, 'Cohort')]
    
    all_Data  = na.omit(all_Data)
    
  } else if(combined_variables == 'clinical+TMB'){     ### comblined Clinical and TMB
    
   TrainData = combinedClinicalandTMB22(TMBData, clinicalDataomit)
    
   TMB_variables_names = TMB_variables_names1
   
   all_Data = TrainData[c("response", clinical_variables , TMB_variables_names, 'Cohort')]
   
   all_Data  = na.omit(all_Data)
   
   ## TMB categories
   all_Data$`TMB 2 classes` = as.factor(all_Data$`TMB 2 classes`)
   all_Data$`TCGA mutation subtype` = as.factor(all_Data$`TCGA mutation subtype`)
   all_Data$`TMB 3 classes` = as.factor(all_Data$`TMB 3 classes`)
   
  }else if(combined_variables == 'clinical+rna+TMB'){   ## combined RNA , clinical and TMB
    

    if  ( VariablesearchAlgorithm == 'GA parsimony'){
      
      clinicalRNAData = combinedClinicalandRNA(RNaFeaturesAll,  clinicalDataomit)
      
      RNA_variables_names = All_RNA_variables_names
      
    }else{
      
      clinicalRNAData = combinedClinicalandRNA(SelectedRnaFeatures,  clinicalDataomit)
      
      RNA_variables_names = RNA_variables_names1
    }
    
  
    ### comblined Clinical and TMB
    TrainData = combinedClinicalRNAandTMB(TMBData, clinicalRNAData)
    TMB_variables_names = TMB_variables_names1
    
    all_Data = TrainData[c("response", clinical_variables, RNA_variables_names , TMB_variables_names, 'Cohort')]

    all_Data  = na.omit(all_Data)
    
    ## TMB categories
    all_Data$`TMB 2 classes` = as.factor(all_Data$`TMB 2 classes`)
    all_Data$`TCGA mutation subtype` = as.factor(all_Data$`TCGA mutation subtype`)
    all_Data$`TMB 3 classes` = as.factor(all_Data$`TMB 3 classes`)
    
    
  }else if(combined_variables == 'clinical+ihc+TMB'){ ## combined IHC , clinical and TMB
    
    
    if  ( VariablesearchAlgorithm == 'GA parsimony'){
      
      clinicalIHCData = combinedClinicalandIHC(IHCFeaturesAll,  clinicalDataomit)
      
      IHC_variables_names = All_IHC_variables_names
      
    }else{
      
      clinicalIHCData = combinedClinicalandIHC(SelectedIHCFeatures,  clinicalDataomit)
      
      IHC_variables_names = IHC_variables_names1
    }
    
  
    ### comblined ClinicalRNA and TMB
    TrainData = combinedClinicalRNAandTMB(TMBData, clinicalIHCData)
    TMB_variables_names = TMB_variables_names1
    
    all_Data = TrainData[c("response", clinical_variables, IHC_variables_names , TMB_variables_names, 'Cohort')]
    
    all_Data  = na.omit(all_Data)
    
    ## TMB categories
    all_Data$`TMB 2 classes` = as.factor(all_Data$`TMB 2 classes`)
    all_Data$`TCGA mutation subtype` = as.factor(all_Data$`TCGA mutation subtype`)
    all_Data$`TMB 3 classes` = as.factor(all_Data$`TMB 3 classes`)
    
  } else if(combined_variables == 'clinical+rna+ihc'){ ## combined RNA, clinical and IHC
    
    
    ## combined RNA and clinical
    if  ( VariablesearchAlgorithm == 'GA parsimony'){
      
      clinicalRNAData = combinedClinicalandRNA(RNaFeaturesAll,  clinicalDataomit)
      
      RNA_variables_names = All_RNA_variables_names
      
    }else{
      
      clinicalRNAData = combinedClinicalandRNA(SelectedRnaFeatures,  clinicalDataomit)
      
      RNA_variables_names = RNA_variables_names1
    }
    
  
    ### comblined ClinicalIHC and TMB
    
    if  ( VariablesearchAlgorithm == 'GA parsimony'){
      
      TrainData =  combinedClinicalRNAandIHC(IHCFeaturesAll,   clinicalRNAData)
      
      IHC_variables_names = All_IHC_variables_names
      
    }else{
      
      TrainData =  combinedClinicalRNAandIHC(SelectedIHCFeatures,   clinicalRNAData)
      
      IHC_variables_names = IHC_variables_names1
    }
    
    all_Data = TrainData[c("response", clinical_variables, RNA_variables_names , IHC_variables_names, 'Cohort')]
    
    all_Data  = na.omit(all_Data)
    
  }else if(combined_variables == 'clinical+rna+ihc+TMB'){ ## combined RNA, clinical and IHC
    
    ## combined RNA and clinical
    if  ( VariablesearchAlgorithm == 'GA parsimony'){
      
      clinicalRNAData = combinedClinicalandRNA(RNaFeaturesAll,  clinicalDataomit)
      
      RNA_variables_names = All_RNA_variables_names
      
    }else{
      
      clinicalRNAData = combinedClinicalandRNA(SelectedRnaFeatures,  clinicalDataomit)
      
      RNA_variables_names = RNA_variables_names1
    }
    
  
    ### comblined ClinicalIHC and TMB
    if  ( VariablesearchAlgorithm == 'GA parsimony'){
      
      clinicalRNAIHCData =  combinedClinicalRNAandIHC(IHCFeaturesAll,   clinicalRNAData)
      
      IHC_variables_names = All_IHC_variables_names
      
    }else{
      
      clinicalRNAIHCData =  combinedClinicalRNAandIHC(SelectedIHCFeatures,   clinicalRNAData)
      
      IHC_variables_names = IHC_variables_names1
    }
    
    
    TrainData = combinedClinicalRNAandTMB(TMBData, clinicalRNAIHCData)
    TMB_variables_names = TMB_variables_names1
    
    all_Data = TrainData[c("response", clinical_variables, RNA_variables_names , IHC_variables_names, TMB_variables_names, 'Cohort')]
    
    all_Data  = na.omit(all_Data)
    
    
    ## TMB categories
    all_Data$`TMB 2 classes` = as.factor(all_Data$`TMB 2 classes`)
    all_Data$`TCGA mutation subtype` = as.factor(all_Data$`TCGA mutation subtype`)
    all_Data$`TMB 3 classes` = as.factor(all_Data$`TMB 3 classes`)
  }
  

  # set.seed(1234)
  #all_data2 = subset(all_Data, select=-c(Cohort))
  #inTraining <- createDataPartition(all_Data$response, p=0.70, list=FALSE)
  #trainData <- all_data2[ inTraining,]
  #testData  <- all_data2[-inTraining,]
  
  all_Data = na.omit(all_Data)
  
  # #train_data = all_Data[all_Data$Cohort == "Discovery", ]
  # #test_data = all_Data[all_Data$Cohort != "Discovery", ]
  # 
  # 
  # #processed_data =  subset(all_Data, select=-c(Cohort))
  # #trainData =  subset(train_data, select=-c(Cohort))
  # #testData =  subset(test_data, select=-c(Cohort))
  # 
  # 
  # selected_clinical_variables = NULL  
  # 
  # return(list(processed_data=processed_data,
  #             testData = testData,
  #             trainData = trainData,
  #             all_Data = DataImport,  
  #             selected_clinical_variables=selected_clinical_variables))
  
  return(all_Data)
}






# test_clinicaldata_func222 = function(train_clinicaldata_func_object, test_data){
#   
#   test_data[test_data==""] <- NA
#  
#   names(test_data) = c("response",
#                        "Sex",
#                        "Age",
#                        "Treatment" ,
#                        "Line" ,
#                        "Primary melanoma" ,
#                        "Cutaneous primary" ,
#                        "Previous MAPKi",
#                        "Baseline LDH" ,
#                        "M Stage" ,
#                        "ECOG",
#                        "Hemoglobin",
#                        "Neutro Lympho ratio",
#                        "Brain metastasis",
#                        "Lung metastasis",
#                        "Liver metastasis")
#   
#   test_data$Sex = as.factor(test_data$Sex)
#   test_data$Treatment = as.factor(test_data$Treatment)
#   test_data$Line = as.factor(test_data$Line)
#   test_data$`Primary melanoma` = as.factor(test_data$`Primary melanoma`)
#   test_data$`Cutaneous primary` = as.factor(test_data$`Cutaneous primary`)
#   test_data$`Previous MAPKi` = as.factor(test_data$`Previous MAPKi`)
#   test_data$`Baseline LDH` = as.factor(test_data$`Baseline LDH`)
#   # m stage at entry
#   test_data$`M Stage` = as.factor(test_data$`M Stage`)
#   test_data$ECOG = as.factor(test_data$ECOG) 
#   test_data$`Baseline LDH` = as.factor(test_data$`Baseline LDH`)
#   test_data$`Lung metastasis` = as.factor(test_data$`Lung metastasis`)
#   test_data$`Liver metastasis` = as.factor(test_data$`Liver metastasis`)
#   
#  test_data$`Neutro Lympho ratio`[test_data$`Neutro Lympho ratio`<0.1] = 0.1 # round ratio < 0.1 to 0.1
#   test_data$`Neutro Lympho ratio`[test_data$`Neutro Lympho ratio`>11] = 11   # round ratio < 0.1 to 0.1
#   
#   
#   
#   test_data = test_data[names(train_clinicaldata_func_object)]
#   
#   test_data$response = factor(test_data$response, levels = c("Responder", 'Non-responder'))
#   levels(test_data$response) <- c("R", "NR")  # NR: Non responder, R: Responder
#   
#   
#   test_data = na.omit(test_data)
#   
#   return(test_data)
# }
# 
# 
# 
# 
# 
# 
# 
# 
# test_clinicaldata_func222PIPtest = function(train_clinicaldata_func_object, 
#                                             
#                                             test_data, treatment = 'IPI+PD1'){
#   
#   test_data[test_data==""] <- NA
#   
#   test_data =  subset(test_data, select=-c(Neutro.Lympho.ratio, Hemoglobin))
#   
#   names(test_data) = c("Sex",
#                        "Age",
#                        "Line" ,
#                        "Primary melanoma" ,
#                        "Cutaneous primary" ,
#                        "Previous MAPKi",
#                        "Baseline LDH" ,
#                        "M Stage" ,
#                        "ECOG",
#                        #"Hemoglobin",
#                        #"Neutro Lympho ratio",
#                        "Brain metastasis",
#                        "Lung metastasis",
#                        "Liver metastasis")
#   
#   
#   
#   
#   
#   test_data$Sex = as.factor(test_data$Sex)
#   test_data$Line = as.factor(test_data$Line)
#   test_data$`Primary melanoma` = as.factor(test_data$`Primary melanoma`)
#   test_data$`Cutaneous primary` = as.factor(test_data$`Cutaneous primary`)
#   test_data$`Previous MAPKi` = as.factor(test_data$`Previous MAPKi`)
#   test_data$`Baseline LDH` = as.factor(test_data$`Baseline LDH`)
#   # m stage at entry
#   test_data$`M Stage` = as.factor(test_data$`M Stage`)
#   test_data$ECOG = as.factor(test_data$ECOG) 
#   test_data$`Baseline LDH` = as.factor(test_data$`Baseline LDH`)
#   test_data$`Lung metastasis` = as.factor(test_data$`Lung metastasis`)
#   test_data$`Liver metastasis` = as.factor(test_data$`Liver metastasis`)
#   
#   #test_data$`Neutro Lympho ratio`[test_data$`Neutro Lympho ratio`<0.1] = 0.1 # round ratio < 0.1 to 0.1
#   #test_data$`Neutro Lympho ratio`[test_data$`Neutro Lympho ratio`>11] = 11   # round ratio < 0.1 to 0.1
#   
#   
#   
#   test_data$Treatment = as.factor(rep(treatment, length=length(test_data$Sex)))
#   
#   #test_data = test_data[names(train_clinicaldata_func_object)]
#   
#   test_data = na.omit(test_data)
#   
#   return(test_data)
# }
# 
# 
# 
# 
# 
# 
# 
# 
# prospectiveCombinedClincalRNA = function(train_clinicaldata_func_object, combined_variables){
# 
#   prospective_clinical = read.csv('data preprocessed/validation_clinicalCombined.csv')
#   names(prospective_clinical) <- gsub("\\.", " ", names(prospective_clinical))
#   prospective_clinical =  prospective_clinical[!duplicated(prospective_clinical$PIN),]
#   names(prospective_clinical )[1] = "Pathology ID"
#  
#   
#   ################### RNA prospective data
#   prospective_rna = read.csv('data preprocessed/New RNA training sample.csv')
#   names(prospective_rna) <- gsub("\\.", " ", names(prospective_rna))
#   #prospective_rnaProspective = dplyr::filter(prospective_rna, `Reasons for extraction` == 'PIP Prospective' | `Reasons for extraction` == 'PIP validation')
#   
#   # retain one melpins out of duplicates
#   #RNaduplicateRemoved = prospective_rnaProspective[!duplicated(prospective_rnaProspective$`Path ID`),]
#   
#   RNaduplicateRemoved = prospective_rna
#   names(RNaduplicateRemoved)[1] = "PIN"
#   names(RNaduplicateRemoved)[2] = "Pathology ID"
#   ############################################
#   
#   
#   ############################# TMB prospective data
#   TMBData = read.csv('DataTMB.csv', header = TRUE)
#   names(TMBData) <- gsub("\\.", " ", names(TMBData))
#   TMBData$`Pathology ID` = gsub("\\_", "-", TMBData$`Pathology ID`)
#   
#   
#   TMBDataProspective = TMBData %>% dplyr::filter(Cohort == "PIP_Prospective" | Cohort ==  'PIP Validation cohort')
#   TMBDataProspectiveNonNA =  TMBDataProspective[!is.na( TMBDataProspective$`Pathology ID`), ]
#   
#   # retain one melpins out of duplicates
#   TMBduplicateRemoved = TMBDataProspectiveNonNA[!duplicated(TMBDataProspectiveNonNA$`Pathology ID`),]
#     #################################
#   
#   
#   if(combined_variables == 'clinical'){
#     
#     dataS <-     prospective_clinical
#     data1_NAomit =  dataS[!is.na(dataS$Response_PIP), ]
#     
#     row.names( data1_NAomit) =  data1_NAomit$`Pathology ID`
#     
#     
#   }else if(combined_variables == 'clinical+rna'){
#     
#     dataS <-     prospective_clinical %>% right_join(RNaduplicateRemoved, by=c("Pathology ID"))
#     data1_NAomit =  dataS[!is.na(dataS$Response_PIP), ]
#     
#     row.names( data1_NAomit) =  data1_NAomit$`Pathology ID`
#     
#     
#   }else if(combined_variables == 'clinical+TMB'){
#     
#     dataS <-  prospective_clinical %>% right_join(TMBduplicateRemoved , by=c("Pathology ID"))
#     data1_NAomit =  dataS[!is.na( dataS$Response_PIP), ]
#     data1_NAomit = data1_NAomit[!is.na(data1_NAomit$`TMB  mutations per Mb  all somatic`), ]
#     
#     ## TMB categories
#     data1_NAomit$`TMB 2 classes` = as.factor(data1_NAomit$`TMB 2 classes`)
#     data1_NAomit$`TCGA mutation subtype` = as.factor(data1_NAomit$`TCGA mutation subtype`)
#     data1_NAomit$`TMB 3 classes` = as.factor(data1_NAomit$`TMB 3 classes`)
#     
#     row.names( data1_NAomit) =  data1_NAomit$`Pathology ID`
#     
#   }else if(combined_variables == 'clinical+rna+TMB'){
#     
#     clinicalRNAdata <-     prospective_clinical %>% right_join(RNaduplicateRemoved, by=c("Pathology ID"))
#     clinicalRNAdataomit =  clinicalRNAdata[!is.na(clinicalRNAdata$Response_PIP), ]
#     
#     
#     dataS <-  clinicalRNAdataomit %>% right_join(TMBduplicateRemoved , by=c("Pathology ID"))
#     data1_NAomit =  dataS[!is.na( dataS$Response_PIP), ]
#     data1_NAomit = data1_NAomit[!is.na(data1_NAomit$`TMB  mutations per Mb  all somatic`), ]
#     
#     
#     ## TMB categories
#     data1_NAomit$`TMB 2 classes` = as.factor(data1_NAomit$`TMB 2 classes`)
#     data1_NAomit$`TCGA mutation subtype` = as.factor(data1_NAomit$`TCGA mutation subtype`)
#     data1_NAomit$`TMB 3 classes` = as.factor(data1_NAomit$`TMB 3 classes`)
#    
#     row.names( data1_NAomit) =  data1_NAomit$`Pathology ID` 
#   }
# 
#   
#   
#   # new factor
#   data1_NAomit$Response_PIP = factor( data1_NAomit$Response_PIP,  levels = c("Responder", 'Non-responder'))
#   levels( data1_NAomit$Response_PIP) <- c("R", "NR")  # NR: Non responder, R: Responder
#   
#   
#   names(data1_NAomit)[c(3:18)]  = c("response",
#                                     "Sex",
#                                     "Age",
#                                     "Treatment" ,
#                                     "Line" ,
#                                     "Primary melanoma" ,
#                                     "Cutaneous primary" ,
#                                     "Previous MAPKi",
#                                     "Baseline LDH" ,
#                                     "M Stage" ,
#                                     "ECOG",
#                                     "Hemoglobin",
#                                     "Neutro Lympho ratio",
#                                     "Brain metastasis",
#                                     "Lung metastasis",
#                                     "Liver metastasis")
#   
#   data1_NAomit$Sex = as.factor(data1_NAomit$Sex)
#   data1_NAomit$Treatment = as.factor(data1_NAomit$Treatment)
#   data1_NAomit$Line = as.factor(data1_NAomit$Line)
#   data1_NAomit$`Primary melanoma` = as.factor(data1_NAomit$`Primary melanoma`)
#   data1_NAomit$`Cutaneous primary` = as.factor(data1_NAomit$`Cutaneous primary`)
#   data1_NAomit$`Previous MAPKi` = as.factor(data1_NAomit$`Previous MAPKi`)
#   data1_NAomit$`Baseline LDH` = as.factor(data1_NAomit$`Baseline LDH`)
#   # m stage at entry
#   data1_NAomit$`M Stage` = as.factor(data1_NAomit$`M Stage`)
#   data1_NAomit$ECOG = as.factor(data1_NAomit$ECOG) 
#   data1_NAomit$`Brain metastasis` = as.factor(data1_NAomit$`Brain metastasis`)
#   data1_NAomit$`Lung metastasis` = as.factor(data1_NAomit$`Lung metastasis`)
#   data1_NAomit$`Liver metastasis` = as.factor(data1_NAomit$`Liver metastasis`)
#   
# 
#   data1_NAomit$`Neutro Lympho ratio`[ data1_NAomit$`Neutro Lympho ratio`<0.1] = 0.1 # round ratio < 0.1 to 0.1
#   data1_NAomit$`Neutro Lympho ratio`[ data1_NAomit$`Neutro Lympho ratio`>11] = 11   # round ratio < 0.1 to 0.1
#   
#   names(data1_NAomit) <- gsub("\\.", " ", names(data1_NAomit))
#   
#   test_data = data1_NAomit[names(train_clinicaldata_func_object)]
# 
#   return(test_data)
# }
# 





# data_preprocesedTrain = function(train_clinicaldata_func_object,
#                               combined_variables= 'clinical+pathology'){
# 
# 
# 
#   ############################################# clinical dataset
# 
#   clinicalDataImport <- read.csv("DataclinicalPlusBloodCounts - discoveryOnly.csv")
#   names(clinicalDataImport) <- gsub("\\.", " ", names(clinicalDataImport))
# 
#   clinicalDataImport = clinicalDataImport[!duplicated(clinicalDataImport$PIN),]
# 
# 
#   ##############################################
# 
# 
#   ############################## read RNA-seq data
# 
#   RNAdata <- read.csv("DataRNA - DiscoveryOnly.csv",  header = TRUE)
# 
#   ## RNA selected Features
#   SelectedRnaFeaturesResultsList = corThreshold(RNAdata,
#                                                 clinicalDataImport,
#                                                 kbestthreshol=1,
#                                                 th.corr = 1,
#                                                 omics_variables = 'rna')
# 
#   # selected RNA Features
#   SelectedRnaFeatures = SelectedRnaFeaturesResultsList$selectedFeaturesDFFinal
#   ############################
# 
# 
#   ############################# TMB data
#   TMBData = read.csv('DataTMB.csv', header = TRUE)
#   names(TMBData) <- gsub("\\.", " ", names(TMBData))
#   TMBData$`Pathology ID`  <- gsub("\\_", "-", TMBData$`Pathology ID`)
# 
#   TMBDataDiscovery = TMBData %>% filter(Cohort == "1_PIP_Retrospective")
#   #################################
# 
# 
# 
# 
#   ############################# IHC data
#   IHCData = read.csv('DataIHC.csv', header = TRUE)
# 
#   # -ve : Negative
#   # +ve : Positive
#   # - :
#   # %: Percentage
# 
#   names(IHCData) <- gsub("\\.", "_", names(IHCData))
# 
#   IHCData$Sample_ID = gsub("\\.tif_object_Data_output", "", IHCData$Sample_ID)
#   names(IHCData)[1] = "Pathology Reference No"
#   IHCData = IHCData[!duplicated(IHCData$`Pathology Reference No`),]
# 
# 
#   ## RNA selected Features
#   SelectedIHCFeaturesResultsList = corThreshold(IHCData,
#                                                 clinicalDataImport,
#                                                 kbestthreshol=0.001,
#                                                 th.corr = 0.85,
#                                                 omics_variables = 'ihc')
# 
#   # selected RNA Features
#   SelectedIHCFeatures = SelectedIHCFeaturesResultsList$selectedFeaturesDFFinal
#   ############################
# 
# 
# 
# 
# 
# 
# 
#   # excluded rows from clinical data
#   data1 = clinicalDataImport[(clinicalDataImport$Comment != 'Exclude') ,]
# 
# 
#   # hemoglobin
#   ## clinicalDataomit =  data1[!is.na(data1$Hemoglobin), ]
# 
# 
#   clinicalDataomit =  subset(data1, select=-c(`Neutro Lympho ratio`, Hemoglobin))
# 
#   # new factor
#   clinicalDataomit$Response_PIP = factor( clinicalDataomit$Response_PIP,  levels = c("Responder", 'Non-responder'))
#   levels( clinicalDataomit$Response_PIP) <- c("R", "NR")  # NR: Non responder, R: Responder
# 
# 
#   # Neutro.Lympho.ratio
#   ##clinicalDataomit$`Neutro Lympho ratio`[ clinicalDataomit$`Neutro Lympho ratio`<0.1] = 0.1 # round ratio < 0.1 to 0.1
#   ##clinicalDataomit$`Neutro Lympho ratio`[ clinicalDataomit$`Neutro Lympho ratio`>11] = 11   # round ratio < 0.1 to 0.1
# 
#   names(clinicalDataomit)[c(6:10, 13:21)]  = c("response",
#                                                "Sex",
#                                                "Age",
#                                                "Treatment" ,
#                                                "Line" ,
#                                                "Primary melanoma" ,
#                                                "Cutaneous primary" ,
#                                                "Previous MAPKi",
#                                                "Baseline LDH" ,
#                                                "M Stage" ,
#                                                "ECOG",
#                                                #"Hemoglobin",
#                                                #"Neutro Lympho ratio",
#                                                "Brain metastasis",
#                                                "Lung metastasis",
#                                                "Liver metastasis")
# 
#   clinicalDataomit$Sex = as.factor(clinicalDataomit$Sex)
#   clinicalDataomit$Treatment = as.factor(clinicalDataomit$Treatment)
#   clinicalDataomit$Line = as.factor(clinicalDataomit$Line)
#   clinicalDataomit$`Primary melanoma` = as.factor(clinicalDataomit$`Primary melanoma`)
#   clinicalDataomit$`Cutaneous primary` = as.factor(clinicalDataomit$`Cutaneous primary`)
#   clinicalDataomit$`Previous MAPKi` = as.factor(clinicalDataomit$`Previous MAPKi`)
#   clinicalDataomit$`Baseline LDH` = as.factor(clinicalDataomit$`Baseline LDH`)
#   # m stage at entry
#   clinicalDataomit$`M Stage` = as.factor(clinicalDataomit$`M Stage`)
#   clinicalDataomit$ECOG = as.factor(clinicalDataomit$ECOG)
#   clinicalDataomit$`Brain metastasis` = as.factor(clinicalDataomit$`Brain metastasis`)
#   clinicalDataomit$`Lung metastasis` = as.factor(clinicalDataomit$`Lung metastasis`)
#   clinicalDataomit$`Liver metastasis` = as.factor(clinicalDataomit$`Liver metastasis`)
# 
# 
#   if (combined_variables == 'clinical'){  # clinical
# 
#     TrainData = clinicalDataomit
# 
#     clinical_variables  =  names(TrainData)[c(7:10, 13:21)]
# 
#     processed_data = clinicalDataomit[c("response", clinical_variables)]
# 
#     processed_data = na.omit(processed_data)
# 
#   } else if (combined_variables == 'clinical+rna'){ # combined RNA and clinical Features
# 
# 
#     TrainData = combinedClinicalandRNA(SelectedRnaFeatures,  clinicalDataomit)
# 
#     clinical_variables  =  names(TrainData)[c(5:8, 11:19)]
#     RNA_variables_names = names(SelectedRnaFeatures)[-c(length(names(SelectedRnaFeatures)):(length(names(SelectedRnaFeatures))-1) )]
# 
#     processed_data = TrainData[c("response", clinical_variables , RNA_variables_names)]
# 
#     processed_data = na.omit(processed_data)
# 
#   }else if (combined_variables == 'clinical+ihc'){   # combined IHC and clinical Features
# 
#     TrainData = combinedClinicalandIHC(SelectedIHCFeatures,  clinicalDataomit)
# 
#     clinical_variables  =  names(TrainData)[c(6:9, 12:20)]
#     IHC_variables_names = names(SelectedIHCFeatures)[-c( (length(names(SelectedIHCFeatures)) -1):length(names(SelectedIHCFeatures)))]
# 
#     processed_data = TrainData[c("response", clinical_variables , IHC_variables_names)]
# 
#     processed_data = na.omit(processed_data)
# 
# 
#   } else if(combined_variables == 'clinical+TMB'){     ### comblined Clinical and TMB
# 
#     TrainData = combinedClinicalandTMB22(TMBDataDiscovery, clinicalDataomit)
# 
#     clinical_variables  =  names(TrainData)[c(6:9, 12:20)]
#     TMB_variables_names = names(TMBDataDiscovery)[c(5, 6, 14)]
# 
#     ## TMB categories
#     TrainData$`TMB 2 classes` = as.factor(TrainData$`TMB 2 classes`)
#     TrainData$`TCGA mutation subtype` = as.factor(TrainData$`TCGA mutation subtype`)
#     TrainData$`TMB 3 classes` = as.factor(TrainData$`TMB 3 classes`)
# 
#     processed_data = TrainData[c("response", clinical_variables , TMB_variables_names)]
#     processed_data = na.omit(processed_data)
# 
#   }else if(combined_variables == 'clinical+rna+TMB'){   ## combined RNA , clinical and TMB
# 
#     ## combined RNA and clinical
#     clinicalRNAData = combinedClinicalandRNA(SelectedRnaFeatures,  clinicalDataomit)
# 
#     ### comblined Clinical and TMB
#     TrainData = combinedClinicalRNAandTMB(TMBDataDiscovery, clinicalRNAData)
# 
#     clinical_variables  =  names(TrainData)[c(6:9, 12:20)]
#     RNA_variables_names =  names(SelectedRnaFeatures)[-c( (length(names(SelectedRnaFeatures))-1):length(names(SelectedRnaFeatures)))]
#     TMB_variables_names = names(TMBDataDiscovery)[c(5, 6, 14)]
# 
#     ## TMB categories
#     TrainData$`TMB 2 classes` = as.factor(TrainData$`TMB 2 classes`)
#     TrainData$`TCGA mutation subtype` = as.factor(TrainData$`TCGA mutation subtype`)
#     TrainData$`TMB 3 classes` = as.factor(TrainData$`TMB 3 classes`)
# 
#     processed_data = TrainData[c("response", clinical_variables, RNA_variables_names , TMB_variables_names)]
#     processed_data = na.omit(processed_data)
# 
#   }else if(combined_variables == 'clinical+ihc+TMB'){ ## combined IHC , clinical and TMB
# 
#     clinicalIHCData =  combinedClinicalandIHC(SelectedIHCFeatures,  clinicalDataomit)
# 
#     ### comblined ClinicalRNA and TMB
#     TrainData = combinedClinicalRNAandTMB(TMBDataDiscovery, clinicalIHCData)
# 
#     clinical_variables  =  names(TrainData)[c(6:9, 12:20)]
#     IHC_variables_names = names(SelectedIHCFeatures)[-c( (length(names(SelectedIHCFeatures)) -1):length(names(SelectedIHCFeatures)))]
#     TMB_variables_names = names(TMBDataDiscovery)[c(5, 6, 14)]
# 
#     ## TMB categories
#     TrainData$`TMB 2 classes` = as.factor(TrainData$`TMB 2 classes`)
#     TrainData$`TCGA mutation subtype` = as.factor(TrainData$`TCGA mutation subtype`)
#     TrainData$`TMB 3 classes` = as.factor(TrainData$`TMB 3 classes`)
# 
#     processed_data = TrainData[c("response", clinical_variables, IHC_variables_names , TMB_variables_names)]
# 
#     processed_data = na.omit(processed_data)
# 
#   } else if(combined_variables == 'clinical+rna+ihc'){ ## combined RNA, clinical and IHC
# 
#     ## combined RNA and clinical
#     clinicalRNAData = combinedClinicalandRNA(SelectedRnaFeatures,  clinicalDataomit)
# 
#     ### comblined ClinicalIHC and TMB
#     TrainData =  combinedClinicalRNAandIHC(SelectedIHCFeatures,   clinicalRNAData)
# 
#     clinical_variables  =  names(TrainData)[c(6:9, 12:20)]
#     RNA_variables_names = names(SelectedRnaFeatures)[-c( (length(names(SelectedRnaFeatures))-1):length(names(SelectedRnaFeatures)))]
#     IHC_variables_names = names(SelectedIHCFeatures)[-c( (length(names(SelectedIHCFeatures)) -1):length(names(SelectedIHCFeatures)))]
# 
#     processed_data = TrainData[c("response", clinical_variables, RNA_variables_names , IHC_variables_names)]
# 
#     processed_data = na.omit(processed_data)
# 
#   }else if(combined_variables == 'clinical+rna+ihc+TMB'){ ## combined RNA, clinical and IHC
# 
#     ## combined RNA and clinical
#     clinicalRNAData = combinedClinicalandRNA(SelectedRnaFeatures,  clinicalDataomit)
# 
#     ### comblined ClinicalIHC and TMB
#     clinicalRNAIHCData =  combinedClinicalRNAandIHC(SelectedIHCFeatures,   clinicalRNAData)
# 
# 
#     TrainData = combinedClinicalRNAandTMB(TMBDataDiscovery, clinicalRNAIHCData)
# 
#     clinical_variables  =  names(TrainData)[c(6:9, 12:20)]
#     RNA_variables_names =  names(SelectedRnaFeatures)[-c( (length(names(SelectedRnaFeatures))-1):length(names(SelectedRnaFeatures)))]
#     IHC_variables_names = names(SelectedIHCFeatures)[-c( (length(names(SelectedIHCFeatures)) -1):length(names(SelectedIHCFeatures)))]
#     TMB_variables_names = names(TMBDataDiscovery)[c(5, 6, 14)]
# 
#     ## TMB categories
#     TrainData$`TMB 2 classes` = as.factor(TrainData$`TMB 2 classes`)
#     TrainData$`TCGA mutation subtype` = as.factor(TrainData$`TCGA mutation subtype`)
#     TrainData$`TMB 3 classes` = as.factor(TrainData$`TMB 3 classes`)
# 
# 
#     processed_data = TrainData[c("response", clinical_variables, RNA_variables_names , IHC_variables_names, TMB_variables_names )]
# 
#     processed_data = na.omit(processed_data)
#   }
# 
# 
#   processed_data = na.omit(processed_data)
#   selected_clinical_variables = NULL
# 
#   processed_data = processed_data[names(train_clinicaldata_func_object)]
# 
#   return(processed_data)
# }



