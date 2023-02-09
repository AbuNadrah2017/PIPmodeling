
library(dplyr)


## unlist and change mode to numeric
makenumcols<-function(df){
  
  variables_names= names(df)
  for(i in 1:length(variables_names)){
    df[i] = as.numeric(as.matrix(unlist(df[i])))
  }
  return(df)
}




## This function first signature with   missing genes greater than 40%. It also check if there is duplicate melpins 
## it retains only one melpin when there are duplicates
RNA_preprocessed = function(data, cellNumber = 'all'){
  
  ## extract RNA features
  rnaDF = subset(data, select=-c(Melpin, `Path.ID`, Cohort))

  rnanames = names(rnaDF)
  
  
  # # if (cellNumber == 'all'){
  # #   
  # #   # remove RNA with more than 30 missing cells
  # #   rnanamesgreater03 = as.vector(rnanames[as.numeric(as.matrix(unname(rnaDF [1,]))) >= 0])
  # #   
  # # }else if (cellNumber == 'fraction'){
  # #   
  # #   # remove RNA with more than 30 missing cells
  # #   rnanamesgreater03 = as.vector(rnanames[as.numeric(as.matrix(unname(rnaDF [1,]))) <= 0.3])# <= 0.3
  # #   
  # # }
  # 
  # ## remove genes with missing genes greater than 30% 
  # rnaDF = rnaDF[,rnanamesgreater03] 
  
  rnaDF['PIN'] = data$Melpin
  rnaDF['Path ID'] = data$Path.ID
  rnaDF['Cohort'] = data$Cohort
  
  RNADF = rnaDF    #[-c(1:3),]
  
  #RNADF$PIN[duplicated(RNADF$PIN)]  # duplicate melpins
  
  # retain one melpins out of duplicates
  RNaduplicateRemoved = RNADF[!duplicated(RNADF$`Path ID`),]
  
  #RNaduplicateRemoved$PIN[duplicated(RNaduplicateRemoved$PIN)]  #  ## check if there is still duplicate melpins
  
  ### Turn the columns to numerical values
  RNAnumeric = makenumcols(subset(RNaduplicateRemoved, select=-c(PIN, `Path ID`,  Cohort)))
  
  RNAnumeric['PIN'] =  RNaduplicateRemoved['PIN']
  RNAnumeric['Path ID'] =  RNaduplicateRemoved['Path ID']
  RNAnumeric['Cohort'] =  RNaduplicateRemoved['Cohort']
  
  
  names(RNAnumeric) <- gsub("\\.", " ", names(RNAnumeric))
  
  return(RNAnumeric)
}





# Select Kbest using Wilcoxon test 
selectKBest <- function(RNaProcessed, clinicaData, kbestthreshol= 0.001, omics_variables = 'ihc'){
  

  
  if(omics_variables == 'rna'){
    
    ## select response class, pin, hemoglin, commnet 
    clinicaDataSubset  = clinicaData[c('PIN', 'Response_PIP')]
    
    # mergeddatasets
    mergedDataset = merge(clinicaDataSubset, RNaProcessed,  by.x = "PIN",  by.y = "PIN", all.x = TRUE)
    names(mergedDataset) <- gsub("\\.", " ", names(mergedDataset))
    
  }else if (omics_variables == 'ihc'){
    
    ## select response class, pin, hemoglin, commnet 
    clinicaDataSubset  = clinicaData[c('PIN', 'Response_PIP', "Pathology Reference No")]
    
    mergedDataset = merge(clinicaDataSubset, RNaProcessed,  by.x = "Pathology Reference No",  by.y = "Pathology Reference No", all.x = TRUE)
    
    names(mergedDataset) <- gsub("\\.", " ", names(mergedDataset))
  }
  
  
  if(omics_variables == 'rna'){
    
    mergedDatasetomit = mergedDataset[mergedDataset$Cohort == 'Discovery', ]
    
    processed_data = na.omit(subset(mergedDatasetomit, select=-c(PIN,  `Path ID`,  Cohort)))

  }else if (omics_variables == 'ihc'){
    
    mergedDatasetomit = mergedDataset
    
    processed_data = na.omit(subset(mergedDatasetomit, select=-c(PIN,  `Pathology Reference No`)))
    
  }
 
  names(processed_data)[1] = 'response'

  class = names(processed_data)[1];  features= names(processed_data)[-1]
  
  # Get values of the features and the class
  features <- unlist(features)
  feature.data <- processed_data[, features, drop = FALSE]
  
  # classes <- levels( factor(unlist(processed_data[class]),
  #                           levels = c("Responder", "Non-responder") ) )


  classes = levels(as.factor(processed_data$response))
  
  # nrow of data
  nDataItems = nrow(feature.data)
  # n features
  nFeatures  = ncol(feature.data)
  
  nonresponse.class  = classes[1];  response.class  = classes[2] 
  
  
  
  data.response  = feature.data[processed_data$response == response.class,]
  data.nonresponse  = feature.data[processed_data$response ==nonresponse.class,]
  
  pValues  = numeric(nFeatures)
  
  cat("Wilcoxon testing", nFeatures, "features", fill=TRUE)
  for (i in 1:nFeatures){
    
    current.response = as.numeric(data.response[, i])
    current.nonresponse = as.numeric(data.nonresponse[, i])
    
    outputObject    = wilcox.test(current.response, current.nonresponse, exact=FALSE)
    pValues[i]      = outputObject$p.value
    
  }
  
  selected_features = c()
  for(j in 1:length(pValues)){
    
    if(pValues[j] < kbestthreshol){
      
      selected_features[j] = features[j]
      
    }else{
      selected_features[j] = NA
    }
  }
  
  selected_features = as.vector(na.omit(selected_features))
  
  return(selected_features)
}





# check for highly correlated features
corThreshold = function(RNAdata, 
                        clinicaData,
                        kbestthreshol=0.001, 
                        th.corr = 0.85,
                        omics_variables = 'ihc'){

  if(omics_variables == 'rna'){
    # run rna preporcessed using  RNA_preprocessed function
    RNaProcessed =  RNA_preprocessed(RNAdata, cellNumber = 'fraction')
    RNaProcessed2 = RNaProcessed[RNaProcessed$Cohort == 'Discovery', ]
    
  }else if(omics_variables == 'ihc'){
    
    RNaProcessed =  RNAdata
    
    RNaProcessed2 = RNaProcessed
    
    ## select response class, pin, hemoglin, commnet 
    clinicaDataSubset  = clinicaData[c('PIN', "Pathology Reference No")]
    
    mergedDataset = merge(clinicaDataSubset,
                          RNaProcessed,
                          by.x = "Pathology Reference No",
                          by.y = "Pathology Reference No", 
                          all.x = TRUE)
  }
  
  ## apply selectKBest
  KBestselectedFeatures = selectKBest(RNaProcessed, 
                                      clinicaData,
                                      kbestthreshol=kbestthreshol,
                                      omics_variables = omics_variables)
  
  
  # select features from kbest method: RNaProcessed2 is for discovery cohort alone.
  selectedFeaturesKBest = RNaProcessed2[KBestselectedFeatures]
  
  # remove reducdant variables
  cor_matrix <- cor(makenumcols(selectedFeaturesKBest), method = 'spearman')
  
  cor_matrix[upper.tri(cor_matrix)] <- 0
  
  diag(cor_matrix) <- 0
  
  # SelectedFeactures using correaltion threshold
  SelectedFeactureCOR =  KBestselectedFeatures[!apply(cor_matrix, 2, function(x) any(abs(x) > th.corr, na.rm = TRUE))]
  
  if(omics_variables == 'rna'){
  ### add this if dropped
  toadd = c("IFNg expanded immune 18",
            "Antigen Presentation",
            "Cytotoxic cells")

  toaddcheck = c()
  for(i in 1:length(toadd)){
    if ( !(toadd[i]  %in%  SelectedFeactureCOR) ){
      toaddcheck[i]  = toadd[i]
    }
  }

  toaddcheck = toaddcheck[!is.na(toaddcheck)]
  
  #toaddcheck = NULL
  
  
  }else if(omics_variables == 'ihc'){
    
    toaddcheck = NULL
  }
  
  ## selected features based on correlation threshold
  SelectedFeactureCOR = unique(c(SelectedFeactureCOR, toaddcheck))
  

  if(omics_variables == 'rna'){
  
    ### selected cells based on kbest and correlation approach
    selectedFeaturesDFFinal <- RNaProcessed[SelectedFeactureCOR]
    
    
    selectedFeaturesDFFinal['PIN'] = RNaProcessed$PIN
    selectedFeaturesDFFinal['Path ID'] = RNaProcessed$`Path ID`
    selectedFeaturesDFFinal['Cohort']  = RNaProcessed$Cohort
    
  }else if(omics_variables == 'ihc'){
    
    ### selected cells based on kbest and correlation approach
    selectedFeaturesDFFinal <- mergedDataset[SelectedFeactureCOR]
    
    selectedFeaturesDFFinal['PIN'] = mergedDataset$PIN
    selectedFeaturesDFFinal['Pathology Reference No'] = mergedDataset$`Pathology Reference No`
    
  }
  
  return(list(KBestselectedFeatures = selectedFeaturesKBest,
              
              cor_matrix = cor_matrix,
              
              selectedFeaturesDFFinal = selectedFeaturesDFFinal))
}





# combine clincal and RNA
combinedClinicalandRNA = function(otherData, 
                                     clinicaData){
  
  names(otherData)[dim(otherData)[2]-1] =  "Pathology Reference No"
  
  #otherData['Pathology ID'] = otherData["Pathology Reference No"]
  
  # mergeddatasets
  mergedDataset = merge(clinicaData, 
                        otherData, 
                        by.x = "Pathology Reference No", 
                        by.y = "Pathology Reference No", 
                        all.x = TRUE)
  
  
  row.names(mergedDataset) =  mergedDataset$PIN.x
  mergedDataset = subset(mergedDataset, select=-c(PIN.x, PIN.y,  Cohort.y))
  names(mergedDataset) <- gsub("\\.", " ", names(mergedDataset))

  #names(mergedDataset)[dim(mergedDataset)[2]] = "Pathology Reference No"
  names(mergedDataset)[2] = "Cohort"
  
  return(mergedDataset)
}


# combine clincal and IHC
combinedClinicalandIHC = function(IHC, 
                                  clinicaData){
  
  # mergeddatasets
  
  IHC = subset(IHC, select=-c(`Pathology Reference No`))
  
  mergedDataset = merge(clinicaData, 
                        IHC, 
                        by.x = "PIN", 
                        by.y = "PIN", 
                        all.x = TRUE)
  
  
  row.names(mergedDataset) =  mergedDataset$PIN
  mergedDataset = subset(mergedDataset, select=-c(PIN))
  names(mergedDataset) <- gsub("\\.", " ", names(mergedDataset))
  
  
  return(mergedDataset)
}




# combine clinical and TMB
combinedClinicalandTMB22 = function(TmbData, clinicalData){
  
   names(TmbData)[2] = "Pathology Reference No"
   TmbData$`Pathology Reference No` = gsub("\\_", "-", TmbData$`Pathology Reference No`)
   
   names(clinicalData)[2] = "Pathology Reference No"
     
     
   ## 2 repeated samples: "17H06302"   "71576-15SK"
   #TmbData$`Pathology Reference No`[duplicated(TmbData$`Pathology Reference No`)] 
   
   
   # retain one duplicate
   TmbDataduplicateRemoved = TmbData[!duplicated(TmbData$`Pathology Reference No`),]
   ## no repeats TMB 
   # TmbDataduplicateRemoved$`Pathology Reference No`[duplicated( TmbDataduplicateRemoved$`Pathology Reference No`)] 
   
   
   
   ## sort TmbDataduplicateRemoved by `Pathology Reference No`
   TmbDataduplicateRemoved =    TmbDataduplicateRemoved[order(TmbDataduplicateRemoved$`Pathology Reference No`),]  

  
   ## sort clinicalData by `Pathology Reference No`
   clinicalDataSort =    clinicalData[order(clinicalData$`Pathology Reference No`),]  

  
  # mergeddatasets
  mergedDataset = merge(clinicalDataSort, 
                        TmbDataduplicateRemoved, 
                        by.x = "Pathology Reference No", 
                        by.y = "Pathology Reference No", 
                        all.x = TRUE)
  
  

  row.names(mergedDataset) =  clinicalDataSort$PIN
  mergedDataset = subset(mergedDataset, select=-c(PIN, Cohort.y))
  
  names(mergedDataset) <- gsub("\\.", " ", names(mergedDataset))
  
  names(mergedDataset)[2] = "Cohort"
  
  return(mergedDataset)
}




# combine clinical, RNA and TMB
combinedClinicalRNAandTMB = function(TmbData, clinicalRNAData){
  
  names(TmbData)[2] = "Pathology Reference No"
  TmbData$`Pathology Reference No` = gsub("\\_", "-", TmbData$`Pathology Reference No`)
  
  
 
  
  ## 2 repeated samples: "17H06302"   "71576-15SK"
  #TmbData$`Pathology Reference No`[duplicated(TmbData$`Pathology Reference No`)] 
  
  
  # retain one duplicate
  TmbDataduplicateRemoved = TmbData[!duplicated(TmbData$`Pathology Reference No`),]
  ## no repeats TMB 
  # TmbDataduplicateRemoved$`Pathology Reference No`[duplicated( TmbDataduplicateRemoved$`Pathology Reference No`)] 
  
  
  
  ## sort TmbDataduplicateRemoved by `Pathology Reference No`
  TmbDataduplicateRemoved =    TmbDataduplicateRemoved[order(TmbDataduplicateRemoved$`Pathology Reference No`),]  
  
  
  ## append MELPIN to clinical dataset
  clinicalRNAData['PIN'] = row.names(clinicalRNAData)
  
  
  ## sort clinicalRNAData by `Pathology Reference No`
  clinicalRNADataSort =    clinicalRNAData[order(clinicalRNAData$`Pathology Reference No`),]  
  
  
  
  # mergeddatasets
  mergedDataset = merge(clinicalRNADataSort, 
                        TmbDataduplicateRemoved, 
                        by.x = "Pathology Reference No", 
                        by.y = "Pathology Reference No", 
                        all.x = TRUE)
  
  
  
  row.names(mergedDataset) =  clinicalRNADataSort$PIN
  mergedDataset = subset(mergedDataset, select=-c(PIN, Cohort.y))
  
  names(mergedDataset) <- gsub("\\.", " ", names(mergedDataset))
  names(mergedDataset)[2] = 'Cohort'
  
  return(mergedDataset)
}





# combine clincal and IHC
combinedClinicalRNAandIHC = function(IHC, 
                                  clinicaRNAData){
  
  # mergeddatasets
  
  clinicaRNAData['PIN'] = row.names(clinicaRNAData)
  
  mergedDataset = merge(clinicaRNAData, 
                        IHC, 
                        by.x = "Pathology Reference No", 
                        by.y = "Pathology Reference No", 
                        all.x = TRUE)
  
  
  row.names(mergedDataset) =   mergedDataset$PIN.x
  mergedDataset = subset(mergedDataset, select=-c(PIN.x, PIN.y))
  names(mergedDataset) <- gsub("\\.", " ", names(mergedDataset))
  
  
  return(mergedDataset)
}



############################################# clinical dataset

clinicalDataImport <- read.csv("New clinical samples combined.csv")
names(clinicalDataImport) <- gsub("\\.", " ", names(clinicalDataImport))

clinicalDataImport = clinicalDataImport[!duplicated(clinicalDataImport$PIN),]


##############################################
############################## read RNA-seq data
RNAdata <- read.csv("New RNA combined.csv",  header = TRUE)

## RNA selected Features
SelectedRnaFeaturesResultsList = corThreshold(RNAdata,
                                              clinicalDataImport,
                                              kbestthreshol=0.01, 
                                              th.corr = 0.85,
                                              omics_variables = 'rna')

# selected RNA Features
SelectedRnaFeatures = SelectedRnaFeaturesResultsList$selectedFeaturesDFFinal

## all RNA features
RNaFeaturesAll =  RNA_preprocessed(RNAdata, cellNumber = 'all')
############################


############################# TMB data
TMBData = read.csv('New TMB samples combined.csv', header = TRUE)
names(TMBData) <- gsub("\\.", " ", names(TMBData))
TMBData$`Pathology ID`  <- gsub("\\_", "-", TMBData$`Pathology ID`)
TMBData = TMBData %>% filter(`TMB VAF  5 ` == "PASS")
TMBData = TMBData %>% filter(`Overall TMB data QC` == 'PASS')
#################################


############################# IHC data
IHCData = read.csv('New IHC samples.csv', header = TRUE)

# -ve : Neg
# +ve : Pos
# - : 
# %: Per

IHCData = na.omit(IHCData)

names(IHCData) <- gsub("\\.", " ", names(IHCData))

IHCData$`Pathology Reference No` = gsub("\\_", "-", IHCData$`Pathology Reference No`)
IHCData = IHCData[!duplicated(IHCData$`Pathology Reference No`),]

## make column numeric
for (i in 1:ncol(IHCData)){
  if(i == 1){
    IHCData[, i] = unlist(IHCData[, 1])
  }else{
    
    IHCData[, i] = as.numeric(unlist(IHCData[, i]))
  }
}


## IHC selected Features
SelectedIHCFeaturesResultsList = corThreshold(IHCData,
                                              clinicalDataImport,
                                              kbestthreshol=0.01,
                                              th.corr = 0.85,
                                              omics_variables = 'ihc')

# selected IHC Features
SelectedIHCFeatures = SelectedIHCFeaturesResultsList$selectedFeaturesDFFinal


### All IHC features
clinicalDataImportcopy  = clinicalDataImport[c('PIN', "Pathology Reference No")]

IHCmergedDataset = merge( clinicalDataImportcopy,
                        IHCData,
                         by.x = "Pathology Reference No",
                         by.y = "Pathology Reference No", 
                         all.x = TRUE)

IHCFeaturesAll = IHCmergedDataset
#######################








################################################################################
###################### appended pathology reference number to internal validation
# data_intervalclinical = read.csv('data preprocessed/internalValidation Set PIP_Clinical Data_JCO.csv')
# pin_intervalclinical = read.csv('data preprocessed/internalvalidation_PIN_pathology_references_JCO.csv')
# mergedDataset = merge(data_intervalclinical, 
#                       pin_intervalclinical, 
#                       by.x = "PIN", 
#                       by.y = "PIN", 
#                       all.x = TRUE)
# names(mergedDataset) <- gsub("\\.", " ", names(mergedDataset))
# mergedDataset =  mergedDataset[!is.na(mergedDataset$`Pathology Reference No`), ]
# write.csv(mergedDataset, 'data preprocessed/internalselected_validationAppendedPathologyreference_JCO.csv')
# 