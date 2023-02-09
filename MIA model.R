
## Full Model (Clinical Model by Ines P. 2021)
MIAOutcomeRisk <- function(dataN){
  ## This program is to estimate individual's risk of death, progression and response rate wihtin 1 year for diagnosis date
  ## 9 input parameters are required
  ##1. ECOG (0,1,2,3)
  ##2. Hemoglobin - continuous from -- to --
  ##3. LDH ('Normal','Elevated')
  ##4  Neutrophil count (x 10^9/L)
  ##5. Lymphocyte count (x 10^9/L)
  ##6. Neutro-Lympho ratio - continuous from -- to --
  ##7. Brain metastases ('No', 'Yes') 
  ##8. Lung metastases ('No', 'Yes')
  ##9. Liver metastases ('No', 'Yes')
  ##10. Treatment ('PD1', 'IPI+PD1')  
  ##11. First line treatment ('Yes','No')
  # SL 25-01-22: changed risk of outcome to probability of outcome-free (IS query)
  # SL 27-01-22: changed risk of outcome to probability of outcome-free (IS query)
  
  ## Get Risk of Objective Response:
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  ECOG = dataN[1]
  LDH = dataN[2]
  NLR = dataN[3]
  Lung  = dataN[4]
  Liver = dataN[5]
  Treatment = dataN[6]
  FirstLine = dataN[7]
  
  ## ECOG
  {if (ECOG == 0) {c1 = 0}
    else {c1 = -0.2982}}
  ## LDH
  {if (LDH=="Normal") {c2 = 0}
    else if (LDH == "Elevated"){c2= -0.6714}}
  ## NLR
  {c3 = -0.1545*NLR}
  ## Lung metastases
  {if (Lung=="No") {c4 = 0}
    else if (Lung == "Yes"){c4= 0.4636}}
  ## Liver metastases
  {if (Liver=="No") {c5 = 0}
    else if (Liver == "Yes"){c5 = -0.4714}}
  ## Treatment
  {if (Treatment=="PD1") {c6 = 0}
    else if(Treatment == "IPI+PD1"){c6 = 0.1684}}
  ## Line of Treatment
  {if (FirstLine=="1st Line") {c7 = 0}
    else {c7 = -0.6579}}
  
  ### Get Prognosis index (Beta*X)
  BX_ORR = 0.8887 + c1 + c2 + c3 + c4 + c5 + c6 + c7
  ### Absolute probability of overall response (CR or PR)
  Risk_ORR = exp(BX_ORR)/(1+exp(BX_ORR))
  
  return(Risk_ORR)
}


MIAPrediction = function(dataN){
  
  predictedprob = prediction= c()
  for(i in 1:nrow(dataN)){
    ##Example ECOG, Hb,LDH,Neutrophil,Lymphocyte,Brain,Lung,Liver,Treatment,FirstLine
   predictedprob[i] = MIAOutcomeRisk(dataN[i,])
   prediction[i] = ifelse(predictedprob[i] > 0.5, "R", "NR")
  }
  
  prediction = factor(as.vector(prediction), levels = c('R', 'NR'))
  predictedprob = as.data.frame(as.numeric(predictedprob))
  colnames(predictedprob) = 'R'
  predictedprob['NR'] = 1- as.numeric(unlist(predictedprob))
  
  
  return(list(predicted =  prediction,
              predictedprob=  predictedprob ))
}


