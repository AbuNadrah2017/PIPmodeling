library(survival)
library(survminer)


plot_quantiles=function (score, splits = 10, fill = "deepskyblue", model_name = NA, subtitle = NA, 
                         table = FALSE, save = FALSE, subdir = NA, file_name = "viz_ncuts.png"){
  if (splits > 25) 
    stop("You should try with less splits!")
  deciles <- quantile(score, probs = seq((1/splits), 1, length = splits), 
                      names = TRUE)
  deciles <- data.frame(range = row.names(as.data.frame(deciles)), 
                        cuts = as.vector(signif(deciles, 6)))
  rownames(deciles) <- NULL
  p <- deciles %>% ggplot(aes(x = reorder(.data$range, .data$cuts), 
                              y = .data$cuts * 100)) + geom_col(fill = fill, 
                                                                colour = "transparent") +
    xlab("Cumulative volume") + 
    ylab("Score") +
    geom_text(aes(label = round(100 * .data$cuts, 1), vjust = ifelse(.data$cuts * 100 < 50, -0.3, 1.3)), size = 3, colour = "black", inherit.aes = TRUE, check_overlap = TRUE) + guides(colour = "none") + 
    labs(title = sprintf("Score cuts (%s quantiles)", splits)) + theme_lares() +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.title.x = element_text(hjust = 0.5)) +
    theme(plot.subtitle = element_text(hjust = 0.5))
  
  if (!is.na(subtitle)) 
    p <- p + labs(subtitle = subtitle)
  if (!is.na(model_name)) 
    p <- p + labs(caption = model_name)
  if (save) 
    export_plot(p, file_name, subdir = subdir, width = 6, height = 6)
  if (table) {
    return(deciles)
  }
  else {
    return(p)
  }
}





PFS_and_response_prop_sumr = function(aa, 
                                           percentiles = c(0.25, 0.75), 
                                           Kaplan_mair_title,
                                      type,
                                      legendPosition = "bottom",
                                      title = NULL, 
                                      x_lim = NULL,
                                      y_lim = NULL,
                                      xlab = "Predicted proportion",
                                      ylab = "Observed proportion",
                                      points_col_list = NULL,
                                      data_summary = FALSE){
  
  train_data_used = aa$cal_trainData
  # Load the data into R
  prediction_prob <- as.numeric(unlist(unlist(train_data_used$pred)))
  
  # Calculate the quantiles
  bottom_percentile <- quantile(prediction_prob, probs = percentiles[1])
  middle_percentile <- quantile(prediction_prob, probs = c(percentiles[1], percentiles[2]))
  top_percentile <- quantile(prediction_prob, probs = percentiles[2])
  
  # Create a new variable to store the values
  predict_class_names <- numeric(length(prediction_prob))
  
  # Assign 1, 2, and 3 based on the percentiles
  predict_class_names[prediction_prob < bottom_percentile] <- paste0(' Bottom ', percentiles[1]*100, 'th percentile')
  predict_class_names[prediction_prob >= bottom_percentile &  prediction_prob < middle_percentile[2]] <- paste0(percentiles[1]*100, 'th', '-', percentiles[2]*100, 'th percentile')
  predict_class_names[ prediction_prob >= top_percentile] <- paste0('Top ', percentiles[2]*100, 'th percentile')
  
  responder_class = train_data_used$response
  predict_class2 = as.data.frame(cbind(predict_class_names, prediction_prob))
  colnames(predict_class2) = c('predict_class_names', 'prediction_prob')
  row.names(predict_class2)  = row.names(aa$prediction_r$train_predictionsProb)
  predict_class2['response class'] = responder_class
  
  res.prop.table = prop.table(predict_class2 , type='R')
  nonres.prop.table = prop.table(predict_class2 , type='NR')  

  
  #### extract dataset
  clincal_data = read.csv('New clinical samples combined.csv')
  clincal_data = clincal_data[clincal_data$Cohort == "Discovery", ]
  row.names(clincal_data) = clincal_data$PIN
  
  
  train_data_used['PIN'] = row.names(train_data_used)
  train_data_used['predict_class'] = predict_class_names
  train_data_used['predict_prob'] = prediction_prob
  train_data_used['train_predicted'] = train_data_used$train_predicted
  
  
  # mergeddatasets
  comBTrainData=merge(train_data_used,
                      clincal_data, 
                      by.x = "PIN",
                      by.y = "PIN",
                      all.x = TRUE)
  
  row.names(comBTrainData) = comBTrainData$PIN
  
  extract_data = comBTrainData[c("response",
                                 "predict_class",
                                 "predict_prob",
                                 "Progressed.Yes.or.No",
                                 "Progression.free.survival.in.day",
                                 'PIN',
                                 'train_predicted')]
  
  ## convert PFS to month
  extract_data["PFS (Months)"] = (extract_data$Progression.free.survival.in.day/30.44)
  
  
  ### export recist dataset
  recist_data = read.csv("PIP_Retrospective_RECIST.csv")
  
  # Rearrange the columns of the data frame
  comBTrainData_Recist=merge(extract_data,
                             recist_data, 
                             by.x = "PIN",
                             by.y = "PIN",
                             all.x = TRUE)
  row.names(comBTrainData_Recist) = comBTrainData_Recist$PIN
  
  ## conver Progress to 1 and 0
  comBTrainData_Recist$progressed[comBTrainData_Recist["Progressed.Yes.or.No"] == "Yes"] = 1   #(recurrence : 1)
  comBTrainData_Recist$progressed[comBTrainData_Recist["Progressed.Yes.or.No"] == "No"] = 0   #(censoring or non-recurrence: 0)
  
  
  # #######
  # comBTrainDatalower = comBTrainData_Recist[comBTrainData_Recist$predict_class == levels(as.factor(predict_class_names))[1], ]
  # comBTrainDatamiddle = comBTrainData_Recist[comBTrainData_Recist$predict_class == levels(as.factor(predict_class_names))[2], ]
  # comBTrainDataupper = comBTrainData_Recist[comBTrainData_Recist$predict_class == levels(as.factor(predict_class_names))[3], ]
  # 
  # ################# find ORR 
  # number_of_lower = nrow(comBTrainDatalower)
  # number_of_middle = nrow(comBTrainDatamiddle )
  # number_of_upper = nrow(comBTrainDataupper)
  # 
  # ## number of number_of CR and PR
  # number_of_CR_PR_lower = nrow(subset(comBTrainDatalower, RECIST == "CR" | RECIST == "PR"))
  # number_of_CR_PR_middle = nrow(subset(comBTrainDatamiddle, RECIST == "CR" | RECIST == "PR"))
  # number_of_CR_PR_upper = nrow(subset(comBTrainDataupper, RECIST == "CR" | RECIST == "PR"))
  #  
  # ## ORR
  # ORRlower <-  number_of_CR_PR_lower/ number_of_lower
  # ORRmiddle <-   number_of_CR_PR_middle/ number_of_middle
  # ORRupper <-   number_of_CR_PR_upper/ number_of_upper
  # 
  # 
  # 
  # fit1 <- survfit(Surv(`PFS (Months)`, progressed) ~ 1, data = comBTrainDatalower)
  # fit2 <-  survfit(Surv(`PFS (Months)`, progressed) ~ 1,  data = comBTrainDatamiddle)
  # fit3 <-  survfit(Surv(`PFS (Months)`, progressed) ~ 1,  data =comBTrainDataupper)
  # 
  # summary(fit1, times = 6)
  # summary(fit2, times = 6)
  # summary(fit3, times = 6)
  
  fit <- survfit(Surv(`PFS (Months)`, progressed) ~ predict_class, data = comBTrainData_Recist)
  summary(fit, times = 6)
  
  labelS =  c(  paste0(' Bottom ', percentiles[1]*100, 'th'),
                paste0(percentiles[1]*100, 'th', ' - ', percentiles[2]*100, 'th'),
                paste0('Top ', percentiles[2]*100, 'th'))
  
  pfs_plot =   PFS_plots(complete_dataset=comBTrainData_Recist,
                         main_title=Kaplan_mair_title,
                         percentiles=percentiles,
                         fit=fit,
                         labelS =   labelS,
                         legend.title = 'percentile')  
 
   
  
  
   fit_response1 <- survfit(Surv(`PFS (Months)`, progressed) ~ train_predicted, 
                          data = comBTrainData_Recist)
  #summary(fit_response, times = 6)
  labelS = c('R', 'NR')
  pfs_plot_response1 =   PFS_plots(complete_dataset=comBTrainData_Recist,
                         main_title=Kaplan_mair_title,
                         percentiles=percentiles,
                         fit =  fit_response1,
                         labelS =    labelS,
                         legend.title = 'predicted response class')  
  
  
  fit_response <- survfit(Surv(`PFS (Months)`, progressed) ~ response + train_predicted, 
                          data = comBTrainData_Recist)
  #summary(fit_response, times = 6)
  labelS = NULL# c('NR', 'R', 'a', 'c')
  pfs_plot_response =   PFS_plots(complete_dataset=comBTrainData_Recist,
                                  main_title=Kaplan_mair_title,
                                  percentiles=percentiles,
                                  fit =  fit_response,
                                  labelS =    labelS,
                                  legend.title = 'response class')  
  
 
  dataDec_mods = rbind(res.prop.table, nonres.prop.table)

  calibPlot_obj <- ggplot(data = dataDec_mods, aes(y = obsRate, x = predRate, color = as.factor(type_name))) + 
    geom_point(size = 4) + 
    lims(x = ifelse(rep(is.null(x_lim), 2), c(min(dataDec_mods$predRate), max(dataDec_mods$predRate)), x_lim), 
         y = ifelse(rep(is.null(y_lim), 2), c(min(dataDec_mods$obsRate_LCL), max(dataDec_mods$obsRate_UCL)), y_lim)) + 
    geom_errorbar(aes(ymax = obsRate_UCL, ymin = obsRate_LCL), width = 0.1) + 
    geom_abline(intercept = 0, slope = 1) + 
    scale_color_manual(values = c("red", "blue"), labels = c("Non responder", "Responder")) + 
    labs(x = ifelse(is.null(xlab), pred, xlab), y = ifelse(is.null(ylab), obs, ylab), title = title) + 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"),
          axis.text = element_text(colour = "black", size = 12), 
          legend.position = legendPosition) + 
    theme(plot.title = element_text(hjust = 0.5)) +
    guides(colour = guide_legend(title = "Response class")) + 
    geom_text(data = subset(dataDec_mods, type_name == 'R' & decile == levels(as.factor(predict_class_names))[1]), 
              aes(label = levels(as.factor(predict_class_names))[1]), x = dataDec_mods$predRate[1], y = dataDec_mods$obsRate[1], size = 3) + 
    geom_text(data = subset(dataDec_mods, type_name == 'R' & decile == levels(as.factor(predict_class_names))[2]), 
              aes(label = levels(as.factor(predict_class_names))[2]), x = (dataDec_mods$predRate)[2], y = (dataDec_mods$obsRate)[2], size = 3) + 
    geom_text(data = subset(dataDec_mods, type_name == 'R' & decile == levels(as.factor(predict_class_names))[3]), 
              aes(label = levels(as.factor(predict_class_names))[3]), x = (dataDec_mods$predRate)[3], y = (dataDec_mods$obsRate)[3], size = 3) + 
    geom_text(data = subset(dataDec_mods, type_name == 'NR' & decile == levels(as.factor(predict_class_names))[1]), 
              aes(label = levels(as.factor(predict_class_names))[1]), x = (dataDec_mods$predRate)[4]+0.01, y = (dataDec_mods$obsRate)[4], size = 3) + 
    geom_text(data = subset(dataDec_mods, type_name == 'NR' & decile == levels(as.factor(predict_class_names))[2]), 
              aes(label = levels(as.factor(predict_class_names))[2]), x = (dataDec_mods$predRate)[5], y = (dataDec_mods$obsRate)[5], size = 3) +
  geom_text(data = subset(dataDec_mods, type_name == 'NR' & decile == levels(as.factor(predict_class_names))[3]), 
            aes(label = levels(as.factor(predict_class_names))[3]), x = (dataDec_mods$predRate)[6], y = (dataDec_mods$obsRate)[6], size = 3) 
  
  
  results = list(KM_plot = pfs_plot,
                fit = fit,
                calibPlot_obj =  calibPlot_obj,
                dataDec_mods = dataDec_mods,
                res.prop.table = res.prop.table,
                nonres.prop.table =  nonres.prop.table,
                PFS_fit_6months = summary(fit, times = 6),
                PFS_fit_1year = summary(fit, times = 12),
                fit_response1 =  fit_response1,
                fit_response =  fit_response, 
                pfs_plot_response1  = pfs_plot_response1,
                pfs_plot_response = pfs_plot_response)
  
  return(results)
}



prop.table= function(predict_class2 , type){
  ## get proportion and numbers of responders with each percentiles
  obsRate = obsRate_SE =  number_of_responders =  number_of_nonresponders=  predRate = c()
  
  type_name =  decile = c()
  
  for (i in 1:length(levels(as.factor(predict_class2$predict_class_names)))){
    
    if(type == 'R'){
      
      aaan=predict_class2[predict_class2$predict_class_names == levels(as.factor(predict_class2$predict_class_names))[i], ]
      
      aaan_01 = ifelse(aaan$`response class` == 'R', 1, 0)
      
      obsRate[i]= mean(aaan_01)#table(aaan$`response class`)[1]/sum(table(aaan$`response class`))
      
      obsRate_SE[i]= sd(aaan_01)/sqrt(length(aaan_01))
      
      number_of_responders[i] = table(aaan$`response class`)[1]
      
      number_of_nonresponders[i] = table(aaan$`response class`)[2]
      
      
      predRate[i]  = mean(as.numeric(aaan[aaan$`response class` == 'R',]$prediction_prob))
      
      type_name[i] = 'R'
      
      decile[i] = levels(as.factor(predict_class2$predict_class_names))[i]
    }else{
      
      aaan=predict_class2[predict_class2$predict_class_names == levels(as.factor(predict_class2$predict_class_names))[i], ]
      
      aaan_01 = ifelse(aaan$`response class` == 'NR', 1, 0)
      
      obsRate[i]= mean(aaan_01)#table(aaan$`response class`)[2]/sum(table(aaan$`response class`))
      
      obsRate_SE[i]= sd(aaan_01)/sqrt(length(aaan_01))
      
      number_of_responders[i] = table(aaan$`response class`)[1]
      
      number_of_nonresponders[i] = table(aaan$`response class`)[2]
      
      
      predRate[i]  = mean(as.numeric(aaan[aaan$`response class` == 'NR',]$prediction_prob))
      
      type_name[i] = 'NR'
      
      decile[i] = levels(as.factor(predict_class2$predict_class_names))[i]
    }
  }
  
  results.prop = as.data.frame(number_of_responders) 
  names(results.prop) = c('no. of R')
  
  row.names(results.prop) = levels(as.factor(predict_class2$predict_class_names))
  results.prop['no. of R'] =  number_of_responders
  results.prop['no. of NR'] =  number_of_nonresponders
  results.prop['total patients'] =  (number_of_nonresponders + number_of_responders)
  
  results.prop['obsRate'] = obsRate 
  results.prop['obsRate_SE'] =  obsRate_SE

  
  #results.prop.ordered = results.prop[order(results.prop$`Prop. of R`),]
  results.prop$obsRate_UCL <- results.prop$obsRate + 1.96 * results.prop$obsRate_SE
  results.prop$obsRate_LCL <- results.prop$obsRate - 1.96 * results.prop$obsRate_SE
  
  results.prop['predRate'] = predRate
  
  results.prop['type_name'] = type_name
  results.prop['decile'] = decile
  
  return(results.prop)
}



PFS_plots = function(complete_dataset, main_title, percentiles,  fit,  labelS, legend.title){
  

  #fit <- survfit(Surv(`PFS (Months)`, progressed) ~ predict_class, data = complete_dataset)
  
    
  #labelS =  c(  paste0(' Bottom ', percentiles[1]*100, 'th'),
  #              paste0(percentiles[1]*100, 'th', ' - ', percentiles[2]*100, 'th'),
  #  paste0('Top ', percentiles[2]*100, 'th'))

  #labelS =  c(paste0(' Bottom ', percentiles[1]*100, 'th percentile'),
  #            paste0(percentiles[1]*100, 'th', '-', percentiles[2]*100, 'th percentile'),
  #           paste0('Top ', percentiles[2]*100, 'th percentile'))

    
  ggtheme = theme_minimal() + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  surv_plot1 <- ggsurvplot(
    fit,                     # survfit object with calculated statistics.
    data = complete_dataset,             # data used to fit survival curves.
    risk.table = TRUE,       # show risk table.
    pval = TRUE,             # show p-value of log-rank test.
    pval.method = TRUE,
    surv.scale = "percent",
    # conf.int = TRUE,         # show confidence intervals for 
    # point estimates of survival curves.
    #  palette = c("#E7B800", "#2E9FDF"),
    #  xlim = c(0, round(max(complete_dataset$`Progression-free survival (Months)`))),         # present narrower X axis, but not affect
    xlim = c(0, 72),         # present narrower X axis, but not affect
    title = main_title,
    legend.title=legend.title,
    pval.coord = c(55, 0.82),
    pval.method.coord = c(55, 0.92),
    # survival estimates.
    xlab = "Follow-up time (Months)",   # customize X axis label.
    ylab = "Progression-free survival (PFS)",   # customize X axis label.
    break.time.by = 6,     # break X axis in time intervals by 500.
    risk.table.y.text.col = T,# colour risk table text annotations.
    risk.table.height = 0.25, # the height of the risk table
    risk.table.y.text = FALSE,# show bars instead of names in text annotations
    risk.table.col = "strata",# Risk table color by groups
    # in legend of risk table.
    # ncensor.plot = TRUE,      # plot the number of censored subjects at time t
    #  ncensor.plot.height = 0.25,
    #conf.int.style = "step",  # customize style of confidence intervals
    surv.median.line = "hv",  # add the median survival pointer.
    legend.labs =    labelS,   # change legend labels.
    ggtheme =  theme_light() + theme(plot.title = element_text(colour = "black", face = "bold", size = 15, hjust = 0.5),
                                     legend.text = element_text(colour = "black", face = "bold", size = 12),
                                     legend.title = element_text(colour = "black", face = "bold", size = 12),
                                     axis.text = element_text(colour = "black", face='bold', size = 12),
                                     axis.text.y = element_text(colour = "black", face='bold', size = 12),
                                     axis.text.x = element_text(colour = "black", face='bold', size = 12),
                                     axis.title.y = element_text(colour = "black", face='bold', size = 12),
                                     axis.title.x = element_text(colour = "black", face='bold', size = 12))  # customize plot and risk table with a theme.
  )
  
  
  return(surv_plot1)  
}

