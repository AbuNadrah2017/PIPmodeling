

calibrationPlot = function (data, obs, follow_up = NULL, pred, group = NULL, nTiles = 10, 
                            legendPosition = "right", title = NULL, x_lim = NULL, 
                            y_lim = NULL, xlab = "Prediction", ylab = "Observation", 
                            points_col_list = NULL, data_summary = FALSE) {
 
  library(dplyr)
  library(magrittr)

  
   if (!exists("obs") | !exists("pred")) 
    stop("obs and pred can not be null.")

  n_groups <- length(unique(data[, group]))
  data %<>% mutate(decile = ntile(!!sym(pred), nTiles))
  
  if (is.null(follow_up)) 
    data$follow_up <- 1
  
  if (!is.null(group)) {
    
    dataDec_mods <- data %>% group_by(.data$decile, !!sym(group)) %>% 
      summarise(obsRate = mean(!!sym(obs)/follow_up, na.rm = T), 
                obsRate_SE = sd(!!sym(obs)/follow_up, na.rm = T)/sqrt(n()), 
                obsNo = n(), 
                predRate = mean(!!sym(pred), na.rm = T))
    colnames(dataDec_mods)[colnames(dataDec_mods) == "group"] <- group
    
  } else {
    dataDec_mods <- data %>% group_by(.data$decile) %>% summarise(obsRate = mean(!!sym(obs)/follow_up, 
                                                                                 na.rm = T), obsRate_SE = sd(!!sym(obs)/follow_up, 
                                                                                                             na.rm = T)/sqrt(n()), obsNo = n(), predRate = mean(!!sym(pred), 
                                                                                                                                                                na.rm = T))
  }
  
  
  dataDec_mods$obsRate_UCL <- dataDec_mods$obsRate + 1.96 * dataDec_mods$obsRate_SE
  dataDec_mods$obsRate_LCL <- dataDec_mods$obsRate - 1.96 * dataDec_mods$obsRate_SE
  dataDec_mods <- as.data.frame(dataDec_mods)
  
  
  if (!is.null(group)) {
    dataDec_mods[, group] <- factor(dataDec_mods[, group])
    calibPlot_obj <- ggplot(data = dataDec_mods, aes(y = .data$obsRate, 
                                                     x = .data$predRate, group = !!sym(group), color = !!sym(group))) + 
      geom_point() + lims(x = ifelse(rep(is.null(x_lim), 
                                         2), c(min(dataDec_mods$predRate), max(dataDec_mods$predRate)), 
                                     x_lim), y = ifelse(rep(is.null(y_lim), 2), c(min(dataDec_mods$obsRate_LCL), 
                                                                                  max(dataDec_mods$obsRate_UCL)), y_lim)) + geom_errorbar(aes(ymax = .data$obsRate_UCL, 
                                                                                                                                              ymin = .data$obsRate_LCL)) + geom_abline(intercept = 0, 
                                                                                                                                                                                       slope = 1) + scale_color_manual(values = ifelse(rep(is.null(points_col_list), 
                                                                                                                                                                                                                                           n_groups), (ggplot2::scale_colour_brewer(palette = "Set3")$palette(8)[c(4:8, 
                                                                                                                                                                                                                                                                                                                   1:3)])[c(1:n_groups)], points_col_list)) + labs(x = ifelse(is.null(xlab), 
                                                                                                                                                                                                                                                                                                                                                                              pred, xlab), y = ifelse(is.null(ylab), obs, ylab), 
                                                                                                                                                                                                                                                                                                                                                                   title = title) + theme(panel.grid.major = element_blank(), 
                                                                                                                                                                                                                                                                                                                                                                                          panel.grid.minor = element_blank(), panel.background = element_blank(), 
                                                                                                                                                                                                                                                                                                                                                                                          axis.line = element_line(colour = "black"), 
                                                                                                                                                                                                                                                                                                                                                                                          legend.key = element_rect(fill = "white"), 
                                                                                                                                                                                                                                                                                                                                                                                          axis.text = element_text(colour = "black", 
                                                                                                                                                                                                                                                                                                                                                                                                                   size = 12), legend.position = legendPosition)
  }
  else {
    calibPlot_obj <- ggplot(data = dataDec_mods, aes(y = .data$obsRate, 
                                                     x = .data$predRate)) + 
      geom_point(color = ggplot2::scale_colour_brewer(palette = "Set3")$palette(8)[5]) + 
     
      lims(x = ifelse(rep(is.null(x_lim), 2), c(min(dataDec_mods$predRate), 
                                                max(dataDec_mods$predRate)), x_lim), y = ifelse(rep(is.null(y_lim), 
                                                                                                    2), c(min(dataDec_mods$obsRate_LCL), max(dataDec_mods$obsRate_UCL)), 
                                                                                                y_lim)) + geom_errorbar(aes(ymax = .data$obsRate_UCL, 
                                                                                                                            ymin = .data$obsRate_LCL), col = ifelse(is.null(points_col_list), 
                                                                                                                                                                    ggplot2::scale_colour_brewer(palette = "Set3")$palette(8)[5], 
                                                                                                                                                                    points_col_list)) + geom_abline(intercept = 0, slope = 1) + 
      scale_color_manual(values = ifelse(is.null(points_col_list), 
                                         ggplot2::scale_colour_brewer(palette = "Set3")$palette(8)[5], 
                                         points_col_list)) + labs(x = ifelse(is.null(xlab), 
                                                                             pred, xlab), y = ifelse(is.null(ylab), obs, ylab), 
                                                                  title = title) + theme(panel.grid.major = element_blank(), 
                                                                                         panel.grid.minor = element_blank(), panel.background = element_blank(), 
                                                                                         axis.line = element_line(colour = "black"), 
                                                                                         axis.text = element_text(colour = "black", 
                                                                                                                  size = 12), legend.position = legendPosition)
  }
  res_list <- list(calibration_plot = calibPlot_obj)
  if (data_summary) 
    res_list$data_summary <- dataDec_mods
  return(res_list)
}



# 
# 
# library(predtools)
# library(magrittr)
# library(dplyr)
# #> 
# #> Attaching package: 'dplyr'
# #> The following objects are masked from 'package:stats':
# #> 
# #>     filter, lag
# #> The following objects are masked from 'package:base':
# #> 
# #>     intersect, setdiff, setequal, union
# library(ggplot2)
# 
# data(dev_data)
# data(val_data)
# 
# 
# reg <- glm(y~sex+age+severity+comorbidity,data=dev_data,family=binomial(link="logit"))
# summary(reg)
# 
# 
# 
# dev_data$pred <- predict.glm(reg, type = 'response')
# val_data$pred <- predict.glm(reg, newdata = val_data, type = 'response')
# 
# cal_result = calibrationPlot(data = dev_data, 
#                              obs = "y", 
#                              nTiles = 4,
#                              pred = "pred",
#                              title = "Calibration plot for development data",
#                              data_summary = TRUE)
# 
# cal_result
# 
# 
# #> $calibration_plot
# #> 
# #> 
# 
# a=calibration_plot(data = val_data,
#                    obs = "y", 
#                    pred = "pred",
#                    y_lim = c(0, 0.6),
#                    title = "Calibration plot for validation data", group = "sex")
# 



# 
# calibPlot_obj <- ggplot(data = dataDec_mods, aes(y = obsRate, x = predRate)) + 
#   geom_point(color = scale_colour_brewer(palette = "Set3")$palette(8)[5]) + 
#   lims(x = ifelse(rep(is.null(x_lim), 2), c(min(dataDec_mods$predRate), max(dataDec_mods$predRate)), x_lim), 
#        y = ifelse(rep(is.null(y_lim), 2), c(min(dataDec_mods$obsRate_LCL), max(dataDec_mods$obsRate_UCL)), y_lim)) + 
#   geom_errorbar(aes(ymax = obsRate_UCL, ymin = obsRate_LCL), 
#                 col = ifelse(is.null(points_col_list), scale_colour_brewer(palette = "Set3")$palette(8)[5], points_col_list)) + 
#   geom_abline(intercept = 0, slope = 1) + 
#   scale_color_manual(values = ifelse(is.null(points_col_list), scale_colour_brewer(palette = "Set3")$palette(8)[5], points_col_list)) + 
#   labs(x = ifelse(is.null(xlab), pred, xlab), y = ifelse(is.null(ylab), obs, ylab), title = title) + 
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
#         axis.line = element_line(colour = "black"), axis.text = element_text(colour = "black", size = 12), 
#         legend.position = legendPosition)
# 
# 
# 
# calibPlot_obj <- ggplot(data = dataDec_mods, aes(y = obsRate, x = predRate, color = as.factor(type))) + 
#   geom_point() + 
#   lims(x = ifelse(rep(is.null(x_lim), 2), c(min(dataDec_mods$predRate), max(dataDec_mods$predRate)), x_lim), 
#        y = ifelse(rep(is.null(y_lim), 2), c(min(dataDec_mods$obsRate_LCL), max(dataDec_mods$obsRate_UCL)), y_lim)) + 
#   geom_errorbar(aes(ymax = obsRate_UCL, ymin = obsRate_LCL)) + 
#   geom_abline(intercept = 0, slope = 1) + 
#   scale_color_manual(values = c("red", "blue")) + 
#   labs(x = ifelse(is.null(xlab), "predRate", xlab), 
#        y = ifelse(is.null(ylab), "obsRate", ylab), 
#        title = ifelse(is.null(title), "Calibration Plot", title)) + 
#   theme(panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(), 
#         panel.background = element_blank(), 
#         axis.line = element_line(colour = "black"), 
#         axis.text = element_text(colour = "black", size = 12), 
#         legend.position = "top")
