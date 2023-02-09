library(caret)
library(boot)
library(boot.pval)
library(doParallel)
library(foreach)
library(parallel)



ODD_ratio_estimates = function(print_result, 
                               learner = 'glmnet', 
                               boot.size = 1000L){

  cl <- makeCluster(detectCores())
  doParallel::registerDoParallel(cl)
  
  a = print_result  
  best_model = a$trainedModel_r$modelResults$best_model$`best model`
  new_train_data = best_model$call$data
  class = 'response'
  learner = best_model$call$method
  
  coeff_ = coef(best_model$finalModel, best_model$bestTune$lambda)
  predictors <- row.names(coef(best_model$finalModel, best_model$bestTune$lambda))
  
  
  fittingParams <- list(metric = best_model$call$metric,
                        family = "binomial",
                        importance = TRUE,
                        standardize = FALSE,
                        tuneGrid = best_model$bestTune, 
                        tuneLength = best_model$call$tuneLength)
  
  ctrl2 <- best_model$call$trControl

  # refit model to see if the estimates are correct.
  # set.seed(1234)
  # final_fit <-  do.call(caret::train, append(list(
  #   form = as.formula(paste( class, ".", sep = "~")),
  #   data = new_train_data,
  #   method = learner,
  #   trControl = ctrl2),
  #   fittingParams))
  
  #coef(final_fit$finalModel, final_fit$bestTune$lambda)
  #coef(best_model$finalModel, best_model$bestTune$lambda)
  

  bootSamples <- boot(new_train_data, function(data, idx) {
    
    bootstrapData <- new_train_data[idx, ]
    
    set.seed(1234)
    bootstrapMod <- do.call(caret::train, append(list(
      form = as.formula(paste( class, ".", sep = "~")),
      data = bootstrapData,
      method = learner,
      trControl = ctrl2),
      fittingParams))
    
    as.vector(exp(coef(bootstrapMod$finalModel, best_model$bestTune$lambda)))
  }, R=boot.size)
  
  foreach::registerDoSEQ()
  stopCluster(cl)
  
  
  
  # Extract the coefficients and standard errors from the results
  coefs <- bootSamples$t
  colnames(coefs) = predictors
  #coefs
  
  
  # Calculate the rate risk: model estimate
  rate_risk <- apply(log(coefs), 2, mean)
  se_rr <- apply(log(coefs), 2, sd)
  
  ## odd ratio
  odd_ration_mean = apply(coefs, 2, mean)
  #odd_ration_mean
  
  odd_ration_sd = apply(coefs, 2, sd)
  #odd_ration_sd
  
  bootstrapSE <- apply(coefs, 2, function(x) sd(x)/sqrt(length(x)))

  # confidence interval
  #apply(coefs, 2, quantile, probs = c(0.025, 0.5, 0.975))
  
  conff = loweroddratio = upperoddratio = pval = list()
  for(i in 1:length(predictors)){
    #  print('odd ratio')
    conff[[i]] = boot.ci(bootSamples, index = c(i), type = c("perc"))
    loweroddratio[[i]] = conff[[i]]$percent[4]
    upperoddratio[[i]] = conff[[i]]$percent[5]
    # pvalues
    pval[[i]] = boot.pval(bootSamples,
                          type = "perc", 
                          theta_null = coeff_[i], 
                          index = c(i))
    
  }

  boot_results = cbind(round(as.numeric(coeff_), 3),
                       round(as.numeric(exp(coeff_)), 3),
                       round(odd_ration_mean, 3),
                       round(bootstrapSE, 3),
                       round(unlist(loweroddratio), 3),
                       round(unlist(upperoddratio), 3),
                       unlist(pval))
  
  boot_results = as.data.frame(boot_results)
  names(boot_results) = c('model coefs',
                          'observed OR',
                          'boot. OR',
                          'boot. SE',
                          'boot. lower limit',
                          'boot. upper limit',
                          'p-value')
  
  return(boot_results)
}





