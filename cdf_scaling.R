

# marginal fit
preprocess_CDF = function(dat=train_dat,  criterion = "AICc"){
  
  vars = as.vector(get_typesMine(dat, coarse = TRUE))
  orig_vars <- vars
  vars <- vars %in% c("integer", "numeric", "double")
  names(vars) <- names(orig_vars)
  
  marginal_distr <- check_distributions(dat[vars], criterion = criterion)$distr
  
  Y_continuous_marginals <- marginal_parameters_splits_AND_distri(dat[vars], marginal_distr)$marginal_distributions
  
  dat[names(dat[vars])] = Y_continuous_marginals
  
  return( list(numerical_varibles = names(dat[vars]),
               
               train_data = dat,
               
               marginal_distr = marginal_distr, 
               
               Y_all_marginals = Y_continuous_marginals ))
}




predict_CDF = function(preprocess_CDF_object,
                       test_data){
  
  test_dat_cotinous = test_data[preprocess_CDF_object$numerical_varibles]
  
  test_continuous_marginals <- marginal_parameters_splits_AND_distri(test_dat_cotinous, 
                                                                     
                                                                     preprocess_CDF_object$marginal_distr)$marginal_distributions
  
  test_data[names(test_data[preprocess_CDF_object$numerical_varibles])] = test_continuous_marginals
  
  return(test_data)
}



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




check_distributions = function(Y, criterion = "BICc"){
  
  unbiased_mean <- apply(Y, MARGIN = 2, FUN = mean)
  
  st_deviation <- apply(Y, MARGIN = 2, FUN = sd)
  
  n <- dim(Y)[1] ;   p <- dim(Y)[2]
  
  detail <- as.list(1:p); names(detail) <- colnames(Y)
  
  irank <- function (criterion) {

    switch(criterion,
           loglik = rank(loglik),
           AIC = rank(AIC),
           AICc = rank(AICc),
           BIC = rank(BIC),
           BICc = rank(BICc))
  }
  
  dist <- rep(0,p) 
  n1 <- rep(NA, p)
  location_param <- rep(NA,p) 
  scale_param <- rep(NA,p) 
  t_df <- rep(NA, p)
  nu_param <- rep(NA,p) 
  tau_param <- rep(NA,p)

  for (k in 1:p) {
    
    y <- Y[,k] ; n1[k] <- length(y[y > 0])
    
    # Fit models
    models <-     c( "LO", "GU",  "NO",  "TF",  'LQNO')#,   'GT'
    nparms <- c(2, 2, 2, 3,  2)#  4
    mdl1 <- suppressWarnings(gamlss::gamlss(y~1, family=LO, control = gamlss.control(trace = FALSE)))  # logistics distribution
    mdl2 <- suppressWarnings(gamlss::gamlss(y~1,sigma.formula=~1,  family=GU, control = gamlss.control(trace = FALSE))) # gumbell distribution
    mdl3 <- suppressWarnings(gamlss::gamlss(y~1,sigma.formula=~1,  family=NO, control = gamlss.control(trace = FALSE))) # Normal
    mdl4 <- suppressWarnings(gamlss::gamlss(y~1,sigma.formula=~1,  family=TF, control = gamlss.control(trace = FALSE))) # t family distribution for fitting a GAMLSS
    mdl5 <- suppressWarnings(gamlss::gamlss(y~1,sigma.formula=~1,  family=LQNO, control = gamlss.control(trace = FALSE))) # t family distribution for fitting a GAMLSS
    #mdl6 <- suppressWarnings(gamlss(y~1,sigma.formula=~1,  family=GT, control = gamlss.control(trace = FALSE))) # generalized t distribution for fitting a GAMLSS

    
    
    # Put these models in a list
    mdls <- list(mdl1, mdl2,  mdl3, mdl4,   mdl5)#,  mdl6
    
    if(min(y) > 0){

      models2 <-  c("LNO", "EXP",  "GA", "WEI", 'IGAMMA')
      nparms2 <- c(2, 1, 2, 2, 2)
      mdl7 <- suppressWarnings(gamlss::gamlss(y~1,sigma.formula=~1,  family=LNO, control = gamlss.control(trace = FALSE))) # Log normal distribution
      mdl8 <- suppressWarnings(gamlss::gamlss(y~1,sigma.formula=~1,  family=EXP, control = gamlss.control(trace = FALSE))) # Exponential
      mdl9 <- suppressWarnings(gamlss::gamlss(y~1,sigma.formula=~1,  family=GA, control = gamlss.control(trace = FALSE))) # Gamma distribution for fitting a GAMLSS
      mdl10 <- suppressWarnings(gamlss::gamlss(y~1,sigma.formula=~1,  family=WEI, control = gamlss.control(trace = FALSE))) # Weibull distribution  for fitting a GAMLSS
      mdl11 <- suppressWarnings(gamlss::gamlss(y~1,sigma.formula=~1,  family=IGAMMA, control = gamlss.control(trace = FALSE))) # Weibull distribution  for fitting a GAMLSS


      mdls2 <- list(mdl7, mdl8,  mdl9, mdl10, mdl11)

      models <- c(models, models2)
      nparms <- c(nparms, nparms2)
      mdls <- c(mdls, mdls2)
    }else{

      models <- models
      nparms <- nparms
      mdls <- mdls
    }
    
    names(mdls) <- models
    
    loglik <- rep(0, length(nparms))
    
    for (i in 1:length(nparms)) {
      
      loglik[i] = logLik(mdls[[i]])
      
    }
    
    AIC <- 2*nparms - 2*loglik
    AICc <- AIC + 2*nparms*(nparms+1)/(n-nparms-1)
    BIC <- nparms*log(n) - 2*loglik
    BICc <-  - 2*loglik + (nparms * log(n) * n) / (n - nparms - 1)
    min.AIC <- min(AIC)
    
    weights <- exp((min.AIC - AIC)/2)
    results <- cbind(nparms, loglik, AIC,AICc,BIC, BICc,weights,irank(criterion))
    rownames(results) <- models
    colnames(results) <- c("nparms","loglik","AIC","AICc","BIC", "BICc", "weight",paste("rank",criterion,sep="."))
    
    results.ord <- results[order(irank(criterion)),]

    detail[[k]] <- results.ord
    
    dist[k] <- rownames(results.ord)[1]
    
    if(dist[k] == "LO") {location_param[k] <- fitted(mdl1,"mu")[1]; scale_param[k] <- fitted(mdl1,"sigma")[1]}
    if(dist[k] == "GU") {location_param[k] <- fitted(mdl2,"mu")[1] ; scale_param[k] <- fitted(mdl2,"sigma")[1]}
    if(dist[k] == "NO") {location_param[k] <- fitted(mdl3,"mu")[1] ; scale_param[k] <- fitted(mdl3,"sigma")[1]}
    if(dist[k] == "TF") {location_param[k] <- fitted(mdl4,"mu")[1] ; scale_param[k] <- fitted(mdl4,"sigma")[1]; nu_param[k] = fitted(mdl4,"nu")[1]}
    #if(dist[k] == "LQNO") {location_param[k] <- fitted(mdl5,"mu")[1] ; scale_param[k] <- fitted(mdl5,"sigma")[1]}
   
    if(dist[k] == "GT") {location_param[k] <- fitted(mdl6,"mu")[1] ; scale_param[k] <- fitted(mdl6,"sigma")[1]; nu_param[k] = fitted(mdl6,"nu")[1]; tau_param[k] = fitted(mdl6,"tau")[1]}
    if( min(y) > 0){
      if(dist[k] == "LNO") {location_param[k] <- fitted(mdl7,"mu")[1] ; scale_param[k] <- fitted(mdl7,"sigma")[1]}
      if(dist[k] == "EXP") {location_param[k] <- fitted(mdl8,"mu")[1]}
      if(dist[k] == "GA") {location_param[k] <- fitted(mdl9,"mu")[1] ; scale_param[k] <- fitted(mdl9,"sigma")[1]}
      if(dist[k] == "WEI") {location_param[k] <- fitted(mdl10,"mu")[1] ; scale_param[k] <- fitted(mdl10,"sigma")[1]}
      if(dist[k] == "IGAMMA") {location_param[k] <- fitted(mdl11,"mu")[1] ; scale_param[k] <- fitted(mdl11,"sigma")[1]}


    }
    
  }
  
  distr <- data.frame(colnames(Y), dist, location_param, scale_param, nu_param, tau_param)
  
  colnames(distr) <- c("Variable", "distrib","mu","sigma", "nu", 'tau')
  
  output <- list(distr, detail)
  
  names(output) <- c("distrib", "detail")
  
  return(output)
}


marginal_parameters <- function(name, parameters) {
  
  if (toupper(name) == "LO") {
    return(list(distr = "LO", parameters = list(mu = parameters[[1]],
                                                sigma = parameters[[2]])))
  }
  else if (toupper(name) == "LNO") {
    return(list(distr = "LNO", parameters = list(mu = parameters[[1]],
                                                  sigma = parameters[[2]])))
  }else if (toupper(name) == "EXP") {
    return(list(distr = "EXP", parameters = list(mu = parameters[[1]])))
  }else if (toupper(name) == "GU") {
    return(list(distr = "GU", parameters = list(mu = parameters[[1]],
                                                 sigma = parameters[[2]])))
  }else if (toupper(name) == "GA") {
    return(list(distr = "GA", parameters = list(mu = parameters[[1]],
                                                 sigma = parameters[[2]])))
  }else if (toupper(name) == "IGAMMA") {
    return(list(distr = "IGAMMA", parameters = list(mu = parameters[[1]],
                                                sigma = parameters[[2]])))
  }else if (toupper(name) == "NO") {
    return(list(distr = "NO", parameters = list(mu = parameters[[1]],
                                                 sigma = parameters[[2]])))
  }else if (toupper(name) == "TF") {
    return(list(distr = "TF", parameters = list(mu = parameters[[1]],
                                                 sigma = parameters[[2]],
                                                nu = parameters[[3]])))
  }
  else if (toupper(name) == "WEI") {
    return(list(distr = "WEI", parameters = list(mu = parameters[[1]],
                                                 sigma = parameters[[2]])))
  }else if (toupper(name) == "LQNO") {
    return(list(distr = "LQNO", parameters = list(mu = parameters[[1]],
                                                 sigma = parameters[[2]])))
  }
  # else if (toupper(name) == "GT") {
  #   return(list(distr = "GT", parameters = list(mu = parameters[[1]],
  #                                                sigma = parameters[[2]],
  #                                             nu = parameters[[3]],
  #                                             tau = parameters[[4]])))
  # }
  else {
    stop(sprintf("%s distribution type is not yet implemented.", name))
  }
}



 
 
marginal_parameters_splits_AND_distri <- function(Y,  marginal_distr) {

  n_species <- dim(Y)[2]

  marginals <- as.list(rep(NA, n_species))

  marginal_distributions = matrix(NA, nrow = nrow(Y), ncol = n_species)

  for (j in 1:n_species) {

    marginals[[j]] <- marginal_parameters(marginal_distr[j, "distrib"],

                                          marginal_distr[j, c("mu","sigma", "nu", 'tau')])

    if(marginals[[j]]$distr == "LO"){

      marginal_distributions[,j] = pLO(Y[,j],

                                         mu = marginals[[j]]$parameters$mu,

                                       sigma = marginals[[j]]$parameters$sigma)

    }
    else if(marginals[[j]]$distr == "LNO"){
      marginal_distributions[,j] = pLNO(Y[,j],

                                         mu = marginals[[j]]$parameters$mu,

                                         sigma = marginals[[j]]$parameters$sigma)

    }else if(marginals[[j]]$distr == "EXP"){
      marginal_distributions[,j] = pEXP(Y[,j],

                                        mu = marginals[[j]]$parameters$mu)

    }else if(marginals[[j]]$distr == "GU"){
      marginal_distributions[,j] = pGU(Y[,j],

                                        mu = marginals[[j]]$parameters$mu,

                                        sigma = marginals[[j]]$parameters$sigma)

    }
    else if(marginals[[j]]$distr == "GA"){
      marginal_distributions[,j] = pGA(Y[,j],

                                        mu = marginals[[j]]$parameters$mu,

                                        sigma = marginals[[j]]$parameters$sigma)

    }else if(marginals[[j]]$distr == "IGAMMA"){
      marginal_distributions[,j] = pIGAMMA(Y[,j],
                                       
                                       mu = marginals[[j]]$parameters$mu,
                                       
                                       sigma = marginals[[j]]$parameters$sigma)
      
    }
    else if(marginals[[j]]$distr == "NO"){
      marginal_distributions[,j] = pNO(Y[,j],

                                        mu = marginals[[j]]$parameters$mu,

                                         sigma = marginals[[j]]$parameters$sigma)

    }
    else if(marginals[[j]]$distr == "TF"){
      marginal_distributions[,j] = pTF(Y[,j],

                                        mu = marginals[[j]]$parameters$mu,

                                        sigma = marginals[[j]]$parameters$sigma,
                                        nu = marginals[[j]]$parameters$nu)

    }

    
    
    
    else if(marginals[[j]]$distr == "WEI"){
      marginal_distributions[,j] = pWEI(Y[,j],

                                        mu = marginals[[j]]$parameters$mu,

                                        sigma = marginals[[j]]$parameters$sigma)

    }
    else if(marginals[[j]]$distr == "LQNO"){
      marginal_distributions[,j] = pLQNO(Y[,j],

                                        mu = marginals[[j]]$parameters$mu,

                                        sigma = marginals[[j]]$parameters$sigma)

    }
    # else if(marginals[[j]]$distr == "GT"){
    #   marginal_distributions[,j] = pGT(Y[,j],
    # 
    #                                     mu = marginals[[j]]$parameters$mu,
    # 
    #                                     sigma = marginals[[j]]$parameters$sigma,
    #                                     nu = marginals[[j]]$parameters$nu,
    #                                     tau = marginals[[j]]$parameters$tau)
    # 
    # }
    else {
      stop(sprintf("%s distribution type is not yet implemented.", name))
    }
  }

  margins = c( ) ;   paramMargins = list()

  for(j in 1:n_species){

    margins[j] = marginals[[j]]$distr

    paramMargins[[j]] = marginals[[j]]$parameters

  }

  colnames(marginal_distributions) = colnames(Y)

  return(list(margins =  margins,

              paramMargins = paramMargins,

              marginal_distributions = marginal_distributions))
}








QDistr <- function(Y,  marginal_distr) {
  
  n_species <- dim(Y)[2]
  
  marginals <- as.list(rep(NA, n_species))
  
  marginal_distributions = matrix(NA, nrow = nrow(Y), ncol = n_species)
  
  for (j in 1:n_species) {
    
    marginals[[j]] <- marginal_parameters(marginal_distr[j, "distrib"],
                                          
                                          marginal_distr[j, c("mu","sigma", "nu", "tau")])
    
    if(marginals[[j]]$distr == "LO"){
      
      marginal_distributions[,j] = qLO(Y[,j],
                                       
                                       mu = marginals[[j]]$parameters$mu,
                                       
                                       sigma = marginals[[j]]$parameters$sigma)
      
    }
    else if(marginals[[j]]$distr == "LNO"){
      marginal_distributions[,j] = qLNO(Y[,j],
                                        
                                        mu = marginals[[j]]$parameters$mu,
                                        
                                        sigma = marginals[[j]]$parameters$sigma)
      
    }else if(marginals[[j]]$distr == "EXP"){
      marginal_distributions[,j] = qEXP(Y[,j],
                                        
                                        mu = marginals[[j]]$parameters$mu)
      
    }else if(marginals[[j]]$distr == "GU"){
      marginal_distributions[,j] = qGU(Y[,j],
                                       
                                       mu = marginals[[j]]$parameters$mu,
                                       
                                       sigma = marginals[[j]]$parameters$sigma)
      
    }
    else if(marginals[[j]]$distr == "GA"){
      marginal_distributions[,j] = qGA(Y[,j],
                                       
                                       mu = marginals[[j]]$parameters$mu,
                                       
                                       sigma = marginals[[j]]$parameters$sigma)
      
    }else if(marginals[[j]]$distr == "IGAMMA"){
      marginal_distributions[,j] = qIGAMMA(Y[,j],
                                           
                                           mu = marginals[[j]]$parameters$mu,
                                           
                                           sigma = marginals[[j]]$parameters$sigma)
      
    }
    else if(marginals[[j]]$distr == "NO"){
      marginal_distributions[,j] = qNO(Y[,j],
                                       
                                       mu = marginals[[j]]$parameters$mu,
                                       
                                       sigma = marginals[[j]]$parameters$sigma)
      
    }
    else if(marginals[[j]]$distr == "TF"){
      marginal_distributions[,j] = qTF(Y[,j],

                                       mu = marginals[[j]]$parameters$mu,

                                       sigma = marginals[[j]]$parameters$sigma,
                                       nu = marginals[[j]]$parameters$nu)

    }
    else if(marginals[[j]]$distr == "ST1"){
      marginal_distributions[,j] = qST1(Y[,j],
                                       
                                       mu = marginals[[j]]$parameters$mu,
                                       
                                       sigma = marginals[[j]]$parameters$sigma,
                                       nu = marginals[[j]]$parameters$nu,
                                       tau = marginals[[j]]$parameters$tau)
      
    }else if(marginals[[j]]$distr == "ST2"){
      marginal_distributions[,j] = qST2(Y[,j],
                                        
                                        mu = marginals[[j]]$parameters$mu,
                                        
                                        sigma = marginals[[j]]$parameters$sigma,
                                        nu = marginals[[j]]$parameters$nu,
                                        tau = marginals[[j]]$parameters$tau)
      
    }
    else if(marginals[[j]]$distr == "ST3"){
      marginal_distributions[,j] = qST3(Y[,j],
                                        
                                        mu = marginals[[j]]$parameters$mu,
                                        
                                        sigma = marginals[[j]]$parameters$sigma,
                                        nu = marginals[[j]]$parameters$nu,
                                        tau = marginals[[j]]$parameters$tau)
      
    }
    else if(marginals[[j]]$distr == "ST4"){
      marginal_distributions[,j] = qST4(Y[,j],
                                        
                                        mu = marginals[[j]]$parameters$mu,
                                        
                                        sigma = marginals[[j]]$parameters$sigma,
                                        nu = marginals[[j]]$parameters$nu,
                                        tau = marginals[[j]]$parameters$tau)
      
    }
    else if(marginals[[j]]$distr == "ST5"){
      marginal_distributions[,j] = qST5(Y[,j],
                                        
                                        mu = marginals[[j]]$parameters$mu,
                                        
                                        sigma = marginals[[j]]$parameters$sigma,
                                        nu = marginals[[j]]$parameters$nu,
                                        tau = marginals[[j]]$parameters$tau)
      
    } else if(marginals[[j]]$distr == "WEI"){
      marginal_distributions[,j] = qWEI(Y[,j],
                                        
                                        mu = marginals[[j]]$parameters$mu,
                                        
                                        sigma = marginals[[j]]$parameters$sigma)
      
    }
    else if(marginals[[j]]$distr == "LQNO"){
      marginal_distributions[,j] = qLQNO(Y[,j],
                                         
                                         mu = marginals[[j]]$parameters$mu,
                                         
                                         sigma = marginals[[j]]$parameters$sigma)
      
    }
    else if(marginals[[j]]$distr == "GT"){
      marginal_distributions[,j] = qGT(Y[,j],

                                        mu = marginals[[j]]$parameters$mu,

                                        sigma = marginals[[j]]$parameters$sigma,
                                        nu = marginals[[j]]$parameters$nu,
                                        tau = marginals[[j]]$parameters$tau)

    }
    else {
      stop(sprintf("%s distribution type is not yet implemented.", name))
    }
  }
  
  margins = c( ) ;   paramMargins = list()
  
  for(j in 1:n_species){
    
    margins[j] = marginals[[j]]$distr
    
    paramMargins[[j]] = marginals[[j]]$parameters
    
  }
  
  colnames(marginal_distributions) = colnames(Y)
  
  return(list(margins =  margins,
              
              paramMargins = paramMargins,
              
              marginal_distributions = marginal_distributions))
}

