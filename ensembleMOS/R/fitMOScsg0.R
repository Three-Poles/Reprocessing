fitMOScsg0 <-
function(ensembleData, control = controlMOScsg0(), exchangeable = NULL)
{
  if(!is.integer(control$maxIter)){
      control$maxIter <- as.integer(1e7)
      warning("maxIter not an integer")
  }
  na.rm <- FALSE
  ### Functions for imputing missing values
  "rmNArows" <- function(x)
  {
   forecasts <- ensembleForecasts(x)
   M <- apply(forecasts,1,function(z) any(is.na(z)))
   x <- x[!M,]
   return(x)
  }
  "impute" <- function(x)
  {
    M <- apply(x, 1, anyNA)
    for(row in which(M)){
      x[row,which(is.na(x[row,]))  ] <- mean(x[row,], na.rm = TRUE)
    }
    return(x)
  }
  "meanNA" <- function(x){
      v <- mean(x, na.rm = TRUE)
  }
  
  # remove instances missing all forecasts or obs
  M <- apply(ensembleForecasts(ensembleData), 1, function(z) all(is.na(z)))
  M <- M | is.na(ensembleVerifObs(ensembleData))
  ensembleData <- ensembleData[!M,]

  if (is.null(obs <- ensembleVerifObs(ensembleData)))
      stop("verification observations required")

  nObs <- ensembleNobs(ensembleData)
  ensMemNames <- ensembleMemberLabels(ensembleData)
  nForecasts <- length(ensMemNames)
  forecasts <- ensembleForecasts(ensembleData)

  if (is.null(exchangeable)) exchangeable <- ensembleGroups(ensembleData)

  if (length(unique(exchangeable)) == length(exchangeable))
      exchangeable <- NULL
  nEX <- 2
  if (!(nullX <- is.null(exchangeable))) {
      ### Helper function to name exchangeable Groups
      "catStrings" <- function(strings) {
          if(length(strings) == 1) return(strings)
          if(is.null(strings[1]) || !is.character(strings)) {
          stop("Not a valid string")
      }
      newstrings <- strings[1]
      for(i in 2:length(strings)) {
          newstrings <- paste(newstrings, strings[i], sep = ".")
      }
      newstrings
    }

    namEX <- as.character(exchangeable)
    uniqueEX <- unique(namEX)
    nEX <- length(uniqueEX)
    splitEX <- split(seq(along = exchangeable), exchangeable)
    groupNames <- NULL
    for(i in 1:nEX) {
        groupNames <- c(groupNames,
                        catStrings(ensMemNames[splitEX[[i]]]))
    }

    if(!na.rm) {
     for(i in 1:nEX) {
         if(length(splitEX[[i]]) > 1) {
             M <- apply(forecasts[,splitEX[[i]]],1,function(z) all(is.na(z)))
             forecasts[!M,splitEX[[i]]] <- impute(forecasts[!M,splitEX[[i]]])
         }
     }
     D <- apply(forecasts, 1, function(z) !any(is.na(z)))
    }else D <- 1:nObs

    nObs <- dim(forecasts[D,])[1]
    matEX <- matrix(NA, nObs , ncol = nEX)
    dimnames(matEX) <- list(row.names(forecasts[D,]), groupNames)
    for(k in 1:nEX) {
      if(length(splitEX[[k]]) > 1) {
       matEX[,k] <- apply(forecasts[D,splitEX[[k]]], 1, mean)
      }else{
       matEX[,k] <- forecasts[D, splitEX[[k]]] }
    }

    ensembleEX <- cbind(matEX,
                        dates = ensembleData$dates[D],
                        observations = ensembleData$obs[D])


  }
  if(nullX) {
    fitData <- rmNArows(ensembleData)
  }else{

    fcstEX <- data.frame(matEX)
    fcstHour <- ensembleFhour(ensembleData)
    initTime <- ensembleItime(ensembleData)
    fitData <- ensembleData(forecasts = fcstEX,
                            dates = ensembleData$dates[D],
                            observations = ensembleData$obs[D],
                            forecastHour = fcstHour,
                            initializationTime = initTime)
  }
  fitForecasts <- ensembleForecasts(fitData)
  obs <- ensembleVerifObs(fitData)
  xTrain <- cbind(rep(1,length(obs)),fitForecasts)
  if(nullX){		
      var <- apply(fitForecasts,1,meanNA)
  }else{
      var <- apply(ensembleForecasts(ensembleData)[D,],1,meanNA)  
  }

  if(any(is.null(c(control$start$q, control$start$c, control$start$d)))) {
      control$start$c <- 5
      control$start$d <- 1
      control$start$q <- 1
  }
  
  pozobs <- (obs>0)
  
  if (any(is.null(c(control$start$B,control$start$a)))) {
    if (sum(pozobs)> 2*(nEX+1)) {
    olsCoefs <- lm(obs[pozobs]~ensembleForecasts(fitData[pozobs,]))$coef
    control$start$a <- olsCoefs[1]
    B <- olsCoefs[-1]    
    if(control$coefRule == "square") {
    control$start$B <- sqrt(abs(B))
    } else control$start$B <- B} 
    else 
      {control$start$a <- 1
      control$start$B <- rep(1,nEX)}
  }
  names(control$start$B) <- if(nullX) ensMemNames else groupNames
  c <- control$start$c
  if(control$varRule == "square") c <- sqrt(c)
  pars <- c(sqrt(control$start$q), c, sqrt(control$start$d), control$start$a, control$start$B)

# optimize pars on chosen scoring rule

  
  switch(control$scoringRule,
         crps = {
  objectiveFUN <- function(pars) {
      if(control$coefRule == "square") {  
        pars[-c(1:4)] <- pars[-c(1:4)]^2
      }
      if(control$varRule == "square") {
        pars[2] <- pars[2]^2
      }
      mu <- xTrain%*%pars[-c(1:3)]  #mean of gamma
      sig.Sq <- pars[2] + pars[3]^2*var  #variance of gamma
      Shp <- mu^2/sig.Sq  #shape of gamma
      Scl <- sig.Sq/mu    #scale of gamma
      
      Z <- (obs + pars[1]^2)/Scl
      C <- pars[1]^2/Scl 
      
      crps <- Scl*Z*(2*pgamma(Z,Shp,1)-1)-Scl*C*(pgamma(C,Shp,1))^2 + mu*(1+2*pgamma(C,Shp,1)*pgamma(C,Shp+1,1)
	      -(pgamma(C,Shp,1))^2-2*pgamma(Z,Shp+1,1)) - mu*(1-pgamma(2*C,2*Shp,1))*beta(.5,Shp+.5)/pi

      ###Formally "mean" but sum and mean produce different results.
      sum(crps)
  }
                  },
         log = {
  objectiveFUN <- function(pars) { 
      if(control$coefRule == "square") {  
        pars[-c(1:4)] <- pars[-c(1:4)]^2
      }
      if(control$varRule == "square") {
        pars[2] <- pars[2]^2
      }
      mu <- xTrain%*%pars[-c(1:3)]  #mean of gamma
      sig.Sq <- pars[2] + pars[3]^2*var  #variance of gamma
      Shp <- mu^2/sig.Sq  #shape of gamma
      Scl <- sig.Sq/mu    #scale of gamma
      
      Z <- (obs + pars[1]^2)/Scl
      C <- pars[1]^2/Scl     
      
      ign <- -log(pgamma(C,Shp,1)) * (obs==0) -log(dgamma(Z,Shp,1)/Scl) * (obs>0)
     
      return(sum(ign,na.rm=T))
      
  }
                  }
  )
  
 
  switch(control$scoringRule,
         crps = {
         opt <- optim(pars, fn=objectiveFUN, method="L-BFGS-B", control=list(maxit=control$maxIter),lower=rep(0.0001,(4+nEX)),upper=rep(20,(4+nEX)))
         },
         log = {
            if (control$optimRule == 'L-BFGS-B'){
              opt <- optim(pars, fn=objectiveFUN, method="L-BFGS-B", control=list(maxit=control$maxIter),lower=rep(0.0001,(4+nEX)),upper=rep(20,(4+nEX)))
            } else {
              opt <- optim(pars, fn=objectiveFUN, method=control$optimRule, control=list(maxit=control$maxIter))
            }
          }
         )

  optB <- opt$par[-c(1:4)]
  B <- rep(NA,nForecasts)
  if(nullX) {
      B <- optB
  } else{
   for(i in 1:nEX) {
      B[splitEX[[i]]] <- optB[i]
   }
  }
  names(B) <- ensMemNames
  if (control$coefRule == "square") {
      B <- B^2
      if(!nullX){
          for(i in 1:nEX){
              B[splitEX[[i]]] <- B[splitEX[[i]]]/length(splitEX[[i]])
          }
      }
      if (control$varRule == "square") {
          opt$par[2] <- opt$par[2]^2
      }
      fit <- structure(list(a = opt$par[4], B = B,
                            c = opt$par[2], d = opt$par[3]^2, q = opt$par[1]^2,
                            exhangeable = exchangeable),
                       class = "fitMOScsg0")
      return(fit)
  } else if (control$coefRule == "positive") {
      if(control$varRule == "positive") {
          opt$par[2] <- opt$par[2]^2
      }
      fit <- structure( list(a = opt$par[4], B = B, 
                             c = opt$par[2], d = opt$par[3]^2, q = opt$par[1]^2,
                             exhangeable = exchangeable),
                       class = "fitMOScsg0")
      posfit <- fit
      posEnsNames <- ensMemNames[posfit$B>0]
      control$coefRule <- "none"
      posControl <- control
      posControl$start$B <- control$start$B[posfit$B>0]
      while(any(posfit$B<0)) {
          posExchangeable <- exchangeable[posEnsNames]
          posEnsembleData <- ensembleData(forecasts = ensembleForecasts(ensembleData)[,posEnsNames], 
                                          dates = ensembleData$dates, 
                                          observations = ensembleData$observations,
                                          station = ensembleData$station,
                                          forecastHour = ensembleData$forecastHour, 
                                          initializationTime = ensembleData$initializationTime,
                                          exchangeable = posExchangeable)
          posfit <- fitMOScsg0(posEnsembleData,
                                 control = posControl,
                                 exchangeable = posExchangeable)
          posEnsNames <- posEnsNames[posfit$B>0]
          if (any(is.na(control$start$B[posEnsNames]))) browser()
          posControl$start$B <- control$start$B[posEnsNames]
      }
      fit$B <- rep(0,nForecasts)
      names(fit$B) <- ensMemNames
      fit$B[posEnsNames] <- posfit$B
      fit$a <- posfit$a
      fit$c <- posfit$c
      fit$d <- posfit$d
      fit$q <- posfit$q
      return(fit)
  }else{
      if (control$varRule == "square") {
          opt$par[2] <- opt$par[2]^2
      }
      fit <- structure( list(a = opt$par[4], B = B, 
                             c = opt$par[2], d = opt$par[3]^2, q = opt$par[1]^2,
                             exhangeable = exchangeable),
                       class = "fitMOScsg0")
      return(fit)
  }
}

