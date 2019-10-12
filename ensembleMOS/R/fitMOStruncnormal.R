fitMOStruncnormal <-
function(ensembleData, control = controlMOStruncnormal(), exchangeable = NULL)
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
  "varNA" <- function(x){
      v <- var(x, na.rm = TRUE)
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
  
  if (nullX) {
    var <- apply(fitForecasts,1,varNA)
  }else{
    var <- apply(ensembleForecasts(ensembleData)[D,],1,varNA)
  }
  
 
  if(any(is.null(c(control$start$c, control$start$d)))) {
      control$start$c <- 5
      control$start$d <- 1
  }
  if (any(is.null(c(control$start$B,control$start$a)))) {
    olsCoefs <- lm(obs~ensembleForecasts(fitData))$coef
    control$start$a <- olsCoefs[1]
    B <- olsCoefs[-1]
    if(control$coefRule == "square") {
    control$start$B <- sqrt(abs(B))
    } else control$start$B <- B
  }
  names(control$start$B) <- if(nullX) ensMemNames else groupNames
  c <- control$start$c
  if(control$varRule == "square") c <- sqrt(c)
  pars <- c(c, sqrt(control$start$d), control$start$a, control$start$B)

# optimize pars on chosen scoring rule
  
  switch(control$scoringRule,
         crps = {
  objectiveFUN <- function(pars) {  
      if(control$coefRule == "square") {  
        pars[-c(1:3)] <- pars[-c(1:3)]^2
      }
      if(control$varRule == "square") {
        pars[1] <- pars[1]^2
      }
      mu <- xTrain%*%pars[-c(1:2)]
      sig <- sqrt(pars[1] + pars[2]^2*var)
      z <- (obs - mu)/sig
      w <- mu/sig
      crps <- sig*(z*pnorm(w)*(2*pnorm(z)+pnorm(w)-2) + 2*dnorm(z)*pnorm(w) - pnorm(w*sqrt(2))/sqrt(pi))/(pnorm(w))^2

      ###Formally "mean" but sum and mean produce different results.
      sum(crps)
  }
                  },

         log = {
  objectiveFUN <- function(pars) {
      if(control$coefRule == "square") {
       pars[-c(1:3)] <- pars[-c(1:3)]^2
      }
      if(control$varRule == "square") {
       pars[1] <- pars[1]^2
      }
      mu <- xTrain%*%c(pars[3],pars[-c(1:3)])
      sig.sq <- pars[1] + pars[2]^2*var
      ign <- 0.5*log(2*pi*sig.sq) + ((obs - mu)^2)/(2*sig.sq)+log(pnorm(mu/sqrt(sig.sq)))
      mean(ign)
  }
                  }
  )

  opt <- optim(pars, fn=objectiveFUN, method=control$optimRule, control=list(maxit=control$maxIter))

  optB <- opt$par[-c(1:3)]
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
          opt$par[1] <- opt$par[1]^2
      }
      fit <- structure(list(
                            a = opt$par[3], B = B,
                            c = opt$par[1], d = opt$par[2]^2,
                            exhangeable = exchangeable),
                       class = "fitMOStruncnormal")
      return(fit)
  } else if (control$coefRule == "positive") {
      if(control$varRule == "positive") {
          opt$par[1] <- opt$par[1]^2
      }
      fit <- structure( list(
                             a = opt$par[3], B = B,
                             c = opt$par[1], d = opt$par[2]^2,
                             exhangeable = exchangeable),
                       class = "fitMOStruncnormal")
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
          posfit <- fitMOStruncnormal(posEnsembleData,
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
      return(fit)
  }else{
      if (control$varRule == "square") {
          opt$par[1] <- opt$par[1]^2
      }
      fit <- structure( list(
                             a = opt$par[3], B = B,
                             c = opt$par[1], d = opt$par[2]^2,
                             exhangeable = exchangeable),
                       class = "fitMOStruncnormal")
      return(fit)
  }
}

