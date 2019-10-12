fitMOSgev0 <-
function(ensembleData, control = controlMOSgev0(), exchangeable = NULL)
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
  
  "gini.md" <- function(x,na.rm=FALSE)  {     ## Michael Scheuerer's code
   if(na.rm & any(is.na(x)))  x <- x[!is.na(x)] 
   n <-length(x)
   return(4*sum((1:n)*sort(x,na.last=TRUE))/(n^2)-2*mean(x)*(n+1)/n)
    }
    
  "crps.GEVneq0" <- function(pars, obs, xTrain, p0, md)  { ## Michael Scheuerer's code
  
      
     
      MEAN <- xTrain%*%pars[-c(1:4)]+pars[2]*p0  #mean of GEV
      SCALE <- pars[3] + pars[4]^2*md  #scale of GEV
      SHAPE <- pars[1]
      
      LOC <- MEAN - SCALE*(gamma(1-SHAPE)-1)/SHAPE #location of GEV

      SCdSH <- SCALE/SHAPE
      Gam1mSH <- gamma(1-SHAPE)
      prob0 <- pgev(0, loc=LOC, scale=SCALE, shape=SHAPE)
      probY <- pgev(obs, loc=LOC, scale=SCALE, shape=SHAPE)

      T1 <- (obs-LOC)*(2*probY-1) + LOC*prob0^2
      T2 <- SCdSH * ( 1-prob0^2 - 2^SHAPE*Gam1mSH*pgamma(-2*log(prob0),1-SHAPE) )
      T3 <- -2*SCdSH * ( 1-probY - Gam1mSH*pgamma(-log(probY),1-SHAPE) )
      
      return( mean(T1+T2+T3) )
}  

  eps <- 1e-5
  
  
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
      var <- apply(fitForecasts,1,varNA)
      md <- apply(fitForecasts,1,gini.md,na.rm = TRUE)
      p0 <- rowMeans(fitForecasts==0,na.rm = TRUE)
  }else{
      var <- apply(ensembleForecasts(ensembleData)[D,],1,varNA) 
      md <- apply(ensembleForecasts(ensembleData)[D,],1,gini.md,na.rm=TRUE)
      p0 <- rowMeans(ensembleForecasts(ensembleData)[D,]==0, na.rm = TRUE)      
  }
 
  
  if(any(is.null(c(control$start$q, control$start$s, control$start$c, control$start$d)))) {
      control$start$c <- .1
      control$start$d <- 1
      control$start$q <- .1
      control$start$s <- 1
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
  pars <- c(control$start$q, control$start$s, c, sqrt(control$start$d), control$start$a, control$start$B)

# optimize pars on crps

  
  objectiveFUN <- function(param) { 
      if(control$coefRule == "square") {  
        param[-c(1:5)] <- param[-c(1:5)]^2
      }
      if(control$varRule == "square") {
        param[3] <- param[3]^2
      }
      
     if( abs(param[1]) < eps )  {
      crps.eps.m <- crps.GEVneq0 (c(-eps,param[-1]), obs, xTrain, p0, md)
      crps.eps.p <- crps.GEVneq0 (c(eps,param[-1]), obs, xTrain, p0, md)
      w.m <- (eps-param[1])/(2*eps)
      w.p <- (eps+param[1])/(2*eps)
      res <- w.m*crps.eps.m + w.p*crps.eps.p
   } else  {
      res <- crps.GEVneq0 (param, obs, xTrain, p0, md)
   }
  
   return(res) 
     
  }
 
 
  if (control$optimRule == 'L-BFGS-B'){
    opt <- optim(pars, fn=objectiveFUN, method="L-BFGS-B", control=list(maxit=control$maxIter),lower=c(-0.278,-10,rep(0.0001,(3+nEX))), upper=c(0.99999,10,rep(10,(3+nEX))))
    } else {
    opt <- optim(pars, fn=objectiveFUN, method=control$optimRule, control=list(maxit=control$maxIter))
  }

  optB <- opt$par[-c(1:5)]
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
          opt$par[3] <- opt$par[3]^2
      }
      fit <- structure(list(a = opt$par[5], B = B, s = opt$par[2],
                            c = opt$par[3], d = opt$par[4]^2, q = opt$par[1], 
                            exhangeable = exchangeable),
                       class = "fitMOSgev0")
      return(fit)
  } else if (control$coefRule == "positive") {
      if(control$varRule == "positive") {
          opt$par[3] <- opt$par[3]^2
      }
      fit <- structure( list(a = opt$par[5], B = B, s = opt$par[2],
                             c = opt$par[3], d = opt$par[4]^2, q = opt$par[1],
                             exhangeable = exchangeable),
                       class = "fitMOSgev0")
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
          posfit <- fitMOSgev0(posEnsembleData,
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
      fit$s <- posfit$s
      return(fit)
  }else{
      if (control$varRule == "square") {
          opt$par[3] <- opt$par[3]^2
      }
      fit <- structure( list(a = opt$par[5], B = B, s = opt$par[2],
                             c = opt$par[3], d = opt$par[4]^2, q = opt$par[1],
                             exhangeable = exchangeable),
                       class = "fitMOSgev0")
      return(fit)
  }
}

