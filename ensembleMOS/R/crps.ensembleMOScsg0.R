crps.ensembleMOScsg0 <-
function(fit, ensembleData, dates=NULL, ...)
{

 crpsFunc <- function(mu, sig.Sq,qval,y)
  {
   
      Shp <- mu^2/sig.Sq  #shape of gamma
      Scl <- sig.Sq/mu    #scale of gamma
      
      if (!is.na(Scl) & (Scl <= 0))  stop("scale of gamma distribution is not positive")
      
      Z <- (y + qval)/Scl
      C <- qval/Scl 
 
      crps <- Scl*Z*(2*pgamma(Z,Shp,1)-1)-Scl*C*(pgamma(C,Shp,1))^2 + Scl*Shp*(1+2*pgamma(C,Shp,1)*pgamma(C,Shp+1,1)
	      -(pgamma(C,Shp,1))^2-2*pgamma(Z,Shp+1,1)) - Scl*Shp*(1-pgamma(2*C,2*Shp,1))*beta(.5,Shp+.5)/pi
   
    crps
  }
 
 matchITandFH(fit,ensembleData)

 M <- matchEnsembleMembers(fit,ensembleData)
 nForecasts <- ensembleSize(ensembleData)
 if (!all(M == 1:nForecasts)) ensembleData <- ensembleData[,M]

## remove instances missing all forecasts

 M <- apply(ensembleForecasts(ensembleData), 1, function(z) all(is.na(z)))
 M <- M | is.na(ensembleVerifObs(ensembleData))
 ensembleData <- ensembleData[!M,]

## match specified dates with dateTable in fit

 dateTable <- dimnames(fit$B)[[2]]

 if (!is.null(dates)) {

   dates <- sort(unique(as.character(dates)))

   if (length(dates) > length(dateTable))
     stop("parameters not available for some dates")

   K <- match( dates, dateTable, nomatch=0)

   if (any(!K) || !length(K))
     stop("parameters not available for some dates")

 }
 else {

   dates <- dateTable
   K <- 1:length(dateTable)

  }

 ensDates <- ensembleValidDates(ensembleData)

## match dates in data with dateTable
 if (is.null(ensDates) || all(is.na(ensDates))) {
   if (length(dates) > 1) stop("date ambiguity")
   nObs <- nrow(ensembleData)
   Dates <- rep( dates, nObs)
 }
 else {
## remove instances missing dates
   if (any(M <- is.na(ensDates))) {
     ensembleData <- ensembleData[!M,]
     ensDates <- ensembleValidDates(ensembleData)
   }
   Dates <- as.character(ensDates)
   L <- as.logical(match( Dates, dates, nomatch=0))
   if (all(!L) || !length(L))
     stop("model fit dates incompatible with ensemble data")
   Dates <- Dates[L]
   ensembleData <- ensembleData[L,]
   nObs <- length(Dates)
 }


 obs <- ensembleVerifObs(ensembleData)
 nForecasts <- ensembleSize(ensembleData)

 CRPS <- rep(NA, nObs)

 ensembleData <- ensembleForecasts(ensembleData)
 
 crpsEns1 <- apply(abs(sweep(ensembleData, MARGIN=1,FUN ="-",STATS=obs))
                   ,1,mean,na.rm=TRUE)
 
 if (nrow(ensembleData) > 1) {
   crpsEns2 <- apply(apply(ensembleData, 2, function(z,Z) 
     apply(abs(sweep(Z, MARGIN = 1, FUN = "-", STATS = z)),1,sum,na.rm=TRUE),
     Z = ensembleData),1,sum, na.rm = TRUE)
 }
 else {
   crpsEns2 <- sum(sapply(as.vector(ensembleData), 
                          function(z,Z) sum( Z-z, na.rm = TRUE),
                          Z = as.vector(ensembleData)), na.rm = TRUE)
 }
 
 crpsEns <- crpsEns1 - crpsEns2/(2*(nForecasts*nForecasts))
 

 l <- 0
 for (d in dates) {

    l <- l + 1
    k <- K[l]

    B <- fit$B[,k]
    if (all(Bmiss <- is.na(B))) next

    A <- fit$a[,k]
    C <- fit$c[,k]
    D <- fit$d[,k]
    Q <- fit$q[,k]

    I <- which(as.logical(match(Dates, d, nomatch = 0)))

    for (i in I) {

       f <- ensembleData[i,]
       S.sq <- mean(f)
       f <- c(1,f)
       Mu <- c(A,B)%*%f
       Sig <- C + D*S.sq              
       CRPS[i] <- crpsFunc(Mu, Sig, Q, obs[i])
    }
}
 if (any(is.na(c(crpsEns,CRPS)))) warning("NAs in crps values") 
 
 cbind(ensemble = crpsEns, EMOS = CRPS)
}

