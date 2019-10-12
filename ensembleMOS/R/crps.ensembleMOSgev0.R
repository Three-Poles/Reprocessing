crps.ensembleMOSgev0 <-
function(fit, ensembleData, dates=NULL, ...)
{

  gini.md <- function(x,na.rm=FALSE)  {     ## Michael Scheuerer's code
   if(na.rm & any(is.na(x)))  x <- x[!is.na(x)] 
   n <-length(x)
   return(4*sum((1:n)*sort(x,na.last=TRUE))/(n^2)-2*mean(x)*(n+1)/n)
    }
    
  crps.GEVneq0 <- function(MEAN,SCALE,SHAPE, obs)  { ## Michael Scheuerer's code
  
      LOC <- MEAN - SCALE*(gamma(1-SHAPE)-1)/SHAPE

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
    S <- fit$s[,k]
    
    

    I <- which(as.logical(match(Dates, d, nomatch = 0)))

    for (i in I) {

       f <- ensembleData[i,]
       
       
      MEAN <- c(A,B)%*%c(1,f)+S*mean(f==0, na.rm = TRUE)  #mean of GEV
      SCALE <- C + D*gini.md(f,na.rm = TRUE)  #scale of GEV
      
      if( abs(Q) < eps )  {
      crps.eps.m <- crps.GEVneq0 (MEAN,SCALE,-eps, obs[i])
      crps.eps.p <- crps.GEVneq0 (MEAN,SCALE,eps, obs[i])
      w.m <- (eps-Q)/(2*eps)
      w.p <- (eps+Q)/(2*eps)
      CRPS[i] <- w.m*crps.eps.m + w.p*crps.eps.p
   } else  {
      CRPS[i] <- crps.GEVneq0 (MEAN,SCALE,Q, obs[i])
   }

}

}
 if (any(is.na(c(crpsEns,CRPS)))) warning("NAs in crps values") 
 
 cbind(ensemble = crpsEns, EMOS = CRPS)
 
}

