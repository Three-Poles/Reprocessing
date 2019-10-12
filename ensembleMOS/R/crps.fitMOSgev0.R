crps.fitMOSgev0 <-
function(fit, ensembleData, dates=NULL, ...)
{

 if(!is.null(dates)) warning("dates ignored")
 
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

 M <- matchEnsembleMembers(fit,ensembleData)
 nForecasts <- ensembleSize(ensembleData)
 if (!all(M == 1:nForecasts)) ensembleData <- ensembleData[,M]

# remove instances missing all forecasts or obs

 M <- apply(ensembleForecasts(ensembleData), 1, function(z) all(is.na(z)))
 M <- M | is.na(ensembleVerifObs(ensembleData))
 ensembleData <- ensembleData[!M,]

 if (is.null(obs <- ensembleVerifObs(ensembleData)))
   stop("verification observations required")

 nObs <- ensembleNobs(ensembleData)

 nForecasts <- ensembleSize(ensembleData)

 CRPS <- rep(NA, nObs)

 obs <- ensembleVerifObs(ensembleData)
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
 
 B <- fit$B

 if (!all(Bmiss <- is.na(B))) {

    A <- fit$a
    C <- fit$c
    D <- fit$d
    Q <- fit$q
    S <- fit$s

    for (i in 1:nObs) {

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
 cbind(ensemble = crpsEns, EMOS = CRPS)
}

