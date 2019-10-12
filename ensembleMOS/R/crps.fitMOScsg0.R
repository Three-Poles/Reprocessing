crps.fitMOScsg0 <-
function(fit, ensembleData, dates=NULL, ...)
{

 if(!is.null(dates)) warning("dates ignored")
 crpsFunc <- function(mu, sig.Sq,qval,y)
  {
   
      Shp <- mu^2/sig.Sq  #shape of gamma
      Scl <- sig.Sq/mu    #scale of gamma
      
      
      if (!is.na(Scl) & (Scl <= 0))  stop("scale of gamma distribution is not positive")
      
      Z <- (y + qval)/Scl
      C <- qval/Scl 
      
      crps <- Scl*Z*(2*pgamma(Z,Shp,1)-1)-Scl*C*(pgamma(C,Shp,1))^2 + mu*(1+2*pgamma(C,Shp,1)*pgamma(C,Shp+1,1)
	      -(pgamma(C,Shp,1))^2-2*pgamma(Z,Shp+1,1)) - mu*(1-pgamma(2*C,2*Shp,1))*beta(.5,Shp+.5)/pi
   
    crps
  }
 

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

    for (i in 1:nObs) {

       f <- ensembleData[i,]
       S.sq <- mean(f)
       f <- c(1,f)
      
       Mu <- c(A,B)%*%f
       Sig <- C + D*S.sq       
       
       CRPS[i] <- crpsFunc(Mu, Sig, Q, obs[i])

   }

}
 cbind(ensemble = crpsEns, EMOS = CRPS)
}

