crps.fitMOStruncnormal <-
function(fit, ensembleData, dates=NULL, ...)
{

 if(!is.null(dates)) warning("dates ignored")
 crpsFunc <- function(mu, sig, y)
  {
    z <- (y - mu)/sig
    w <- mu/sig
    crps <- sig*(z*pnorm(w)*(2*pnorm(z)+pnorm(w)-2) + 2*dnorm(z)*pnorm(w) - pnorm(w*sqrt(2))/sqrt(pi))/(pnorm(w))^2
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

    for (i in 1:nObs) {

       f <- ensembleData[i,]
       S.sq <- var(f)
       f <- c(1,f)
       Mu <- c(A,B)%*%f
       Sig <- sqrt(C + D*S.sq)

       CRPS[i] <- crpsFunc(Mu, Sig, obs[i])

   }

}
 cbind(ensemble = crpsEns, EMOS = CRPS)
}

