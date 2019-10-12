pars.fitMOSlognormal <-
function(fit, ensembleData, dates=NULL, ...)
{

 if(!is.null(dates)) warning("dates ignored")


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

 obsLabels <- ensembleObsLabels(ensembleData)
 nForecasts <- ensembleSize(ensembleData)

 MEAN <- STD <- rep(NA, nObs)
 
 ensembleData <- ensembleForecasts(ensembleData)
 

 B <- fit$B

 if (!all(Bmiss <- is.na(B))) {

    A <- fit$a
    C <- fit$c
    D <- fit$d

    for (i in 1:nObs) {

       f <- ensembleData[i,]
       S.sq <- var(f)
       f <- c(1,f)
       Mval <- c(A,B)%*%f
       Vval <- C + D*S.sq       
       MEAN[i] <- log(Mval^2/sqrt(Vval + Mval^2))  #mean of normal
       STD[i] <- sqrt(log(1 +(Vval/(Mval^2))))  #std dev of normal
       
   }

}
 parValues <- cbind(meanlog = MEAN, sdlog = STD)
 row.names(parValues) <- obsLabels
 
 parValues
}

