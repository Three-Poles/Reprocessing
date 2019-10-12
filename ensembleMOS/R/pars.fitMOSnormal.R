pars.fitMOSnormal <-
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
       MEAN[i] <- c(A,B)%*%f
       STD[i] <- sqrt(C + D*S.sq)


   }

 }
 parValues <- cbind(mean = MEAN, sd = STD)
 row.names(parValues) <- obsLabels
 
 parValues
}

