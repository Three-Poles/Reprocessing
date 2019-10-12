pars.fitMOScsg0 <-
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

 SHAPE <- SCALE <- SHIFT <- rep(NA, nObs)
 
 ensembleData <- ensembleForecasts(ensembleData)
 

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
       
       SHAPE[i] <- Mu^2/Sig  #shape of gamma
       SCALE[i] <- Sig/Mu    #scale of gamma
       
       SHIFT[i] <- Q
       


   }

 }
 parValues <- cbind(shape = SHAPE, scale = SCALE, shift = SHIFT)
 row.names(parValues) <- obsLabels
 
 parValues
}

