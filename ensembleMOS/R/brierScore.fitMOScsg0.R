brierScore.fitMOScsg0 <-
function(fit, ensembleData, thresholds, dates = NULL, ...)
{
 
 M <- matchEnsembleMembers(fit,ensembleData)
 nForecasts <- ensembleSize(ensembleData)
 if (!all(M == 1:nForecasts)) ensembleData <- ensembleData[,M]

## remove instances missing all forecasts or dates

 M <- apply(ensembleForecasts(ensembleData), 1, function(z) all(is.na(z)))

 ensembleData <- ensembleData[!M,]

 nObs <- nrow(ensembleData)

 if (!is.null(dates)) warning("dates ignored")

 BScores <- matrix(NA, nObs, length(thresholds))
 dimnames(BScores) <- list(ensembleObsLabels(ensembleData),as.character(thresholds))
 Mu <- rep(NA,nObs)
 Sig.sq <- rep(NA,nObs)

 obs <- ensembleVerifObs(ensembleData)
 ensembleData <- ensembleForecasts(ensembleData)
 x <- c(fit$a,fit$B)
 A <- cbind(rep(1,nObs),ensembleData)
 Q <- fit$q

 S.sq <- apply(ensembleData,1,mean)

 Mu <- A%*%x
 Sig.sq <- rep(fit$c,nObs) + rep(fit$d,nObs)*S.sq
 
 Shp <- Mu^2/Sig.sq  #shape of gamma
 Scl <- Sig.sq/Mu    #scale of gamma


 for (i in 1:length(thresholds)) 
    {BScores[,i] <- (pgamma(thresholds[i]+Q, shape = Shp, scale = Scl)* (thresholds[i]>=0) - (obs <= thresholds[i]))^2}
        
 BScores
}

