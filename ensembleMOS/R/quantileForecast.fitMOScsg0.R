quantileForecast.fitMOScsg0 <-
function(fit, ensembleData, quantiles = 0.5, dates = NULL, ...)
{
 
 M <- matchEnsembleMembers(fit,ensembleData)
 nForecasts <- ensembleSize(ensembleData)
 if (!all(M == 1:nForecasts)) ensembleData <- ensembleData[,M]

## remove instances missing all forecasts or dates

 M <- apply(ensembleForecasts(ensembleData), 1, function(z) all(is.na(z)))

 ensembleData <- ensembleData[!M,]

 nObs <- nrow(ensembleData)

 if (!is.null(dates)) warning("dates ignored")

 Quants <- matrix(NA, nObs, length(quantiles))
 dimnames(Quants) <- list(ensembleObsLabels(ensembleData),as.character(quantiles))
 Mu <- rep(NA,nObs)
 Sig.sq <- rep(NA,nObs)

 ensembleData <- ensembleForecasts(ensembleData)
 x <- c(fit$a,fit$B)
 A <- cbind(rep(1,nObs),ensembleData)
 Q <- fit$q

 S.sq <- apply(ensembleData,1,mean)

 Mu <- A%*%x
 Sig.sq <- rep(fit$c,nObs) + rep(fit$d,nObs)*S.sq
 
 Shp <- Mu^2/Sig.sq  #shape of gamma
 Scl <- Sig.sq/Mu    #scale of gamma


 for (i in 1:length(quantiles)) 
        {zero.prec <- pgamma(Q,shape=Shp,scale=Scl)>=quantiles[i]
         Quants[is.na(zero.prec),i]<-NA
         Quants[(!is.na(zero.prec))&zero.prec,i]<-0
         Quants[(!is.na(zero.prec))&(!zero.prec),i] <- qgamma(quantiles[i], shape=Shp[(!is.na(zero.prec))&(!zero.prec)], scale=Scl[(!is.na(zero.prec))&(!zero.prec)])-fit$q}
 Quants
}

