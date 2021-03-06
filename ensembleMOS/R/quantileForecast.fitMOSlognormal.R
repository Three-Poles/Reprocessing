quantileForecast.fitMOSlognormal <-
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
 Sig <- rep(NA,nObs)

 ensembleData <- ensembleForecasts(ensembleData)
 x <- c(fit$a,fit$B)
 A <- cbind(rep(1,nObs),ensembleData)

 S.sq <- apply(ensembleData,1,var)

 Mval <- A%*%x
 Vval <- rep(fit$c,nObs) + rep(fit$d,nObs)*S.sq
 Mu <- log(Mval^2/sqrt(Vval + Mval^2))  #mean of normal
 Sig <- sqrt(log(1 +Vval/Mval^2))  #std dev of normal

 for(i in 1:length(quantiles)) Quants[,i] <- qlnorm(quantiles[i], Mu, Sig)
 Quants
}

