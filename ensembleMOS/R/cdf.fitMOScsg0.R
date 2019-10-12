cdf.fitMOScsg0 <-
function(fit, ensembleData, values, dates = NULL, randomizeATzero = FALSE, ...)
{

 M <- matchEnsembleMembers(fit,ensembleData)
 nForecasts <- ensembleSize(ensembleData)
 if (!all(M == 1:nForecasts)) ensembleData <- ensembleData[,M]

## remove instances missing all forecasts or dates

 M <- apply(ensembleForecasts(ensembleData), 1, function(z) all(is.na(z)))

 ensembleData <- ensembleData[!M,]

 nObs <- nrow(ensembleData)

 if (!is.null(dates)) warning("dates ignored")

 CDF <- matrix(NA, nObs, length(values))
 dimnames(CDF) <- list(ensembleObsLabels(ensembleData),as.character(values))
 Mu <- rep(NA,nObs)
 Sig.sq <- rep(NA,nObs)

 ensembleData <- ensembleForecasts(ensembleData)
 x <- c(fit$a,fit$B)
 A <- cbind(rep(1,nObs),ensembleData)

 S.sq <- apply(ensembleData,1,mean)
 Mu <- A%*%x
 Sig.sq <- rep(fit$c,nObs) + rep(fit$d,nObs)*S.sq
 
 Shp <- Mu^2/Sig.sq  #shape of gamma
 Scl <- Sig.sq/Mu    #scale of gamma
 
 for (i in 1:length(values)){
   if (randomizeATzero & (values[i]==0)){
     cdfval <- pgamma(fit$q, shape = Shp, scale = Scl)
     CDF[,i] <- runif(nObs,0,cdfval)
   }
   else {
     CDF[,i] <- pgamma(values[i]+fit$q, shape = Shp, scale = Scl)
   }
 }
 CDF
}

