quantileForecast.fitMOSgev0 <-
function(fit, ensembleData, quantiles = 0.5, dates = NULL, ...)
{
 
  gini.md <- function(x,na.rm=FALSE)  {     ## Michael Scheuerer's code
    if(na.rm & any(is.na(x)))  x <- x[!is.na(x)] 
    n <-length(x)
    return(4*sum((1:n)*sort(x,na.last=TRUE))/(n^2)-2*mean(x)*(n+1)/n)
  }
  
  
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

 ensembleData <- ensembleForecasts(ensembleData)
 x <- c(fit$a,fit$B)
 A <- cbind(rep(1,nObs),ensembleData)
 SHAPE <- fit$q
 S <- fit$s

 S.sq <- apply(ensembleData,1,gini.md, na.rm = TRUE)


 MEAN <- A%*%x + S*rowMeans(ensembleData==0, na.rm = TRUE)
 SCALE <- rep(fit$c,nObs) + rep(fit$d,nObs)*S.sq
 LOC <- MEAN - SCALE*(gamma(1-SHAPE)-1)/SHAPE
 

 for (i in 1:length(quantiles)) 
        {zero.prec <- as.vector(pgev(0, loc=LOC, scale=SCALE, shape=SHAPE))>=quantiles[i]
         Quants[is.na(zero.prec),i] <-NA
         Quants[(!is.na(zero.prec))&zero.prec,i]<-0
         Quants[(!is.na(zero.prec))&(!zero.prec),i] <- qgev(quantiles[i], loc=LOC[(!is.na(zero.prec))&(!zero.prec)], scale=SCALE[(!is.na(zero.prec))&(!zero.prec)], shape=SHAPE)
 }
 Quants
}

