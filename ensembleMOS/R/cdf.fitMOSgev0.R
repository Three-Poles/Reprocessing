cdf.fitMOSgev0 <-
function(fit, ensembleData, values, dates = NULL, randomizeATzero = FALSE, ...)
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

 CDF <- matrix(NA, nObs, length(values))
 dimnames(CDF) <- list(ensembleObsLabels(ensembleData),as.character(values))
 
 ensembleData <- ensembleForecasts(ensembleData)
 x <- c(fit$a,fit$B)
 A <- cbind(rep(1,nObs),ensembleData)
 SHAPE <- fit$q
 S <- fit$s
 
 S.sq <- apply(ensembleData,1,gini.md, na.rm = TRUE)
 
 
 MEAN <- A%*%x + S*rowMeans(ensembleData==0, na.rm = TRUE)
 SCALE <- rep(fit$c,nObs) + rep(fit$d,nObs)*S.sq
 LOC <- MEAN - SCALE*(gamma(1-SHAPE)-1)/SHAPE
 
 for (i in 1:length(values)){
   if (randomizeATzero & (values[i]==0)){
     cdfval <- pgev(0, loc=LOC, scale=SCALE, shape=SHAPE)
     CDF[,i] <- runif(nObs,0,cdfval)
   }
   else {
     CDF[,i] <- pgev(values[i], loc=LOC, scale=SCALE, shape=SHAPE)
   }
 }
 CDF
}

