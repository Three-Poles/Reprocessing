pars.fitMOSgev0 <-
function(fit, ensembleData, dates=NULL, ...)
{

  gini.md <- function(x,na.rm=FALSE)  {     ## Michael Scheuerer's code
    if(na.rm & any(is.na(x)))  x <- x[!is.na(x)] 
    n <-length(x)
    return(4*sum((1:n)*sort(x,na.last=TRUE))/(n^2)-2*mean(x)*(n+1)/n)
  }
  
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

 LOC <- SCALE <- SHAPE <- rep(NA, nObs)
 
 ensembleData <- ensembleForecasts(ensembleData)
 

 B <- fit$B

 if (!all(Bmiss <- is.na(B))) {

    A <- fit$a
    C <- fit$c
    D <- fit$d
    Q <- fit$q
    S <- fit$s

    for (i in 1:nObs) {

       f <- ensembleData[i,]
       MEAN <- as.numeric(c(A,B)%*%c(1,f)+S*mean(f==0))  #location of GEV
       SCALE[i] <- C + D*gini.md(f)  #scale of GEV
       LOC[i] <- as.numeric(MEAN - SCALE[i]*(gamma(1-Q)-1)/Q)
       SHAPE[i] <- Q
   }

 }
 parValues <- cbind(loc = LOC, scale = SCALE, shape = SHAPE)
 row.names(parValues) <- obsLabels
 
 parValues
}

