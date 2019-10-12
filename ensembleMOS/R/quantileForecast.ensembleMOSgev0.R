quantileForecast.ensembleMOSgev0 <-
function(fit, ensembleData, quantiles = 0.5, dates = NULL, ...)
{
 
 gini.md <- function(x,na.rm=FALSE)  {     ## Michael Scheuerer's code
   if(na.rm & any(is.na(x)))  x <- x[!is.na(x)] 
   n <-length(x)
   return(4*sum((1:n)*sort(x,na.last=TRUE))/(n^2)-2*mean(x)*(n+1)/n)
    }

 matchITandFH(fit,ensembleData)

 M <- matchEnsembleMembers(fit,ensembleData)
 nForecasts <- ensembleSize(ensembleData)
 if (!all(M == 1:nForecasts)) ensembleData <- ensembleData[,M]

## remove instances missing all forecasts

 M <- apply(ensembleForecasts(ensembleData), 1, function(z) all(is.na(z)))
 ensembleData <- ensembleData[!M,]

## match specified dates with dateTable in fit

 dateTable <- dimnames(fit$B)[[2]]

 if (!is.null(dates)) {

   dates <- sort(unique(as.character(dates)))

   if (length(dates) > length(dateTable))
     stop("parameters not available for some dates")

   K <- match( dates, dateTable, nomatch=0)

   if (any(!K) || !length(K))
     stop("parameters not available for some dates")

 }
 else {

   dates <- dateTable
   K <- 1:length(dateTable)

  }

 ensDates <- ensembleValidDates(ensembleData)

## match dates in data with dateTable
 if (is.null(ensDates) || all(is.na(ensDates))) {
   if (length(dates) > 1) stop("date ambiguity")
   nObs <- nrow(ensembleData)
   Dates <- rep( dates, nObs)
 }
 else {
## remove instances missing dates
   if (any(M <- is.na(ensDates))) {
     ensembleData <- ensembleData[!M,]
     ensDates <- ensembleValidDates(ensembleData)
   }
   Dates <- as.character(ensDates)
   L <- as.logical(match( Dates, dates, nomatch=0))
   if (all(!L) || !length(L))
     stop("model fit dates incompatible with ensemble data")
   Dates <- Dates[L]
   ensembleData <- ensembleData[L,]
   nObs <- length(Dates)
 }

 nForecasts <- ensembleSize(ensembleData)

 Quants <- matrix(NA, nObs, length(quantiles))
 dimnames(Quants) <- list(ensembleObsLabels(ensembleData),
                          as.character(quantiles))

 ensembleData <- ensembleForecasts(ensembleData)

 l <- 0
 for (d in dates) {
   

    l <- l + 1
    k <- K[l]

    B <- fit$B[,k]
    if (all(Bmiss <- is.na(B))) next
    A <- fit$a[,k]
    C <- fit$c[,k]
    D <- fit$d[,k]
    SHAPE <- fit$q[,k]
    S <- fit$s[,k]

    I <- which(as.logical(match(Dates, d, nomatch = 0)))

    for (i in I) {
        f <- ensembleData[i,]

         
        MEAN <- as.numeric(c(A,B)%*%c(1,f)+S*mean(f==0, na.rm = TRUE))  #location of GEV
        SCALE <- C + D*gini.md(f, na.rm = TRUE)  #scale of GEV
        LOC <- as.numeric(MEAN - SCALE*(gamma(1-SHAPE)-1)/SHAPE)
        
        zero.prec <- as.vector(pgev(0, loc=LOC, scale=SCALE, shape=SHAPE))>=quantiles
        
        Quants[i,is.na(zero.prec)] <-NA
        
      	Quants[i,(!is.na(zero.prec))&zero.prec] <- 0

        Quants[i,(!is.na(zero.prec))&(!zero.prec)] <- qgev(quantiles, loc=LOC, scale=SCALE, shape=SHAPE)[(!is.na(zero.prec))&(!zero.prec)]
        
    }
 }

 if (any(is.na(Quants))) warning("NAs in quantiles values") 
 Quants
}

