pars.ensembleMOScsg0 <-
function(fit, ensembleData, dates=NULL, ...)
{
 
 matchITandFH(fit,ensembleData)

 M <- matchEnsembleMembers(fit,ensembleData)
 nForecasts <- ensembleSize(ensembleData)
 if (!all(M == 1:nForecasts)) ensembleData <- ensembleData[,M]

## remove instances missing all forecasts

 M <- apply(ensembleForecasts(ensembleData), 1, function(z) all(is.na(z)))
 M <- M | is.na(ensembleVerifObs(ensembleData))
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

 
 obsLabels <- ensembleObsLabels(ensembleData)
 nForecasts <- ensembleSize(ensembleData)
 

 SHAPE <- SCALE <- SHIFT <- rep(NA, nObs)
 
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
    Q <- fit$q[,k]

    I <- which(as.logical(match(Dates, d, nomatch = 0)))

    for (i in I) {

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
 
 if (any(is.na(c(SHAPE,SCALE,SHIFT)))) warning("NAs in parameters") 
 
 
 parValues <- cbind(shape = SHAPE, scale = SCALE, shift = SHIFT)
 row.names(parValues) <- obsLabels
 
 parValues
 
}

