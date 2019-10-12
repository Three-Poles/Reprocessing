`getExchangeable` <-
function (argument, attribute, nForecasts) 
{
 if (is.null(argument)) {
   return(attribute)
 } else {
   if (length(argument) != nForecasts) {
     stop("exchangeable specification not consistent with forecasts")
     }
   if (is.null(attribute)){
     namArg <- names(argument)
     if (is.null(namArg)) {
       stop("no member names associated with exchangeable specification")
     }
     return(argument)
   } else {
     namAtr <- names(attribute)
     if (is.null(namAtr))
       stop("exchangeable data attribute has no member names")
     namArg <- names(argument)
     if (is.null(namArg)) {
       warning("no member names associated with exchangeable specification")
       names(argument) <- names(attribute)
       return(argument)
     }
     m <- match( namArg, namAtr, nomatch = 0)
     if (any(!m) || (length(unique(m)) != length(m))) {
       stop("exchangeable argument names inconsistent with data")
     }
     return(argument[order(m)])
   }
 }

} 

