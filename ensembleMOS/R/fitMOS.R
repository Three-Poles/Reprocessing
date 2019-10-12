fitMOS <-
function(ensembleData, control = NULL, model = NULL,
         exchangeable = NULL)
{
 if (!inherits(ensembleData,"ensembleData")) stop("not an ensembleData object")
 mc <- match.call()
 mc$model <- NULL

 if (!is.null(model)) {
   switch( model,
           "normal" = {
             mc[[1]] <- as.name("fitMOSnormal")
           },
           "truncnormal" = {
             mc[[1]] <- as.name("fitMOStruncnormal")
           },
           "lognormal" = {
             mc[[1]] <- as.name("fitMOSlognormal")
           },
           "csg0" = {
             mc[[1]] <- as.name("fitMOScsg0")
           },
           "gev0" = {
             mc[[1]] <- as.name("fitMOSgev0")
           },
           stop("unrecognized model")
   )
 } else stop("unspecified model")
 
 if (length(attr(ensembleData, "class")) > 2) {
   attr(ensembleData, "class") <- attr(ensembleData, "class")[-1]
   mc$ensembleData <- ensembleData
 }

 eval(mc, parent.frame())
}

