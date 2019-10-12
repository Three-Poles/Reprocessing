ensembleMOS <-
function(ensembleData, trainingDays, consecutive = FALSE, dates = NULL,
         control = NULL,  warmStart = FALSE,
         model = NULL, exchangeable = NULL)
{
 if (!inherits(ensembleData,"ensembleData")) stop("not an ensembleData object")
 mc <- match.call()
 mc$model <- NULL

 if (!is.null(model)) {
   switch( model,
        "normal" = {
             mc[[1]] <- as.name("ensembleMOSnormal")
          },
        "truncnormal" = {
             mc[[1]] <- as.name("ensembleMOStruncnormal")
          },
        "lognormal" = {
          mc[[1]] <- as.name("ensembleMOSlognormal")
        },
        "csg0" = {
          mc[[1]] <- as.name("ensembleMOScsg0")
        },
        "gev0" = {
          mc[[1]] <- as.name("ensembleMOSgev0")
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

