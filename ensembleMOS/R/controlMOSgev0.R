controlMOSgev0 <-
function(optimRule = c("Nelder-Mead","L-BFGS-B","BFGS"),
         coefRule = c("square", "none", "positive"),
         varRule = c("square","none"),
         start = list(a = NULL, B = NULL, s = NULL, c = NULL, d = NULL, q = NULL),
         maxIter = Inf)
{
 if (is.infinite(maxIter) && maxIter > 0) maxIter <- .Machine$integer.max
 list(optimRule = optimRule[1], coefRule = coefRule[1],
      varRule = varRule[1], start = start, maxIter = maxIter)
}

