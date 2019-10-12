controlMOSnormal <-
function(scoringRule = c("crps","log"),
         optimRule = c("BFGS","Nelder-Mead"),
         coefRule = c("square", "none", "positive"),
         varRule = c("square","none"),
         start = list(a = NULL, B = NULL, c = NULL, d = NULL),
         maxIter = Inf)
{
 if (is.infinite(maxIter) && maxIter > 0) maxIter <- .Machine$integer.max
 list(scoringRule = scoringRule[1], optimRule = optimRule[1], coefRule = coefRule[1],
      varRule = varRule[1], start = start, maxIter = maxIter)
}