getSummary.crch <- function(obj, alpha = 0.05, ...)
{
  ## internally force setting of the summary template (so that
  ## memisc does not have to be registered in the NAMESPACE already)
  memisc::setSummaryTemplate("crch" = c(
    "Log-likelihood" = "($logLik:f#)",
    "AIC" = "($AIC:f#)",
    "BIC" = "($BIC:f#)",
    "N" = "($N:d)"
  ))

  ## extract coefficient summary
  cf <- summary(obj)$coefficients
  ## augment with confidence intervals
  cval <- qnorm(1 - alpha/2)
  for(i in seq_along(cf)) cf[[i]] <- cbind(cf[[i]],
    cf[[i]][, 1] - cval * cf[[i]][, 2],
    cf[[i]][, 1] + cval * cf[[i]][, 2])
  ## collect in array
  nam <- unique(unlist(lapply(cf, rownames)))
  acf <- array(dim = c(length(nam), 6, length(cf)),
    dimnames = list(nam, c("est", "se", "stat", "p", "lwr", "upr"), names(cf)))
  for(i in seq_along(cf)) acf[rownames(cf[[i]]), , i] <- cf[[i]]

  ## contrasts (omitting duplicates between location and scale part) and factor levels
  ctr <- c(obj$contrasts$location, obj$contrasts$scale)
  ctr <- ctr[!duplicated(names(ctr))]
  xlev <- obj$levels$full
  
  ## return everything
  return(list(
    coef = acf,
    sumstat = c(
      "N" = nobs(obj),
      "logLik" = as.vector(logLik(obj)),
      "AIC" = AIC(obj),
      "BIC" = BIC(obj)
    ),
    contrasts = ctr,
    xlevels = xlev,
    call = obj$call
  ))
}

