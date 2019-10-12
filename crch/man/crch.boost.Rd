\name{crch.boost}
\alias{crch.boost}
\alias{crch.boost.fit}

\title{Auxiliary functions to fit \code{crch} models via boosting.}

\description{
  Auxiliary functions to fit \code{crch} models via boosting
}
\usage{
crch.boost(maxit = 100, nu = 0.1, start = NULL, dot = "separate", 
  mstop = c("max", "aic", "bic", "cv"),  nfolds = 10, foldid = NULL, 
  maxvar = NULL)

crch.boost.fit(x, z, y, left, right, truncated = FALSE, dist = "gaussian",
  df = NULL, link.scale = "log", type = "ml", weights = NULL, offset = NULL, 
  control = crch.boost())
}

\arguments{
  \item{maxit}{the maximum number of boosting iterations.}
  \item{nu}{boosting step size. Default is 0.1.}
  \item{start}{a previously boosted but not converged \code{"crch.boost"} 
    object to continue.}
  \item{dot}{character specifying how to process formula parts with a dot
    (\code{.}) on the right-hand side. This can either be \code{"separate"}
    so that each formula part is expanded separately or \code{"sequential"}
    so that the parts are expanded sequentially conditional on all prior parts.
    Default is \code{"separate"}}
  \item{mstop}{method to find optimum stopping iteration. Default is \code{"max"}
    which is \code{maxit}. Alternatives are \code{"aic"} and \code{"bic"} for 
    AIC and BIC
    optimization and \code{"cv"} for cross validation.
    \code{mstop} can also be a positive integer to set the number of boosting
    iterations. Then \code{maxit} is set to \code{mstop} and \code{mstop="max"}.}
  \item{nfolds}{if \code{mstopopt = "cv"}, number of folds in cross validation.}
  \item{foldid}{if \code{mstopopt = "cv"}, an optional vector of values between 1 
    and \code{nfold} identifying the fold each observation is in. If supplied,
    \code{nfolds} can be missing.}
  \item{maxvar}{Positive \code{numeric}. Maximum number of parameters to be selected
    during each iteration (not including intercepts). Used for stability selection.}
  \item{x, z, y, left, right, truncated, dist, df, link.scale, type, weights, offset, control}{see \code{\link{crch.fit}} for details.}
}

\details{
\code{crch.boost} extends \code{crch} to fit censored (tobit) or 
truncated regression models with conditional heteroscedasticy by
boosting. If \code{crch.boost()} is supplied as \code{control} in 
\code{crch} then \code{crch.boost.fit} is used as lower level fitting
function. Note that \code{\link{crch.control}()} with \code{method=boosting}
is equivalent to \code{crch.boost()}. Thus, boosting can more 
conveniently be called with \code{crch(\ldots, method = "boosting")}.
}

\value{
  For \code{crch.boost}: A list with components named as the arguments.
  For \code{crch.boost.fit}: An object of class \code{"crch.boost"}, 
    i.e., a list with the following elements.
  \item{coefficients}{list of coefficients for location and scale. Scale
  coefficients are in log-scale. Coefficients are of optimum stopping
  stopping iteration specified by \code{mstop}.}
  \item{df}{if \code{dist = "student"}: degrees of freedom of student-t
  distribution. else \code{NULL}.}
  \item{residuals}{the residuals, that is response minus fitted values.}
  \item{fitted.values}{list of fitted location and scale parameters at
  optimum stopping iteration specified by \code{mstop}.}
  \item{dist}{assumed distribution for the dependent variable \code{y}.}
  \item{cens}{list of censoring points.}
  \item{control}{list of control parameters.}
  \item{weights}{case weights used for fitting.}
  \item{offset}{list of offsets for location and scale.}
  \item{n}{number of observations.}
  \item{nobs}{number of observations with non-zero weights.}
  \item{loglik}{log-likelihood.}
  \item{link}{a list with element \code{"scale"} containing
          the link objects for the scale model.}
  \item{truncated}{logical indicating wheter a truncated model has been 
  fitted.}
  \item{iterations}{number of boosting iterations.}
  \item{stepsize}{boosting stepsize \code{nu}.}
  \item{mstop}{criterion used to find optimum stopping iteration.}
  \item{mstopopt}{optimum stopping iterations for different criteria.}
  \item{standardize}{list of center and scale values used to standardize 
  response and regressors.}
}

\references{
  Messner JW, Mayr GJ, Zeileis A (2016). Non-Homogeneous Boosting for Predictor 
    Selection in Ensemble Post-Processing. \emph{Working Papers}, Faculty of 
    Economics and Statistics, University of Innsbruck, url:http://econpapers.repec.org/paper/innwpaper/2016-04.htm.
}

\seealso{\code{\link{crch}}, \code{\link{crch.control}}}

\examples{
# generate data
suppressWarnings(RNGversion("3.5.0"))
set.seed(5)
x <- matrix(rnorm(1000*20),1000,20)
y <- rnorm(1000, 1 + x[,1] - 1.5 * x[,2], exp(-1 + 0.3*x[,3]))
y <- pmax(0, y)
data <- data.frame(cbind(y, x))

# fit model with maximum likelihood
CRCH <- crch(y ~ .|., data = data, dist = "gaussian", left = 0)

# fit model with boosting
boost <- crch(y ~ .|.,  data = data, dist = "gaussian", left = 0,
  control = crch.boost(mstop = "aic"))

# more conveniently, the same model can also be fit through
# boost <- crch(y ~ .|.,  data = data, dist = "gaussian", left = 0,
#   method = "boosting", mstop = "aic")

# AIC comparison
AIC(CRCH, boost)

# summary
summary(boost)

# plot
plot(boost)
}

\keyword{regression, boosting}
