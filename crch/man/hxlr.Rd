\name{hxlr}
\alias{hxlr}

\title{Heteroscedastic Extended Logistic Regression}

\description{
  This is a wrapper function for \code{\link[ordinal]{clm}} (from package
  \pkg{ordinal}) to fit (heteroscedastic) extended logistic regression (HXLR)
  models (Messner et al. 2013).
}
\usage{
hxlr(formula, data, subset, na.action, weights, thresholds, link, control, \dots)
}
\arguments{
  \item{formula}{a formula expression of the form \code{y ~ x | z} where
  \code{y} is the response and \code{x} and \code{z} are regressor variables
  for the location and the scale of the latend distribution respectively.
  Response can either be a continuous variable or a factor.}
  \item{data}{an optional data frame containing the variables occurring in the
  formulas.}
  \item{subset}{an optional vector specifying a subset of observations to be
  used for fitting.}
  \item{na.action}{a function which indicates what should happen when the data
  contain \code{NA}s. Default is na.omit}
  \item{weights}{optional case weights in fitting.}
  \item{thresholds}{vector of (transformed) thresholds that are used to cut the
  continuous response into categories. Data frames or matrices with multiple
  columns are allowed as well. Then each column is used as separate predictor
  variable for the intercept model.}
  \item{link}{link function, i.e., the type of location-scale distribution
  assumed for the latent distribution. Default is \code{logit}.}
  \item{control}{a list of control parameters passed to \code{\link{optim}}.
  Default is \code{\link{hxlr.control}}}
  \item{\dots}{arguments to be used to form the default \code{control} argument
  if it is not supplied directly.}
}

\details{
Extended logistic regression (Wilks 2009) extends binary logistic regression to
multi-category responses by including the thresholds, that are used to cut a
continuous variable into categories, in the regression equation.
Heteroscedastic extended logistic regression (Messner et al. 2013) extends this
model further and allows to add additional predictor variables that are used to
predict the scale of the latent logistic distribution.
}

\value{
  An object of class \code{"hxlr"}, i.e., a list with the following elements.
  \item{coefficients}{list of CLM coefficients for intercept, location,
  and scale model.}
  \item{fitted.values}{list of fitted location and scale parameters.}
  \item{optim}{output from optimization from \code{\link{optim}}.}
  \item{method}{Optimization method used for \code{\link{optim}}.}
  \item{control}{list of control parameters passed to \code{\link{optim}}}
  \item{start}{starting values of coefficients used in the optimization.}
  \item{weights}{case weights used for fitting.}
  \item{n}{number of observations.}
  \item{nobs}{number of observations with non-zero weights.}
  \item{loglik}{log-likelihood.}
  \item{vcov}{covariance matrix.}
  \item{converged}{logical variable whether optimization has converged or not.}
  \item{iterations}{number of iterations in optimization.}
  \item{call}{function call.}
  \item{scale}{the formula supplied.}
  \item{terms}{the \code{terms} objects used.}
  \item{levels}{list of levels of the factors used in fitting for location and
  scale respectively.}
  \item{thresholds}{the thresholds supplied.}
}

\references{
  Messner JW, Mayr GJ, Zeileis A, Wilks DS (2014). Extending Extended Logistic
    Regression to Effectively Utilize the Ensemble Spread. \emph{Monthly Weather
    Review}, \bold{142}, 448--456. \doi{10.1175/MWR-D-13-00271.1}.

  Wilks DS (2009).  Extending Logistic Regression to Provide
    Full-Probability-Distribution MOS Forecasts.
    \emph{Meteorological Applications}, \bold{368}, 361--368.
}

\seealso{\code{\link{predict.hxlr}}, \code{\link[ordinal]{clm}}}

\examples{
data("RainIbk")
## mean and standard deviation of square root transformed ensemble forecasts
RainIbk$sqrtensmean <- 
  apply(sqrt(RainIbk[,grep('^rainfc',names(RainIbk))]), 1, mean)
RainIbk$sqrtenssd <- 
  apply(sqrt(RainIbk[,grep('^rainfc',names(RainIbk))]), 1, sd)

## climatological deciles
q <- unique(quantile(RainIbk$rain, seq(0.1, 0.9, 0.1)))

## fit ordinary extended logistic regression with ensemble mean as 
## predictor variable
XLR <- hxlr(sqrt(rain) ~ sqrtensmean, data = RainIbk, thresholds = sqrt(q))
## print
XLR
## summary
summary(XLR)


## fit ordinary extended logistic regression with ensemble mean 
## and standard deviation as predictor variables
XLRS <- hxlr(sqrt(rain) ~ sqrtensmean + sqrtenssd, data = RainIbk, 
  thresholds = sqrt(q))
## fit heteroscedastic extended logistic regression with ensemble 
## standard deviation as predictor for the scale
HXLR <- hxlr(sqrt(rain) ~ sqrtensmean | sqrtenssd, data = RainIbk, 
  thresholds = sqrt(q))

## compare AIC of different models
AIC(XLR, XLRS, HXLR)

## XLRS and HXLR are nested in XLR -> likelihood-ratio-tests
if(require("lmtest")) {
  lrtest(XLR, XLRS)
  lrtest(XLR, HXLR)
}

\dontrun{
###################################################################
## Cross-validation and bootstrapping RPS for different models 
## (like in Messner 2013). 
N <- NROW(RainIbk)
## function that returns model fits
fits <- function(data, weights = rep(1, N)) {
  list(
    "XLR"    = hxlr(sqrt(rain) ~ sqrtensmean, data = data, 
      weights = weights, thresholds = sqrt(q)),
    "XLR:S"  = hxlr(sqrt(rain) ~ sqrtensmean + sqrtenssd, data = data, 
      weights = weights, thresholds = sqrt(q)),
    "XLR:SM" = hxlr(sqrt(rain) ~ sqrtensmean + I(sqrtensmean*sqrtenssd), 
      data = data, weights = weights, thresholds = sqrt(q)),
    "HXLR"   = hxlr(sqrt(rain) ~ sqrtensmean | sqrtenssd, data = data, 
      weights = weights, thresholds = sqrt(q)),
    "HXLR:S" = hxlr(sqrt(rain) ~ sqrtensmean + sqrtenssd | sqrtenssd, 
      data = data, weights = weights, thresholds = sqrt(q))
  )
}


## cross validation
id <- sample(1:10, N, replace = TRUE)
obs <- NULL
pred <- list(NULL)
for(i in 1:10) {
  ## splitting into test and training data set
  trainIndex <- which(id != i)     
  testIndex <- which(id == i)			     
  ## weights that are used for fitting the models
  weights <- as.numeric(table(factor(trainIndex, levels = c(1:N))))
  ## testdata
  testdata <- RainIbk[testIndex,]
  ## observations    
  obs <- c(obs, RainIbk$rain[testIndex])
  ## estimation
  modelfits <- fits(RainIbk, weights)
  ## Prediction
  pred2 <- lapply(modelfits, predict, newdata = testdata, type = "cumprob")
  pred <- mapply(rbind, pred, pred2, SIMPLIFY = FALSE)
}
names(pred) <- c(names(modelfits))

## function to compute RPS
rps <- function(pred, obs) {
  OBS <- NULL
  for(i in 1:N) 
    OBS <- rbind(OBS, rep(0:1, c(obs[i] - 1, length(q) - obs[i] + 1)))
  apply((OBS-pred)^2, 1, sum)
}

## compute rps
RPS <- lapply(pred, rps, obs = as.numeric(cut(obs, c(-Inf, q, Inf))))

## bootstrapping mean rps 
rpsall <- NULL
for(i in 1:250) {
  index <- sample(length(obs), replace = TRUE)
  rpsall <- rbind(rpsall, sapply(RPS, function(x) mean(x[index])))
}
  
rpssall <- 1 - rpsall/rpsall[,1]
boxplot(rpssall[,-1], ylab = "RPSS", main = "RPSS relative to XLR")
abline(h = 0, lty = 2)
}
}

\keyword{regression}
