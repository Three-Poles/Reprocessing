\name{controlMOSgev0}
\alias{controlMOSgev0}
\title{
  Control parameters for censored generalized extreme value distribution EMOS models
}
\description{
  Specifies a list of values controling the censored generalized extreme value distribution
  EMOS fit of ensemble forecasts.
}
\usage{
controlMOSgev0(optimRule = c("Nelder-Mead", "L-BFGS-B", "BFGS"),
               coefRule = c("square", "none", "positive"),
               varRule = c("square", "none"),
               start = list(a = NULL, B = NULL, 
                            s = NULL, c = NULL, 
                            d = NULL, q = NULL),
               maxIter = Inf)
}
\arguments{
  \item{optimRule}{
    Numerical optimization method to be supplied to \code{optim}. 
    Options are "BFGS" for the Broyden-Fletcher-Goldfarb-Shanno algorithm,
    "L-BFGS-B" for a constrained version thereof, and "Nelder-Mead" for the
    Nelder-Mead method, see \link{optim} for details.
  }
  \item{coefRule}{
    Method to control non-negativity of regression
    estimates. Options are:
    \itemize{
      \item \code{"square"} EMOS coefficients are
      parameterized as squares and thus gauranteed to be non-negative.
      \item \code{"positive"} finds non-negative coefficents
      iteratively by setting negative estimates at the current iteration
      to zero. 
      \item \code{"none"} no restriction on the coefficient
      estimates.
    }
  }
  \item{varRule}{
    Method to control non-negativity of the scale parameters.
    Options \code{"square"} and \code{"none"} are the same as in
    \code{coefRule}.
  }
  \item{start}{
    A list of starting parameters, \code{a}, \code{B}, \code{s}, 
    \code{c}, \code{d} and \code{q} specifying initial values for the 
    location, scale and shape coefficients supplied to \code{optim}. See details.
  }
  \item{maxIter}{
    An integer specifying the upper limit of the number of iterations
    used to fit the model.
  }
}
\details{
  Note that only minimum CRPS estimation is available and
  chosen by default.

  If no value is assigned to an argument, the first entry of 
  the list of possibly choices will be used by default.
  
  Given an ensemble of size \eqn{m}: \eqn{X_1, \ldots , X_m}, the
  following generalized extreme value distribution EMOS
  model left-censored at 0 is fit by \code{ensembleMOSgev0}: 
      
  \deqn{Y ~ GEV_0(\mu,\sigma,q)}
      
  where \eqn{GEV_0} denotes the generalized extreme value distribution
  left-censored at zero,
  with location \eqn{\mu}, scale \eqn{\sigma} and shape \eqn{q}. The model is 
  parametrized such that the mean \eqn{m} is a linear function
  \eqn{a + b_1 X_1 + \ldots + b_m X_m + s p_0}
  of the ensemble forecats, where \eqn{p_0} denotes the ratio of ensemble forecasts
  that are exactly 0, and the shape parameter \eqn{\sigma} is a linear
  function of the ensemble variance \eqn{c + d MD(X_1,\ldots,X_m)}, where
  \eqn{MD(X_1,\ldots,X_m)} is Gini's mean difference. 
  See \link{ensembleMOSgev0} for details.
}
\value{
  A list whose components are the input arguments and their assigned
  values.
}

\references{
  M. Scheuerer, Probabilistic quantitative precipitation forecasting using ensemble
  model output statistics. \emph{Quarterly Journal of the Royal Meteorological
  Society} 140:1086--1096, 2014.
}

\seealso{
  \code{\link{ensembleMOScsg0}},
  \code{\link{fitMOScsg0}}
}
\examples{
data("ensBMAtest", package = "ensembleBMA")

ensMemNames <- c("gfs","cmcg","eta","gasp","jma","ngps","tcwb","ukmo")

obs <- paste("PCP24","obs", sep = ".")
ens <- paste("PCP24", ensMemNames, sep = ".")
prcpTestData <- ensembleData(forecasts = ensBMAtest[,ens],
                             dates = ensBMAtest[,"vdate"],
                             observations = ensBMAtest[,obs],
                             station = ensBMAtest[,"station"],
                             forecastHour = 48,
                             initializationTime = "00")
                             
prcpTestFitGEV0 <- ensembleMOSgev0(prcpTestData, trainingDays = 25,
                        dates = "2008010100",
                        control = controlMOSgev0(maxIter = as.integer(100),
                           optimRule = "Nelder-Mead",
                           coefRule= "none", 
                           varRule = "square"))
}
\keyword{models}