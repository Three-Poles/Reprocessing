\name{controlMOSlognormal}
\alias{controlMOSlognormal}
\title{
  Control parameters for log-normal EMOS models
}
\description{
  Specifies a list of values controling the log-normal EMOS fit
  of ensemble forecasts.
}
\usage{
controlMOSlognormal(scoringRule = c("crps", "log"),
                    optimRule = c("BFGS","Nelder-Mead"),
                    coefRule = c("square", "none", "positive"),
                    varRule = c("square", "none"),
                    start = list(a = NULL, B = NULL,
                                 c = NULL, d = NULL),
                    maxIter = Inf)
}
\arguments{
  \item{scoringRule}{
    The scoring rule to be used in optimum score estimation.  Options
    are "crps" for the continuous ranked probability score and "log" for
    the logarithmic score.
  }
  \item{optimRule}{
    Numerical optimization method to be supplied to \code{optim}. 
    Options are "BFGS" for the
    Broyden-Fletcher-Goldfarb-Shanno algorithm and "Nelder-Mead" for the
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
    Method to control non-negativity of the variance parameters.
    Options \code{"square"} and \code{"none"} are the same as in
    \code{coefRule}.
  }
  \item{start}{
    A list of starting parameters, \code{a}, \code{B}, \code{c} and
    \code{d} specifying initial values for the intercept coefficient
    and scale parameters supplied to \code{optim}. See details.
  }
  \item{maxIter}{
    An integer specifying the upper limit of the number of iterations
    used to fit the model.
  }
}
\details{
  If no value is assigned to an argument, the first entry of 
  the list of possibly choices will be used by default.
  Given an ensemble of size \eqn{m}: \eqn{X_1, \ldots , X_m},  the
  following log-normal model is fit by \code{ensembleMOSlognormal}: 

  \deqn{Y ~ LN(\mu, \sigma)}

  where \eqn{LN} denotes the log-normal distrbution with \code{meanlog} 
  parameter \eqn{\mu} and \code{scalelog} parameter \eqn{\sigma}, see
  \link{Lognormal}. The model is parametrized such that the mean value of 
  the log-normal distribution is a linear function \eqn{a + b_1 X_1 + \ldots + b_m X_m}
  of the ensemble forecats, and the variance is a linear function
  \eqn{c + d S^2}. For transformations between \eqn{\mu, \sigma} and mean
  and variance of the log-normal distribution, see Baran and Lerch (2015).
  See \link{ensembleMOSlognormal} for details.
  
  Note that in case of \code{scoringRule = "log"}, forecast cases in the
  training period with observation values of 0 are ignored in the model 
  estimation as 0 is not included in the support of the log-normal 
  distribution.
}
\value{
  A list whose components are the input arguments and their assigned
  values.
}

\references{
  S. Baran and S. Lerch, Log-normal distribution based Ensemble Model
  Output Statistics models for probabilistic wind-speed forecasting. 
  \emph{Quarterly Journal of the Royal Meteorological Society} 141:2289--2299,
  2015.
}

\seealso{
  \code{\link{ensembleMOSlognormal}},
  \code{\link{fitMOSlognormal}}
}
\examples{
data("ensBMAtest", package = "ensembleBMA")

ensMemNames <- c("gfs","cmcg","eta","gasp","jma","ngps","tcwb","ukmo")

obs <- paste("MAXWSP10","obs", sep = ".")
ens <- paste("MAXWSP10", ensMemNames, sep = ".")
windTestData <- ensembleData(forecasts = ensBMAtest[,ens],
                             dates = ensBMAtest[,"vdate"],
                             observations = ensBMAtest[,obs],
                             station = ensBMAtest[,"station"],
                             forecastHour = 48,
                             initializationTime = "00")

windTestFitLN <- ensembleMOSlognormal(windTestData, trainingDays = 25,
                                      dates = "2008010100",
                                      control = controlMOSlognormal(maxIter = as.integer(100),
                                                                    scoringRule = "log",
                                                                    optimRule = "BFGS",
                                                                    coefRule= "none", 
                                                                    varRule = "square"))
}
\keyword{models}