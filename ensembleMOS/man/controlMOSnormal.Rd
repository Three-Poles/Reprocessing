\name{controlMOSnormal}
\alias{controlMOSnormal}
\title{
  Control parameters for Gaussian (normal) EMOS models
}
\description{
  Specifies a list of values controling the Gaussian (normal) EMOS fit
  of ensemble forecasts.
}
\usage{
controlMOSnormal(scoringRule = c("crps", "log"),
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
    and variance parameters supplied to \code{optim}. See details.
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
  following Gaussian model is fit by \code{ensembleMOSnormal}: 

  \deqn{ Y ~ N(a + b_1 X_1 + ... + b_m X_m , c + dS^2).}

  \code{B} is the array of fitted regression coefficients \eqn{b_1,
    \ldots, b_m} for each date.  See \link{ensembleMOSnormal} for details.
}
\value{
  A list whose components are the input arguments and their assigned
  values.
}

\references{
T. Gneiting, A. E. Raftery, A. H. Westveld and T. Goldman,
calibrated probabilistic forecasting using ensemble model output
statistics and minimum CRPS estimation.
\emph{Monthly Weather Review} 133:1098--1118, 2005.
}

\seealso{
  \code{\link{ensembleMOSnormal}},
  \code{\link{fitMOSnormal}}
}
\examples{
data("ensBMAtest", package = "ensembleBMA")

ensMemNames <- c("gfs","cmcg","eta","gasp","jma","ngps","tcwb","ukmo")

obs <- paste("T2", "obs", sep = ".")
ens <- paste("T2", ensMemNames, sep = ".")
tempTestData <- ensembleData(forecasts = ensBMAtest[,ens],
                             dates = ensBMAtest[,"vdate"],
                             observations = ensBMAtest[,obs],
                             station = ensBMAtest[,"station"],
                             forecastHour = 48,
                             initializationTime = "00")

tempTestFit <- ensembleMOSnormal(tempTestData, trainingDays = 25, 
                                 dates = "2008010100",
                                 control = controlMOSnormal(maxIter = as.integer(100),
                                                            scoringRule = "log",
                                                            optimRule = "BFGS",
                                                            coefRule= "none", 
                                                            varRule = "square"))
}
\keyword{models}