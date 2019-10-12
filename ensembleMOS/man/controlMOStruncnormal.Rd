\name{controlMOStruncnormal}
\alias{controlMOStruncnormal}
\title{
  Control parameters for truncated normal EMOS models
}
\description{
  Specifies a list of values controling the truncated normal EMOS fit
  of ensemble forecasts.
}
\usage{
controlMOStruncnormal(scoringRule = c("crps", "log"),
                      optimRule = c("BFGS","Nelder-Mead"),
                      coefRule = c("square", "none", "positive"),
                      varRule = c("square", "none"),
                      start = list(a = NULL, B = NULL,
                                   c = NULL, d = NULL),
                      maxIter = Inf)
}
\arguments{
  \item{scoringRule}{
    The scoring rule to be used in optimum score estimation. Options
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
    Method to control non-negativity of the scale parameters.
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
  
  
  Given an ensemble of size \eqn{m}: \eqn{X_1, \ldots , X_m}, the
  following truncated normal model is fit by \code{ensembleMOStruncnormal}: 
      
  \deqn{Y ~ N_0(a + b_1 X_1 + ... + b_m X_m, c + dS^2),}
      
  where \eqn{N_0} denotes the normal distribution truncated at zero,
  with location \eqn{a + b_1 X_1 + ... + b_m X_m} and squared scale
  \eqn{c + dS^2}.
  \code{B} is a vector of fitted regression coefficients \eqn{b_1,
  \ldots, b_m}. See \link{ensembleMOStruncnormal} for details.  
}
\value{
  A list whose components are the input arguments and their assigned
  values.
}

\references{
T. L. Thorarinsdottir and T. Gneiting,
Probabilistic forecasts of wind speed: 
Ensemble model output statistics by using 
heteroscedastic censored regression.
\emph{Journal of the Royal Statistical Society Series A} 
173:371--388, 2010.
}

\seealso{
  \code{\link{ensembleMOStruncnormal}},
  \code{\link{fitMOStruncnormal}}
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

windTestFitTN <- ensembleMOStruncnormal(windTestData, trainingDays = 25,
                                        dates = "2008010100",
                                        control = controlMOStruncnormal(maxIter = as.integer(100),
                                                                        scoringRule = "log",
                                                                        optimRule = "BFGS",
                                                                        coefRule= "none", 
                                                                        varRule = "square"))
}
\keyword{models}