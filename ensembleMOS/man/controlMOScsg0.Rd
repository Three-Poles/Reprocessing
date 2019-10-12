\name{controlMOScsg0}
\alias{controlMOScsg0}
\title{
  Control parameters for censored and shifted gamma EMOS models
}
\description{
  Specifies a list of values controling the censored and shifted gamma EMOS fit
  of ensemble forecasts.
}
\usage{
controlMOScsg0(scoringRule = c("crps", "log"),
               optimRule = c("Nelder-Mead", "BFGS", "L-BFGS-B"),
               coefRule = c("square", "none", "positive"),
               varRule = c("square", "none"),
               start = list(a = NULL, B = NULL,
                            c = NULL, d = NULL, q = NULL),
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
    Nelder-Mead method, see \link{optim} for details. Note that these 
    options are only available for \code{scoringRule = "log"}. 
    In case of \code{scoringRule = "crps"}, the optimization method
    is set to "L-BFGS-B" by default.
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
    A list of starting parameters, \code{a}, \code{B}, \code{c}, \code{d} 
    and \code{q} specifying initial values for the intercept coefficient
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
  Note that optimMethod options are only available for \code{scoringRule = "log"}. 
  In case of \code{scoringRule = "crps"}, the optimization method
  is set to "L-BFGS-B" by default.
  
 Given an ensemble of size \eqn{m}: \eqn{X_1, \ldots , X_m}, the
  following shifted gamma model left-censored at 0
  is fit by \code{ensembleMOScsg0}: 
      
  \deqn{Y ~ Gamma_0(\kappa,\theta,q)}
      
  where \eqn{Gamma_0} denotes the shifted gamma distribution left-censored at zero,
  with shape \eqn{\kappa}, scale \eqn{\theta} and shift \eqn{q}. The model is 
  parametrized such that the mean \eqn{\kappa\theta} is a linear function
  \eqn{a + b_1 X_1 + \ldots + b_m X_m}
  of the ensemble forecats, and the variance \eqn{\kappa\theta^2} is a linear
  function of the ensemble mean \eqn{c+d \overline{f}}, see \link{ensembleMOScsg0}
  for details.
}
\value{
  A list whose components are the input arguments and their assigned
  values.
}

\references{
  M. Scheuerer and T. M. Hamill, Statistical post-processing of ensemble precipitation
  forecasts by fitting censored, shifted gamma distributions. 
  \emph{Monthly Weather Review} 143:4578--4596, 2015.

  S. Baran and D. Nemoda, Censored and shifted gamma distribution based EMOS
  model for probabilistic quantitative precipitation forecasting. 
  \emph{Environmetrics} 27:280--292, 2016.
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

prcpTestFitCSG0 <- ensembleMOScsg0(prcpTestData, trainingDays = 25,
                                   dates = "2008010100",
                                   control = controlMOScsg0(maxIter = as.integer(100),
                                                            scoringRule = "log",
                                                            optimRule = "Nelder-Mead",
                                                            coefRule= "none", 
                                                            varRule = "square"))

}
\keyword{models}