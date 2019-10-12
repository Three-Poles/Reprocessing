\name{fitMOScsg0}
\alias{fitMOScsg0}
\title{
  Censored and shifted gamma EMOS modeling
}
\description{
  Fits a censored and shifted gamma EMOS model to a given training set.
}
\usage{
fitMOScsg0(ensembleData, control = controlMOScsg0(),
           exchangeable = NULL)
}
\arguments{
  \item{ensembleData}{
     An \code{ensembleData} object including ensemble forecasts and
     verification observations.
     Missing values (indicated by \code{NA}) are allowed. Dates are ignored
     if they are included.  This is the training set for the model.
}
  \item{control}{
    A list of control values for the fitting functions specified via the
    function  \link{controlMOScsg0}. For details
    and default values, see \link{controlMOScsg0}.
}
  \item{exchangeable}{
    An optional numeric or character vector or factor indicating groups of
    ensemble members that are exchangeable (indistinguishable).
    The models have equal EMOS coefficients within each group.
    If supplied, this argument will override any specification of
    exchangeability in \code{ensembleData}.
}
}
\details{
  Given an ensemble of size \eqn{m}: \eqn{X_1, \ldots , X_m}, the
  following shifted gamma model left-censored at 0
  is fit by \code{ensembleMOScsg0}: 
      
  \deqn{Y ~ Gamma_0(\kappa,\theta,q)}
      
  where \eqn{Gamma_0} denotes the shifted gamma distribution left-censored at zero,
  with shape \eqn{\kappa}, scale \eqn{\theta} and shift \eqn{q}. The model is 
  parametrized such that the mean \eqn{\kappa\theta} is a linear function
  \eqn{a + b_1 X_1 + \ldots + b_m X_m}
  of the ensemble forecats, and the variance \eqn{\kappa\theta^2} is a linear
  function of the ensemble mean \eqn{c+d \overline{f}}, see Baran and Nemoda (2016) 
  for details.

  \code{B} is a vector of fitted regression coefficients: \eqn{b_1,
  \ldots, b_m}.  Specifically, \eqn{a, b_1,\ldots, b_m, c, d} are
  fitted to optimize
  \code{control$scoringRule} over the specified training period using
  \code{optim} with \code{method = control$optimRule}.
}
\value{
A list with the following output components:
\item{training}{
    A list containing information on the training length and lag and
    the number of instances used for training for each modeling date.
    }
    \item{a}{
    A vector of fitted EMOS intercept parameters for each date.
    }
    \item{B}{
    A matrix of fitted EMOS coefficients for each date.
    }
    \item{c,d}{
    The fitted parameters for the variance, see details.
    }
    \item{q}{
    Fitted shift parameter, see details.
    }
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

  \code{\link{controlMOScsg0}},
  \code{\link{ensembleMOScsg0}},
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

prcpTrain <- trainingData(prcpTestData, trainingDays = 30,
                             date = "2008010100")

prcpTestFit <- fitMOScsg0(prcpTrain)
}
\keyword{models}