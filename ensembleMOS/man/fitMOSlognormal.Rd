\name{fitMOSlognormal}
\alias{fitMOSlognormal}
\title{
  Log-normal EMOS model fit to a training set
}
\description{
  Fits a log-normal EMOS model to a given training set.
}
\usage{
fitMOSlognormal(ensembleData, control = controlMOSlognormal(),
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
    function  \link{controlMOSlognormal}. For details
    and default values, see \link{controlMOSlognormal}.
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
}
\value{

A list with the following output components:
  \item{a}{
    The fitted intercept.
  }
  \item{B}{
    The fitted EMOS coefficients.
  }
  \item{c,d}{
    The fitted parameters for the variance, see details.
  }
}

\references{
  S. Baran and S. Lerch, Log-normal distribution based Ensemble Model
  Output Statistics models for probabilistic wind-speed forecasting. 
  \emph{Quarterly Journal of the Royal Meteorological Society} 141:2289--2299,
  2015.
}
\seealso{

  \code{\link{controlMOSlognormal}},
  \code{\link{ensembleMOSlognormal}},
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
                             
windTrain <- trainingData(windTestData, trainingDays = 30,
                          date = "2008010100")

windTestFit <- fitMOSlognormal(windTrain)
}
\keyword{models}