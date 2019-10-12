\name{pars}
\alias{pars}
\alias{pars.ensembleMOSnormal}
\alias{pars.ensembleMOStruncnormal}
\alias{pars.ensembleMOSlognormal}
\alias{pars.ensembleMOScsg0}
\alias{pars.ensembleMOSgev0}
\alias{pars.fitMOSnormal}
\alias{pars.fitMOStruncnormal}
\alias{pars.fitMOSlognormal}
\alias{pars.fitMOScsg0}
\alias{pars.fitMOSgev0}
\title{
  Forecast distribution parameters
}
\description{
  Computes the parameters of the forecast distribution
  resulting for univariate ensemble forecasting models.
}
\usage{
  pars(fit, ensembleData, dates=NULL, \dots)
}
\arguments{
  \item{fit}{
    A model fit to ensemble forecasting data, obtained using
    \code{\link{fitMOS}} or \cr \code{\link{ensembleMOS}}.
  }
  \item{ensembleData}{
    An \code{ensembleData} object that includes ensemble forecasts,
    verification observations and possibly dates.
    Missing values (indicated by \code{NA}) are allowed.
    This need not be the data used for the model \code{fit},
    although it must include the same ensemble members.
  }
  \item{dates}{
    The dates for which the parameters will be computed.
    These dates must be consistent with \code{fit} and \code{ensembleData}.
    The default is to use all of the dates in \code{fit}.
    The dates are ignored if \code{fit} originates from \code{fitMOS},
    which also ignores date information.
  }
  \item{\dots}{
    Included for generic function compatibility.
  }
}
\value{
  \code{pars} is a matrix giving the distribution parameter values
  for each instance in the data. The returned parameters denoted in
  the column names depend on  the model used to generate the 
  \code{fit} object.
}
\details{
  These methods are generic, and can be applied to all ensemble
  forecasting models. 
}
\seealso{
  \code{\link{ensembleMOS}},
  \code{\link{fitMOS}}
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

tempTestFit <- ensembleMOS(tempTestData, trainingDays = 25, 
                           dates = c("2008010100", "2008010200"),
                           model = "normal")

parValues <- pars(tempTestFit, tempTestData)
parValues
}
\keyword{models}