\name{cdf}
\alias{cdf}
\alias{cdf.ensembleMOSnormal}
\alias{cdf.ensembleMOStruncnormal}
\alias{cdf.ensembleMOSlognormal}
\alias{cdf.ensembleMOScsg0}
\alias{cdf.ensembleMOSgev0}
\alias{cdf.fitMOSnormal}
\alias{cdf.fitMOStruncnormal}
\alias{cdf.fitMOSlognormal}
\alias{cdf.fitMOScsg0}
\alias{cdf.fitMOSgev0}
\title{
   Cummulative distribution function for ensemble forcasting models
}
\description{
   Computes the cumulative distribution function (CDF)
   of an ensemble forecasting model at observation locations.
}
\usage{
cdf(fit, ensembleData, values, dates = NULL, \dots)
}
\arguments{
   \item{fit}{
     A model fit to ensemble forecasting data, obtained using
     \code{\link{fitMOS}} or \code{\link{ensembleMOS}}.
  }
  \item{ensembleData}{
     An \code{ensembleData} object that includes ensemble forecasts,
     verification observations and possibly dates.
     Missing values (indicated by \code{NA}) are allowed.
     This need not be the data used for the model \code{fit},
     although it must include the same ensemble members.
  }
  \item{values}{
    The vector of desired values at which the CDF of the ensemble
    forecasting model is to be evaluated.
  }
 \item{dates}{
     The dates for which the CDF will be computed.
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
  A matrix of probabilities corresponding to the CDF at the desired values.
  Useful for determining propability of freezing, precipitation, etc.
}
\details{
   This method is generic, and can be applied to any ensemble forecasting
   model obtained using \code{\link{fitMOS}} or \code{\link{ensembleMOS}}. 
   
   For the EMOS models that allow for point masses at 0, i.e. the censored and 
   shifted gamma, and the censored generalized extreme value distribution EMOS model,
   the function contains an addition logical argument \code{randomizeATzero} that specifies
   whether the value of the CDF at zero should be chosen randomly from the interval between
   0 and the value of the CDF at zero. The default choice if \code{FALSE}, setting
   \code{randomizeATzero = TRUE} is practical for computing randomized PIT values.
}
\seealso{
  \code{\link{ensembleMOS}},
  \code{\link{fitMOS}},
  \code{\link{quantileForecast}}
}
\references{
T. Gneiting, A. E. Raftery, A. H. Westveld and T. Goldman,
Calibrated probabilistic forecasting using ensemble model output
statistics and minimum CRPS estimation.
\emph{Monthly Weather Review} 133:1098--1118, 2005.
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
                           model = "normal",
                           dates = "2008010100")

tempTestCDF <- cdf(tempTestFit, tempTestData,
                   values = seq(from=277, to=282, by = 1))

tempTestCDF
}
\keyword{models}