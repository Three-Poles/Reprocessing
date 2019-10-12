\name{quantileForecast}
\alias{quantileForecast}
\alias{quantileForecast.ensembleMOSnormal}
\alias{quantileForecast.fitMOSnormal}
\alias{quantileForecast.ensembleMOStruncnormal}
\alias{quantileForecast.fitMOStruncnormal}
\alias{quantileForecast.ensembleMOSlognormal}
\alias{quantileForecast.fitMOSlognormal}
\alias{quantileForecast.ensembleMOScsg0}
\alias{quantileForecast.fitMOScsg0}
\alias{quantileForecast.ensembleMOSgev0}
\alias{quantileForecast.fitMOSgev0}
\title{
   Quantile forecasts at observation locations
}
\description{
   Computes quantiles for the probability distribution function
   (PDF) for ensemble forecasting models.
}
\usage{
quantileForecast(fit, ensembleData, quantiles = 0.5, dates = NULL, \dots)
}
\arguments{
  \item{fit}{
    A model fit to ensemble forecasting data.
  }
  \item{ensembleData}{
     An \code{ensembleData} object that includes ensemble
     forecasts, verification observations and possibly dates.
     Missing values (indicated by \code{NA}) are allowed. 
     This need not be the data used for the model \code{fit},
     although it must include the same ensemble members.
  }
  \item{quantiles}{
     The vector of desired quantiles for the PDF of the EMOS model.
  }
 \item{dates}{
     The dates for which the quantile forecasts will be computed.
     These dates must be consistent with \code{fit} and
     \code{ensembleData}.
     The default is to use all of the dates in \code{fit}.
     If \code{ensembleData} does not include dates, they will
     be inferred from \code{fit} and \code{dates}.
  }
 \item{\dots}{
     Included for generic function compatibility.
 }
}
\value{
A matrix of forecasts corresponding to the desired quantiles.
}
\details{
   This method is generic, and can be applied to any ensemble forecasting
   model.  This can be used to compute prediction intervals for the PDF. 
}
\seealso{
  \code{\link{ensembleMOS}},
  \code{\link{fitMOS}},
  \code{\link{cdf}}
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
                           dates = c("2008010100", "2008010200"),
                           model = "normal")

tempTestForc <- quantileForecast(tempTestFit, tempTestData)
}
\keyword{models}