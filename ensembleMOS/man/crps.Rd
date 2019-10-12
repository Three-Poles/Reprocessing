\name{crps}
\alias{crps}
\alias{crps.ensembleMOSnormal}
\alias{crps.ensembleMOStruncnormal}
\alias{crps.ensembleMOSlognormal}
\alias{crps.ensembleMOScsg0}
\alias{crps.ensembleMOSgev0}
\alias{crps.fitMOSnormal}
\alias{crps.fitMOStruncnormal}
\alias{crps.fitMOSlognormal}
\alias{crps.fitMOScsg0}
\alias{crps.fitMOSgev0}
\title{
   Continuous Ranked Probability Score
}
\description{
   Computes the continuous ranked probability score (CRPS)
   for univariate ensemble forecasting models.
}
\usage{
crps(fit, ensembleData, dates=NULL, \dots)
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
     The dates for which the CRPS will be computed.
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
  \code{crps} is a matrix giving the CRPS for each instance in the data for 
  both the raw ensemble and the probabilistic forecast. 
}
\details{
   These methods are generic, and can be applied to all ensemble
   forecasting models. Missing values in forecasts and/or observations
   result in \code{NA} values in the CRPS vector.
}
\seealso{
  \code{\link{ensembleMOS}},
  \code{\link{fitMOS}}
}
\references{
T. Gneiting and A. E. Raftery,
Strictly proper scoring rules, prediction and estimation,
\emph{Journal of the American Statistical Association} 102:359--378, 2007.
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
                           dates = "2008010100",
                           model = "normal")

crpsValues <- crps(tempTestFit, tempTestData)
colMeans(crpsValues)
}
\keyword{models}