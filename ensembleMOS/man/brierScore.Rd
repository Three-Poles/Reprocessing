\name{brierScore}
\alias{brierScore}
\alias{brierScore.ensembleMOScsg0}
\alias{brierScore.ensembleMOSgev0}
\alias{brierScore.fitMOScsg0}
\alias{brierScore.fitMOSgev0}
\title{
   Brier Score
}
\description{
   Computes the Brier score for the probability of exceedance of precipitation threshold values 
   for univariate ensemble forecasting models.
}
\usage{
brierScore(fit, ensembleData, thresholds, dates=NULL, \dots)
}
\arguments{
  \item{fit}{
     A model fit to ensemble forecasting data, obtained using
     \code{\link{fitMOS}} or \cr \code{\link{ensembleMOS}}. Only available
     for the censored and shifted gamma, and the censored generalized
     extreme value distribution model.
  }
  \item{ensembleData}{
     An \code{ensembleData} object that includes ensemble forecasts,
     verification observations and possibly dates.
     Missing values (indicated by \code{NA}) are allowed.
     This need not be the data used for the model \code{fit},
     although it must include the same ensemble members.
  }
  \item{thresholds}{
     Threshold values for which the probability of exceedance is evaluated, set to 0 
     to evaluate probability of precipitation forecasts.
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
  \code{BScores} is a vector giving the Brier scores for each instance in the
  data. 
}
\details{
   Note that the Brier scores are only available for EMOS models suitable for
   precipitation accumulation, i.e. the censored and shifted gamma, and the censored generalized
   extreme value distribution EMOS model.
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

obs <- paste("PCP24","obs", sep = ".")
ens <- paste("PCP24", ensMemNames, sep = ".")
prcpTestData <- ensembleData(forecasts = ensBMAtest[,ens],
                             dates = ensBMAtest[,"vdate"],
                             observations = ensBMAtest[,obs],
                             station = ensBMAtest[,"station"],
                             forecastHour = 48,
                             initializationTime = "00")

prcpTestFitCSG0 <- ensembleMOScsg0(prcpTestData, trainingDays = 25,
                                   dates = "2008010100")

brierScore(prcpTestFitCSG0, ensembleData = prcpTestData, 
           thresholds = 0)
}
\keyword{models}