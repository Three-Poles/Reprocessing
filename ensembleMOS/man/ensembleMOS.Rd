\name{ensembleMOS}
\alias{ensembleMOS}
\alias{ensembleObsLabels}
\alias{ensembleSize}
\alias{ensembleNobs}
\alias{ensembleValidDates}
\alias{ensembleVerifObs}
\alias{ensembleObsLabels.ensembleData}
\alias{ensembleSize.ensembleData}
\alias{ensembleNobs.ensembleData}
\alias{ensembleValidDates.ensembleData}
\alias{ensembleVerifObs.ensembleData}
\alias{matchEnsembleMembers}
\alias{matchEnsembleMembers.ensembleMOSnormal}
\alias{matchEnsembleMembers.fitMOSnormal}
\alias{ensembleMemberLabels}
\alias{ensembleMemberLabels.ensembleData}
\alias{matchITandFH.default}
\alias{matchITandFH}
\alias{getExchangeable}
\alias{getHH}
\title{
  EMOS modeling
}
\description{
  Fits a EMOS model to ensemble forecasts.
   Allows specification of a model, training rule, and forecasting dates.
}
\usage{
ensembleMOS(ensembleData, trainingDays, consecutive = FALSE,
            dates = NULL, control = NULL, warmStart = FALSE,
            model = NULL, exchangeable = NULL)
}
\arguments{
  \item{ensembleData}{
     An \code{ensembleData} object including ensemble forecasts with the
     corresponding verifying observations and their dates.
     Missing values (indicated by \code{NA}) are allowed.
  }
  \item{trainingDays}{
    An integer giving the number of time steps (e.g. days)
    in the training period. There is no default.
  } 
  \item{consecutive}{
    If \code{TRUE} then the sequence of dates in the training set are
    treated as consecutive, i.e. date gaps are ignored
  }
  \item{dates}{
    The dates for which EMOS forecasting models are desired.
    By default, this will be all dates in \code{ensembleData}
    for which modeling is allowed given the training rule.
  }
  \item{control}{
     A list of control values for the fitting functions. The corresponding
     control function has to be chosen in accordance with the selected 
     \code{model}. For the Gaussian (normal) EMOS model see \link{controlMOSnormal},
     for the truncated normal model see \link{controlMOStruncnormal},
     for the log-normal model see \link{controlMOSlognormal},
     for the censored and shifted gamma model see \link{controlMOScsg0},
     and for the censored generalized extreme value distribution model see \link{controlMOSgev0}.
  }
  \item{warmStart}{
    If \code{TRUE}, then starting values for parameters in optimization
    are set to the estimates of the preceding date's fit.
  } 
  \item{model}{
      A character string describing the EMOS model to be fit.
      Current choices are \code{"normal"} (typically used for temperature
      or pressure data), \code{"truncnormal"} (typically used for wind 
      speed data), \code{"lognormal"} (typically used for wind 
      speed data), \code{"csg0"} (typically used for precipitation accumulation data),
      and \code{"gev0"} (typically used for precipitation accumulation data).
      For specific details on model fitting see \link{ensembleMOSnormal},
      \link{ensembleMOStruncnormal}, \link{ensembleMOSlognormal}, 
      \link{ensembleMOScsg0}, or \link{ensembleMOSgev0}.
  }
  \item{exchangeable}{
     A numeric or character vector or factor indicating groups of
     ensemble members that are exchangeable (indistinguishable).
     The model fit will have equal weights and parameters within each
     group.
     The default determines exchangeability from \code{ensembleData}.
}
}
\details{
   If dates are specified in \code{dates} that cannot be forecast
   with the training rule, the corresponding EMOS model parameter
   outputs will be missing (\code{NA}) but not \code{NULL}. \cr
   The training rule uses the number of days corresponding to its
   \code{length} regardless of whether or not the dates are consecutive.
   \cr

}
\value{
 A list containing information on the training (length, lag and the 
 number of instances used for training for each modeling date), the
 exchangeability, and vectors and/or matrics containing the estimated
 regression and variance coefficient values depending on the specified
 \code{model}.
}
\references{
  \emph{Gaussian (normal) EMOS model}: \cr
  T. Gneiting, A. E. Raftery, A. H. Westveld and T. Goldman,
  Calibrated probabilistic forecasting using ensemble model output
  statistics and minimum CRPS estimation.
  \emph{Monthly Weather Review} 133:1098--1118, 2005.
  \cr\cr
  \emph{Truncated normal EMOS model}: \cr
  T. L. Thorarinsdottir and T. Gneiting,
  Probabilistic forecasts of wind speed: 
  Ensemble model output statistics by using 
  heteroscedastic censored regression.
  \emph{Journal of the Royal Statistical Society Series A} 
  173:371--388, 2010.
  \cr\cr
  \emph{Log-normal EMOS model}: \cr
  S. Baran and S. Lerch, Log-normal distribution based Ensemble Model
  Output Statistics models for probabilistic wind-speed forecasting. 
  \emph{Quarterly Journal of the Royal Meteorological Society} 141:2289--2299,
  2015.
  \cr\cr
  \emph{Censored and shifted gamma EMOS model}: \cr
  M. Scheuerer and T. M. Hamill, Statistical post-processing of ensemble precipitation
  forecasts by fitting censored, shifted gamma distributions. 
  \emph{Monthly Weather Review} 143:4578--4596, 2015.\cr
  S. Baran and D. Nemoda, Censored and shifted gamma distribution based EMOS
  model for probabilistic quantitative precipitation forecasting. 
  \emph{Environmetrics} 27:280--292, 2016.  
  \cr\cr
  \emph{Censored generalized extreme value distribution EMOS model}: \cr
  M. Scheuerer, Probabilistic quantitative precipitation forecasting using ensemble
  model output statistics. \emph{Quarterly Journal of the Royal Meteorological
  Society} 140:1086--1096, 2014.
}
\seealso{
  \code{\link{trainingData}},
  \code{\link{ensembleMOSnormal}},
  \code{\link{ensembleMOStruncnormal}},
  \code{\link{ensembleMOSlognormal}},
  \code{\link{ensembleMOScsg0}},
  \code{\link{ensembleMOSgev0}},
  \code{\link{controlMOSnormal}},
  \code{\link{controlMOStruncnormal}},
  \code{\link{controlMOSlognormal}},
  \code{\link{controlMOScsg0}},
  \code{\link{controlMOSgev0}},
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
                           model = "normal")

 ## Same as
 ## tempTestFit <- ensembleMOSnormal(tempTestData, trainingDays = 25)
}
\keyword{models}