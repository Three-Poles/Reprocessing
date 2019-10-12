\name{fitMOS}
\alias{fitMOS}
\title{
    EMOS model fit to a training set
}
\description{
   Fits an EMOS model to a given training set.
}
\usage{
fitMOS(ensembleData, control = NULL, model = NULL,
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
     A list of control values for the fitting functions. The corresponding
     control function has to be chosen in accordance with the selected 
     \code{model}. For the Gaussian (normal) EMOS model see \link{controlMOSnormal},
     for the truncated normal model see \link{controlMOStruncnormal},
     for the log-normal model see \link{controlMOSlognormal},
     for the censored and shifted gamma model see \link{controlMOScsg0},
     and for the censored generalized extreme value distribution model see \link{controlMOSgev0}.
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
     The model fit will have equal weights and parameters
     within each group.
     The default determines exchangeability from \code{ensembleData}.
  }
}
\value{
  A list with estimated coefficient values. The specific content depends on 
  the chosen \code{model}.
}
\references{
  \emph{Gaussian (normal) EMOS model}: \cr
  T. Gneiting, A. E. Raftery, A. H. Westveld and T. Goldman,
  Calibrated probabilistic forecasting using ensemble model output
  statistics and minimum CRPS estimation.
  \emph{Monthly Weather Review} 133:1098--1118, 2005.
  \cr\cr
  \emph{Truncated (normal) EMOS model}: \cr
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
  \link{fitMOSnormal}
  \link{fitMOStruncnormal}
  \link{fitMOSlognormal}
  \link{fitMOScsg0}
  \link{fitMOSgev0}
  \link{controlMOSnormal}
  \link{controlMOStruncnormal}
  \link{controlMOSlognormal}
  \link{controlMOScsg0}
  \link{controlMOSgev0}
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

tempTrain <- trainingData(tempTestData, trainingDays = 30,
                          date  = "2008010100")

tempTrainFit <- fitMOS(tempTrain, model = "normal")

## equivalent to
##    tempTrainFit <- fitMOSnormal(tempTrain)
}
\keyword{models}