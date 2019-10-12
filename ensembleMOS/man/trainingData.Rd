\name{trainingData}
\alias{trainingData}
\title{
    Extract Training Data
}
\description{
  Extracts a subset of an \code{ensembleData} object corresponding
  to a given date and number of  training days.
}
\usage{
trainingData(ensembleData, trainingDays, consecutive = FALSE, date)
}
\arguments{
  \item{ensembleData}{
     An \code{ensembleData} object that includes ensemble
     forecasts, observations and dates.
    }
  \item{trainingDays}{
    An integer specifying the number of days in the training period.

  }
  \item{consecutive}{
    If \code{TRUE} then dates in training set are treated as consecutive,
  i.e. date gaps are ignored.
  }
  \item{date}{
     The date for which the training data is desired.
  }
}
\value{
An \code{ensembleData} object corresponding to the training data for
the given date relative to \cr \code{ensembleData}.
}

\seealso{
  \code{\link{ensembleMOSnormal}},
  \code{\link{fitMOSnormal}}
}

\references{
A. E. Raftery, T. Gneiting, F. Balabdaoui and M. Polakowski,
Using Bayesian model averaging to calibrate forecast ensembles,
\emph{Monthly Weather Review} 133:1155-1174, 2005.

J. M. Sloughter, A. E. Raftery, T. Gneiting and C. Fraley,
Probabilistic quantitative precipitation forecasting
using Bayesian model averaging,
\emph{Monthly Weather Review} 135:3309--3320, 2007.

C. Fraley, A. E. Raftery, T. Gneiting and J. M. Sloughter,
\code{ensembleBMA}: An \code{R} Package for Probabilistic Forecasting
using Ensembles and Bayesian Model Averaging,
Technical Report No. 516R, Department of Statistics, University of
Washington, December 2008. \cr
Available at: http://www.stat.washington.edu/research/reports/

C. Fraley, A. E. Raftery and T. Gneiting,
Calibrating multi-model forecast ensembles with exchangeable and missing
members using Bayesian model averaging, \emph{Monthly Weather Review}
138:190-202, 2010.
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
                             date = "2008010100")

tempTrainFit <- fitMOSnormal(tempTrain)
}
\keyword{models}