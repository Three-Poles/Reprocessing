\name{fitMOSnormal}
\alias{fitMOSnormal}
\title{
  Gaussian (normal) EMOS model fit to a training set
}
\description{
  Fits a Gaussian (normal) EMOS model to a given training set.
}
\usage{
fitMOSnormal(ensembleData, control = controlMOSnormal(),
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
    function  \link{controlMOSnormal}. For details
    and default values, see \link{controlMOSnormal}.
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
  following Gaussian model is fit by \code{ensembleMOSnormal}: 

  \deqn{ Y ~ N(a + b_1 X_1 + ... + b_m X_m , c + dS^2).}

  \code{B} is a vector of fitted regression coefficients: \eqn{b_1,
    \ldots ,b_m}.  Specifically, \eqn{a,b_1,\ldots , b_m, c, d} are
    fitted to optimize
    \code{control$scoringRule} over the specified training period using
    \code{optim} with \code{method = control$optimRule}.
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
    The fitted variance parameters, see details.
  }
}

\references{
T. Gneiting, A. E. Raftery, A. H. Westveld and T. Goldman,
Calibrated probabilistic forecasting using ensemble model output
statistics and minimum CRPS estimation.
\emph{Monthly Weather Review} 133:1098--1118, 2005.
}
\seealso{

  \code{\link{controlMOSnormal}},
  \code{\link{ensembleMOSnormal}},
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