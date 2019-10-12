\name{fitMOStruncnormal}
\alias{fitMOStruncnormal}
\title{
  Truncated normal EMOS model fit to a training set
}
\description{
  Fits a truncated normal EMOS model to a given training set.
}
\usage{
fitMOStruncnormal(ensembleData, control = controlMOStruncnormal(),
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
    function  \link{controlMOStruncnormal}. For details
    and default values, see \link{controlMOStruncnormal}.
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
  Given an ensemble of size \eqn{m}: \eqn{X_1, \ldots , X_m}, the
  following truncated normal model is fit by \code{ensembleMOStruncnormal}: 
      
  \deqn{Y ~ N_0(a + b_1 X_1 + ... + b_m X_m, c + dS^2),}
      
  where \eqn{N_0} denotes the normal distribution truncated at zero,
  with location \eqn{a + b_1 X_1 + ... + b_m X_m} and squared scale
  \eqn{c + dS^2}.
  \code{B} is a vector of fitted regression coefficients: \eqn{b_1,
  \ldots, b_m}.  Specifically, \eqn{a, b_1,\ldots, b_m, c, d} are
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
    The fitted parameters for the squared scale, see details.
  }
}

\references{
  T. L. Thorarinsdottir and T. Gneiting,
  Probabilistic forecasts of wind speed: 
  Ensemble model output statistics by using 
  heteroscedastic censored regression.
  \emph{Journal of the Royal Statistical Society Series A} 
  173:371--388, 2010.
}
\seealso{

  \code{\link{controlMOStruncnormal}},
  \code{\link{ensembleMOStruncnormal}},
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

windTestFit <- fitMOStruncnormal(windTrain)
}
\keyword{models}