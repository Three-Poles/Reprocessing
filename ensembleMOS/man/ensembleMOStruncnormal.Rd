\name{ensembleMOStruncnormal}
\alias{ensembleMOStruncnormal}
\title{
  Truncated normal EMOS modeling
}
\description{
  Fits a truncated normal EMOS model to ensemble forecasts for
  specified dates.
}
\usage{
ensembleMOStruncnormal(ensembleData, trainingDays, consecutive = FALSE,
                       dates = NULL, control = controlMOStruncnormal(),
                       warmStart = FALSE, exchangeable = NULL)
}
\arguments{
  \item{ensembleData}{
    An \code{ensembleData} object including ensemble forecasts with
    the corresponding verifying observations and their dates.
    Missing values (indicated by \code{NA}) are allowed.
  }
  \item{trainingDays}{
    An integer giving the number of time steps (e.g. days)
    in the training period. There is no default.
  }
  \item{consecutive}{
    If \code{TRUE} then the sequence of dates in the training set are
    treated as consecutive, i.e. date gaps are ignored.
  }
  \item{dates}{
    The dates for which EMOS forecasting models are desired.
    By default, this will be all dates in \code{ensembleData}
    for which modeling is allowed given the training rule.
  }
  \item{control}{
    A list of control values for the fitting functions specified via the
    function  \link{controlMOStruncnormal}. For details
    and default values, see \link{controlMOStruncnormal}.
  }
  \item{warmStart}{
    If \code{TRUE}, then starting values for parameters in optimization
    are set to the estimates of the preceding date's fit.
  }
  \item{exchangeable}{
    A numeric or character vector or factor indicating groups of
    ensemble members that are exchangeable (indistinguishable).
    The modeling will have equal parameters within each group.
    The default determines exchangeability from \code{ensembleData}.
    }
  }
    
\value{
A list with the following output components:
\item{training}{
    A list containing information on the training length and lag and
    the number of instances used for training for each modeling date.
    }
    \item{a}{
    A vector of fitted EMOS intercept parameters for each date.
    }
    \item{B}{
    A matrix of fitted EMOS coefficients for each date.
    }
    \item{c,d}{
    The fitted parameters for the squared scale, see details.
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
    \code{\link{fitMOStruncnormal}}
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

windTestFitTN <- ensembleMOStruncnormal(windTestData, trainingDays = 25)
}
\keyword{models}