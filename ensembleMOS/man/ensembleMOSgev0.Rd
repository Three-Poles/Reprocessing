\name{ensembleMOSgev0}
\alias{ensembleMOSgev0}
\title{
  Censored generalized extreme value distribution EMOS modeling
}
\description{
  Fits a Censored generalized extreme value distribution EMOS model to ensemble forecasts for
  specified dates.
}
\usage{
ensembleMOSgev0(ensembleData, trainingDays, consecutive = FALSE,
                dates = NULL, control = controlMOSgev0(),
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
    function  \link{controlMOSgev0}. For details
    and default values, see \link{controlMOSgev0}.
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
    \item{s}{
    A vector of fitted EMOS coefficients for \eqn{p_0} for each date, 
    see details.
    }
    \item{c,d}{
    The fitted coefficients for the shape parameter, see details.
    }
    \item{q}{
    Fitted shape parameter, see details.
    }
}
    
\details{
  Given an ensemble of size \eqn{m}: \eqn{X_1, \ldots , X_m}, the
  following generalized extreme value distribution EMOS
  model left-censored at 0 is fit by \code{ensembleMOSgev0}: 
      
  \deqn{ Y ~ GEV_0(\mu,\sigma,q)}
      
  where \eqn{GEV_0} denotes the generalized extreme value distribution
  left-censored at zero,
  with location \eqn{\mu}, scale \eqn{\sigma} and shape \eqn{q}. The model is 
  parametrized such that the mean \eqn{m} is a linear function
  \eqn{a + b_1 X_1 + \ldots + b_m X_m + s p_0}
  of the ensemble forecats, where \eqn{p_0} denotes the ratio of ensemble forecasts
  that are exactly 0, and the shape parameter \eqn{\sigma} is a linear
  function of the ensemble variance \eqn{c + d MD(X_1,\ldots,X_m)}, where
  \eqn{MD(X_1,\ldots,X_m)} is Gini's mean difference. 
  See \link{ensembleMOSgev0} for details.

  \code{B} is a vector of fitted regression coefficients: \eqn{b_1,
  \ldots, b_m}.  Specifically, \eqn{a, b_1,\ldots, b_m, s, c, d, q} are
  fitted to optimize
  the mean CRPS over the specified training period using
  \code{optim}.
}
\references{
  M. Scheuerer, Probabilistic quantitative precipitation forecasting using ensemble
  model output statistics. \emph{Quarterly Journal of the Royal Meteorological
  Society} 140:1086--1096, 2014.
}
\seealso{
    \code{\link{controlMOSgev0}},
    \code{\link{fitMOSgev0}}
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
                             
prcpTestFitGEV0 <- ensembleMOSgev0(prcpTestData, trainingDays = 25,
                                   dates = "2008010100")
}
\keyword{models}