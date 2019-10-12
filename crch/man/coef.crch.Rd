\name{coef.crch}
\alias{coef.crch}
\alias{vcov.crch}
\alias{print.crch}
\alias{summary.crch}
\alias{print.summary.crch}
\alias{logLik.crch}
\alias{terms.crch}
\alias{model.frame.crch}
\alias{model.matrix.crch}
\alias{residuals.crch}
\alias{fitted.crch}
\alias{getSummary.crch}
\alias{estfun.crch}
\alias{crps.crch}

\title{Methods for CRCH Objects}

\description{
Methods for extracting information from fitted \code{crch} objects.
}

\usage{
\method{coef}{crch}(object, model = c("full", "location", "scale", "df"), \dots)
\method{vcov}{crch}(object, model = c("full", "location", "scale", "df"), \dots)
\method{terms}{crch}(x, model = c("location", "scale", "full"), \dots)
\method{fitted}{crch}(object, type = c("location", "scale"), \dots)
}

\arguments{
  \item{object, x}{an object of class \code{"crch"}.}
  \item{model}{model for which coefficients shall be returned.}
  \item{type}{type of fitted values.}
  \item{\dots}{further arguments passed to or from other methods.}
}

\details{
  In addition to the methods above, a set of standard extractor functions for 
  \code{"crch"} objects is available, including methods to the generic
  functions \code{\link[base]{print}}, \code{\link[base]{summary}}, 
  \code{\link[stats]{logLik}}, and  \code{\link[stats]{residuals}}. 
}

\seealso{\code{\link{crch}}}

\keyword{regression}
