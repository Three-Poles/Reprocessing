\name{coef.crch.boost}
\alias{coef.crch.boost}
\alias{print.crch.boost}
\alias{summary.crch.boost}
\alias{print.summary.crch.boost}
\alias{logLik.crch.boost}

\title{Methods for boosted CRCH Objects}

\description{
Methods for extracting information from fitted \code{crch.boost} objects.
}

\usage{
\method{coef}{crch.boost}(object, model = c("full", "location", "scale", "df"), 
  mstop = NULL, zero.coefficients = FALSE, standardize = FALSE, \dots)
\method{print}{crch.boost}(x, digits = max(3, getOption("digits") - 3),
  mstop = NULL, zero.coefficients = FALSE, \dots)
\method{summary}{crch.boost}(object, mstop = NULL, zero.coefficients = FALSE, \dots)
\method{logLik}{crch.boost}(object, mstop = NULL, \dots)
}

\arguments{
  \item{object, x}{an object of class \code{"crch.boost"}.}
  \item{model}{model for which coefficients shall be returned.}
  \item{mstop}{stopping iteration for which coefficients shall be returned.
  Can be either a character (\code{"max"}, \code{"aic"}, \code{"bic"}, \code{"cv"}) 
  or a numeric value.}
  \item{zero.coefficients}{logical whether zero coefficients are returned.}
  \item{standardize}{logical whether coefficients shall be standardized.}
  \item{digits}{the number of significant digits to use when printing.}
  \item{\dots}{further arguments passed to or from other methods.}
}

\details{
  In addition to the methods above, the \code{"crch"} methods 
  \code{\link{terms}}, \code{\link{model.frame}}, \code{\link{model.matrix}},
  \code{\link{residuals}}, and \code{\link{fitted}} can be used also for
  \code{"crch.boost"} objects .
}

\seealso{\code{\link{crch.boost}},\code{\link{coef.crch}}}

\keyword{regression}
