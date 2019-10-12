\name{coef.hxlr}
\alias{coef.hxlr}
\alias{vcov.hxlr}
\alias{print.hxlr}
\alias{summary.hxlr}
\alias{print.summary.hxlr}
\alias{logLik.hxlr}
\alias{terms.hxlr}


\title{Methods for HXLR Objects}

\description{
Methods for extracting information from fitted \code{hxlr} 
objects.
}

\usage{
\method{coef}{hxlr}(object, model = c("full", "intercept", "location", "scale"),
  type = c("CLM", "latent"), \dots)
\method{vcov}{hxlr}(object, model = c("full", "intercept", "location", "scale"), 
  type = c("CLM", "latent"), \dots)
\method{terms}{hxlr}(x, model = c("full", "location", "scale"), \dots)
}

\arguments{
  \item{object, x}{an object of class \code{"hxlr"}.}
  \item{model}{model for which coefficients shall be returned.}
  \item{type}{type of coefficients. Default are CLM type coefficients.  For
  type \code{"latent"} coefficients are converted in coefficients for location
  and scale of the latent distribution (analog to \code{crch} models).}
  \item{\dots}{further arguments passed to or from other methods.}
}

\details{
  In addition to the methods above, a set of standard extractor functions for 
  \code{"hxlr"} objects is available, including methods to the generic
  functions \code{\link[base]{print}}, \code{\link[base]{summary}}, and
  \code{\link[stats]{logLik}}.
}

\seealso{\code{\link{hxlr}}}

\keyword{regression}
