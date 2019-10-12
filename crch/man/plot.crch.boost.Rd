\name{plot.crch.boost}
\alias{plot.crch.boost}

\title{Plot coefficient paths of boosted CRCH objects.}

\description{
Plot paths of coefficients or log-likelihood contributions for \code{crch.boost} models.
}

\usage{
\method{plot}{crch.boost}(x, loglik = FALSE, 
  standardize = TRUE, which = c("both", "location", "scale"), 
  mstop = NULL, coef.label = TRUE, col = NULL, \dots)
}

\arguments{
  \item{x}{an object of class \code{"crch.boost"}.}
  \item{loglik}{logical whether log-likelihood contribution shall be plotted
  instead of coefficient value.}
  \item{standardize}{logical whether coefficients shall be standardized.
  Not used if \code{loglik = TRUE} }
  \item{which}{which coefficients: \code{"location"} and \code{"scale"} plots
  only the coefficients for the location and scale part of the model 
  respectively. \code{"both"} plots the coefficient paths of both parts in one
  graph.}
  \item{mstop}{Stopping iteration for which a vertical line is plotted. 
  Possible choices are \code{"max"}, \code{"aic"}, \code{"bic"}, 
  \code{"cv"}, \code{"all"}, or \code{"no"}. Default is the stopping iteration
  used for fitting.}
  \item{coef.label}{logical whether paths shall be labeled.}
  \item{col}{Color(s) for the paths. If \code{which="both"} a vector of 
  two colors where the paths for the location are plotted in the first
  color and for the scale in the second color.}
  \item{\dots}{further arguments passed to \code{\link{plot.ts}}.}
}

\seealso{\code{\link{crch.boost}},\code{\link{plot.ts}}}

\keyword{regression}
