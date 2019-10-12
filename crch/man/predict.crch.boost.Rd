\name{predict.crch.boost}
\alias{predict.crch.boost}

\title{Predicted/Fitted Values for boosted CRCH Fits}

\description{
Obtains various types of predictions for \code{crch.boost} models.
}

\usage{
\method{predict}{crch.boost}(object, newdata = NULL, mstop = NULL, \dots)
}

\arguments{
  \item{object}{an object of class \code{"crch.boost"}.}
  \item{newdata}{an optional data frame in which to look for variables with
  which to predict.}
  \item{mstop}{stopping iteration. Can be either a character (\code{"max"}, 
  \code{"aic"}, \code{"bic"}, \code{"cv"}) or a numeric value. If not \code{NULL},
  newdata has to be supplied.}
  \item{\dots}{further arguments passed to or from other methods.}
}

\value{
For type \code{"response"}, \code{"location"}, or \code{"scale"} a vector with
either the location or the scale of the predicted distribution.

For type \code{"quantile"} a matrix of predicted quantiles each column
corresponding to an element of \code{at}. 
}

\seealso{\code{\link{crch.boost}},\code{\link{predict.crch}}}

\keyword{regression}
