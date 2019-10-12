\name{predict.crch}
\alias{predict.crch}

\title{Predicted/Fitted Values for CRCH Fits}

\description{
Obtains various types of predictions for \code{crch} models.
}

\usage{
\method{predict}{crch}(object, newdata = NULL, type = c("location", "scale", 
  "response", "parameter", "density", "probability", "quantile", "crps"), 
  na.action = na.pass, at = 0.5, left = NULL, right = NULL, \dots)
}

\arguments{
  \item{object}{an object of class \code{"crch"}.}
  \item{newdata}{an optional data frame in which to look for variables with
  which to predict.}
  \item{type}{type of prediction: \code{"location"} returns the location of 
  the predicted distribution.
  \code{"scale"} returns the scale of the predicted distribution.
  \code{"response"} returns the expected value of the predicted distribution
  (not equal to location for censored and truncated distributions).
  \code{"parameter"} returns a data frame with predicted location and scale 
  parameters.
  \code{"density"} evaluates the predictive density at \code{at}.
  \code{"probability"} evaluates the predictive CDF at \code{at}.
  \code{"quantile"} returns a matrix of predicted quantiles with quantile
  probabilities \code{at}.
  \code{"crps"} returns the CRPS of the predictive distributions at \code{at}.}
  \item{na.action}{a function which indicates what should happen when the data
  contain \code{NA}s. Default is na.pass}
  \item{at}{a vector of values to evaluate the predictive density (\code{type =
  "density"}), probability (\code{type = "probability"}), or CRPS (\code{type =
  "crps"}) or a vector of quantile probabilities used for \code{type =
  "quantile"}.
  Alternatively, with \code{at = "function"} a function is returned that 
  takes \code{at} as an argument.}
  \item{left}{left censoring or truncation point. Only used for \code{type = "quantile"}.
  If NULL, censoring or truncation point is obtained from \code{object}.}
  \item{right}{right censoring or truncation point. Only used for \code{type = "quantile"}.
  If NULL, censoring or truncation point is obtained from \code{object}.}
  \item{\dots}{further arguments passed to or from other methods.}
}

\value{
For type \code{"response"}, \code{"location"}, or \code{"scale"} a vector with
either the location or the scale of the predicted distribution.

For type \code{"quantile"} a matrix of predicted quantiles each column
corresponding to an element of \code{at}. 
}

\seealso{\code{\link{crch}}}

\keyword{regression}
