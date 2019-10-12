\name{clogis}
\alias{dclogis}
\alias{pclogis}
\alias{rclogis}
\alias{qclogis}

\title{The Censored Logistic Distribution}

\description{
  Density, distribution function, quantile function, and random generation
  for the left and/or right censored logistic distribution.  
}
\usage{
dclogis(x, location = 0, scale = 1, left = -Inf, right = Inf, log = FALSE)

pclogis(q, location = 0, scale = 1, left = -Inf, right = Inf, 
  lower.tail = TRUE, log.p = FALSE)

rclogis(n, location = 0, scale = 1, left = -Inf, right = Inf)

qclogis(p, location = 0, scale = 1, left = -Inf, right = Inf,
  lower.tail = TRUE, log.p = FALSE)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. If \code{length(n) > 1}, the length is
  taken to be the number required.}
  \item{location}{location parameter.}
  \item{scale}{scale parameter.}
  \item{left}{left censoring point.}
  \item{right}{right censoring point.}
  \item{log, log.p}{logical; if TRUE, probabilities p are given as log(p).}
  \item{lower.tail}{logical; if TRUE (default), probabilities are P[X <= x]
  otherwise, P[X > x].}
}

\details{
If \code{location} or \code{scale} are not specified they assume the default values
of \code{0} and \code{1}, respectively. \code{left} and \code{right} have the defaults \code{-Inf} and \code{Inf} respectively.

The censored logistic distribution has density \eqn{f(x)}:
\tabular{ll}{
\eqn{\Lambda((left - \mu)/\sigma)} \tab if \eqn{x \le left}\cr
\eqn{1 - \Lambda((right - \mu)/\sigma)} \tab if \eqn{x \ge right}\cr
\eqn{\lambda((x - \mu)/\sigma)/\sigma} \tab if \eqn{left < x < right}
}
     
where \eqn{\Lambda} and \eqn{\lambda} are the cumulative distribution function
and probability density function of the standard logistic distribution
respectively, \eqn{\mu} is the location of the distribution, and \eqn{\sigma} 
the scale.
}

\value{
 \code{dclogis} gives the density, \code{pclogis} gives the distribution
     function, \code{qclogis} gives the quantile function, and \code{rclogis}
     generates random deviates.
}

\seealso{\code{\link{dlogis}}}

\keyword{distribution}
