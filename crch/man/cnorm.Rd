\name{cnorm}
\alias{dcnorm}
\alias{pcnorm}
\alias{rcnorm}
\alias{qcnorm}

\title{The Censored Normal Distribution}

\description{
  Density, distribution function, quantile function, and random generation
  for the left and/or right censored normal distribution.  
}
\usage{
dcnorm(x, mean = 0, sd = 1, left = -Inf, right = Inf, log = FALSE)

pcnorm(q, mean = 0, sd = 1, left = -Inf, right = Inf, 
  lower.tail = TRUE, log.p = FALSE)

rcnorm(n, mean = 0, sd = 1, left = -Inf, right = Inf)

qcnorm(p, mean = 0, sd = 1, left = -Inf, right = Inf,
  lower.tail = TRUE, log.p = FALSE)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. If \code{length(n) > 1}, the length is
  taken to be the number required.}
  \item{mean}{vector of means.}
  \item{sd}{vector of standard deviations.}
  \item{left}{left censoring point.}
  \item{right}{right censoring point.}
  \item{log, log.p}{logical; if TRUE, probabilities p are given as log(p).}
  \item{lower.tail}{logical; if TRUE (default), probabilities are P[X <= x]
  otherwise, P[X > x].}
}

\details{
If \code{mean} or \code{sd} are not specified they assume the default values
of \code{0} and \code{1}, respectively. \code{left} and \code{right} have the defaults \code{-Inf} and \code{Inf} respectively.

The censored normal distribution has density \eqn{f(x)}:
\tabular{ll}{
\eqn{\Phi((left - \mu)/\sigma)} \tab if \eqn{x \le left}\cr
\eqn{1 - \Phi((right - \mu)/\sigma)} \tab if \eqn{x \ge right}\cr
\eqn{\phi((x - \mu)/\sigma)/\sigma} \tab if \eqn{left < x < right}
}
     
where \eqn{\Phi} and \eqn{\phi} are the cumulative distribution function
and probability density function of the standard normal distribution
respectively, \eqn{\mu} is the mean of the distribution, and \eqn{\sigma} 
the standard deviation.
}

\value{
 \code{dcnorm} gives the density, \code{pcnorm} gives the distribution
     function, \code{qcnorm} gives the quantile function, and \code{rcnorm}
     generates random deviates.
}

\seealso{\code{\link{dnorm}}}

\keyword{distribution}
