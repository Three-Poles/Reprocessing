\name{tnorm}
\alias{dtnorm}
\alias{ptnorm}
\alias{rtnorm}
\alias{qtnorm}

\title{The Truncated Normal Distribution}

\description{
  Density, distribution function, quantile function, and random generation
  for the left and/or right truncated normal distribution.  
}
\usage{
dtnorm(x, mean = 0, sd = 1, left = -Inf, right = Inf, log = FALSE)

ptnorm(q, mean = 0, sd = 1, left = -Inf, right = Inf, 
  lower.tail = TRUE, log.p = FALSE)

rtnorm(n, mean = 0, sd = 1, left = -Inf, right = Inf)

qtnorm(p, mean = 0, sd = 1, left = -Inf, right = Inf, 
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

The truncated normal distribution has density
\deqn{f(x) = 1/\sigma \phi((x - \mu)/\sigma) / 
  (\Phi((right - \mu)/\sigma) - \Phi((left - \mu)/\sigma))}
for \eqn{left \le x \le right}, and 0 otherwise.
     
\eqn{\Phi} and \eqn{\phi} are the cumulative distribution function
and probability density function of the standard normal distribution
respectively, \eqn{\mu} is the mean of the distribution, and \eqn{\sigma} 
the standard deviation.
}

\value{
 \code{dtnorm} gives the density, \code{ptnorm} gives the distribution
     function, \code{qtnorm} gives the quantile function, and \code{rtnorm}
     generates random deviates.
}

\seealso{\code{\link{dnorm}}}

\keyword{distribution}
