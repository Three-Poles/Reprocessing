\name{ct}
\alias{dct}
\alias{pct}
\alias{rct}
\alias{qct}

\title{The Censored Student-t Distribution}

\description{
  Density, distribution function, quantile function, and random generation
  for the left and/or right censored student-t distribution with \code{df}
  degrees of freedom.  
}
\usage{
dct(x, location = 0, scale = 1, df, left = -Inf, right = Inf, log = FALSE)

pct(q, location = 0, scale = 1, df, left = -Inf, right = Inf, 
  lower.tail = TRUE, log.p = FALSE)

rct(n, location = 0, scale = 1, df, left = -Inf, right = Inf)

qct(p, location = 0, scale = 1, df, left = -Inf, right = Inf, 
  lower.tail = TRUE, log.p = FALSE)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. If \code{length(n) > 1}, the length is
  taken to be the number required.}
  \item{location}{location parameter.}
  \item{scale}{scale parameter.}
  \item{df}{degrees of freedom (> 0, maybe non-integer).  \code{df = Inf} is
  allowed.}
  \item{left}{left censoring point.}
  \item{right}{right censoring point.}
  \item{log, log.p}{logical; if TRUE, probabilities p are given as log(p).}
  \item{lower.tail}{logical; if TRUE (default), probabilities are P[X <= x]
  otherwise, P[X > x].}
}

\details{
If \code{location} or \code{scale} are not specified they assume the default values
of \code{0} and \code{1}, respectively. \code{left} and \code{right} have the defaults \code{-Inf} and \code{Inf} respectively.

The censored student-t distribution has density \eqn{f(x)}:
\tabular{ll}{
\eqn{T((left - \mu)/\sigma)} \tab if \eqn{x \le left}\cr
\eqn{1 - T((right - \mu)/\sigma)} \tab if \eqn{x \ge right}\cr
\eqn{\tau((x - \mu)/\sigma)/\sigma} \tab if \eqn{left < x < right}
}
     
where \eqn{T} and \eqn{\tau} are the cumulative distribution function
and probability density function of the student-t distribution with
\code{df} degrees of freedom respectively, \eqn{\mu} is the location of the
distribution, and \eqn{\sigma} the scale.
}

\value{
 \code{dct} gives the density, \code{pct} gives the distribution
     function, \code{qct} gives the quantile function, and \code{rct}
     generates random deviates.
}

\seealso{\code{\link{dt}}}

\keyword{distribution}
