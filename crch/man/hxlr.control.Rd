\name{hxlr.control}
\alias{hxlr.control}

\title{Auxiliary Function for Controlling HXLR Fitting}

\description{
  Auxiliary function for \code{hxlr} fitting. Specifies a list of values passed
  to \code{\link{optim}}.
}
\usage{
hxlr.control(method = "BFGS", maxit = 5000, hessian = TRUE, 
  trace = FALSE, start = NULL, \dots)}
\arguments{
  \item{method}{optimization method used in \code{\link{optim}}}
  \item{maxit}{the maximum number of iterations.}
  \item{hessian}{logical. Should a numerically differentiated Hessian matrix be
  returned?}
  \item{trace}{non-negative integer. If positive, tracing information on the
  progress of the optimization is produced.}
  \item{start}{initial values for the parameters to be optimized over.}
  \item{\dots}{Additional parameters passed to \code{\link{optim}}. }
}

\value{
  A list with components named as the arguments.
}

\seealso{\code{\link{hxlr}}, \code{\link{optim}}}

\keyword{regression}
