\name{crch.control}
\alias{crch.control}

\title{Auxiliary Function for Controlling crch Fitting}

\description{
  Auxiliary function for \code{crch} fitting. Specifies a list of values passed
  to \code{\link{optim}}.
}
\usage{
crch.control(method = "BFGS", maxit = NULL, hessian = NULL,
  trace = FALSE, start = NULL, dot = "separate",
  lower = -Inf, upper = Inf, \dots)}
\arguments{
  \item{method}{optimization method passed to \code{\link{optim}}}
  \item{maxit}{the maximum number of iterations. Default is 5000 except for 
    \code{method="boosting"} where the default is 100.}
  \item{hessian}{logical or NULL. If TRUE the numerical Hessian matrix from the 
    \code{\link{optim}} output is used for estimation of the covariance matrix.
    If FALSE no covariance matrix is computed.
    If NULL (the default) the Hessian matrix is computed analytically for 
    \code{dist="gaussian"}, \code{dist="logistic"}, and \code{dist="student"}
    with predefined \code{df}. 
    For \code{dist="student"} without prespecified \code{df}, no analytical 
    solution is available and a numerical Hessian matrix is forced.}
  \item{trace}{non-negative integer. If positive, tracing information on the
  progress of the optimization is produced.}
  \item{start}{initial values for the parameters to be optimized over.}
  \item{dot}{character specifying how to process formula parts with a dot
    (\code{.}) on the right-hand side. This can either be \code{"separate"}
    so that each formula part is expanded separately or \code{"sequential"}
    so that the parts are expanded sequentially conditional on all prior parts.
    Default is \code{"separate"}}
  \item{lower, upper}{bounds on the variables for the \code{"L-BFGS-B"} method,
    or bounds in which to search for method \code{"Brent"}.}
  \item{\dots}{additional parameters passed to \code{\link{optim}}. }
}

\value{
  A list with components named as the arguments.
}

\seealso{\code{\link{crch}}, \code{\link{optim}}}

\keyword{regression}
