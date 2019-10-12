\name{crch.stabsel}
\alias{crch.stabsel}

\title{Auxiliary functions to perform stability selection using boosting.}

\description{
   Auxilirary function which allows to do stability selection on heteroscedastic
   \code{\link{crch}} models based on \code{\link{crch.boost}}.
}
\usage{
crch.stabsel(formula, data, ..., nu = 0.1, q, B = 100, thr = 0.9, 
  maxit = 2000, data_percentage = 0.5)
}
\arguments{
  \item{formula}{a formula expression of the form \code{y ~ x | z} where
  \code{y} is the response and \code{x} and \code{z} are regressor variables
  for the location and the scale of the fitted distribution respectively.}
  \item{data}{an optional data frame containing the variables occurring in the
  formulas.}
  \item{\dots}{Additional attributes to control the \code{\link{crch}} model.
      Note that \code{control} is \emph{not} allowed; \code{crch.stabsel}
      uses \code{crch.boost} by default.}
  \item{nu}{Boosting step size (see \code{\link{crch.boost}}) default is \code{0.1}
      as for \code{\link{crch.boost}} while lower values might yield better results
      frequently and should be considered.}
  \item{q}{Positive \code{numeric}. Maximum number of parameters to be selected
      during each iteration (not including intercepts).}
  \item{B}{\code{numeric}, total number of iterations.}
  \item{thr}{\code{numeric} threshold (\code{(0.5-1.0)}). Used to generate the new
      formula and the computation of the per-family error rate.}
  \item{maxit}{Positive \code{numeric} value. Maximum number for the boosting algorithm.
      If \code{q} is not reached before \code{maxit} the algorithm will stop.}
  \item{data_percentage}{Percentage of data which should be sampled in each of the
      iterations. Default (and suggested) is \code{0.5}.}  
}

\details{
\code{\link{crch.boost}} allows to perform gradient boosting on heteroscedastic
additive models. \code{crch.stabsel} is a wrapper around the core \code{\link{crch.boost}}
algorithm to perform stability selection (see references).

Half of the data set (\code{data}) is sampled \code{B} times to perform boosting
(based on \code{\link{crch.boost}}). Rather than perform the boosting iterations
until a certain stopping criterion is reached (e.g., maximum number of iterations
\code{maxit}) the algorithm stops as soon as \code{q} parameters have been selected.
The number of parameters is computed across both parameters location and scale.
Intercepts are not counted.
}

\value{
  Returns an object of class \code{"stabsel.crch"} containing the stability
  selection summary and the new formula based on the stability selection.
  \item{table}{A table object containing the parameters which have been
   selected and the corresponding frequency of selection.}
  \item{formula.org}{Original formula used to perform the stability selection.}
  \item{formula.new}{New formula based including the coefficients selected during
   stability selection.}
  \item{family}{A list object which contains the distribution-specification from
   the \code{crch.stabsel} call including: \code{dist}, \code{cens}, and \code{truncated}.}
  \item{parameter}{List with the parameters used to perform the stability selection
   including \code{q}, \code{B}, \code{thr}, \code{p}, and \code{PFER} (per-family error rate).}
}

\references{
  Meinhausen N, Buehlmann P (2010). Stability selection.
    \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)}, \bold{72}(4),
    417--473. \doi{10.1111/j.1467-9868.2010.00740.x}.
}

\seealso{\code{\link{crch}}, \code{\link{crch.boost}}}

\examples{
# generate data
suppressWarnings(RNGversion("3.5.0"))
set.seed(5)
x <- matrix(rnorm(1000*20),1000,20)
y <- rnorm(1000, 1 + x[,1] - 1.5 * x[,2], exp(-1 + 0.3*x[,3]))
y <- pmax(0, y)
data <- data.frame(cbind(y, x))

# fit model with maximum likelihood
CRCH1 <- crch(y ~ .|., data = data, dist = "gaussian", left = 0)

# Perform stability selection
stabsel <- crch.stabsel(y ~ .|.,  data = data, dist = "gaussian", left = 0,
           q = 8, B = 5)

# Show stability selection summary
print(stabsel); plot(stabsel)

CRCH2 <- crch(stabsel$formula.new, data = data, dist = "gaussian", left = 0 )
BOOST <- crch(stabsel$formula.new, data = data, dist = "gaussian", left = 0,
              control = crch.boost() )

### AIC comparison
sapply( list(CRCH1,CRCH2,BOOST), logLik )
}
\keyword{regression}
