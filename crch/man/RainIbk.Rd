\name{RainIbk}
\alias{RainIbk}

\title{Precipitation Observations and Forecasts for Innsbruck}

\description{
Accumulated 5-8 days precipitation amount for Innsbruck. Data includes GEFS
reforecasts (Hamill et al. 2013) and observations from SYNOP station Innsbruck
Airport (11120) from 2000-01-01 to 2013-09-17.
}

\usage{data("RainIbk")}

\format{
  A data frame with 4977 rows. The first column (\code{rain}) are 3 days
  accumulated precipitation amount observations, Columns 2-12 (\code{rainfc})
  are 5-8 days accumulated precipitation amount forecasts from the individual
  ensemble members.
}

\source{
Observations: http://www.ogimet.com/synops.phtml.en

Reforecasts: http://www.esrl.noaa.gov/psd/forecasts/reforecast2/
}

\references{
  Hamill TM, Bates GT, Whitaker JS, Murray DR, Fiorino M, Galarneau Jr TJ,
    Zhu Y, Lapenta W (2013). NOAA's Second-Generation Global Medium-Range
    Ensemble Reforecast Data Set. \emph{Bulletin of the American Meteorological
    Society}, 94(10), 1553-1565.
}

\examples{
## Spread skill relationship ##

## load and prepare data
data(RainIbk)

## mean and standard deviation of square root transformed ensemble forecasts
RainIbk$sqrtensmean <- 
  apply(sqrt(RainIbk[,grep('^rainfc',names(RainIbk))]), 1, mean)
RainIbk$sqrtenssd <- 
  apply(sqrt(RainIbk[,grep('^rainfc',names(RainIbk))]),  1, sd)

## quintiles of sqrtenssd
sdcat <- cut(RainIbk$sqrtenssd, c(-Inf, quantile(RainIbk$sqrtenssd, 
  seq(0.2,0.8,0.2)), Inf), labels = c(1:5))

## mean forecast errors for each quintile
m <- NULL
for(i in levels(sdcat)) {
  m <- c(m, mean((sqrt(RainIbk$rain)[sdcat == i] -
  RainIbk$sqrtensmean[sdcat == i])^2, na.rm = TRUE))
}

## plot
boxplot((sqrt(rain) - sqrtensmean)^2~sdcat, RainIbk, 
  xlab = "Quintile of ensemble standard deviation", 
  ylab = "mean squared error", main = "Spread skill relationship")

}

\keyword{datasets}
