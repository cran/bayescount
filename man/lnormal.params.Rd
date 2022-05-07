\name{lnormal.params}
\alias{lnormal.params}
\title{Calculate the Log-Normal Mean and Standard Deviation Using the Normal Mean and Standard Deviation}
\description{
   Function to calculate the equivalent values for the mean and standard deviation of a log-normal distribution from the mean and standard deviation of the normal distribution.  Outputs from this function can be used with the dlnorm() function, and with the lognormal distribution in JAGS.
}
\usage{
lnormal.params(mean, sd, coeff.variation)}
\arguments{
   \item{mean}{either a single value or vector of values for the mean of the normal distribution.}
   \item{sd}{either a single value or vector of values for the standard deviation of the normal distribution.  Ignored if values are supplied for coeff.variation.}
   \item{coeff.variation}{either a single value or vector of values for the coefficient of dispersion.}
}
\value{A list with elements representing the mean of the lognormal distribution, the standard deviation of the lognormal distribution, and the coefficient of variation.
}
\seealso{
   \code{\link{normal.params}}
}

\examples{

mean <- 10
sd <- 2
lmean <- lnormal.params(mean,sd)[[1]]
lsd <- lnormal.params(mean,sd)[[2]]

curve(dnorm(x, mean, sd), from=0, to=20, col="blue")
curve(dlnorm(x, lmean, lsd), from=0, to=20, add=TRUE, col="red")
}
\keyword{methods}
