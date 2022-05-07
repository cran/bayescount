\name{normal.params}
\alias{normal.params}
\title{Calculate the Normal Mean and Standard Deviation Using the Log-Normal Mean and Standard Deviation}
\description{
   Function to calculate the equivalent values for the mean and standard deviation of a normal distribution from the mean and standard deviation of the log-normal distribution.  Outputs from this function can be used with the dnorm() function, and with the normal distribution in JAGS.
}
\usage{
normal.params(log.mean, log.sd, coeff.variation=sqrt(exp(log.sd^2)-1))}
\arguments{
   \item{log.mean}{either a single value or vector of values for the mean of the lognormal distribution.}
   \item{log.sd}{either a single value or vector of values for the standard deviation of the lognormal distribution.  Ignored if values are supplied for coeff.variation.}
   \item{coeff.variation}{either a single value or vector of values for the coefficient of dispersion.}
}
\value{A list with elements representing the mean of the normal distribution, the standard deviation of the normal distribution, and the coefficient of variation.
}
\seealso{
   \code{\link{lnormal.params}}
}

\examples{

lmean <- 2.5
lsd <- 0.2
mean <- normal.params(lmean,lsd)[[1]]
sd <- normal.params(lmean,lsd)[[2]]

curve(dlnorm(x, lmean, lsd), from=0, to=25, col="blue")
curve(dnorm(x, mean, sd), from=0, to=25, add=TRUE, col="red")

}
\keyword{methods}
