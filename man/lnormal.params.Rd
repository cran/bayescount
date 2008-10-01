\name{lnormal.params}
\alias{lnormal.params}
\title{CALCULATE THE LOG-NORMAL MEAN AND STANDARD DEVIATION USING THE NORMAL MEAN AND STANDARD DEVIATION}
\description{
   Function to calculate the equivalent values for the mean and standard deviation of a log-normal distribution from the mean and standard deviation of the normal distribution.  Outputs from this function can be used with the dlnorm() function, and with the lognormal distribution in JAGS.
}
\usage{
lnormal.params(mean, sd)}
\arguments{
   \item{mean}{either a single value or vector of values for the mean of the normal distribution}
   \item{sd}{either a single value or vector of values for the standard deviation of the normal distribution}
}
\value{Either a vector of length 2 representing the mean of the lognormal distribution as the first element and the standard deviation of the lognormal distribution as the second element, or a matrix with column 1 as the mean and column 2 as the standard deviation.
}
\seealso{
   \code{\link{bayescount}}
   \code{\link{normal.params}}
}

\author{Matthew Denwood \email{m.denwood@vet.gla.ac.uk} funded as part of the DEFRA VTRI project 0101.}

\examples{

\dontrun{
mean <- 10
sd <- 2
lmean <- lnormal.params(mean,sd)[1]
lsd <- lnormal.params(mean,sd)[2]

curve(dnorm(x, mean, sd), from=0, to=20)
x11()
curve(dlnorm(x, lmean, lsd), from=0, to=20)
}
}
\keyword{methods}