\name{normal.params}
\alias{normal.params}
\title{CALCULATE THE NORMAL MEAN AND STANDARD DEVIATION USING THE LOG-NORMAL MEAN AND STANDARD DEVIATION}
\description{
   Function to calculate the equivalent values for the mean and standard deviation of a normal distribution from the mean and standard deviation of the log-normal distribution.  Outputs from this function can be used with the dnorm() function, and with the normal distribution in JAGS.
}
\usage{
normal.params(log.mean, log.sd)}
\arguments{
   \item{log.mean}{either a single value or vector of values for the mean of the log-normal distribution}
   \item{log.sd}{either a single value or vector of values for the standard deviation of the log-normal distribution}
}
\value{Either a vector of length 2 representing the mean of the normal distribution as the first element and the standard deviation of the normal distribution as the second element, or a matrix with column 1 as the mean and column 2 as the standard deviation.
}
\seealso{
   \code{\link{bayescount}}
   \code{\link{lnormal.params}}
}

\author{Matthew Denwood \email{m.denwood@vet.gla.ac.uk} funded as part of the DEFRA VTRI project 0101.}

\examples{

\dontrun{
lmean <- 2.5
lsd <- 0.2
mean <- normal.params(lmean,lsd)[1]
sd <- normal.params(lmean,lsd)[2]

curve(dlnorm(x, lmean, lsd), from=0, to=25)
x11()
curve(dnorm(x, mean, sd), from=0, to=25)
}

}
\keyword{methods}