\name{maximise.likelihood}
\alias{maximise.likelihood}
\title{CALCULATE THE MAXIMUM LIKELIHOOD PARAMETERS OF A CONTINUOUS OR COUNT DISTRIBUTION}
\description{
   Crude function to maximise the likelihood of one of the following distributions:  Poisson (P), gamma (G), lognormal (L), Weibull (W), gamma Poisson (GP), lognormal Poisson (LP), Weibull Poisson (WP), all with or without zero-inflation (ZI).  Uses the likelihood() function to calculate the likelihood at each iteration.  For mixture models, the likelihood is calculated for the data by integrating over each possible value of lambda for each data point, which may take some time for large datasets.  Starting values for each parameter are optional, but may improve the speed and reliability of the function if appropriate values are provided.  If missing, starting values will be calculated from the data.  This function is provided for interest only, and is vastly inferior as a data analysis tool to the Bayesian MCMC methods used by bayescount().  *THIS SOFTWARE IS INTENDED FOR EDUCATIONAL PURPOSES ONLY AND SHOULD NOT BE RELIED UPON FOR REAL WORLD APPLICATIONS*
}
\usage{
maximise.likelihood(data=stop("Data must be specified"), model=stop("Please specify a distribution"), mean=NA, variance=NA, 
   zi=NA, shape=NA, scale=NA, silent=FALSE)
}
\arguments{
   \item{data}{a vector of data to fit the distribution to.  Count data for the count models, continuous data for the continuous distributions.}
   \item{model}{the distribution to fit to the data.  Choices are:  "P", "ZIP", "G", "ZIG", "L", "ZIL", "W", "ZIW", "GP", "ZIGP", "LP", "ZILP", "WP", "ZIWP" (case insensitive).  No default.}
   \item{mean}{the starting value for mean.  Optional.}
   \item{variance}{the starting value for variance.  Optional.}
   \item{zi}{the starting value for zero-inflation.  Optional.}
   \item{shape}{the starting value for the shape parameter.  Optional.}
   \item{scale}{the starting value for the scale parameter.  Optional.}
   \item{silent}{should warning messages and progress indicators be supressed?  Default FALSE.}
}
\value{The values for each parameter at the maximum likelihood are output.  No standard error is given (Bayesian MCMC using bayescount() gives a much better analysis)
}
\seealso{
   \code{\link{bayescount}}
   \code{\link{bayescount.single}}
   \code{\link{likelihood}}
}

\author{Matthew Denwood \email{m.denwood@vet.gla.ac.uk} funded as part of the DEFRA VTRI project 0101.}

\examples{

# obtain values for mean and zero-inflation of a zero-inflated gamma Poisson model:

data <- rpois(100, rgamma(100, shape=1, scale=8))
data[1:15] <- 0

maximise.likelihood(data, "ZIGP")

}
\keyword{models}