\name{likelihood}
\alias{likelihood}
\title{Calculate the (Log) Likelihood of Obtaining Data from a Distribution}
\description{
   Function to calculate the (log) likelihood of obtaining a set of data from one of the following distributions:  Poisson (P), gamma (G), lognormal (L), Weibull (W), gamma Poisson (GP), lognormal Poisson (LP), Weibull Poisson (WP), all with or without zero-inflation (ZI).  For mixture models, the likelihood is calculated for the data by integrating over each possible value of lambda for each data point, which may take some time for large datasets.  
   
   * NOTE: THIS FUNCTION IS DEPRECATED AND WILL BE REMOVED FROM BAYESCOUNT VERSION 1 *
   
}
\usage{
likelihood(model, data, mean=NA, variance=NA, 
   zi=NA, shape=NA, scale=NA, 
   iterations=min(1000, (length(as.matrix(mean)[,1])*length(shape))),
   log=TRUE, silent=FALSE, raw.output=FALSE)
}
\arguments{
   \item{model}{the distribution to fit to the data.  Choices are:  'IP', 'P', 'ZIP', 'G', 'ZIG', 'L', 'ZIL', 'W', 'ZIW', 'GP', 'ZIGP', 'LP', 'ZILP', 'WP', 'ZIWP' (case insensitive).  No default.}
   \item{data}{a vector of data to fit the distribution to.  Count data for the count models, continuous data for the continuous distributions.}
   \item{mean}{a vector of mean values over which to calculate the likelihood.  Must be a matrix with number of columns equal to the number of datapoints for the IP model only, since using the IP model each datapoint has an independant mean (this is effectively a shortcut for repeating the likelihood function for each datapoint).  Mean is required for the IP, (ZI)P, (ZI)L, and (ZI)L models.  Otherwise ignored (with a warning if silent=FALSE).}
   \item{variance}{a vector of values for variance over which to calculate the likelihood, with length equal to length of mean.  Required for the (ZI)L, and (ZI)L models.  Otherwise ignored (with a warning if silent=FALSE).}
   \item{zi}{a vector of values for zero-inflation, with length equal to either scale/shape or mean/variance as appropriate.  Required for zero-inflated models only; if supplied for other models the appropriate zero-inflated model will be used instead (with a warning if silent=FALSE).}
   \item{shape}{a vector of values for the shape parameter.  Must be of length equal to that of the scale parameter.  Required for the (ZI)W, (ZI)WP, (ZI)G, and (ZI)GP models.  Otherwise ignored (with a warning if silent=FALSE).}
   \item{scale}{a vector of values for the scale parameter.  Must be of length equal to that of the shape parameter.  Required for the (ZI)W, (ZI)WP, (ZI)G, and (ZI)GP models.  Otherwise ignored (with a warning if silent=FALSE).}
   \item{iterations}{the total number of iterative samples for which to calculate the likelihood.  If this is smaller than the number of iterations provided for mean/variance/shape/scale, then the chain will be thinned to the number of iterations required (this is useful for the mixture models where integrating the likelihood for each datapoint may take some time).  Default is 1000 iterations, or the length of the parameters supplied if this is shorter than 1000.}
   \item{log}{choose to output the log likelihood (TRUE) or the likelihood (FALSE).  Default TRUE.}
   \item{silent}{should warning messages and progress indicators be supressed?  Default FALSE.}
   \item{raw.output}{output a vector of length equal to the number of iterations representing the likelihood at each iteration if TRUE, or a summary of the results (median and 95 percent confidence interval) if FALSE.  Default FALSE.}
}
\value{The (log) likelihood is returned either as a vector of length equal to the number of iterations representing the likelihood at each iteration if raw.output==TRUE, or a summary of the results (median, maximum and 95\% highest posterior density interval (see \code{\link[coda]{HPDinterval}}) if raw.output==FALSE.  It is possible for the functions that use integrate to produce incalculable probabilities for some iterations, in which case the likelihood for these iterations are returned as missing data (and therefore the summary statistics for the whole simulation are also returned as missing data).  If only 1 iteration of values for mean/variance/shape/scale are provided, a single value for likelihood is output.
}
\seealso{
   \code{\link{count.analysis}}
}

\examples{

# calculate the likelihood of obtaining a set of count data from 
# a zero-inflated Poisson distribution with set mean and zero-inflation 
# values

data <- rpois(100, 10)
data[1:15] <- 0

likelihood('ZIP', data, mean=10, zi=15)

# now calculate the likelihood for the same data using an MCMC object 
# to provide the values for mean and zero-inflation

\dontrun{
values <- bayescount.single(data, model='ZISP', raw.output=TRUE)$mcmc
means <- c(values[,'mean'][[1]], values[,'mean'][[2]])
zis <- (1-c(values[,'prob'][[1]], values[,'prob'][[2]]))*100 
# The function outputs the prevalence of disease when raw.ouput is 
# TRUE, so zero-inflation must be calculated from this

likes <- likelihood('ZIP', data, mean=means, zi=zis, 
raw.output=TRUE)$likelihood
hist(likes, breaks='fd', col='red')
}

}
\keyword{models}