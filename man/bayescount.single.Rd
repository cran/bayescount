\name{bayescount.single}
\alias{bayescount.single}
\title{Analyse Count data using MCMC}
\description{
   Apply a Bayesian [zero-inflated] gamma / Weibull / lognormal / independant / simple Poisson model to count data to return possible values for mean count, coefficient of variation, and zero-inflation, as either summary statistics or mcmc objects.  Convergence is assessed for each dataset by calculating the Gelman-Rubin statistic for each parameter, see \code{\link[runjags]{autorun.jags}}.  Optionally, the log likelihood for the model fit is also calculated.  The time taken to complete each analysis (not including calculation of the likelihood) is also recorded.  The lower level functions in the runjags package are used for calling JAGS.  

Note: The GUI interface for R in Windows may not continually refresh the output window, making it difficult to track the progress of the simulation (if silent.jags is FALSE).  To avoid this, you can run the function from the terminal version of R (located in the Program Files/R/bin/ folder).  

*THIS SOFTWARE IS INTENDED FOR EDUCATIONAL PURPOSES ONLY AND SHOULD NOT BE RELIED UPON FOR REAL WORLD APPLICATIONS*
}
\usage{

bayescount.single(data = stop("Data must be specified"), model="ZILP", 
   alt.prior = FALSE, adjust.zi.mean = FALSE, raw.output = FALSE, 
   likelihood=FALSE, ...)
}
\arguments{
   \item{data}{an existing R object containing the data.  Data can either be specified as a numeric vector of single counts for each sample, or as a matrix of repeated counts (columns) for each sample (rows).  Repeated counts are modelled as part of the same Poisson process.  Data can also be specified as a (ragged) list of repeated counts, with each element of the list representing a seperate sample.  Finally, data can be specified as a list containing the elements 'totals' representing the sum of the repeated McMasters counts, and 'repeats' representing the number of McMasters counts performed per sample.  Missing data (or unused elements of non-ragged arrays) may be represented using NA, which will be removed from the data before analysis.  The likelihood calculation is not available for ragged arrays and missing data, and will print a warning.  No default.}
   \item{model}{the model to use.  Choices are "GP" (gamma Poisson = negative binomial), "ZIGP" (zero-inflated gamma Poisson = zero-inflated negative binomial), "LP" (lognormal Poisson), "ZILP" (zero-inflated lognormal Poisson), "WP" (Wiebull Poisson), "ZIWP" (zero-inflated Weibull Poisson), "SP" (simple Poisson), "ZISP" (zero-inflated simple Poisson) or "IP" (independant Poisson).  Case insensitive.  The simple Poisson model forces each count to have the same mean, wheras the independant Poisson process allows each count to have an unrelated mean (therefore a zero-inflated version is not possible).  Default "ZILP".}
   \item{alt.prior}{should the model run the [ZI] [WP|GP|LP] models using the standard or the alternative prior distribution for variance? (logical)  Can also be a character value of a user-specified prior distribution.  Default FALSE.  Where information concerning overdispersion in the data is sparse, the choice of prior distribution will have an affect on the posterior distribution for ALL parameters.  It is recommended to run a simulation using both types of prior when working with small datasets, to make sure results are consistent.}
   \item{adjust.zi.mean}{should the mean count parameter of the zero-inflated models be adjusted to reflect the mean of the whole population? (logical)  If FALSE the mean count of the zero-inflated models reflects the mean of the gamma or Poisson distribution only, if TRUE the mean includes extra zeros.  Used for comparing results between zero-inflated and non zero-inflated models.  Default FALSE.}
   \item{raw.output}{the function can return either a summary of the results (as with bayescount), or an MCMC object representing the estimates at each iteration for both chains (including the likelihood estimates where appropriate).  If TRUE, the latter is output.  (logical)  Default FALSE.}
   \item{likelihood}{should the (log) likelihood for the fit of the model to the dataset be calculated? (logical)  The likelihood for the [ZI] WP, LP and GP models are calculated using a likelihood function integrated over all possible values for lambda, which can take some time.  The likelihood is calculated using a thinned chain of 1000 values to reduce the time taken.  Default FALSE.}
   \item{...}{additional arguments to be passed directly to \code{\link[runjags]{autorun.jags}}.}
}

\value{Either a vector similar to that obtained from \code{\link{bayescount}} containing an indication of the error/crash/convergence status, the number of sampled updates used, and a lower/upper 95\% highest posterior density interval (see \code{\link[coda]{HPDinterval}}), and median estimate for each relevant parameter (optionally including the likelihood), or an MCMC object representing the estimates at each iteration for both chains (optionally including the likelihood).
}

\seealso{
   \code{\link{bayescount}}, \code{\link{run.model}}, \code{\link{likelihood}}
}

\examples{

# use a zero-inflated lognormal Poisson model to analyse some count 
# data, and suppressing JAGS output:

\dontrun{
bayescount.single(data=c(0,5,3,7,0,4,3,8,0), model="ZILP", 
silent.jags=TRUE)
}

}
\keyword{models}
