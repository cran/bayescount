\name{bayescount.single}
\alias{bayescount.single}
\title{ANALYSE COUNT DATA USING JAGS}
\description{
   Apply a Bayesian (zero-inflated) (gamma / Weibull / lognormal / independant / simple) Poisson model to count data to return possible values for mean count, variance, shape paramater, scale parameter (overdispersion or 'k') and zero-infaltion where appropriate to the model selected.  Convergence is assessed for each dataset by calculating the Gelman-Rubin statistic for each parameter.  Optionally, the (log) likelihood for the model fit is also calculated.  Requires Just Another Gibbs Sampler (JAGS).  *THIS SOFTWARE IS INTENDED FOR EDUCATIONAL PURPOSES ONLY AND SHOULD NOT BE RELIED UPON FOR REAL WORLD APPLICATIONS*  When used with the GUI interface for R in Windows, this function spawns lots of terminal windows in which JAGS is run (if silent.jags is FALSE), which makes it difficult to track the progress of the simulations.  To avoid this, set silent.jags to TRUE, or run the function from the terminal version of R.
}
\usage{
bayescount.single(data = stop("Data must be specified"), model = "ZILP",
   burnin = 5000, updates = c(10000,100000,500000), jags = "jags", 
   alt.prior = FALSE, adjust.mean = FALSE, silent.jags = FALSE, 
   raw.output = FALSE, likelihood = FALSE)
}
\arguments{
   \item{data}{an existing R vector containing the data (integer vector).  No default.}
   \item{model}{model to use.  Choices are "GP" (gamma Poisson = negative binomial), "ZIGP" (zero-inflated gamma Poisson = zero-inflated negative binomial), "LP" (lognormal Poisson), "ZILP" (zero-inflated lognormal Poisson), "WP" (Wiebull Poisson), "ZIWP" (zero-inflated Weibull Poisson), "SP" (simple Poisson), "ZISP" (zero-inflated simple Poisson) or "IP" (independant Poisson).  The simple Poisson model forces each count to have the same mean, wheras the independant Poisson process allows each count to have an unrelated mean (therefore a zero-inflated version is not possible).  Default "ZILP".}
   \item{burnin}{the number of burnin iterations (not sampled) to use (numeric).  Default 5000 iterations.}
   \item{updates}{the number of sampling iterations to use (numeric).  Can be a single number or a vector of numbers (sorted from low to high by the function).  If a vector is supplied, the model run is extended to each successive number of iterations until convergence is achieved or the maximum number of iterations is reached.  Default c(10000, 100000, 500000).}
   \item{jags}{the system call or path for activating JAGS.  Default for Linux is 'jags', in Windows try 'C:/JAGS-0.90/jags.exe' etc.}
   \item{alt.prior}{should the model run the [ZI] [WP|GP|LP] models using the standard or the alternative prior distribution for variance? (logical)  Can also be a character value of a user-specified prior distribution.  Default FALSE.  Where information concerning overdispersion in the data is sparse, the choice of prior distribution will have an affect on the posterior distribution for ALL parameters.  It is recommended to run a simulation using both types of prior when working with small datasets, to make sure results are consistent.}
   \item{adjust.mean}{should the mean count parameter of the zero-inflated models be adjusted to reflect the mean of the whole population? (logical)  If FALSE the mean count of the zero-inflated models reflects the mean of the gamma or Poisson distribution only, if TRUE the mean includes extra zeros.  Used for comparing results between zero-inflated and non zero-inflated models.  Default FALSE.}
   \item{silent.jags}{should the JAGS output be suppressed? (logical)  If TRUE, no indication of the progress of individual models is supplied.  Also applies to the likelihood calculation.  Default FALSE.}
   \item{raw.output}{the function can return either a summary of the results (as with bayescount), or an MCMC object representing the estimates at each iteration for both chains (including the likelihood estimates where appropriate).  If TRUE, the latter is output.  (logical)  Default FALSE.}
   \item{likelihood}{should the (log) likelihood for the fit of the model to the dataset be calculated? (logical)  The likelihood for the [ZI] WP, LP and GP models are calculated using a likelihood function integrated over all possible values for lambda, which can take some time.  The likelihood is calculated using a thinned chain of 1000 values to reduce the time taken.  Default FALSE.}
}

\value{Either a vector containing an indication of the error/crash/convergence status, the number of sampled updates used, and a lower/upper 95% confidence interval estimate and median estimate for each relevant parameter (optionally including the likelihood), or an MCMC object representing the estimates at each iteration for both chains (optionally including the likelihood).
}

\seealso{
   \code{\link{bayescount}}
   \code{\link{run.model}}
   \code{\link{run.jags}}
   \code{\link{likelihood}}
}

\author{Matthew Denwood \email{m.denwood@vet.gla.ac.uk} funded as part of the DEFRA VTRI project 0101.}

\examples{

# use a zero-inflated lognormal Poisson model to analyse some count data using only 10000 updates and suppressing JAGS output:

\dontrun{
bayescount.single(data=c(0,5,3,7,0,4,3,8,0), model="ZILP", updates=10000, silent.jags=TRUE)
}

}
\keyword{models}
