\name{run.model}
\alias{run.model}
\title{ANALYSE COUNT DATA USING JAGS}
\description{
  Apply a Bayesian (zero-inflated) (gamma / Weibull / lognormal / independant / simple) Poisson model to count data to return possible values for mean count, variance, shape paramater, scale parameter (overdispersion or 'k') and zero-infaltion where appropriate to the model selected.  This function generates the model specifications and starting values, and is used by the higher level functions bayescount and bayescount.single, but can also be called directly.  Requires Just Another Gibbs Sampler (JAGS).  *THIS SOFTWARE IS INTENDED FOR EDUCATIONAL PURPOSES ONLY AND SHOULD NOT BE RELIED UPON FOR REAL WORLD APPLICATIONS*  When used with the GUI interface for R in Windows, this function spawns lots of terminal windows in which JAGS is run (if silent.jags is FALSE), which makes it difficult to track the progress of the simulations.  To avoid this, set silent.jags to TRUE, or run the function from the terminal version of R.
}
\usage{
run.model(data=stop("No data supplied"), 
   model=stop("No model specified"), burnin = 5000, 
   updates = 10000, call.jags = TRUE, alt.prior=FALSE, 
   jags="jags", silent.jags = FALSE, check.conv=TRUE, 
   monitor.lambda=FALSE)}
\arguments{
   \item{data}{an existing R vector containing the data (integer vector).  No default.}
   \item{model}{model to use.  Choices are "GP" (gamma Poisson = negative binomial), "ZIGP" (zero-inflated gamma Poisson = zero-inflated negative binomial), "LP" (lognormal Poisson), "ZILP" (zero-inflated lognormal Poisson), "WP" (Wiebull Poisson), "ZIWP" (zero-inflated Weibull Poisson), "SP" (simple Poisson), "ZISP" (zero-inflated simple Poisson) or "IP" (independant Poisson).  The simple Poisson model forces each count to have the same mean, wheras the independant Poisson process allows each count to have an unrelated mean (therefore a zero-inflated version is not possible).  Default "ZILP".}
   \item{burnin}{the number of burnin iterations (not sampled) to use (numeric).  Default 5000 iterations.}
   \item{updates}{the number of sampling iterations to use (numeric).  Can be a single number or a vector of numbers (sorted from low to high by the function).  If a vector is supplied, the model run is extended to each successive number of iterations until convergence is achieved or the maximum number of iterations is reached.  Default c(10000, 100000, 500000).}
   \item{call.jags}{should the function run the model using run.jags?  If not, the model specification, initial values and data in the correct format are returned so that the model can be run using run.jags with no modifications.  If TRUE, the model is run and the results are returned.}
   \item{alt.prior}{should the model run the [ZI] [WP|GP|LP] models using the standard or the alternative prior distribution for variance? (logical)  Can also be a character value of a user-specified prior distribution.  Default FALSE.  Where information concerning overdispersion in the data is sparse, the choice of prior distribution will have an affect on the posterior distribution for ALL parameters.  It is recommended to run a simulation using both types of prior when working with small datasets, to make sure results are consistent.}
   \item{jags}{the system call or path for activating JAGS.  Default for Linux is 'jags', in Windows try 'C:/JAGS-0.90/jags.exe' etc.}
   \item{silent.jags}{should the JAGS output be suppressed? (logical)  If TRUE, no indication of the progress of individual models is supplied.  Default FALSE.}
   \item{check.conv}{should the convergence be assessed after the model has been run using run.jags?  If TRUE, each monitored variable will be assessed for a potential scale reduction factor of less than 1.05, which indicates convergence.  2 or more chains are required to assess convergence.  Ignored if call.jags is FALSE.  Default TRUE.}
   \item{monitor.lambda}{should the model or model specification monitor the mean of the Poisson process for each datapoint?  This is required to calculate the likelihood for the Independant Poisson model only, but may be useful for other purposes.  Default FALSE.}
  }

\value{Either a list of character strings representing [[1]] the model specification, [[2]] the data string, [[3]] the initial value string(s) and [[4]] the monitors required (if call.jags is FALSE) or an MCMC object returned from run.jags (if call.jags is TRUE).
}

\seealso{
   \code{\link{bayescount}}
   \code{\link{bayescount.single}}
   \code{\link{run.jags}}
}

\author{Matthew Denwood \email{m.denwood@vet.gla.ac.uk} funded as part of the DEFRA VTRI project 0101.}

\examples{

#  Return the model specification and starting values for a lognormal Poisson, then run the model using run.jags:

\dontrun{
data <- rpois(100, rlnorm(3, 0.2))

strings <- run.model(model="LP", data=data, call.jags=FALSE)
modelstring <- strings[[1]]
datastring <- strings[[2]]
initial.values <- strings[[3]]
monitors <- strings[[4]]

run.jags(model=modelstring, inits=initial.values, check.conv=TRUE, data=datastring, monitor=monitors)

#  This is equivalent to using:

run.model(model="LP", data=data, call.jags=TRUE, check.conv=TRUE)
}

}
\keyword{models}
