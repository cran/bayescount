\name{run.model}
\alias{run.model}
\title{Analyse Count Data Using Jags}
\description{
  Apply a Bayesian (zero-inflated) (gamma / Weibull / lognormal / independant / simple) Poisson model to count data to return possible values for mean count, variance, shape paramater, scale parameter (overdispersion or 'k') and zero-infaltion where appropriate to the model selected.  This function generates the model specifications and starting values, and is used by the higher level functions bayescount and bayescount.single, but can also be called directly.  Requires Just Another Gibbs Sampler (JAGS).  *THIS SOFTWARE IS INTENDED FOR EDUCATIONAL PURPOSES ONLY AND SHOULD NOT BE RELIED UPON FOR REAL WORLD APPLICATIONS*  The GUI interface for R in Windows may not continually refresh the output window, making it difficult to track the progress of the simulation (if silent.jags is FALSE).  To avoid this, you can run the function from the terminal version of R (located in the Program Files/R/bin/ folder).
}
\usage{
run.model(data=stop("No data supplied"), 
   model=stop("No model specified"), call.jags = TRUE, 
   alt.prior=FALSE, monitor.lambda=FALSE, ...)}
\arguments{
   \item{data}{an existing R vector containing the data (integer vector).  No default.}
   \item{model}{model to use.  Choices are "GP" (gamma Poisson = negative binomial), "ZIGP" (zero-inflated gamma Poisson = zero-inflated negative binomial), "LP" (lognormal Poisson), "ZILP" (zero-inflated lognormal Poisson), "WP" (Wiebull Poisson), "ZIWP" (zero-inflated Weibull Poisson), "SP" (simple Poisson), "ZISP" (zero-inflated simple Poisson) or "IP" (independant Poisson).  The simple Poisson model forces each count to have the same mean, wheras the independant Poisson process allows each count to have an unrelated mean (therefore a zero-inflated version is not possible).  Default "ZILP".}
   \item{call.jags}{should the function run the model using run.jags?  If not, the model specification, initial values and data in the correct format are returned so that the model can be run using run.jags with no modifications.  If TRUE, the model is run and the results are returned.}
   \item{alt.prior}{should the model run the [ZI] [WP|GP|LP] models using the standard or the alternative prior distribution for variance? (logical)  Can also be a character value of a user-specified prior distribution.  Default FALSE.  Where information concerning overdispersion in the data is sparse, the choice of prior distribution will have an affect on the posterior distribution for ALL parameters.  It is recommended to run a simulation using both types of prior when working with small datasets, to make sure results are consistent.}
   \item{monitor.lambda}{should the model or model specification monitor the mean of the Poisson process for each datapoint?  This is required to calculate the likelihood for the Independant Poisson model only, but may be useful for other purposes.  Default FALSE.}
   \item{...}{additional arguments to be passed directly to \code{\link[runjags]{autorun.jags}}.}
}

\value{Either a list of character strings representing [[1]] the model specification, [[2]] the data string, [[3]] the initial value string(s) and [[4]] the monitors required (if call.jags is FALSE) or an MCMC object returned from run.jags (if call.jags is TRUE).
}

\seealso{
   \code{\link{bayescount}}, \code{\link{bayescount.single}}
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
