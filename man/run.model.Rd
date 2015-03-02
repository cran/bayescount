\name{count.model}
\alias{count.model}
\alias{fec.model}
\alias{FEC.model}
\alias{run.model}
\title{Analyse Count Data Using Jags}
\description{
  Apply a Bayesian (zero-inflated) (gamma / Weibull / lognormal / independant / simple) Poisson model to count data to return possible values for mean count, variance, shape paramater, scale parameter (overdispersion or 'k') and zero-infaltion where appropriate to the model selected.  This function generates the model specifications and starting values, and is used by the higher level functions bayescount and bayescount.single, but can also be called directly.
  
}
\usage{
count.model(data=stop("No data supplied"), 
   model=stop("No model specified"), call.jags = FALSE, 
   alt.prior=FALSE, monitor.lambda=FALSE, 
   monitor.deviance=FALSE, ...)}

\arguments{
   \item{data}{an existing R vector containing the data (integer vector).  No default.}

   \item{model}{model to use.  Choices are "GP" (gamma Poisson = negative binomial), "ZIGP" (zero-inflated gamma Poisson = zero-inflated negative binomial), "LP" (lognormal Poisson), "ZILP" (zero-inflated lognormal Poisson), "WP" (Wiebull Poisson), "ZIWP" (zero-inflated Weibull Poisson), "SP" (simple Poisson), "ZISP" (zero-inflated simple Poisson) or "IP" (independant Poisson).  The simple Poisson model forces each count to have the same mean, wheras the independant Poisson process allows each count to have an unrelated mean (therefore a zero-inflated version is not possible).  Default "ZILP".}

   \item{call.jags}{should the function run the model using run.jags?  If TRUE, the model is run and the results are returned, otherwise a compiled but not updated object of class \code{\link[runjags]{runjags-class}} is returned. NOTE: This argument is deprecated - in future only call.jags=FALSE will be permitted.}

   \item{alt.prior}{should the model run the [ZI] [WP|GP|LP] models using the standard or the alternative prior distribution for variance? (logical)  Can also be a character value of a user-specified prior distribution.  Default FALSE.  Where information concerning overdispersion in the data is sparse, the choice of prior distribution will have an affect on the posterior distribution for ALL parameters.  It is recommended to run a simulation using both types of prior when working with small datasets, to make sure results are consistent.}

   \item{monitor.lambda}{should the model or model specification monitor the mean of the Poisson process for each datapoint?  This is required to calculate the likelihood for the Independant Poisson model only, but may be useful for other purposes.  Default FALSE.}

   \item{monitor.deviance}{option to monitor the deviance of the model (using the built-in 'deviance' monitor).  Default FALSE.}

   \item{...}{additional arguments to be passed directly to \code{\link[runjags]{autorun.jags}}.}
}

\value{An object of class \code{\link[runjags]{runjags-class}}
}

\seealso{
   \code{\link{count.analysis}}
}

\examples{

#  Return the model specification and starting values for a 
# lognormal Poisson, then run the model using run.jags:

\dontrun{
data <- rpois(100, rlnorm(3, 0.2))

model <- run.model(model="LP", data=data, call.jags=FALSE)
library('runjags')
results <- extend.jags(model, burnin=5000, sample=10000)

}

}
\keyword{models}
