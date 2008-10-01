\name{run.jags}
\alias{run.jags}
\alias{run.JAGS}
\title{RUN A USER SPECIFIED MODEL IN JAGS FROM WITHIN R}
\description{
   Allows any user specified model to be run in JAGS, with the MCMC output returned as a list of MCMC objects for each chain and the end state of each chain.  Data and initial values must be supplied in the R dump format (see dump.format() for an easy way to do this).  A character vector of variables to monitor must also be supplied.  Requires Just Another Gibbs Sampler (JAGS).  When used with the GUI interface for R in Windows, this function spawns lots of terminal windows in which JAGS is run (if silent.jags is FALSE), which makes it difficult to track the progress of the simulations.  To avoid this, set silent.jags to TRUE, or run the function from the terminal version of R.
}
\usage{
run.jags(data=stop("No data supplied"), model=stop("No model supplied"), 
   inits = stop("No inital values supplied"), 
   monitor = stop("No monitored variables supplied"), 
   burnin = 5000, updates = 10000, jags="jags", 
   silent.jags = FALSE, check.conv=TRUE)}
\arguments{
   \item{data}{a character string in the R dump format containing the data.  No default.}
   \item{model}{a character string of the model in the JAGS language.  No default.}
   \item{inits}{a character vector with length equal to the number of chains the model will be run using.  Each element of the vector must be a character string in the R dump format representing the initial values for that chain.  No default.}
   \item{monitor}{a character vector of the names of variables to monitor.  No default.}
   \item{burnin}{the number of burnin iterations (not sampled) to use (numeric).  Default 5000 iterations.}
   \item{updates}{the number of sampling iterations to use (numeric).  Default 10000 iterations.}
   \item{jags}{the system call or path for activating JAGS.  Default for Linux is 'jags', in Windows try 'C:/JAGS-0.90/jags.exe' etc.}
   \item{silent.jags}{should the JAGS output be suppressed? (logical)  If TRUE, no indication of the progress of individual models is supplied.  Default FALSE.}
   \item{check.conv}{should the convergence be assessed after the model has completed?  If TRUE, each monitored variable will be assessed for a potential scale reduction factor of less than 1.05, which indicates convergence.  2 or more chains are required to assess convergence.  Default TRUE.}
}

\value{The results of the simulation are returned as a list of MCMC objects, followed by a description of the model state in the R dump format at the last iteration.  An MCMC object and model state description are returned for each chain.  Model state descriptions can be used as an initial value string for a continuing simulation.
}

\seealso{
   \code{\link{bayescount}}
   \code{\link{bayescount.single}}
   \code{\link{dump.format}}
   \code{\link{run.model}}
   \code{\link{testjags}}
}

\author{Matthew Denwood \email{m.denwood@vet.gla.ac.uk} funded as part of the DEFRA VTRI project 0101.}

\examples{
# run a model to calculate the intercept and slope of the expression y = m x + c, assuming normal observation errors for y:

\dontrun{
# Simulate the data
x <- 1:100
y <- rnorm(length(true.x), 2*true.x + 10, 1)

# Model in the JAGS format
model <- "model {
for(i in 1 : n){
y[i] ~ dnorm(true.y[i], precision);
true.y[i] <- (m * x[i]) + c;
}
m ~ dunif(-1000,1000);
c ~ dunif(-1000,1000);
precision ~ dexp(1);
}"

# Use dump.format to convert the data and initial values files into the R dump format
data <- dump.format(c("x", "y", "n"), list(x, y, length(x)))
inits1 <- dump.format(c("m", "c", "precision"), list(1, 1, 1))
inits2 <- dump.format(c("m", "c", "precision"), list(0.1, 10, 1))

# Run the model
results <- run.jags(data=data, model=model, inits=c(inits1,
inits2), monitor=c("m", "c", "precision"))

# Analyse the results
m <- summary(c(results[[1]][,1], results[[2]][,1]))
c <- summary(c(results[[1]][,2], results[[2]][,2]))
p <- summary(c(results[[1]][,3], results[[2]][,3]))
}

}
\keyword{models}