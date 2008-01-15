\name{bayescount}
\alias{bayescount}
\title{ANALYSE COUNT DATA USING JAGS}
\description{
   Apply a Bayesian (zero-inflated) (gamma / Weibull / lognormal / independant / simple) Poisson model to count data to return possible values for mean count, variance, shape paramater, scale parameter (overdispersion or 'k') and zero-infaltion where appropriate to the model selected.  A text .csv file named *name*.*model*.csv with the results is optionally written to the working directory (before checking if the file already exists), and an object *name*.*model*.results is copied to the Global environment within R.  Where more than 1 model is used, a results file and object is created for each model.  Convergence is assessed for each dataset by calculating the Gelman-Rubin statistic for each parameter.  Optionally, the (log) likelihood for the model fit is also calculated.  This function is a wrapper for bayescount.single(), allowing extra automation.  Requires Just Another Gibbs Sampler (JAGS).  *THIS SOFTWARE IS INTENDED FOR EDUCATIONAL PURPOSES ONLY AND SHOULD NOT BE RELIED UPON FOR REAL WORLD APPLICATIONS*  The GUI interface for R in Windows may not continually refresh the output window, making it difficult to track the progress of the simulation (if silent.jags is FALSE).  To avoid this, you can run the function from the terminal version of R (located in the Program Files/R/bin/ folder).
}
\usage{
bayescount(name = NA, data = NA, setnames = NA, 
   div = 1, model = c("ZILP"), burnin = 5000, 
   updates = c(10000,100000,500000), jags = findjags(), 
   rownames = FALSE, remove.zeros = TRUE, remove.missing = TRUE, 
   test = TRUE, alt.prior = FALSE, write.file = TRUE,
   adjust.mean = FALSE, crash.retry = 1,
   silent.jags = FALSE, likelihood = FALSE)
}
\arguments{
   \item{name}{a name for the analysis (character).  Missing by default (function will require it to be input).}
   \item{data}{either a path to a comma delimited csv file, or an existing R object data frame containing the data.  Missing by default (function will require a path to the data to be input).}
   \item{setnames}{either a character vector of names for each dataset, a logical value indicating if the data contains column labels in the first row, or an 'NA'.  If a character vector, the function will quit if the length does not match the length of the data.  If 'NA', then the function will assume the first row is dataset names if the values cannot be coerced into numeric values, or if the class of the dimnames[[2]] attribute of the matrix is 'character'.  Otherwise, generic names are used.  If the function mistakes data for labels (or vice versa), specify setnames as TRUE or FALSE to prevent it.  Default 'NA'.}
   \item{div}{count division factor to allow egg count data in eggs per gram to be used raw (numeric).  Default 1 (no transformation to data).}
   \item{model}{vector of models to use.  Choices are "GP" (gamma Poisson = negative binomial), "ZIGP" (zero-inflated gamma Poisson = zero-inflated negative binomial), "LP" (lognormal Poisson), "ZILP" (zero-inflated lognormal Poisson), "WP" (Wiebull Poisson), "ZIWP" (zero-inflated Weibull Poisson), "SP" (simple Poisson), "ZISP" (zero-inflated simple Poisson) or "IP" (independant Poisson), or "all" for all of these models (case insensitive).  The simple Poisson model forces each count to have the same mean, wheras the independant Poisson process allows each count to have an unrelated mean (therefore a zero-inflated version is not possible).  Default "ZILP".}
   \item{burnin}{the number of burnin iterations (not sampled) to use (numeric).  Default 5000 iterations.}
   \item{updates}{the number of sampling iterations to use (numeric).  Can be a single number or a vector of numbers (sorted from low to high by the function).  If a vector is supplied, the model run is extended to each successive number of iterations until convergence is achieved or the maximum number of iterations is reached.  Default c(10000, 100000, 500000).}
   \item{jags}{the system call or path for activating JAGS.  Default calls findjags() to attempt to locate JAGS on your system.}
   \item{rownames}{does the data contain row labels in the first column?  (logical)  Default FALSE.}
   \item{remove.zeros}{remove any datasets where the total number of counts is 0, since it is not appropriate to use a count model to analyse these data (logical).  Default TRUE.}
   \item{remove.missing}{remove missing data before passing the data to JAGS?  (logical)  If FALSE, missing data are informed from the posteriors.  Default TRUE.}
   \item{test}{should the function briefly test the model with the first column of data before running the simulation?  (logical)  Affords extra 'user-proofing'.  If set to FALSE and valid values are supplied for 'name', 'data' and 'setnames', the function will not require input at any point (useful for automated data analysis).  Default TRUE.}
   \item{alt.prior}{should the model run the [ZI] [WP|GP|LP] models using the standard or the alternative prior distribution for variance? (logical)  Can also be a character value of a user-specified prior distribution.  Default FALSE.  Where information concerning overdispersion in the data is sparse, the choice of prior distribution will have an affect on the posterior distribution for ALL parameters.  It is recommended to run a simulation using both types of prior when working with small datasets, to make sure results are consistent.}
   \item{write.file}{should the function write a text file to the current directory containing the results? (logical)  Default TRUE.  If FALSE, the text file is written during analysis and then deleted on completion of the dataset.}
   \item{adjust.mean}{should the mean count parameter of the zero-inflated models be adjusted to reflect the mean of the whole population? (logical)  If FALSE the mean count of the zero-inflated models reflects the mean of the gamma or Poisson distribution only, if TRUE the mean includes extra zeros.  Used for comparing results between zero-inflated and non zero-inflated models.  Default FALSE.}
   \item{crash.retry}{How many times should the model retry datasets that fail because of a crash or an error?  Datasets are restarted from the first iteration.  If 0, failed datasets are not retried.  (integer)  Default 1.}
   \item{silent.jags}{should the JAGS output be suppressed? (logical)  If TRUE, no indication of the progress of individual models is supplied.  Also applies to the likelihood calculation.  Default FALSE.}
   \item{likelihood}{should the (log) likelihood for the fit of each model to each dataset be calculated? (logical)  The likelihood for the [ZI] WP, LP and GP models are calculated using a likelihood function integrated over all possible values for lambda, which can take some time.  The likelihood is calculated using a thinned chain of 1000 values to reduce the time taken.  Default FALSE.}
}

\value{No value is returned by this function.  Instead, a text .csv file named *name*.*model*.csv with the results is optionally written to the working directory (before checking if the file already exists), and an object *name*.*model*.results is copied to the Global environment within R.  Where more than 1 model is used, a results file and object is created for each model.  The results files contain the dataset names, an indication of the error/crash/convergence status of each dataset, the number of sampled updates used, and a lower/upper 95% confidence interval estimate and median estimate for each relevant parameter (optionally including the likelihood).  The final column indicates the total time in seconds that the dataset took to be analysed (including calculating the likelihood).
}

\seealso{
   \code{\link{bayescount.single}}
   \code{\link{likelihood}}
}

\author{Matthew Denwood \email{m.denwood@vet.gla.ac.uk} funded as part of the DEFRA VTRI project 0101.}

\examples{

# run the function with all values as default, and 'name', 'data' and 'setnames' to be input by the user when prompted:
\dontrun{
bayescount()
}

# analyse data using zero-inflated gamma Poisson and zero-inflated lognormal Poisson models in 5 text .csv files named 'mydata/data.*numer*.csv' with column labels, using sampling updates increasing in 10000 increments, and calculating the likelihoods:
\dontrun{
for (i in 1:5){
	bayescount(name=paste("Data ", i, sep=""), data=paste("mydata/data.", i, ".csv", sep=""), model=c("ZIGP", "ZILP"), setnames=TRUE, updates = (1:10)*10000, test = FALSE, likelihood = TRUE)
}
}
}
\keyword{models}
