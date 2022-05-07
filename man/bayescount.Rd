\name{bayescount}
\alias{bayescount}
\title{Analyse Count data using MCMC}
\description{
   Apply a Bayesian [zero-inflated] gamma / Weibull / lognormal / independant / simple Poisson model to count data to return possible values for mean count, coefficient of variation, and zero-inflation.  A text .csv file named [name].[model].csv with the results is optionally written to the working directory (before checking if the file already exists), and an object [name].[model].results is copied to the Global environment within R.  Where more than 1 model is used, a results file and object is created for each model.  Convergence is assessed for each dataset by calculating the Gelman-Rubin statistic for each parameter, see \code{\link[runjags]{autorun.jags}}.  Optionally, the log likelihood for the model fit is also calculated.  The time taken to complete each analysis (not including calculation of the likelihood) is also recorded.  This function is a wrapper for bayescount.single(), allowing extra automation.  The lower level functions in the runjags package are used for calling JAGS.

Note: The GUI interface for R in Windows may not continually refresh the output window, making it difficult to track the progress of the simulation (if silent.jags is FALSE).  To avoid this, you can run the function from the terminal version of R (located in the Program Files/R/bin/ folder).

   * NOTE: THIS FUNCTION IS DEPRECATED AND WILL BE REMOVED FROM BAYESCOUNT VERSION 1 - SEE \code{\link{fec.model}} FOR AN ALTERNATIVE *

}
\usage{
bayescount(name = NA, data = NA, setnames = NA,
   model = c("ZILP"), divide.data = 1,
   scale.mean = divide.data, remove.all.zeros = TRUE,
   test = TRUE, alt.prior = FALSE, write.file = TRUE,
   adjust.zi.mean = FALSE, likelihood = FALSE,
   record.chains = FALSE, ...)
}
\arguments{
   \item{name}{a name for the analysis (character).  Missing by default (function will require it to be input).}
   \item{data}{either a path to a comma delimited csv file, or an existing R object containing the data.  Data can be specified in one of the following ways.  If a numeric vector, or as a matrix or array with 1 column only, the data are taken to be a single dataset.  If a matrix, then each column of data is taken to be a dataset.  If an array, then each element of the 3rd dimension represents a dataset, and each column represents repeated McMasters counts from the same sample (each row).  A list containing the elements 'totals' representing the sum of the repeated McMasters counts, and 'repeats' representing the number of McMasters counts performed per sample, is also supported (these may be matrices, in which case each column is a seperate dataset).  If the data is in a comma delimited file, it must be in the format specified for a single matrix.  Missing data, or unused elements of non-ragged arrays, may be represented using NA, which will be removed from the data before analysis.  This argument is 'NA' by default (function will require a path to the data to be input).}
   \item{setnames}{either a character vector of names for each dataset, a logical value indicating if the data contains column labels in the first row, or an NA.  If a character vector, the function will quit if the length does not match the length of the data.  If TRUE, the dimnames[[2]] attribute of the matrix will be used if present, and failing that the first row of data is used.  If NA, then the function will use the dimnames[[2]] attribute of the matrix as setnames if it is not NULL, or otherwise generate them.  If FALSE, generic names are used.  If data is specified as an array or list, setnames must be provided as a character vector or left as FALSE (or NA).  Default NA.}
   \item{model}{vector of models to use.  Choices are "GP" (gamma Poisson = negative binomial), "ZIGP" (zero-inflated gamma Poisson = zero-inflated negative binomial), "LP" (lognormal Poisson), "ZILP" (zero-inflated lognormal Poisson), "WP" (Wiebull Poisson), "ZIWP" (zero-inflated Weibull Poisson), "SP" (simple Poisson), "ZISP" (zero-inflated simple Poisson) or "IP" (independant Poisson), or "all" for all of these models (case insensitive).  The simple Poisson model forces each count to have the same mean, wheras the independant Poisson process allows each count to have an unrelated mean (therefore a zero-inflated version is not possible).  Default "ZILP".}
   \item{divide.data}{count division factor to allow egg count data in eggs per gram to be used raw (numeric).  Default 1 (no transformation to data).}
   \item{scale.mean}{count multiplication factor to allow results to be reported as eggs per gram.  Default equal to divide.data, so that results are reported with the same mean as the input data.}
   \item{remove.all.zeros}{remove any datasets where the total number of counts is 0, since it is not appropriate to use a count model to analyse these data (logical).  Some models may crash repeatedly when trying to analyse a dataset with no observations.  Default TRUE.}
   \item{test}{should the function briefly test the model with the first column of data before running the simulation?  (logical)  Affords extra 'user-proofing'.  If set to FALSE and valid values are supplied for 'name', 'data' and 'setnames', the function will not require input at any point (useful for automated data analysis).  Default TRUE.}
   \item{alt.prior}{should the model run the [ZI] [WP|GP|LP] models using the standard or the alternative prior distribution for variability? (logical)  Can also be a character value of a user-specified prior distribution.  Default FALSE.  Where information concerning the coefficient of variation in the data is sparse, the choice of prior distribution will have an affect on the posterior distribution for ALL parameters.  It is recommended to run a simulation using both types of prior when working with small datasets, to make sure results are consistent.}
   \item{write.file}{should the function write a text file to the current directory containing the results? (logical)  Default TRUE.  If FALSE, the text file is written during analysis and then deleted on completion of the dataset.}
   \item{adjust.zi.mean}{should the mean count parameter of the zero-inflated models be adjusted to reflect the mean of the whole population? (logical)  If FALSE the mean count of the zero-inflated models reflects the mean of the gamma or Poisson distribution only, if TRUE the mean includes extra zeros.  Used for comparing results between zero-inflated and non zero-inflated models.  Default FALSE.}
   \item{likelihood}{should the (log) likelihood for the fit of each model to each dataset be calculated? (logical)  The likelihood for the [ZI] LP and GP models are calculated using a likelihood function integrated over all possible values for lambda, which can take some time.  The likelihood is calculated using a thinned chain of 1000 values to reduce the time taken.  Likelihood calculations for the [ZI] WP model are currently not available for data with repeated observations, or for any model with either missing data or different numbers of repeated observations per sample within a dataset.  Default FALSE.}
   \item{record.chains}{should the final mcmc chains generated for each simulation be saved to disk for future analysis?  If TRUE, an Rsave binary file is saved in the current working directory with a filename made up of the name of the analysis, model and dataset name.  The function \code{\link[runjags]{new_unique}} is used to generate the filename to ensure that no files are ever over-written.  Default FALSE.}
   \item{...}{further arguments to be passed directly to \code{\link[runjags]{autorun.jags}}, including the path to JAGS, number of burnin and sampled iterations, target PSRF, time limit to extend each simulation, etc.  Note that some of the default arguments to autorun.jags are changed by bayescount.single (max.time="1hr", interactive=FALSE, plots=FALSE), but these can be over-ridden.}
}

\value{No value is returned by this function.  Instead, a text .csv file named *name*.*model*.csv with the results is optionally written to the working directory (before checking if the file already exists), and an object *name*.*model*.results is copied to the Global environment within R (this will over-write any existing object of the same name).  Where more than 1 model is used, a results file and object is created for each model.

The results files consist of a table (or matrix) with the following results for each dataset:
\item{name}{the name of the dataset.  This is a dimnames attribute for matrices.}
\item{converged}{a logical flag indicating a PSRF of below the target indicating successful convergence (1), or the opposite (0).  Datasets where the model initially converged but the subsequent final chains had a PSRF of more than the target are classified as 'converged', but a warning message is printed in the function and the datasets will have an error code of 4.}
\item{error.code}{a numeric code indicating the exit status of the simulation.  0 indicates no errors encountered, 1 indicates a repeatedly crashing model, 2 indicates a failure to converge within the specified time limit, 3 indicates a dataset with no observations that was removed from analysis, 4 indicates a model that initially converged but the subsequent final chains had a PSRF of more than the target, and 5 indicates an unknown error in the function not related to the model (i.e. a bug in the runjags/bayescount code...).}
\item{samples}{the number of sampled iterations calculated by \code{\link[coda]{raftery.diag}} to reduce the MC error to the specified level.}
\item{samples.to.conv}{the number of sampled iterations discarded as a result of slow convergence (will be 0 if the pilot chain converged).}
\item{parameter estimates}{the lower 95\% credible interval, median estimate and upper 95\% credible interval for the mean, coefficient of variation, zero-inflation, log mean, log variance, shape parameter, and scale parameter as appropriate to the model.}
\item{multivariate.psrf}{the multivariate PSRF of the final mcmc chains.  If this couldn't be calculated, NA.}
\item{likelihood estimates}{the lower 95\% credible interval, median estimate, upper 95\% credible interval and maximum observed value (this is NOT equivalent to the maximum likelihood) of the likelihood for the model.  For the compound distributions these likelihoods are integrated for all possible values of the Poisson means (calculated using \code{\link{likelihood}}), and may well be different to the deviance calculated using JAGS.}
\item{time.taken}{the total time in seconds taken to run the simulation (not including calculating the likelihood and other summary statistics).}
}

\seealso{
   \code{\link[runjags]{autorun.jags}}, \code{\link{count.analysis}}, and \code{\link{fec.power}} for power analysis methods for FEC.
}

\examples{

# run the function with all values as default, and 'name' and 'data'
# (from a local .csv file) to be input by the user when prompted:
# bayescount()

# analyse local data (2 datasets with 20 animals each with 10 repeat
# samples) using a zero-inflated lognormal Poisson model:
\dontrun{
# Simulate some data:
data <- array(dim=c(20,10,2))
means1 <- rgamma(20, 10, 1)
means2 <- rgamma(20, 5, 1)
for(i in 1:20){
	data[i,,1] <- rpois(10, means1[i])
	data[i,,2] <- rpois(10, means2[i])
}
# Missing data is permissible but means the likelihood cannot be
# calculated - a warning will be printed:
data[sample(1:(20*10*2), 10)] <- NA
try(unlink("analysis.ZILP.csv"), silent=TRUE)
# Run the analysis:
bayescount(name="analysis", data=data, model = "ZILP",
setnames=c("Simulated group A", "Simulated group B"), likelihood=TRUE)
}

}
\keyword{models}
