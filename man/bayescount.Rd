\name{bayescount}
\alias{bayescount}
\alias{BAYESCOUNT}
\title{ANALYSE COUNT DATA USING JAGS}
\description{
   Apply a Bayesian (zero-inflated) (gamma) Poisson model to count data to return possible values for mean count, and overdispersion and/or zero-infaltion where appropriate to the model selected.  A text .csv file named *name*.*model*.csv with the results is optionally written to the working directory (before checking if the file already exists), and an object *name*.*model*.results is copied to the Global environment within R.  Where more than 1 model is used, a results file and object is created for each model.  Convergence is assessed for each dataset by calculating the Gelman-Rubin statistic for each parameter.  Requires Just Another Gibbs Sampler (JAGS).  *THIS SOFTWARE IS INTENDED FOR EDUCATIONAL PURPOSES ONLY AND SHOULD NOT BE RELIED UPON FOR REAL WORLD APPLICATIONS*  When used with the GUI interface for R in Windows, this function spawns lots of terminal windows in which JAGS is run, which makes it difficult to track the progress of the simulations.  To avoid this, I suggest you run the function from the terminal version of R instead.
}
\usage{
bayescount(name = NA, data = NA, setnames = NA, 
   div = 1, p.model = FALSE, zip.model = FALSE, 
   gp.model = FALSE, zigp.model = TRUE, burnin = 5000, 
   updates = c(10000,100000,500000), jags = "jags", 
   rownames = FALSE, remove.zeros = TRUE, remove.missing = TRUE, 
   test = TRUE, log.prior = TRUE, write.file = TRUE, 
   ml.equivalent = FALSE, adjust.mean = FALSE)
}
\arguments{
   \item{name}{a name for the analysis (character string).  Missing by default (function will require it to be input).}
   \item{data}{either a path to a comma delimited csv file, or an existing R object data frame containing the data.  Missing by default (function will require a path to the data to be input).}
   \item{setnames}{either a character vector of names for each dataset, or a logical value indicating if the data contains column labels in the first row.  If a character vector, the function will quit if the length does not match the length of the data.  Missing by default (function will require TRUE or FALSE to be input, a character vector cannot be input manually).}
   \item{div}{count division factor to allow egg count data in eggs per gram to be used raw (numeric).  Default 1 (no transformation to data).}
   \item{p.model}{use the Poisson model? (logical)  Default FALSE.}
   \item{zip.model}{use the zero-inflated Poisson model? (logical)  Default FALSE.}
   \item{gp.model}{use the gamma Poisson (negative binomial) model? (logical)  Default FALSE.}
   \item{zigp.model}{use the zero-inflated gamma Poisson (zero-inflated negative binomial) model? (logical)  Default TRUE.}
   \item{burnin}{the number of burnin iterations (not sampled) to use (numeric).  Default 5000 iterations.}
   \item{updates}{the number of sampling iterations to use (numeric).  Can be a single number or a vector of numbers from low to high.  If a vector is supplied, the function will run the model at each number of iterations until convergence is achieved or until it has tried each number of iterations.  Default c(10000, 100000, 500000).}
   \item{jags}{the system call or path for activating JAGS.  Default for Linux is 'jags', in Windows try 'C:/JAGS-0.90/jags.exe' etc.}
   \item{rownames}{does the data contain row labels in the first column?  (logical)  Default FALSE.}
   \item{remove.zeros}{remove any datasets where the total number of counts is 0, since it is not appropriate to use a count model to analyse these data (logical).  Default TRUE.}
   \item{remove.missing}{remove missing data before passing the data to JAGS?  (logical)  If FALSE, missing data are informed from the posteriors.  Default TRUE.}
   \item{test}{should the function briefly test the model with the first column of data before running the simulation?  (logical)  Affords extra 'user-proofing'.  If set to FALSE and valid values are supplied for 'name', 'data' and 'setnames', the function will not require input at any point (useful for automated data analysis).  Default TRUE.}
   \item{log.prior}{should the model run the (zero-inflated) gamma Poisson models using a log uniform prior distribution for overdispersion? (logical)  If FALSE, a uniform prior is used for overdispersion.  Default TRUE.  Where information concerning overdispersion in the data is sparse, the choice of prior distribution will have an affect on the posterior distribution for ALL parameters.  It is recommended to run a simulation using both types of prior when working with small datasets, to make sure results are consistent.}
   \item{write.file}{should the function write a text file to the current directory containing the results? (logical)  Default TRUE.  If FALSE, the text file is written during analysis and then deleted on completion of the dataset.}
   \item{ml.equivalent}{should the function also produce results that are equivalent to the log likelihood output of many maximum likelihood models (eg. the R package 'zicounts')?  (logical)  If TRUE, an additional 9 columns of data are produced.  Default FALSE.}
   \item{adjust.mean}{should the mean count parameter of the zero-inflated models be adjusted to reflect the mean of the whole population? (logical)  If FALSE the mean count of the zero-inflated models reflects the mean of the gamma or Poisson distribution only, if TRUE the mean includes extra zeros.  Used for comparing results between zero-inflated and non zero-inflated models.  Default FALSE.}
}

\value{No value is returned by this function.  Instead, a text .csv file named *name*.*model*.csv with the results is optionally written to the working directory (before checking if the file already exists), and an object *name*.*model*.results is copied to the Global environment within R.  Where more than 1 model is used, a results file and object is created for each model.  The 14 or 23 (depending on 'ml.equivalent') columns of data represent the following:
   \item{1}{dataset name or number}
   \item{2}{convergence successful (1 or 0)}
   \item{3}{model crashed (1 or 0)}
   \item{4}{error in model (1 or 0)}
   \item{5}{number of sampling iterations achieved}
   \item{6}{lower 95 percent confidence interval estimate for mean count}
   \item{7}{median estimate for mean count}
   \item{8}{upper 95 percent confidence interval estimate for mean count}
   \item{9}{lower 95 percent CI for zero-inflation}
   \item{10}{median for zero-inflation}
   \item{11}{upper 95 percent CI for zero-inflation}
   \item{12}{lower 95 percent CI for overdispersion}
   \item{13}{median for overdispersion}
   \item{14}{upper 95 percent CI for overdispersion}
   \item{15}{lower 95 percent confidence interval estimate for log(mean)}
   \item{16}{median estimate for log(mean)}
   \item{17}{upper 95 percent confidence interval estimate for log(mean)}
   \item{18}{lower 95 percent CI for logit(zero-inflation)}
   \item{19}{median for logit(zero-inflation)}
   \item{20}{upper 95 percent CI for logit(zero-inflation)}
   \item{21}{lower 95 percent CI for log(mean/overdispersion)}
   \item{22}{median for log(mean/overdispersion)}
   \item{23}{upper 95 percent CI for log(mean/overdispersion)}
}

\seealso{
   \code{\link{p.model}}
   \code{\link{gp.model}}
   \code{\link{zip.model}}
   \code{\link{zigp.model}}
}

\author{Matthew Denwood \email{m.denwood.1@research.gla.ac.uk}, funded as part of the DEFRA VTRI project 0101.}

\examples{

# run the function with all values as default, and 'name', 'data' and 'setnames' to be input by the user when prompted:
\dontrun{
bayescount()
}

# analyse data using a zero-inflated Poisson model in 5 text .csv files named 'mydata/data.*numer*.csv' with column labels, using a single tier of sampling updates, and a log uniform prior distribution for overdispersion:
\dontrun{
for (i in 1:5){
	bayescount(name=paste("Data ", i, sep=""), data=paste("mydata/data.", i, ".csv", sep=""), setnames=TRUE, updates = 10000, zip.model=TRUE, zigp.model=FALSE, test = FALSE, log.prior = TRUE)
}
}
}
\keyword{models}
