\name{fecrt}
\alias{fecrt}
\alias{FECRT}
\title{Analyse FECRT Data Using Mcmc to Give a Probability Distribution of Values for the Mean Egg Count Reduction}
\description{
   This function applies a Bayesian [zero-inflated] gamma Poisson (negative binomial) model to faecal egg count reduction test (FECRT) data to return possible values for the mean drug efficacy.  Pre-treatment data are assumed to arise from either a gamma-Poisson or zero-inflated gamma-Poisson distribution, with post-treatment data described by a separate gamma-Poisson or zero-inflated gamma-Poisson distribution.  The change in mean between these distributions is therefore the mean egg count reduction.  MCMC output is also used to provide inference for the drug efficacy in each individual animal, the prevalence of egg shedding (one minus zero-inflation) and the change in variability between egg counts after treatment.  Results are also obtained using non-parametric bootstrapping and the method advocated in the 1992 World Association for the Advancement of Veterinary Parasitology (W.A.A.V.P.) methods for the detection of anthelmintic resistance in nematodes of veterinary importance.  Confidence intervals for the relevant statistics are printed to file if write.file = TRUE, and returned using a custom print method.

Lower level running of JAGS and assessing the simulation for convergence and required run length is done using  \code{\link[runjags]{autorun.jags}}.  Note: The GUI interface for R in Windows may not continually refresh the output window, making it difficult to track the progress of the simulation (if silent.jags is FALSE).  To avoid this, you can run the function from the terminal version of R (located in the Program Files/R/bin/ folder).  

*THIS SOFTWARE IS INTENDED FOR EDUCATIONAL PURPOSES ONLY AND SHOULD NOT BE RELIED UPON FOR REAL WORLD APPLICATIONS*
}
\usage{
fecrt(name = NA, pre.data = NA, post.data = NA, 
   data = list(pre=pre.data, post=post.data), 
   animal.names = FALSE, efficacy=95, restrict.efficacy = TRUE, 
   zero.inflation = FALSE, divide.data = 1, record.chains = FALSE, 
   write.file = FALSE, bootstrap.iters=10000, plot.graph = TRUE, 
   skip.mcmc = FALSE, individual.analysis=TRUE, ...)
}
\arguments{
   \item{name}{a name for the analysis (character).  Missing by default (function will require it to be input).}
   \item{pre.data}{the pre-treatment data.  Either a numeric vector, or a matrix with repeated McMasters counts from the same sample (each row) in different columns.  No default.  Ignored if a value is specified for data.}
   \item{post.data}{the post-treatment data.  Either a numeric vector, or a matrix with repeated McMasters counts from the same sample (each row) in different columns.  No default.  Ignored if a value is specified for data.}
   \item{data}{either a path to a comma delimited csv file, or an existing R object containing the data.  Data can be specified in one of the following ways.  If a list, then the first element is taken as the pre-treatment data, and the second element is taken as the post-treatment data.  A list containing the elements 'totals' representing the sum of the repeated McMasters counts, and 'repeats' representing the number of McMasters counts performed per sample, can be supplied for each element of this list.  If a matrix, the first column is pre-treatment data and the second column is the post-treatment data.  If an array, then the 3rd dimension should be of length 2, with the first element of the 3rd dimension representing pre-treatment data (with repeated McMasters counts from the same sample (each row) in different columns), and the second element of the 3rd dimension representing post-treatment data.en each element of the 3rd dimension represents a dataset.  If the data is in a comma delimited file, it must be in the format specified for a single matrix.  Missing data, or unused elements of non-ragged arrays, may be represented using NA, which will be removed from the data before analysis.  This argument is taken from the specified values of pre.data and post.data by default.  If a value is specified for data then arguments specified for pre.data and post.data are ignored.}
   \item{animal.names}{either a character vector of names to be used for each animal specified in the data, TRUE or FALSE. If TRUE, then animal names are taken from the first column of the matrix supplied for data (a matrix must be supplied for data if animal.names is TRUE).  If FALSE, then animal names are generated automatically.  Default FALSE.}
   \item{efficacy}{the target \% efficacy of the drug used.  Used to calculate the probability of resistance with the MCMC and bootstrap methods.  Default 95\%.}
   \item{restrict.efficacy}{option to restrict the estimate of efficacy from the MCMC method to values of greater than 0 (i.e. the mean faecal egg count cannot increase after treatment).  If TRUE, then the prior for change in mean egg count is uniform between 0 and 1, and if FALSE then the prior for change in mean egg count is uniform between 0 and 10.}
   \item{zero.inflation}{option to use a zero-inflated gamma Poisson in place of the gamma Poisson distribution.  If TRUE, zero inflated distributions are used pre- and post-treatment (with the prevalence fixed between pre- and post-treatment distributions).  To use zero-inflation, the individual.analysis must also be TRUE.  Default FALSE.}
   \item{divide.data}{count division factor to allow egg count data in eggs per gram to be used raw (numeric).  Default 1 (no transformation to data).}
   \item{record.chains}{option to allow the MCMC chains to be recorded for future use.  If TRUE, the function returns the MCMC object as part of the return value (the MCMC object is not printed using the print method).  If write.file==TRUE, the results are also saved using a filename containing the name of the analysis.  Default FALSE.}
   \item{write.file}{option to write the results of the analysis to a text file using a filename containing the name of the analysis.  The contents of the text file are identical to the return value of the function.  Default FALSE.}
   \item{bootstrap.iters}{the number of bootstrap iterations to use for the bootstrap method.  Default 10000.}
   \item{plot.graph}{an option to plot the posterior true egg count reduction from the MCMC method graphically.  If write.file==TRUE then the graph is saved as a PDF, otherwise the graph is plotted on the active device.  Default FALSE.}
   \item{skip.mcmc}{option to omit the MCMC analysis, and return bootstrap and WAAVP method analysis results.  Default FALSE.}
   \item{individual.analysis}{option to make inference on the efficacy on an individual animal basis.  The individual reduction in egg count is monitored, which will increase the amount of time required to run the simulation and may give problems with memory usage with large datasets.  Required to use the zero-inflated distributions.  Default TRUE.}
   \item{...}{other options to be passed directly to \code{\link[runjags]{autorun.jags}}.}
}
\value{Returns a list of the results obtained using each method, printed using a custom print method.
}
\seealso{
   \code{\link{bayescount}}, \code{\link{bayescount.single}}, \code{\link[runjags]{autorun.jags}}
}

\author{Matthew Denwood \email{m.denwood@vet.gla.ac.uk} funded as part of the DEFRA VTRI project 0101.}

\examples{
\dontrun{
fecrt("test", pre.data=rnbinom(10, mu=10, size=1), post.data=rnbinom(10, mu=1, size=1))
}
}
\keyword{models}