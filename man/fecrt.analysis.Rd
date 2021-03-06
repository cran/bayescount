\name{fecrt.analysis}
\alias{fecrt.analysis}
\alias{FECRT.analysis}
\alias{fecrt}
\alias{FECRT}
\title{Analyse FECRT Data Using MCMC to Give a Probability Distribution of Values for the Mean Egg Count Reduction}
\description{
   This function applies a Bayesian [zero-inflated] gamma Poisson (negative binomial) model to faecal egg count reduction test (FECRT) data to return possible values for the mean drug efficacy.  Pre-treatment data are assumed to arise from either a gamma-Poisson or zero-inflated gamma-Poisson distribution, with post-treatment data described by a separate gamma-Poisson or zero-inflated gamma-Poisson distribution.  The change in mean between these distributions is therefore the mean egg count reduction.  If paired.model=TRUE, a slightly different formulation is used whereby the observed count pre-treatment is assumed to follow a compound gamma-gamma-Poisson distribution with the variability within and between animals separated.  The post treatment mean for each animal is derived from the pre-treatment animal mean and FEC reduction.  This formulation allows data with non-random missing post-treatment counts to be analysed correctly, and also allows data with repeat counts from an individual to be analysed - providing a method of increasing the power of the method substantially.  Results are also obtained using non-parametric bootstrapping and the method advocated in the 1992 World Association for the Advancement of Veterinary Parasitology (W.A.A.V.P.) methods for the detection of anthelmintic resistance in nematodes of veterinary importance (unless the data contains repeat values or missing values).  Confidence intervals for the relevant statistics are printed to file if write.file = TRUE, and returned using a custom print method.  Lower level running of JAGS and assessing the simulation for convergence and required run length is done using  \code{\link[runjags]{autorun.jags}}.

** NOTE:  THIS FUNCTION IS CURRENTLY UNDER RE-DEVELOPMENT **

Please feel free to contact the package author if you would like to collaborate on a relevant data analysis

}
\usage{
fecrt.analysis(name = NA, pre.data = NA, post.data = NA,
   data = list(pre=pre.data, post=post.data), animal.names = FALSE, 
   efficacy=95, confidence=0.95, control.animals = FALSE, 
   divide.data = 1, record.chains = FALSE, write.file = FALSE, 
   bootstrap.iters=10000, plot.graph = TRUE, skip.mcmc = FALSE, ...)
}
\arguments{
   \item{name}{a name for the analysis (character).  Missing by default (function will require it to be input).}
   \item{pre.data}{the pre-treatment data.  Either a numeric vector, a matrix with repeated McMasters counts from the same sample (each row) in different columns, or an array with repeat samples from the same animal in dimension 1, repeated McMasters counts from the same sample in dimension 2 (can be of length 1 if only 1 count recorded), and different animals in dimension 3.  If an array, then the use.paired must be TRUE.  No default.  Ignored if a value is specified for data.}
   \item{post.data}{the post-treatment data.  Either a numeric vector, a matrix with repeated McMasters counts from the same sample (each row) in different columns, or an array with repeat samples from the same animal in dimension 1, repeated McMasters counts from the same sample in dimension 2 (can be of length 1 if only 1 count recorded), and different animals in dimension 3.  If an array, then use.paired must be TRUE.  No default.  Ignored if a value is specified for data.}
   \item{data}{either a path to a comma delimited csv file, or an existing R object containing the data.  Data can be specified in one of the following ways.  If a matrix, the first column is pre-treatment data, the second column is the post-treatment data, and the third column (if supplied) indicates control (1) or treatment (0) animal.  Alternatively, the first column is animal names id animal.names = TRUE, with columns 2 and 3 making up the pre- and post-treatment data and optionally column 4 the treatment of control status.  If the data is specified as a list, then the first element is taken as the pre-treatment data, and the second element is taken as the post-treatment data.  These elements of the list can be provided as for pre.data and post.data.  If the data is in a comma delimited file, it must be in the format specified for a single matrix.  Missing data, or unused elements of non-ragged arrays, may be represented using NA, which will be removed from the data before analysis.  This argument is taken from the specified values of pre.data and post.data by default.  If a value is specified for data then arguments specified for pre.data and post.data are ignored.}
   \item{animal.names}{either a character vector of names to be used for each animal specified in the data, TRUE or FALSE. If TRUE, then animal names are taken from the first column of the matrix supplied for data (a matrix must be supplied for data if animal.names is TRUE).  If FALSE, then animal names are generated automatically.  Default FALSE.}
   \item{efficacy}{the target \% efficacy of the drug used.  Used to calculate the probability of resistance with the MCMC and bootstrap methods.  Default 95\%.}
   \item{confidence}{the degree of confidence required with which to report the confidence limits for the true FEC reduction (a probability between 0 and 1).}
   \item{control.animals}{indication of which animals are to be used as controls.  Should be either a vector of TRUE/FALSE (or 1/0) corresponding to whether each animal is a control (TRUE) or treatment (FALSE), or simply FALSE (the default) in which case all animals are assumed to be treated.  Ignored if data is specified as a matrix (or csv file) with 3 columns, in which case the third column should reflect the treatment status.  Default FALSE.}
   \item{divide.data}{count division factor to allow egg count data in eggs per gram to be used raw (numeric).  Default 1 (no transformation to data).}
   \item{record.chains}{option to allow the MCMC chains to be recorded for future use.  If TRUE, the function returns the MCMC object as part of the return value (the MCMC object is not printed using the print method).  If write.file==TRUE, the results are also saved using a filename containing the name of the analysis.  Default FALSE.}
   \item{write.file}{option to write the results of the analysis to a text file using a filename containing the name of the analysis.  The contents of the text file are identical to the return value of the function.  Default FALSE.}
   \item{bootstrap.iters}{the number of bootstrap iterations to use for the bootstrap method.  Default 10000.}
   \item{plot.graph}{an option to plot the posterior true egg count reduction from the MCMC method graphically.  If write.file==TRUE then the graph is saved as a PDF, otherwise the graph is plotted on the active device.  Default FALSE.}
   \item{skip.mcmc}{option to omit the MCMC analysis, and return bootstrap and WAAVP method analysis results alone.  Default FALSE.}
   \item{...}{other options to be passed directly to \code{\link{fecrt.model}} and then \code{\link[runjags]{autoextend.jags}}.}
}
\value{Returns a list of the results obtained using each method, printed using a custom print method.
}
\seealso{
   \code{\link{fecrt.power}} and \code{\link{fecrt.power.limits}} for power analysis methods for the FECRT, and \code{\link{fecrt.model}} and \code{\link[runjags]{autorun.jags}} for relevent lower-level functions to which additional arguments can be passed
}

\references{M. J. Denwood, S. W. J. Reid, S. Love, M. K. Nielsen, L. Matthews, I. J. McKendrick, and G. T. Innocent.
Comparison of three alternative methods for analysis of equine Faecal Egg Count Reduction Test data. Prev. Vet. Med. (2009), doi:10.1016/j.prevetmed.2009.11.009}

\keyword{models}