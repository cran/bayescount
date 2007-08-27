\name{zip.model}
\alias{zip.model}
\alias{ZIP.model}
\title{ANALYSE COUNT DATA USING JAGS}
\description{
   Apply a Bayesian zero-inflated Poisson model to count data to return probable values for mean count and zero-infaltion.  Usually called internally from the bayescount function, but can be called directly by the user (although there is significantly less user-proofing for using this function directly than for bayescount).  Convergence is assessed by calculating the Gelman-Rubin statistic for each parameter.  Requires Just Another Gibbs Sampler (JAGS).  *THIS SOFTWARE IS INTENDED FOR EDUCATIONAL PURPOSES ONLY AND SHOULD NOT BE RELIED UPON FOR REAL WORLD APPLICATIONS*
}
\usage{
zip.model(data, name, jags = "jags", burnin = 5000, 
   updates = 10000, ml.equivalent = FALSE, adjust.mean = FALSE)
}
\arguments{
   \item{data}{an existing integer vector containing the data.  No default.}
   \item{name}{a name for the analysis (character string).  No default.}
   \item{jags}{the system call or path for activating JAGS.  Default for Linux is 'jags', in Windows try 'C:/JAGS-0.90/jags.exe' etc.}
   \item{burnin}{the number of burnin iterations (not sampled) to use (numeric).  Default 5000 iterations.}
   \item{updates}{the number of sampling iterations to use (numeric).  Default 10000 iterations.}
   \item{ml.equivalent}{should the function also produce results that are equivalent to the log likelihood output of many maximum likelihood models (eg. the R package 'zicounts')?  (logical)  If TRUE, an additional 9 columns of data are produced.  Default FALSE.}
   \item{adjust.mean}{should the mean count parameter be adjusted to reflect the mean of the whole population? (logical)  If FALSE the mean count of the reflects the mean of the Poisson distribution only, if TRUE the mean includes extra zeros.  Used for comparing results between zero-inflated and non zero-inflated models.  Default FALSE.}
}
\value{The results of the simulation are returned as a vector of length 14 or 23 (depending on 'ml.equivalent'), containing the following:
   \item{1}{name}
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
   \code{\link{bayescount}}
   \code{\link{p.model}}
   \code{\link{gp.model}}
   \code{\link{zigp.model}}
}

\author{Matthew Denwood \email{m.denwood.1@research.gla.ac.uk}, funded as part of the DEFRA VTRI project 0101.}

\examples{
# run an analysis named 'test' on a dataframe object 'mydata', with all model values as default, and save the results to 'myresults':

mydata <- c(5,3,7,9,5,6,2,1,6,5)
myresults <- zip.model(name="test", data=mydata)

}
\keyword{models}