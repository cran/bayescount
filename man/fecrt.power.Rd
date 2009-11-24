\name{fecrt.power}
\alias{fecrt.power}
\alias{FECRT.power}
\title{FECRT Power Analysis Calculations}
\description{
Finds the power for a faecal egg count reduction test with the given combination of parameters.  This represents the probability that the observed empirical mean reduction will lie between the lower.limit and upper.limit specified.  The power is calculated using Monte Carlo integration and confidence intervals for the true power are produced.  This function can be used to determine the probability that the observed empirical mean reduction will lie within a the range (for example) 90\% to 100\% if the true reduction is 80\% (the false negative rate if considering anthelmintic resistance), or the probability that the observed empirical mean reduction will lie within a the range (for example) 0\% to 90\% if the true reduction is 95\% (the false positive rate if considering anthelmintic resistance)

*THIS SOFTWARE IS INTENDED FOR EDUCATIONAL PURPOSES ONLY AND SHOULD NOT BE RELIED UPON FOR REAL WORLD APPLICATIONS*
}
\usage{
fecrt.power(meanepg=200, reduction = 80, g.faeces=3, 
   sensitivity=1/25, replicates=1, animals=10, pre.coeffvarrep=0.4, 
   pre.coeffvarind=0.3, pre.coeffvargroup=0.7, post.coeffvarrep=0.4, 
   post.coeffvarind=0.3, post.coeffvargroup=0.7, true.sample=FALSE, 
   lower.limit=0, upper.limit=90, maxiterations=1000000, 
   precision=2, confidence = 0.99, feedback=FALSE)
}
\arguments{
   \item{meanepg}{the mean pre-treatment egg count of the group (in EPG).  If this is unknown then use a value likely to be encountered, or iterate over a distribution of values to obtain a distribution of values for power.}
   \item{reduction}{the true mean egg count reduction (in \%).  If this is unknown then use a value likely to be encountered, or iterate over a distribution of values to obtain a distribution of values for power.}
   \item{g.faeces}{the number of grams of faeces used per sample (must be the same for all samples both pre and post treatment).}
   \item{sensitivity}{the minimum egg detection threshold or counting sensitivity of the technique used.  If using the McMasters technique, this is the number of McMasters chambers counted divided by 50.}
   \item{replicates}{the number of different samples (individually analysed) taken from each animal.  Must be the same for all animals, both pre and post treatment.  This would normally be 1, but increasing this provides an effective way of improving power at the cost of additional laboratory work.}
   \item{animals}{the number of animals in the group.  Can be 1, in which case the ability to predict the true mean reduction for the animal (if true.sample=TRUE) or the true mean reduction for a group from which the animal is taken (if true.sample=FALSE) can be calculated.  The same number of animals MUST be sample pre and post treatment.}
   \item{pre.coeffvarrep}{coefficient of variation between sub-samples taken from the same faecal sample, assuming 1g faeces is used per sub-sample (the effective value is automatically adjusted according to g.faeces in the function).  The default value for this is taken from an intensive study of a small group of horses, but may not be accurate in other groups of animals.  This value is applied only to the pre-treatment distributions.}
   \item{pre.coeffvarind}{coefficient of variation between samples taken from different faecal piles from the same animal over a period of days, not including that due to coeffvarrep.  The default value for this is taken from an intensive study of a small group of horses, but may not be accurate in other groups of animals.  This value is applied only to the pre-treatment distributions.}
   \item{pre.coeffvargroup}{coefficient of variation between animals.  This includes the variability between samples taken from different animals within a group, not including that due to within animal variability.  The default value for this is inferred from a combination of several equine datasets and an intensive study of a small group of horses, but may not be accurate in other groups of animals.  This value is applied only to the pre-treatment distributions.}
   \item{post.coeffvarrep}{coefficient of variation between sub-samples taken from the same faecal sample, assuming 1g faeces is used per sub-sample (the effective value is automatically adjusted according to g.faeces in the function).  The default value for this is taken from an intensive study of a small group of horses, but may not be accurate in other groups of animals.  This value is applied only to the post-treatment distributions.}
   \item{post.coeffvarind}{coefficient of variation between samples taken from different faecal piles from the same animal over a period of days, not including that due to coeffvarrep.  The default value for this is taken from an intensive study of a small group of horses, but may not be accurate in other groups of animals.  This value is applied only to the post-treatment distributions.}
   \item{post.coeffvargroup}{coefficient of variation between animals.  This includes the variability between samples taken from different animals within a group, not including that due to within animal variability.  The default value for this is inferred from a combination of several equine datasets and an intensive study of a small group of horses, but may not be accurate in other groups of animals.  This value is applied only to the post-treatment distributions.}
  \item{true.sample}{option to calculate the power for the population mean reduction (if true.sample=FALSE) or the true sample mean reduction (if true.sample=TRUE).  The difference is that conceptually the true mean of a small group of animals may not reflect the mean of the population from which they are derived.  If only one animal is considered, the true sample mean is the true mean of the individual whereas the population mean is the mean of a (theoretical) group of animals from which it is derived.  The power will always be greater (or identical) for the true sample mean.}
   \item{lower.limit}{lower limit for the tolerance (in \%).}
   \item{upper.limit}{lower limit for the tolerance (in \%).}
   \item{maxiterations}{the maximum number of iterations to use for the Monte Carlo integration, or the number to use if precision=NA.  If precision is defined, then only the number of iterations required to estimate the power to the given precision are performed (up to a maximum of maxiterations).}
   \item{precision}{the number of decimal places with which to calculate the power, unless maxiterations is reached first.  A larger precision will give a more precise estimate of the true power but will take longer to calculate.  Specifying this as NA performs the calculation on a fixed (=maxiterations) number of iterations.}
   \item{confidence}{the degree of confidence required with which to report confidence limits for the true power when using Monte Carlo integration to report the power.}
   \item{feedback}{option to display a progress indicator for calculation of the values used for Monte Carlo integration.  Using feedback with some GUI versions of R may slow down the analysis considerably.}
}
\value{Returns a list containing the elements 'roundedci' and 'ci', which specifies the median and confidence limits (as defined by 'confidence') for the true power both rounded by 'precison' and unrounded.  The additional variables 'within' 'without' and 'total' are also returned, and indicate the number of iterations for which the observed mean reduction fell outside and inside the specified limits and the total number of iterations.
}
\seealso{
   \code{\link{fecrt.power.limits}}, \code{\link{fec.power}}, \code{\link{fecrt}}, \code{\link{bayescount}}
}

\author{Matthew Denwood \email{m.denwood@vet.gla.ac.uk} funded as part of the DEFRA VTRI project 0101.}

\keyword{models}
