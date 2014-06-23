\name{fec.power}
\alias{fec.power}
\alias{FEC.power}
\title{FEC Power Analysis Calculations}
\description{
Finds the power for a faecal egg count study with the given combination of parameters.  This represents the probability that the observed empirical mean FEC will lie between the lower.limit and upper.limit specified.  The power is calculated using the negative binomial distribution when considering the true mean of a single individual, or using Monte Carlo integration for more than one animal.  Confidence intervals for the true power are produced for the latter.

*THIS SOFTWARE IS INTENDED FOR EDUCATIONAL PURPOSES ONLY AND SHOULD NOT BE RELIED UPON FOR REAL WORLD APPLICATIONS*
}
\usage{
fec.power(meanepg=200, g.faeces=3, sensitivity=1/25, 
   replicates=1, animals=10, coeffvarrep=0.4, coeffvarind=0.3, 
   coeffvargroup=0.7, true.sample=FALSE, accuracy=0.1, 
   lower.limit=meanepg*(1-accuracy), upper.limit=meanepg*(1+accuracy), 
   maxiterations=1000000, precision=2, confidence = 0.99, 
   feedback=FALSE, forcesim=FALSE)
}
\arguments{
   \item{meanepg}{the mean egg count of the group (in EPG).  If this is unknown then use a value likely to be encountered, or iterate over a distribution of values to obtain a distribution of values for power.}
   \item{g.faeces}{the number of grams of faeces used per sample (must be the same for all samples).}
   \item{sensitivity}{the minimum egg detection threshold or counting sensitivity of the technique used.  If using the McMasters technique, this is the number of McMasters chambers counted divided by 50.}
   \item{replicates}{the number of different samples (individually analysed) taken from each animal.  Must be the same for all animals.  This would normally be 1, but increasing this provides an effective way of improving power at the cost of additional laboratory work.}
   \item{animals}{the number of animals in the group.  Can be 1, in which case the ability to predict the true mean of the animal (if true.sample=TRUE) or the true mean of a group from which the animal is taken (if true.sample=FALSE) can be calculated.}
   \item{coeffvarrep}{coefficient of variation between sub-samples taken from the same faecal sample, assuming 1g faeces is used per sub-sample (the effective value is automatically adjusted according to g.faeces in the function).  The default value for this is taken from an intensive study of a small group of horses, but may not be accurate in other groups of animals.}
   \item{coeffvarind}{coefficient of variation between samples taken from different faecal piles from the same animal over a period of days, not including that due to coeffvarrep.  The default value for this is taken from an intensive study of a small group of horses, but may not be accurate in other groups of animals.}
   \item{coeffvargroup}{coefficient of variation between animals.  This includes the variability between samples taken from different animals within a group, not including that due to within animal variability.  The default value for this is inferred from a combination of several equine datasets and an intensive study of a small group of horses, but may not be accurate in other groups of animals.}
   \item{true.sample}{option to calculate the power for the population mean (if true.sample=FALSE) or the true sample mean (if true.sample=TRUE).  The difference is that conceptually the true mean of a small group of animals may not reflect the mean of the population from which they are derived.  If only one animal is considered, the true sample mean is the true mean of the individual whereas the population mean is the mean of a (theoretical) group of animals from which it is derived.  The power will always be greater (or identical) for the true sample mean.}
   \item{accuracy}{specify the lower and upper limits as the meanepg +/- this fraction of the meanepg.  Ignored if non-default values are specified for lower.limit or upper.limit.}
   \item{lower.limit}{lower limit for the tolerance.  Specifying lower.limit and upper.limit separately allows non-symmetrical tolerance around the mean with which to calculate the power.  The default is to use meanepg*(1-accuracy).}
   \item{upper.limit}{upper limit for the tolerance.  Specifying lower.limit and upper.limit separately allows non-symmetrical tolerance around the mean with which to calculate the power.  The default is to use meanepg*(1+accuracy).}
   \item{maxiterations}{the maximum number of iterations to use for the Monte Carlo integration, or the number to use if precision=NA.  If precision is defined, then only the number of iterations required to estimate the power to the given precision are performed (up to a maximum of maxiterations).}
   \item{precision}{the number of decimal places with which to calculate the power, unless maxiterations is reached first.  A larger precision will give a more precise estimate of the true power but will take longer to calculate.  Specifying this as NA performs the calculation on a fixed (=maxiterations) number of iterations.}
   \item{confidence}{the degree of confidence required with which to report confidence limits for the true power when using Monte Carlo integration to report the power.}
   \item{feedback}{option to display a progress indicator for calculation of the values used for Monte Carlo integration.  Using feedback with some GUI versions of R may slow down the analysis considerably.}
   \item{forcesim}{option to force the function to use the Monte Carlo method of approximating the power when the default would be to use a negative binomial approximation (when true.sample=FALSE and coeffvargroup is very small, or when true.sample=TRUE and animals=1).  Results using the two methods should be very similar, however there may be differences due to the differenct optimisation mechanisms when calculating both lower and upper limits.  In addition, Monte Carlo integration will take longer and gives a confidence interval for the true mean rather than the absolute value.}
}
\value{Returns a list containing the elements 'roundedci' and 'ci', which specifies the median and confidence limits (as defined by 'confidence') for the true power both rounded by 'precison' and unrounded.  For analyses using the Monte Carlo integration method, 'within' 'without' and 'total' are also returned, and indicate the number of iterations for which the observed mean fell outside and inside the specified limits and the total number of iterations.
}
\seealso{
   \code{\link{fec.power.limits}}, \code{\link{fecrt.power}}, \code{\link{fecrt}}, \code{\link{bayescount}}
}

\keyword{models}
