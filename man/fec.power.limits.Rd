\name{fec.power.limits}
\alias{fec.power.limits}
\alias{FEC.power.limits}
\title{FEC Power Analysis Calculations (find tolerance)}
\description{
Finds the appropriate tolerance with which to consider the observed mean for a faecal egg count study with the given combination of power and other parameters.  The power is calculated using the negative binomial distribution when considering the true mean of a single individual, or using Monte Carlo integration for more than one animal.  Confidence intervals for the true power are produced for the latter.  Tolerance can be defined as either a lower limit only if upper.limit is defined, as an upper limit only if lower.limit is defined, or as both (equidistant) limits if neither are defined.

*THIS SOFTWARE IS INTENDED FOR EDUCATIONAL PURPOSES ONLY AND SHOULD NOT BE RELIED UPON FOR REAL WORLD APPLICATIONS*
}
\usage{
fec.power.limits(meanepg=200, g.faeces=3, sensitivity=1/25, 
   replicates=1, animals=10, coeffvarrep=0.4, coeffvarind=0.3, 
   coeffvargroup=0.7, true.sample=FALSE, lower.limit=NA, 
   upper.limit=NA, iterations=100000, power = 0.95, 
   confidence = 0.99, feedback=FALSE, forcesim=FALSE)
}
\arguments{
   \item{meanepg}{the mean egg count of the group (in EPG).  If this is unknown then use a value likely to be encountered, or iterate over a distribution of values to obtain a distribution of values for tolerance.}
   \item{g.faeces}{the number of grams of faeces used per sample (must be the same for all samples).}
   \item{sensitivity}{the minimum egg detection threshold or counting sensitivity of the technique used.  If using the McMasters technique, this is the number of McMasters chambers counted divided by 50.}
   \item{replicates}{the number of different samples (individually analysed) taken from each animal.  Must be the same for all animals.  This would normally be 1, but increasing this provides an effective way of improving power at the cost of additional laboratory work.}
   \item{animals}{the number of animals in the group.  Can be 1, in which case the ability to predict the true mean of the animal (if true.sample=TRUE) or the true mean of a group from which the animal is taken (if true.sample=FALSE) can be calculated.}
   \item{coeffvarrep}{coefficient of variation between sub-samples taken from the same faecal sample, assuming 1g faeces is used per sub-sample (the effective value is automatically adjusted according to g.faeces in the function).  The default value for this is taken from an intensive study of a small group of horses, but may not be accurate in other groups of animals.}
   \item{coeffvarind}{coefficient of variation between samples taken from different faecal piles from the same animal over a period of days, not including that due to coeffvarrep.  The default value for this is taken from an intensive study of a small group of horses, but may not be accurate in other groups of animals.}
   \item{coeffvargroup}{coefficient of variation between animals.  This includes the variability between samples taken from different animals within a group, not including that due to within animal variability.  The default value for this is inferred from a combination of several equine datasets and an intensive study of a small group of horses, but may not be accurate in other groups of animals.}
   \item{true.sample}{option to calculate the power for the population mean (if true.sample=FALSE) or the true sample mean (if true.sample=TRUE).  The difference is that conceptually the true mean of a small group of animals may not reflect the mean of the population from which they are derived.  If only one animal is considered, the true sample mean is the true mean of the individual whereas the population mean is the mean of a (theoretical) group of animals from which it is derived.  The power will always be greater (or identical) for the true sample mean.}
   \item{lower.limit}{(optional) lower limit for the tolerance.  If this is supplied, the required upper limit to maintain the desired power with the specified lower limit will be calculated.  If neither lower limit or upper limit are specified, symmetrical limits about the mean are found that satisfy the power condition.}
   \item{upper.limit}{(optional) upper limit for the tolerance.  If this is supplied, the required lower limit to maintain the desired power with the specified upper limit will be calculated.  If neither lower limit or upper limit are specified, symmetrical limits about the mean are found that satisfy the power condition.}
   \item{iterations}{the number of iterations to use for the Monte Carlo integration.  More iterations will take longer but provide a smaller confidence interval for the true power.}
   \item{power}{the desired power with which to interpret the given tolerances.  Default is 0.95, so that the observed mean FEC will lie within the calculated lower and upper limits 95\% of the time.}
   \item{confidence}{the degree of confidence required with which to report confidence limits for the true power when using Monte Carlo integration to report the power.  Only used for calculating the true power after determining the tolerance limits.}
   \item{feedback}{option to display a progress indicator for calculation of the values used for Monte Carlo integration.  Using feedback with some GUI versions of R may slow down the analysis considerably.}
   \item{forcesim}{option to force the function to use the Monte Carlo method of approximating the power when the default would be to use a negative binomial approximation (when true.sample=FALSE and coeffvargroup is very small, or when true.sample=TRUE and animals=1).  Results using the two methods should be very similar, however there may be differences due to the differenct optimisation mechanisms when calculating both lower and upper limits.  In addition, Monte Carlo integration will take longer and gives a confidence interval for the true mean rather than the absolute value.}
}
\value{Returns a list containing the elements 'limits', which is the calculated lower and upper tolerance level which provides the required power, and 'power' which specifies the median estimate and confidence intervals for the true power when using Monte Carlo integration, or the absolute value (replicated 3 times for consistency) if using the negative binomial approximation method.  The true power returned may not exactly match the required power input due to the integer nature of FEC data.
}
\seealso{
   \code{\link{fec.power}}, \code{\link{fecrt.power.limits}}, \code{\link{fecrt}}, \code{\link{bayescount}}
}

\author{Matthew Denwood \email{m.denwood@vet.gla.ac.uk} funded as part of the DEFRA VTRI project 0101.}

\keyword{models}
