\name{fecrt.power.limits}
\alias{fecrt.power.limits}
\alias{FECRT.power.limits}
\title{FECRT Power Analysis Calculations (find tolerance)}
\description{
Finds the appropriate tolerance with which to consider the observed mean for a faecal egg count reduction test with the given combination of power and other parameters.  The power is calculated using Monte Carlo integration, and confidence intervals for the true power are produced.  Tolerance can be defined as either a lower limit only if upper.limit is defined, as an upper limit only if lower.limit is defined, or as both (equidistant) limits if neither are defined.  This function can be used for example to determine the upper limit for a true reduction of 80\% by setting lower.limit=0, or to determine the lower limit for a true reduction of 99\% by setting upper.limit=100.

*THIS SOFTWARE IS INTENDED FOR EDUCATIONAL PURPOSES ONLY AND SHOULD NOT BE RELIED UPON FOR REAL WORLD APPLICATIONS*
}
\usage{
fecrt.power.limits(meanepg=200, reduction = 95, g.faeces=3, 
   sensitivity=1/25, replicates=1, animals=10, pre.coeffvarrep=0.4, 
   pre.coeffvarind=0.3, pre.coeffvargroup=0.7, post.coeffvarrep=0.4, 
   post.coeffvarind=0.3, post.coeffvargroup=0.7, true.sample=FALSE, 
   lower.limit=NA, upper.limit=NA, iterations=100000, power = 0.95, 
   confidence = 0.99, feedback=FALSE)
}
\arguments{
   \item{meanepg}{the mean pre-treatment egg count of the group (in EPG).  If this is unknown then use a value likely to be encountered, or iterate over a distribution of values to obtain a distribution of values for tolerance.}
   \item{reduction}{the true mean egg count reduction (in \%).  If this is unknown then use a value likely to be encountered, or iterate over a distribution of values to obtain a distribution of values for tolerance.}
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
   \item{true.sample}{option to calculate the power for the population mean (if true.sample=FALSE) or the true sample mean (if true.sample=TRUE).  The difference is that conceptually the true mean of a small group of animals may not reflect the mean of the population from which they are derived.  If only one animal is considered, the true sample mean is the true mean of the individual whereas the population mean is the mean of a (theoretical) group of animals from which it is derived.  The power will always be greater (or identical) for the true sample mean.}
   \item{lower.limit}{(optional) lower limit for the tolerance (in \%).  If this is supplied, the required upper limit to maintain the desired power with the specified lower limit will be calculated.  If neither lower limit or upper limit are specified, symmetrical limits about the true mean reduction are found that satisfy the power condition.}
   \item{upper.limit}{(optional) upper limit for the tolerance (in \%).  If this is supplied, the required lower limit to maintain the desired power with the specified upper limit will be calculated.  If neither lower limit or upper limit are specified, symmetrical limits about the true mean reduction are found that satisfy the power condition.}
   \item{iterations}{the number of iterations to use for the Monte Carlo integration.  More iterations will take longer but provide a smaller confidence interval for the true power.}
   \item{power}{the desired power with which to interpret the given tolerances.  Default is 0.95, so that the observed mean FEC will lie within the calculated lower and upper limits 95\% of the time.}
   \item{confidence}{the degree of confidence required with which to report confidence limits for the true power when using Monte Carlo integration to report the power.  Only used for calculating the true power after determining the tolerance limits.}
   \item{feedback}{option to display a progress indicator for calculation of the values used for Monte Carlo integration.  Using feedback with some GUI versions of R may slow down the analysis considerably.}
}
\value{Returns a list containing the elements 'limits', which is the calculated lower and upper tolerance level which provides the required power, and 'power' which specifies the median estimate and confidence intervals for the true power.  The true power returned may not exactly match the required power input due to the integer nature of FECRT data.
}
\seealso{
   \code{\link{fecrt.power}}, \code{\link{fec.power.limits}}, \code{\link{fecrt}}, \code{\link{bayescount}}
}

\keyword{models}
