\name{assess.variance}
\alias{assess.variance}
\title{ASSESS CERTIAINTY IN VARIANCE ESTIMATE FOR LP/GP/WP MODEL}
\description{
   Function used by bayescount.single() and bayescount() to determine if the variance estimate output by a model is sufficiently vague to warrant a warning.  Uncertain estimates are often heavily prior-driven, so re-running the model using a different prior and comparing the results is recommended
}
\usage{
assess.variance(model, alt.prior, l.95, u.95)
}
\arguments{
   \item{model}{Which model was used to generate the estimate?  (string:  'ZI WP/GP/LP')}
   \item{alt.prior}{Was the model run with the alternative prior?  (logical)}
   \item{l.95}{the lower 95 percent confidence interval for variance (numeric)}
   \item{u.95}{the upper 95 percent confidence interval for variance (numeric)}
}
\value{A logical result; 'TRUE' if the variance estimate is considered vague or 'FALSE' otherwise.
}

\seealso{
   \code{\link{bayescount}}
   \code{\link{bayescount.single}}
}

\author{Matthew Denwood \email{m.denwood@vet.gla.ac.uk} funded as part of the DEFRA VTRI project 0101.}

\keyword{methods}