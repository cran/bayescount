\name{dump.format}
\alias{dump.format}
\title{PRODUCE A CHARACTER STRING IN THE R DUMP FORMAT TO BE USED WITH JAGS}
\description{
   Function used by bayescount.single() and run.model() to convert a numeric vector of data or initial values to a character string in the correct format to be read by JAGS as either data or initial values, using the jags.model function.
}
\usage{
dump.format(variable, value)
}
\arguments{
   \item{variable}{A list of character strings representing the name(s) of the variable(s) as used in the JAGS model}
   \item{value}{A list of numeric vectors, matrices or arrays representing the value(s) of the above variable(s), in the same order}
}
\value{A charcter string in the R dump format.
}

\seealso{
   \code{\link{bayescount.single}}
   \code{\link{run.jags}}
}

\author{Matthew Denwood \email{m.denwood@vet.gla.ac.uk} funded as part of the DEFRA VTRI project 0101.}

\examples{

initial.values.1 <- dump.format(list("mean", "sd", "lambda"), list("1", "0.1", c(4,2,7,0,6,9,1,4,12,1)))
initial.values.2 <- dump.format(list("mean", "sd", "lambda"), list("10", "10", c(4,2,7,0,6,9,1,4,12,1)))
data <- dump.format(list("N", "Count"), list("10", c(4,2,7,0,6,9,1,4,12,1)))

}
\keyword{methods}