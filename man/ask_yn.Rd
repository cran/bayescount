\name{ask_yn}
\alias{ask_yn}
\title{TRUE OR FALSE INPUT FROM USER}
\description{
   A simple function to detect input from the user, and keep prompting until a response matching 'yes' or 'no' is given.  Used by the function 'bayescount'.
}
\usage{
ask_yn(prompt)
}
\arguments{
   \item{prompt}{what text string should be used to prompt the user? (character string)}
}
\seealso{
   \code{\link{bayescount}}
}

\author{Matthew Denwood \email{m.denwood.1@research.gla.ac.uk}, funded as part of the DEFRA VTRI project 0101.}

\examples{
#  Ask the user if they want to proceed
\dontrun{
ask_yn("Do you want to start the program now?  ")
}
}
\keyword{methods}
