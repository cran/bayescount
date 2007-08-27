\name{new_unique}
\alias{new_unique}
\title{CREATE A UNIQUE FILENAME}
\description{
   Search the current working directory for a file or directory matching the input name, and if it exists create a new name by appending a counter to the input name.  Alternatively, the function can ask the user if the existing file should be overwritten.  Used by the bayescount function.
}
\usage{
new_unique(name, suffix = "", ask = FALSE, prompt = "A file or directory with this name already exists.  Overwrite?  ")
}
\arguments{
   \item{name}{the filename to be used (character string).}
   \item{suffix}{the file extension (including '.') to use (character string).  Default none.}
   \item{ask}{if a file exists with the input name, should the function ask to overwrite the file? (logical)  If FALSE, a new filename is used instead and no files will be over-written.  Default FALSE.}
   \item{prompt}{if ask==TRUE, what text string should be used to prompt the user? (charcter string)  A generic default is supplied.}
}
\seealso{
   \code{\link{bayescount}}
}

\author{Matthew Denwood \email{m.denwood.1@research.gla.ac.uk}, funded as part of the DEFRA VTRI project 0101.}

\examples{
#  Create a file name with a .R extension.
new_unique("new_file", ".R", ask=FALSE)
}
\keyword{methods}
