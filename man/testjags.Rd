\name{testjags}
\alias{testjags}
\alias{testJAGS}
\title{ANALYSE THE SYSTEM TO CHECK THAT JAGS IS INSTALLED}
\description{
   Test the users system to determine the operating system, version of R installed, and version of JAGS installed.  Some information is collected from other functions such as .platform and Sys.info.  Used by the bayescount, bayescount.single and call.jags functions.  
}
\usage{
testjags(jags="jags", silent=FALSE)}
\arguments{
   \item{jags}{the system call or path for activating JAGS.  Default for Linux is 'jags', in Windows try 'C:/JAGS-0.90/jags.exe' etc.}
   \item{silent}{should on-screen feedback be suppressed?  Default FALSE.}
}
\value{The following self-explanatory information is returned:  $os, $JAGS.available (logical), $popen.support (internalisation of terminal output), $JAGS.version, $R.version, $R.GUI, $R.package.type, $username.
}
\seealso{
   \code{\link{bayescount}}
   \code{\link{run.jags}}
}

\author{Matthew Denwood \email{m.denwood@vet.gla.ac.uk} funded as part of the DEFRA VTRI project 0101.}

\examples{

testjags()

}
\keyword{methods}