testjags <- function(jags="jags", silent=FALSE){
	tempfile <- paste((new_unique("temp")), ".cmd", sep="")
	write("exit", file=tempfile)
	s.info <- Sys.info()
	p.info <- .Platform
	
	os <- p.info$OS.type
	username <- as.character(s.info["user"])
	rversion <- R.version$version
	gui <- p.info$GUI
	p.type <- p.info$pkgType
	
	if(os=="windows"){
		success <- system(paste(jags, " ", tempfile, sep=""), ignore.stderr = TRUE, wait=TRUE)
		suppressWarnings(try(popen.support <- system(paste(jags, " ", tempfile, sep=""), intern=TRUE, wait=TRUE), silent=TRUE))
		if(class(popen.support)=="try-error"){
			popen <- FALSE
		}else{
			popen <- TRUE
		}
	}
	if(os=="unix"){
		success <- system(paste(jags, " ", tempfile, " > /dev/null", sep=""), ignore.stderr = TRUE, wait=TRUE)
		suppressWarnings(try(popen.support <- system(paste(jags, " ", tempfile, sep=""), intern=TRUE, wait=TRUE), silent=TRUE))
		if(class(popen.support)=="try-error"){
			popen <- FALSE
		}else{
			popen <- TRUE
		}
	}
	if(os != "unix" && os != "windows"){
		stop("Error analysing operating system")
	}
	
	if(silent==FALSE){
		cat("You are currently logged on as ", username, ", on a ", os, " machine\n", sep="")
		cat("You are using ", rversion, ", with the ", gui, " GUI", "\n", sep="")
		if(os=="windows"){
			cat("WARNING:  JAGS will run more slowly under windows than unix, and suppression of JAGS output may not be available\n")
		}
		if(success==0){
			if(popen == TRUE){
				version <- strsplit(popen.support[1], split=" ", fixed=TRUE)[[1]][4]
				cat("JAGS version ", version, " found successfully\n", sep="")
				num.version <- strsplit(version, split=".", fixed=TRUE)
				if(as.numeric(num.version[[1]][1] == 0) && as.numeric(num.version[[1]][2] < 97)){
					cat("This version of JAGS is no longer supported.  Please update JAGS from http://www-fis.iarc.fr/~martyn/software/jags/ and try again\n")
					jags.avail <- FALSE
				}else{
					jags.avail <- TRUE
				}
			}else{
				cat("JAGS was found on your system, but the version cannot be verified due to the absence of popen support.  Please ensure that the latest version of JAGS is installed by visiting http://www-fis.iarc.fr/~martyn/software/jags/\n")
				jags.avail <- TRUE
				num.version <- list("version unknown")
			}
		}else{
			cat("JAGS was not found on your system using the command '", jags, "'.  Please ensure that the command is correct and that the latest version of JAGS from http://www-fis.iarc.fr/~martyn/software/jags/ is installed\n", sep="")
			jags.avail <- FALSE
			num.version <- list("none found")
		}
	}else{
		if(success==0){
			if(popen == TRUE){
				version <- strsplit(popen.support[1], split=" ", fixed=TRUE)[[1]][4]
				num.version <- strsplit(version, split=".", fixed=TRUE)
				if(as.numeric(num.version[[1]][1] == 0) && as.numeric(num.version[[1]][2] < 97)){
					jags.avail <- FALSE
				}else{
					jags.avail <- TRUE
				}
			}else{
				jags.avail <- TRUE
				num.version <- list("version unknown")
			}
		}else{
			jags.avail <- FALSE
			popen <- FALSE
			num.version <- list("none found")
		}
	}
	unlink(tempfile)
	return(c("os"=os, "JAGS.available"=jags.avail, "popen.support"=popen, "JAGS.version"=num.version, "R.version"=rversion, "R.GUI"=gui, "R.package.type"=p.type, "username"=username))
}

testJAGS <- testjags
