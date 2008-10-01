new_unique <- function(name, suffix="", ask=FALSE, prompt="A file or directory with this name already exists.  Overwrite?"){
	temp <- paste(name, sep="")
	exists <- file.exists(paste(temp, suffix, sep=""))
	if(exists==TRUE){
		path.ok <- FALSE
		counter <- 1
		if(ask==TRUE){
			path.ok <- ask(paste("\'", name, suffix, "\'.  ", prompt, "  ", sep=""), type="logical")
			if(path.ok==TRUE){
				unlink(paste(temp, suffix, sep=""), recursive = TRUE)
			}
		}
		while(path.ok == FALSE){
			temp <- paste(name, "_", counter, "", sep="")
			exists <- file.exists(paste(temp, suffix, sep=""))
			if(exists==TRUE){
				counter <- counter + 1
			}else{
				path.ok <- TRUE
				break
			}
		}
	}
	suppressWarnings(try(dir.create(paste(temp, suffix, sep="")), silent=TRUE))
	permissions <- file.exists(paste(temp, suffix, sep=""))
	if(permissions==FALSE){
		cat("Error:  Directory not writable\n")
		return("Directory not writable")
	}else{
		unlink(paste(temp, suffix, sep=""), recursive = TRUE)
	}
	backupforspaces <- file.remove(paste(temp, suffix, sep=""))
	return(paste(temp, suffix, sep=""))
}