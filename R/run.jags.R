run.jags <- function(data=stop("No data supplied"), model=stop("No model supplied"), inits = stop("No inital values supplied"), monitor = stop("No monitored variables supplied"), burnin = 5000, updates = 10000, jags="jags", silent.jags = FALSE, check.conv=TRUE){

if(burnin==0){
	## Having problems with adaptive phase, so:
	burnin <- 100
}

real.runs <- as.integer(updates)
ini.runs <- as.integer(burnin)

if(!require(lattice)){
	stop("The required library 'lattice' is not installed")
}
if(!require(coda)){
	stop("The required library 'coda' is not installed")
}

test <- testjags(jags, silent=TRUE)
if(test[[2]][1]==FALSE){
	cat("Unable to call JAGS using '", jags, "'\n", sep="")
	stop("Unable to call JAGS")
}

if(class(model)!="character" | length(model)!=1){
	stop("The model must be provided in the form of a single character string")
}

if(class(data)!="character" | length(data)!=1){
	stop("The data must be provided in R dump format (see dump.format())")
}

if(class(inits)!="character" | length(inits)==0){
	stop("Initial values must be provided as a character vector in the R dump format (see dump.format()), with length equal to the number of chains required")
}

chains <- length(inits)

if(class(monitor)!="character"){
	stop("Monitored variable(s) must be provided in the form of a character vector")
}

if(check.conv==TRUE & chains==1){
	cat("Warning:  Convergence cannot be assessed with only 1 chain\n")
	check.conv <- FALSE
}

modelstring <- paste(model, "\n", sep="")

monitors <- paste("monitor set <", paste(monitor, collapse=">\nmonitor set <"), ">\n", sep="")
n.params <- length(monitor)
params.names <- monitor
	
initstring <- paste(inits, "\n", sep="")
	
datastring <- paste(data, "\n", sep="")
	
#  Write model file

total.updates <- n.params * updates

save.directory <- getwd()
on.exit(setwd(save.directory))

temp.directory <- new_unique("tempfiles")
if((temp.directory=="Directory not writable")==TRUE){
	cat("Directory not writable\n")
	return(c("Error", "Write permissions"))
}

dir.create(temp.directory)
setwd(temp.directory)


output <- file("model.txt", 'w')
cat(modelstring, file=output,sep="")  
close(output)


## Write data and script and init files

output <- file("data.txt", 'w')
cat(datastring, file=output,sep="")  
close(output)

scriptstring <- paste("model in <\"model.txt\">
data in <\"data.txt\">
compile, nchains(<", as.integer(chains), ">)\n", sep="")
for(i in 1:chains){
scriptstring <- paste(scriptstring, "parameters in <\"inits", i, ".txt\">, chain(<", i, ">)\n", sep="")
}
scriptstring <- paste(scriptstring, "initialize
adapt <", ini.runs, ">
", monitors, "update <", real.runs, ">
coda <*>\n", sep="")
for(i in 1:chains){
scriptstring <- paste(scriptstring, "parameters to <\"out", i, ".Rdump\">, chain(<", i, ">)\n", sep="")
}
scriptstring <- paste(scriptstring, "exit\n", sep="")

output <- file("script.cmd", 'w')
cat(scriptstring, file=output,sep="")  
close(output)

for(i in 1:chains){
	output <- file(paste("inits", i, ".txt", sep=""), 'w')
	cat(initstring[i], file=output,sep="")
	close(output)
}

cat("Calling the simulation... (this may take some time)\n")

jags.status <- testjags(jags, silent=TRUE)

if (jags.status[[1]] == "windows"){
	if(jags.status[[3]] == TRUE){
		success <- system(paste(jags, " script.cmd", sep = ""), intern=TRUE, wait=TRUE)
	}else{
		success <- system(paste(jags, " script.cmd", sep = ""), wait=TRUE)
	}
}else{
	if(silent.jags == FALSE && jags.status[[3]]==TRUE){
		success <- system(paste(jags, "< script.cmd", sep=""), ignore.stderr = FALSE)
	}
	if(silent.jags == TRUE && jags.status[[3]]==TRUE){
		success <- system(paste(jags, "< script.cmd", sep=""), intern=TRUE, ignore.stderr = TRUE)
	}
	if(silent.jags == FALSE && jags.status[[3]]==FALSE){
		success <- system(paste(jags, "< script.cmd", sep=""), ignore.stderr = FALSE)
	}
	if(silent.jags == TRUE && jags.status[[3]]==FALSE){
		success <- system(paste(jags, "< script.cmd > /dev/null", sep=""), ignore.stderr = TRUE)
	}
}

if (file.exists("CODAindex.txt") == FALSE){
  	if (file.exists("JAGS.out") == TRUE){
 		cat("You are using a version of JAGS prior to 0.99.0, which is no longer supported.  Please update JAGS and try again\n")
   		stop("JAGS version not supported")
   	}else{
   		cat("ERROR:  The coda files were not found\n")
   		suppressWarnings(try(cat(success, "\n", sep=""), silent=TRUE))
   		setwd(save.directory)
		unlink(temp.directory, recursive = TRUE)
		results <- c("Unable to load coda files")
		return(results)
   	}
}


if(silent.jags == FALSE){
	for(i in 1:2000000) {
		hold <- exp(100)
	}
	cat("\n")
}

cat("Simulation complete.  Reading coda files...\n")

suppressWarnings(inputsuccess <- try(input.data <- read.openbugs(quiet=TRUE), silent=TRUE))

if((class(inputsuccess)=="try-error")){
	setwd(save.directory)
	unlink(temp.directory, recursive = TRUE)
	cat("Unable to load coda files\n")
	results <- c("Unable to load coda files")
	return(results)
}

if(chains==2){

	new.one <- input.data[[1]]
	new.two <- input.data[[2]]
	
	a <- as.numeric(new.one[,1])
	
	
	
	if(length(new.one) > length(new.two)){
		cat("Chain one was shortened to match the length of chain two\n")
		new.data <- new.two
		for(i in 1:n.params){
			new.data[,i] <- new.one[1:(length(new.two[,1])),i]
		}
		new.one <- as.mcmc(new.data)
	}
	if(length(new.two) > length(new.one)){
		cat("Chain two was shortened to match the length of chain one\n")
		new.data <- new.one
		for(i in 1:n.params){
			new.data[,i] <- new.two[1:(length(new.one[,1])),i]
		}
		new.two <- as.mcmc(new.data)
	}
	
	achieved <- length(new.one[,1])
		
	if(achieved!=updates){
		crashed <- TRUE
		cat("Simulation crashed at ", as.numeric(length(new.one[,1])), " iterations", sep="")
	}
	
	input.data[[1]] <- new.one
	input.data[[2]] <- new.two
	
}else{
	
	lengths = lengths.equal <- numeric(length=chains)
	for(i in 1:chains){
		lengths.equal[i] <- length(input.data[[1]]) != length(input.data[[i]])
		lengths[i] <- length(input.data[[i]])
	}
	if(sum(lengths.equal)!=0){
		pastechains <- paste(chains[lengths[chains]==min(lengths)], collapse=", ")
		cat("Warning:  Chain lengths not equal.  Simulation crashed at ", as.numeric(min(lengths / n.params)), " iterations for chain(s) ", pastechains, ".", sep="")
	}
}

input.end <- character(length=chains)
for(i in 1:chains){
	filename <- paste("out", i, ".Rdump", sep="")
	suppressWarnings(inputsuccess <- try(tempinput <- readLines(filename)))
	if(class(inputsuccess)=="try-error"){
		cat("Error reading output of chain ", i, ".", sep="")
		input.end[i] <- NA
	}else{
		input.end[i] <- ""
		for(j in 1:length(tempinput)){
			input.end[i] <- paste(input.end[i], tempinput[j], "\n", sep="")
		}
	}
}

cat("Coda files loaded successfully\n")

setwd(save.directory)
unlink(temp.directory, recursive = TRUE)

if(check.conv==TRUE){
	success <- try({
	mcmc <- vector('list', length=chains)
	
	for(i in 1:chains){
		mcmc[[i]] <-input.data[[i]][,1:n.params]
	}
	
	convergence <- gelman.diag(mcmc.list(mcmc))
	
	unconverged <- 0
	unconverged.string <- ""
	
	for(j in 1:n.params){
		param.conv <- convergence$psrf[j, 1]
		if(is.na(param.conv)){
			break
		}
		if(param.conv > 1.05){
			unconverged <- unconverged + 1
		}
	}
	
	if(!is.na(param.conv)){
		if(unconverged > 0){
			cat("Convergence failed for this run for ", unconverged, " parameter(s) after ", updates, " iterations (multi-variate psrf = ", round(convergence$mpsrf, digits=3), ")\n", sep="")
		}else{
			cat("Convergence achieved for this run\n")
		}
	}else{
		cat("The potential scale reduction factor could not be calculated for these chains\n")
	}
	
	}, silent=FALSE)
	if(class(success)=="try-error"){
		cat("An error occured when assessing convergence\n")
	}

}

return(c(input.data, input.end))

}

run.JAGS <- run.jags