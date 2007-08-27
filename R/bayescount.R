##########################################################################################################################
##########################################################################################################################

####  Source script to handle/call function to analyse faecal egg count data

####  Created 3rd July 2007 by Matthew Denwood

##### REQUIREMENTS:

##### JAGS, coda, bayesmix.  Data can be a link to a .csv file in COMMA delimited form or a matrix.

##### CSV file is written to current path with results, and 'jobname'.'model'.results is copied to global parameter space


##########################################################################################################################

####  ALSO

####  Function to analyse faecal egg count data with a ZIGP / ZIP / GP / Poisson model then check converegence with 2 chains

####  Created 3rd July 2007 by Matthew Denwood

##### REQUIREMENTS:

##### JAGS 
##### The coda and bayesmix packages for R
##### Data should be a single vector


##########################################################################################################################

##### Put into same source file and P, GP, ZIP models created on 10th July 2007
##### Last modified 23rd July 2007

##########################################################################################################################
##########################################################################################################################


ask_yn <- function (prompt){
	ok <- FALSE
	while(ok==FALSE){
		result <- readline(prompt = prompt)
		if(result=="y"){
			result <- TRUE
			ok <- TRUE
		}
		if(result=="YES"){
			result <- TRUE
			ok <- TRUE
		}
		if(result=="yes"){
			result <- FALSE
			ok <- TRUE
		}
		if(result=="Yes"){
			result <- FALSE
			ok <- TRUE
		}
		if(result=="Y"){
			result <- TRUE
			ok <- TRUE
		}
		if(result=="n"){
			result <- FALSE
			ok <- TRUE
		}
		if(result=="NO"){
			result <- TRUE
			ok <- TRUE
		}
		if(result=="no"){
			result <- TRUE
			ok <- TRUE
		}
		if(result=="No"){
			result <- FALSE
			ok <- TRUE
		}
		if(result=="N"){
			result <- FALSE
			ok <- TRUE
		}
		if(ok==FALSE){
			print("ERROR:  Please enter 'y' or 'n'", quote=FALSE)
		}
	}
	return(result)
}

new_unique <- function(name, suffix="", ask=FALSE, prompt="A file or directory with this name already exists.  Overwrite?  "){
	temp <- paste(name, sep="")
	exists <- file.exists(paste(temp, suffix, sep=""))
	if(exists==TRUE){
		path.ok <- FALSE
		counter <- 1
		if(ask==TRUE){
			path.ok <- ask_yn(paste("[1] '", name, suffix, "'.  ", prompt, sep=""))
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
		print("Error:  Directory not writable", quote=FALSE)
		return("Directory not writable")
	}else{
		unlink(paste(temp, suffix, sep=""), recursive = TRUE)
	}
	backupforspaces <- file.remove(paste(temp, suffix, sep=""))
	return(paste(temp, suffix, sep=""))
}

bayescount <- function (name = NA, data = NA, setnames = NA, div = 1, p.model = FALSE, zip.model = FALSE, gp.model = FALSE, zigp.model = TRUE, burnin = 5000, updates = c(10000,100000,500000), jags = "jags", rownames = FALSE, remove.zeros = TRUE, remove.missing = TRUE, test = TRUE, log.prior = TRUE, write.file=TRUE, ml.equivalent = FALSE, adjust.mean = FALSE)
{

omit.zeros <- remove.zeros
runname <- name
poisson <- p.model
gpoisson <- gp.model
zipoisson <- zip.model
zigpoisson <- zigp.model

if(!require(lattice)){
	stop("The required library 'lattice' is not installed")
}
if(!require(coda)){
	stop("The required library 'coda' is not installed")
}
if(!require(bayesmix)){
	stop("The required library 'bayesmix' is not installed")
}

testjags <- haveJAGS(jags)
if(testjags==FALSE){
	print("", quote=FALSE)
	print(paste("Unable to call JAGS using '", jags, "'", sep=""))
	print("", quote=FALSE)
	stop("Unable to call JAGS")
}

testwritable <- new_unique("test")
if((testwritable=="Directory not writable")==TRUE){
	print("", quote=FALSE)
	print("The working directory is not writable.  Please change the working directory")
	print("", quote=FALSE)
	stop("Directory not writable")
}

datanames <- setnames

data <- as.matrix(data)

print("--- 'Bayescount': Analyse count data using a zero-inflated gamma Poisson model implemented in JAGS ---", quote=FALSE)
print("", quote=FALSE)
print("*PLEASE NOTE:  THIS SOFTWARE IS INTENDED FOR EDUCATIONAL PURPOSES ONLY AND SHOULD NOT BE RELIED UPON FOR REAL WORLD APPLICATIONS*", quote=FALSE)
print("*ANALYSING DATA USING MCMC SAMPLING CAN PRODUCE MISLEADING RESULTS IF YOU DO NOT UNDERSTAND THE PROCESSES INVOLVED*", quote=FALSE)
print("", quote=FALSE)
if(is.na(runname)==TRUE){
	runname <- readline(prompt = "[1] Please enter a name for this analysis:  ")
}

datain <- data
datana <- data
length(datana) <- 1

dataok=FALSE
if(is.na(datana)==FALSE){
	exists <- try(file.exists(datana), silent=TRUE)
	if((class(exists)=="try-error")==FALSE){
		if(exists==TRUE){
			suppressWarnings(data <- try(read.csv(datain, header=FALSE), silent=TRUE))
		}
	}
	suppressWarnings(valid.data <- try((length(as.matrix(data[,1])) > 1), silent=TRUE))
	if((class(valid.data)=="try-error")==TRUE){
		print("ERROR:  The path you have entered does not appear to be valid", quote=FALSE)
	}else{
		if(valid.data==FALSE){
			print("ERROR:  Invalid path / data", quote=FALSE) 
		}else{
			dataok=TRUE
		}
	}
	print("", quote=FALSE)
}
while(dataok==FALSE){
	datain <- readline(prompt = "[1] Please enter the path to a (comma delimited) CSV file containing the data (type 'exit' to quit):  ")
	if((datain=="exit")==TRUE){
		stop("User exited the program")
	}
	exists <- try(file.exists(datain), silent=TRUE)
	if((class(exists)=="try-error")==FALSE){
		if(exists==TRUE){
		data <- try(read.csv(datain, header=FALSE), silent=TRUE)
		if((class(data)=="try-error")==FALSE){
			valid.data <- try(length(data[,1]) > 1, silent=TRUE)
			if((class(valid.data)=="try-error")==TRUE){
				print("ERROR:  The path you have entered does not appear to be valid", quote=FALSE)
			}else{
				if(valid.data==FALSE){
					print("ERROR:  The data you have entered is of length less than 2", quote=FALSE)
				}else{
					dataok=TRUE
				}
			}
		}else{
			print("ERROR:  The path you entered does not appear to be valid; please retry (path may be absolute or relative)", quote=FALSE)
		}
		}else{
			print("ERROR:  The path you entered does not appear to be valid; please retry (path may be absolute or relative)", quote=FALSE)
		}
	}else{
		print("ERROR:  The path you entered does not appear to be valid; please retry (path may be absolute or relative)", quote=FALSE)
	}
	print("", quote=FALSE)
}

data <- as.matrix(data)

names <- NA

if(is.na(datanames[1])==FALSE){
	if(datanames[1]!=TRUE){
		if(datanames[1]!=FALSE){
			if(rownames==FALSE){
				if((length(datanames)==length(data[1,]))==FALSE){
					stop("The length of the character vector setnames does not match the lenth of the data provided")
				}
			}else{
				if((length(datanames)==(length(data[1,])-1))==FALSE){
					stop("The length of the character vector setnames does not match the lenth of the data provided")
				}
			}
			names <- datanames
			datanames <- TRUE
		}
	}
}



if(is.na(datanames[1])==TRUE){
	datanames <- ask_yn(prompt = "[1] Is the first row used for data labels?  ")
	print("", quote=FALSE)
}

if(is.na(names[1])==FALSE){
	if(rownames==TRUE){
		ind.names <- data[ ,1]
		data <- data[ , 2:length(data[1,])]
	}else{
		ind.names <- paste("Row ", 1:length(data[,1]), sep="")
	}
}else{
	if(rownames==TRUE){
		if(datanames==TRUE){
			ind.names <- data[2:length(data[,1]),1]
			names <- data[1,2:length(data[1,])]
			data <- data[2:length(data[,1]), 2:length(data[1,])]
		}else{
			ind.names <- data[ ,1]
			names <- 1:(length(data[1,])-1)
			names <- paste("Dataset ", names, sep="")
			data <- data[ , 2:length(data[1,])]
		}
	}else{
		if(datanames==TRUE){
			ind.names <- paste("Row ", 1:(length(data[,1])-1), sep="")
			names <- data[1,]
			data <- data[2:length(data[,1]), ]
		}else{
			ind.names <- paste("Row ", 1:length(data[,1]), sep="")
			names <- 1:length(data[1,])
			names <- paste("Dataset ", names, sep="")
		}
	}

}


data <- as.matrix(data)

animalsindata1 <- length(data[,1])

names <- as.matrix(names)
suppressWarnings(data[] <- as.numeric(data) / div)

if(rownames==TRUE){
	if(remove.missing==TRUE){
		print("Missing data cannot be removed if rownames are used.  Missing data will not be removed", quote=FALSE)
		print("", quote=FALSE)
		remove.missing <- FALSE
	}
}

tiers <- length(updates)



print("Settings are as follows:", quote=FALSE)
print(paste("Name of analysis:                             ", runname, sep=""), quote=FALSE)
print(paste("First row used for dataset names:             ", datanames, sep=""), quote=FALSE)
print(paste("First column used for individual names:       ", rownames, sep=""), quote=FALSE)
print(paste("Number of datasets:                           ", as.numeric(length(names)), sep=""), quote=FALSE)
print(paste("Number of animals in column one               ", as.numeric(animalsindata1), sep=""), quote=FALSE)
print(paste("Divide data by                                ", div, sep=""), quote=FALSE)
print(paste("Use Poisson model:                            ", poisson, sep=""), quote=FALSE)
print(paste("Use ZI-Poisson model:                         ", zipoisson, sep=""), quote=FALSE)
print(paste("Use Gamma Poisson model:                      ", gpoisson, sep=""), quote=FALSE)
print(paste("Use ZI-Gamma Poisson model:                   ", zigpoisson, sep=""), quote=FALSE)
print(paste("Omit datasets with all zero counts:           ", omit.zeros, sep=""), quote=FALSE)
print(paste("Remove missing data:                          ", remove.missing, sep=""), quote=FALSE)
print(paste("Burn-in updates:                              ", burnin, sep=""), quote=FALSE)
for(i in 1:tiers){
	print(paste("Sample updates (tier ", i, "):                      ", updates[i], sep=""), quote=FALSE)
}
if(log.prior==TRUE){
print(paste("Distribution of overdispersion prior:         Log uniform", sep=""), quote=FALSE)
}else{
print(paste("Distribution of overdispersion prior:         Uniform", sep=""), quote=FALSE)
}
print(paste("Duplicate log likelihood equivalent results:  ", ml.equivalent, sep=""), quote=FALSE)
print(paste("Adjust z-i means to mean of whole population: ", adjust.mean, sep=""), quote=FALSE)
print(paste("System call to activate JAGS:                 ", jags, sep=""), quote=FALSE)
print("", quote=FALSE)

proceed <- as.logical((test-1)*(test-1))
if(test==TRUE){

setdata <- as.integer(na.omit(data[,1]))
print("Testing the model function", quote=FALSE)
output <- ZIGP.model(burnin=1000, updates=1000, name = names[1], jags = jags, data = setdata, log.prior = log.prior)

if(output[3]==1){
	proceed <- ask_yn(prompt = "[1] The test (first) dataset crashed.  Continue with other datasets?  ")
}else{
	if(output[4]==0){
		if(output[2]==0){
			print ("", quote=FALSE)
			print(paste("--- Completed the test dataset ---", sep=""), quote=FALSE)
			print ("", quote=FALSE)
		}
		print("Test completed successfully", quote=FALSE)
		proceed <- TRUE
	}else{
		stop("JAGS call or permissions error.  Aborting")
	}
}

}

if(ml.equivalent==FALSE){
	headers <- paste("dataset", "converged", "crashed", "error", "iterations", "mean.l.95", "mean.median", "mean.u.95", "zi.l.95", "zi.median", "zi.u.95", "od.l.95", "od.median", "od.u.95", "\n", sep = ",")
	zeros <- c("zeros", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA")
	resultscol <- 14
}else{
	headers <- paste("dataset", "converged", "crashed", "error", "iterations", "mean.l.95", "mean.median", "mean.u.95", "zi.l.95", "zi.median", "zi.u.95", "od.l.95", "od.median", "od.u.95", "log.mean.l.95", "log.mean.median", "log.mean.u.95", "logit.zi.l.95", "logit.zi.median", "logit.zi.u.95", "wlog.od.l.95", "wlog.od.median", "wlog.od.u.95", "\n", sep = ",")
	zeros <- c("zeros", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA")
	resultscol <- 14 + 9
}

p.results <- matrix(NA, ncol=resultscol, nrow=length(names))
crashed <- 0
unconverged <- 0
error <- 0

if(poisson==TRUE){
	print("", quote=FALSE)
	p.name <- new_unique(name=paste(runname, ".p", sep=""), suffix=".csv", ask = test, prompt="A results file with this name already exists.  Overwrite?  ")
	p.outfile <- file(p.name, 'w')
	cat(headers, file = p.outfile, sep = "")
	for(i in 1:length(names)){
		tier <- 1
		coda <- 0
		done <- FALSE
		if(remove.missing==TRUE){
			setdata <- as.integer(na.omit(data[,i]))
		}else{
			setdata <- as.integer(data[,i])
		}
		if(omit.zeros==TRUE){
			if((sum(na.omit(setdata)) < 1)==TRUE){
				p.results[i,] <- c(names[i], zeros)
				print(paste("Dataset '", names[i], "' contained all zeros and was therefore omitted", sep=""), quote=FALSE)
				done <- TRUE
			}
		}
		while(done==FALSE){
			output <- P.model(burnin=burnin, updates=updates[tier], name = names[i], jags = jags, data = setdata, ml.equivalent = ml.equivalent)
			if((output[4]=="1")==TRUE){
				print(paste("Dataset '", names[i], "' returned an error", sep=""), quote=FALSE)
				output <- P.model(burnin=burnin, updates=updates[tier], name = names[i], jags = jags, data = setdata, ml.equivalent = ml.equivalent)
			}
			if((output[4]!="1")==TRUE){
				if((output[3]=="1")==TRUE){
					print(paste("Dataset '", names[i], "' crashed", sep=""), quote=FALSE)
					if(is.na(output[2])==FALSE){
						if((output[2]!="1")==TRUE){
							output <- P.model(burnin=burnin, updates=updates[tier], name = names[i], jags = jags, data = setdata, ml.equivalent = ml.equivalent)
						}else{
							p.results[i,] <- output
							done <- TRUE
						}
					}else{
						output <- P.model(burnin=burnin, updates=updates[tier], name = names[i], jags = jags, data = setdata, ml.equivalent = ml.equivalent)
					}
				}
				if((output[4]!="1")==TRUE){
					if((output[3]!="1")==TRUE){
						if((output[2]=="0")==TRUE){
							print(paste("Dataset '", names[i], "' failed to converge at ", updates[tier], " iterations", sep=""), quote=FALSE)
							tier <- tier + 1
							if(tier > tiers){
								unconverged <- unconverged + 1
								p.results[i,] <- output
								done <- TRUE
							}
						}else{
							p.results[i,] <- output
							done <- TRUE
						}
					}else{
						print(paste("Dataset '", names[i], "' crashed", sep=""), quote=FALSE)
						if(is.na(output[2])==FALSE){
							if((output[2]=="0")==TRUE){
								crashed <- crashed + 1
							}
						}else{
							crashed <- crashed + 1
						}
						p.results[i,] <- output
						done <- TRUE
					}
				}else{
					print(paste("Dataset '", names[i], "' returned an error", sep=""), quote=FALSE)
					error <- error + 1
					p.results[i,] <- output
					done <- TRUE
				}
			}else{
				print(paste("Dataset '", names[i], "' returned an error", sep=""), quote=FALSE)
				error <- error + 1
				p.results[i,] <- output
				done <- TRUE
			}
		}
		print(paste("Dataset '", names[i], "' completed", sep=""), quote=FALSE)
		print(paste(i, " of ", as.numeric(length(names)), " datasets completed", sep=""), quote=FALSE)
		print(paste(unconverged, " failed convergence, ", crashed, " crashed, and ", error, " quit with an error", sep=""), quote=FALSE)
		print("", quote=FALSE)
		if(ml.equivalent==FALSE){
			cat(p.results[i,1], p.results[i,2], p.results[i,3], p.results[i,4], p.results[i,5], p.results[i,6], p.results[i,7], p.results[i,8], p.results[i,9], p.results[i,10], p.results[i,11], p.results[i,12], p.results[i,13], p.results[i,14], "\n", file = p.outfile, sep = ",")
		}else{
			cat(p.results[i,1], p.results[i,2], p.results[i,3], p.results[i,4], p.results[i,5], p.results[i,6], p.results[i,7], p.results[i,8], p.results[i,9], p.results[i,10], p.results[i,11], p.results[i,12], p.results[i,13], p.results[i,14], p.results[i,15], p.results[i,16], p.results[i,17], p.results[i,18], p.results[i,19], p.results[i,20], p.results[i,21], p.results[i,22], p.results[i,23], "\n", file = p.outfile, sep = ",")
		}
	}
	close(p.outfile)
	if(write.file==FALSE){
		unlink(p.name)
	}
	assign(paste(runname, ".p.results", sep=""), p.results, pos=".GlobalEnv")
}

gp.results <- matrix(NA, ncol=resultscol, nrow=length(names))
crashed <- 0
unconverged <- 0
error <- 0
gp.largeod <- 0
assign("gp.largeod", gp.largeod, pos=".GlobalEnv")

if(gpoisson==TRUE){
	print("", quote=FALSE)
	gp.name <- new_unique(name=paste(runname, ".gp", sep=""), suffix=".csv", ask = test, prompt="A results file with this name already exists.  Overwrite?  ")
	gp.outfile <- file(gp.name, 'w')
	cat(headers, file = gp.outfile, sep = "")
	for(i in 1:length(names)){
		tier <- 1
		coda <- 0
		done <- FALSE
		if(remove.missing==TRUE){
			setdata <- as.integer(na.omit(data[,i]))
		}else{
			setdata <- as.integer(data[,i])
		}
		if(omit.zeros==TRUE){
			if((sum(na.omit(setdata)) < 1)==TRUE){
				gp.results[i,] <- c(names[i], zeros)
				print(paste("Dataset '", names[i], "' contained all zeros and was therefore omitted", sep=""), quote=FALSE)
				done <- TRUE
			}
		}
		while(done==FALSE){
			output <- GP.model(burnin=burnin, updates=updates[tier], name = names[i], jags = jags, data = setdata, log.prior = log.prior, ml.equivalent = ml.equivalent)
			if((output[4]=="1")==TRUE){
				print(paste("Dataset '", names[i], "' returned an error", sep=""), quote=FALSE)
				output <- GP.model(burnin=burnin, updates=updates[tier], name = names[i], jags = jags, data = setdata, log.prior = log.prior, ml.equivalent = ml.equivalent)
			}
			if((output[4]!="1")==TRUE){
				if((output[3]=="1")==TRUE){
					print(paste("Dataset '", names[i], "' crashed", sep=""), quote=FALSE)
					if(is.na(output[2])==FALSE){
						if((output[2]!="1")==TRUE){
							output <- GP.model(burnin=burnin, updates=updates[tier], name = names[i], jags = jags, data = setdata, log.prior = log.prior, ml.equivalent = ml.equivalent)
						}else{
							gp.results[i,] <- output
							done <- TRUE
						}
					}else{
						output <- GP.model(burnin=burnin, updates=updates[tier], name = names[i], jags = jags, data = setdata, log.prior = log.prior, ml.equivalent = ml.equivalent)
					}
				}
				if((output[4]!="1")==TRUE){
					if((output[3]!="1")==TRUE){
						if((output[2]=="0")==TRUE){
							print(paste("Dataset '", names[i], "' failed to converge at ", updates[tier], " iterations", sep=""), quote=FALSE)
							tier <- tier + 1
							if(tier > tiers){
								unconverged <- unconverged + 1
								gp.results[i,] <- output
								done <- TRUE
							}
						}else{
							gp.results[i,] <- output
							done <- TRUE
						}
					}else{
						print(paste("Dataset '", names[i], "' crashed", sep=""), quote=FALSE)
						if(is.na(output[2])==FALSE){
							if((output[2]=="0")==TRUE){
								crashed <- crashed + 1
							}
						}else{
							crashed <- crashed + 1
						}
						gp.results[i,] <- output
						done <- TRUE
					}
				}else{
					print(paste("Dataset '", names[i], "' returned an error", sep=""), quote=FALSE)
					error <- error + 1
					gp.results[i,] <- output
					done <- TRUE
				}
			}else{
				print(paste("Dataset '", names[i], "' returned an error", sep=""), quote=FALSE)
				error <- error + 1
				gp.results[i,] <- output
				done <- TRUE
			}
		}
		print(paste("Dataset '", names[i], "' completed", sep=""), quote=FALSE)
		print(paste(i, " of ", as.numeric(length(names)), " datasets completed", sep=""), quote=FALSE)
		print(paste(unconverged, " failed convergence, ", crashed, " crashed, and ", error, " quit with an error", sep=""), quote=FALSE)
		print("", quote=FALSE)
		if(ml.equivalent==FALSE){
			cat(gp.results[i,1], gp.results[i,2], gp.results[i,3], gp.results[i,4], gp.results[i,5], gp.results[i,6], gp.results[i,7], gp.results[i,8], gp.results[i,9], gp.results[i,10], gp.results[i,11], gp.results[i,12], gp.results[i,13], gp.results[i,14], "\n", file = gp.outfile, sep = ",")
		}else{
			cat(gp.results[i,1], gp.results[i,2], gp.results[i,3], gp.results[i,4], gp.results[i,5], gp.results[i,6], gp.results[i,7], gp.results[i,8], gp.results[i,9], gp.results[i,10], gp.results[i,11], gp.results[i,12], gp.results[i,13], gp.results[i,14], gp.results[i,15], gp.results[i,16], gp.results[i,17], gp.results[i,18], gp.results[i,19], gp.results[i,20], gp.results[i,21], gp.results[i,22], gp.results[i,23], "\n", file = gp.outfile, sep = ",")
		}
	}
	close(gp.outfile)
	if(write.file==FALSE){
		unlink(gp.name)
	}
	assign(paste(runname, ".gp.results", sep=""), gp.results, pos=".GlobalEnv")
}

zip.results <- matrix(NA, ncol=resultscol, nrow=length(names))
crashed <- 0
unconverged <- 0
error <- 0

if(zipoisson==TRUE){
	print("", quote=FALSE)
	zip.name <- new_unique(name=paste(runname, ".zip", sep=""), suffix=".csv", ask = test, prompt="A results file with this name already exists.  Overwrite?  ")
	zip.outfile <- file(zip.name, 'w')
	cat(headers, file = zip.outfile, sep = "")
	for(i in 1:length(names)){
		tier <- 1
		coda <- 0
		done <- FALSE
		if(remove.missing==TRUE){
			setdata <- as.integer(na.omit(data[,i]))
		}else{
			setdata <- as.integer(data[,i])
		}
		if(omit.zeros==TRUE){
			if((sum(na.omit(setdata)) < 1)==TRUE){
				zip.results[i,] <- c(names[i], zeros)
				print(paste("Dataset '", names[i], "' contained all zeros and was therefore omitted", sep=""), quote=FALSE)
				done <- TRUE
			}
		}
		while(done==FALSE){
			output <- ZIP.model(burnin=burnin, updates=updates[tier], name = names[i], jags = jags, data = setdata, ml.equivalent = ml.equivalent, adjust.mean = adjust.mean)
			if((output[4]=="1")==TRUE){
				print(paste("Dataset '", names[i], "' returned an error", sep=""), quote=FALSE)
				output <- ZIP.model(burnin=burnin, updates=updates[tier], name = names[i], jags = jags, data = setdata, ml.equivalent = ml.equivalent, adjust.mean = adjust.mean)
			}
			if((output[4]!="1")==TRUE){
				if((output[3]=="1")==TRUE){
					print(paste("Dataset '", names[i], "' crashed", sep=""), quote=FALSE)
					if(is.na(output[2])==FALSE){
						if((output[2]!="1")==TRUE){
							output <- ZIP.model(burnin=burnin, updates=updates[tier], name = names[i], jags = jags, data = setdata, ml.equivalent = ml.equivalent, adjust.mean = adjust.mean)
						}else{
							zip.results[i,] <- output
							done <- TRUE
						}
					}else{
						output <- ZIP.model(burnin=burnin, updates=updates[tier], name = names[i], jags = jags, data = setdata, ml.equivalent = ml.equivalent, adjust.mean = adjust.mean)
					}
				}
				if((output[4]!="1")==TRUE){
					if((output[3]!="1")==TRUE){
						if((output[2]=="0")==TRUE){
							print(paste("Dataset '", names[i], "' failed to converge at ", updates[tier], " iterations", sep=""), quote=FALSE)
							tier <- tier + 1
							if(tier > tiers){
								unconverged <- unconverged + 1
								zip.results[i,] <- output
								done <- TRUE
							}
						}else{
							zip.results[i,] <- output
							done <- TRUE
						}
					}else{
						print(paste("Dataset '", names[i], "' crashed", sep=""), quote=FALSE)
						if(is.na(output[2])==FALSE){
							if((output[2]=="0")==TRUE){
								crashed <- crashed + 1
							}
						}else{
							crashed <- crashed + 1
						}
						zip.results[i,] <- output
						done <- TRUE
					}
				}else{
					print(paste("Dataset '", names[i], "' returned an error", sep=""), quote=FALSE)
					error <- error + 1
					zip.results[i,] <- output
					done <- TRUE
				}
			}else{
				print(paste("Dataset '", names[i], "' returned an error", sep=""), quote=FALSE)
				error <- error + 1
				zip.results[i,] <- output
				done <- TRUE
			}
		}
		print(paste("Dataset '", names[i], "' completed", sep=""), quote=FALSE)
		print(paste(i, " of ", as.numeric(length(names)), " datasets completed", sep=""), quote=FALSE)
		print(paste(unconverged, " failed convergence, ", crashed, " crashed, and ", error, " quit with an error", sep=""), quote=FALSE)
		print("", quote=FALSE)
		if(ml.equivalent==FALSE){
			cat(zip.results[i,1], zip.results[i,2], zip.results[i,3], zip.results[i,4], zip.results[i,5], zip.results[i,6], zip.results[i,7], zip.results[i,8], zip.results[i,9], zip.results[i,10], zip.results[i,11], zip.results[i,12], zip.results[i,13], zip.results[i,14], "\n", file = zip.outfile, sep = ",")
		}else{
			cat(zip.results[i,1], zip.results[i,2], zip.results[i,3], zip.results[i,4], zip.results[i,5], zip.results[i,6], zip.results[i,7], zip.results[i,8], zip.results[i,9], zip.results[i,10], zip.results[i,11], zip.results[i,12], zip.results[i,13], zip.results[i,14], zip.results[i,15], zip.results[i,16], zip.results[i,17], zip.results[i,18], zip.results[i,19], zip.results[i,20], zip.results[i,21], zip.results[i,22], zip.results[i,23], "\n", file = zip.outfile, sep = ",")
		}
	}
	close(zip.outfile)
	if(write.file==FALSE){
		unlink(zip.name)
	}
	assign(paste(runname, ".zip.results", sep=""), zip.results, pos=".GlobalEnv")
}

zigp.results <- matrix(NA, ncol=resultscol, nrow=length(names))
crashed <- 0
unconverged <- 0
error <- 0
zigp.largeod <- 0
assign("zigp.largeod", zigp.largeod, pos=".GlobalEnv")

if(zigpoisson==TRUE){
	print("", quote=FALSE)
	zigp.name <- new_unique(name=paste(runname, ".zigp", sep=""), suffix=".csv", ask = test, prompt="A results file with this name already exists.  Overwrite?  ")
	zigp.outfile <- file(zigp.name, 'w')
	cat(headers, file = zigp.outfile, sep = "")
	for(i in 1:length(names)){
		tier <- 1
		coda <- 0
		done <- FALSE
		if(remove.missing==TRUE){
			setdata <- as.integer(na.omit(data[,i]))
		}else{
			setdata <- as.integer(data[,i])
		}
		if(omit.zeros==TRUE){
			if((sum(na.omit(setdata)) < 1)==TRUE){
				zigp.results[i,] <- c(names[i], zeros)
				print(paste("Dataset '", names[i], "' contained all zeros and was therefore omitted", sep=""), quote=FALSE)
				done <- TRUE
			}
		}
		while(done==FALSE){
			output <- ZIGP.model(burnin=burnin, updates=updates[tier], name = names[i], jags = jags, data = setdata, log.prior = log.prior, ml.equivalent = ml.equivalent, adjust.mean = adjust.mean)
			if((output[4]=="1")==TRUE){
				print(paste("Dataset '", names[i], "' returned an error", sep=""), quote=FALSE)
				output <- ZIGP.model(burnin=burnin, updates=updates[tier], name = names[i], jags = jags, data = setdata, log.prior = log.prior, ml.equivalent = ml.equivalent, adjust.mean = adjust.mean)
			}
			if((output[4]!="1")==TRUE){
				if((output[3]=="1")==TRUE){
					print(paste("Dataset '", names[i], "' crashed", sep=""), quote=FALSE)
					if(is.na(output[2])==FALSE){
						if((output[2]!="1")==TRUE){
							output <- ZIGP.model(burnin=burnin, updates=updates[tier], name = names[i], jags = jags, data = setdata, log.prior = log.prior, ml.equivalent = ml.equivalent, adjust.mean = adjust.mean)
						}else{
							zigp.results[i,] <- output
							done <- TRUE
						}
					}else{
						output <- ZIGP.model(burnin=burnin, updates=updates[tier], name = names[i], jags = jags, data = setdata, log.prior = log.prior, ml.equivalent = ml.equivalent, adjust.mean = adjust.mean)
					}
				}
				if((output[4]!="1")==TRUE){
					if((output[3]!="1")==TRUE){
						if((output[2]=="0")==TRUE){
							print(paste("Dataset '", names[i], "' failed to converge at ", updates[tier], " iterations", sep=""), quote=FALSE)
							tier <- tier + 1
							if(tier > tiers){
								unconverged <- unconverged + 1
								zigp.results[i,] <- output
								done <- TRUE
							}
						}else{
							zigp.results[i,] <- output
							done <- TRUE
						}
					}else{
						print(paste("Dataset '", names[i], "' crashed", sep=""), quote=FALSE)
						if(is.na(output[2])==FALSE){
							if((output[2]=="0")==TRUE){
								crashed <- crashed + 1
							}
						}else{
							crashed <- crashed + 1
						}
						zigp.results[i,] <- output
						done <- TRUE
					}
				}else{
					print(paste("Dataset '", names[i], "' returned an error", sep=""), quote=FALSE)
					error <- error + 1
					zigp.results[i,] <- output
					done <- TRUE
				}
			}else{
				print(paste("Dataset '", names[i], "' returned an error", sep=""), quote=FALSE)
				error <- error + 1
				zigp.results[i,] <- output
				done <- TRUE
			}
		}
		print(paste("Dataset '", names[i], "' completed", sep=""), quote=FALSE)
		print(paste(i, " of ", as.numeric(length(names)), " datasets completed", sep=""), quote=FALSE)
		print(paste(unconverged, " failed convergence, ", crashed, " crashed, and ", error, " quit with an error", sep=""), quote=FALSE)
		print("", quote=FALSE)
		if(ml.equivalent==FALSE){
			cat(zigp.results[i,1], zigp.results[i,2], zigp.results[i,3], zigp.results[i,4], zigp.results[i,5], zigp.results[i,6], zigp.results[i,7], zigp.results[i,8], zigp.results[i,9], zigp.results[i,10], zigp.results[i,11], zigp.results[i,12], zigp.results[i,13], zigp.results[i,14], "\n", file = zigp.outfile, sep = ",")
		}else{
			cat(zigp.results[i,1], zigp.results[i,2], zigp.results[i,3], zigp.results[i,4], zigp.results[i,5], zigp.results[i,6], zigp.results[i,7], zigp.results[i,8], zigp.results[i,9], zigp.results[i,10], zigp.results[i,11], zigp.results[i,12], zigp.results[i,13], zigp.results[i,14], zigp.results[i,15], zigp.results[i,16], zigp.results[i,17], zigp.results[i,18], zigp.results[i,19], zigp.results[i,20], zigp.results[i,21], zigp.results[i,22], zigp.results[i,23], "\n", file = zigp.outfile, sep = ",")
		}
	}
	close(zigp.outfile)
	if(write.file==FALSE){
		unlink(zigp.name)
	}
	assign(paste(runname, ".zigp.results", sep=""), zigp.results, pos=".GlobalEnv")
}

final.gp.largeod <- as.integer(getAnywhere(gp.largeod)[1])
final.zigp.largeod <- as.integer(getAnywhere(zigp.largeod)[1])

print("--- End ---", quote=FALSE)


if((final.gp.largeod > 0)==TRUE){
	if(log.prior==FALSE){
		print("", quote=FALSE)
		print(paste("*WARNING*  The 95% confidence interval for overdispersion was very large for a total of ", final.gp.largeod, " dataset(s) using the gamma Poisson model, indicating a lack of information about this parameter in these datasets.  These models should be re-run with a log uniform prior distribution for overdispersion and the two sets of results compared to ensure that the results for all parameters are reliable", sep=""), quote=FALSE)
		
	}else{
		print("", quote=FALSE)
		print(paste("*WARNING*  The 95% confidence interval for overdispersion was very large for a total of ", final.gp.largeod, " dataset(s) using the gamma Poisson model, indicating a lack of information about this parameter in these datasets.  These models should be re-run with a uniform prior distribution for overdispersion and the two sets of results compared to ensure that the results for all parameters are reliable", sep=""), quote=FALSE)
		
	}
}
remove(gp.largeod, pos=".GlobalEnv")
if((final.zigp.largeod > 0)==TRUE){
	if(log.prior==FALSE){
		print("", quote=FALSE)
		print(paste("*WARNING*  The 95% confidence interval for overdispersion was very large for a total of ", final.zigp.largeod, " dataset(s) using the zero-inflated gamma Poisson model, indicating a lack of information about this parameter in these datasets.  These models should be re-run with a log uniform prior distribution for overdispersion and the two sets of results compared to ensure that the results for all parameters are reliable", sep=""), quote=FALSE)
		
	}else{
		print("", quote=FALSE)
		print(paste("*WARNING*  The 95% confidence interval for overdispersion was very large for a total of ", final.zigp.largeod, " dataset(s) using the zero-inflated gamma Poisson model, indicating a lack of information about this parameter in these datasets.  These models should be re-run with a uniform prior distribution for overdispersion and the two sets of results compared to ensure that the results for all parameters are reliable", sep=""), quote=FALSE)
		
	}
}
remove(zigp.largeod, pos=".GlobalEnv")

}



##########################################################################################################################

#######  P MODEL:

##########################################################################################################################


P.model <- function (data, name, jags = "jags", burnin = 5000, updates = 10000, ml.equivalent = FALSE)
{

bugcheck <- FALSE

datana <- data

length(datana) <- 1

if(is.na(name)==TRUE){
	print("This job was not given a name", quote=FALSE)
	return(c("Error", "Name"))
}
if(is.na(datana)==TRUE){
	print("No data specified", quote=FALSE)
	return(c("Error", "Data"))
}


real.runs <- as.integer(updates)
ini.runs <- as.integer(burnin)
jobname <- name


if(!require(lattice)){
	stop("The required library 'lattice' is not installed")
}
if(!require(coda)){
	stop("The required library 'coda' is not installed")
}
if(!require(bayesmix)){
	stop("The required library 'bayesmix' is not installed")
}


test <- haveJAGS(jags)
if(test==FALSE){
	print(paste("Unable to call JAGS using '", jags, "'", sep=""))
	return(c("Error", "JAGS call"))
}

save.directory <- getwd()
on.exit(setwd(save.directory))

temp.directory <- new_unique("tempfiles")
if((temp.directory=="Directory not writable")==TRUE){
	print("Directory not writable", quote=FALSE)
	return(c("Error", "Write permissions"))
}

dir.create(temp.directory)
setwd(temp.directory)

if(ml.equivalent==FALSE){
	errorpadding <- c("NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA")
}else{
	errorpadding <- c("NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA")
}

print ("", quote=FALSE)
print(paste("--- Running Poisson model for dataset '", jobname, "' at ", updates, " iterations ---", sep=""), quote=FALSE)
print ("", quote=FALSE)


counts <- data

N <- length(counts)



##### Write model, data and script file


modelstring <- paste("model {

for(row in 1 : N){
Count[row] ~ dpois(mean);
}

b <- 1;
prob <- 1;

# Priors
mean ~ dunif(0.001,1000);
}")

output <- file("model.txt", 'w')
cat(modelstring, file=output,sep="")  
close(output)

datastring <- paste("\"N\" <- ", N, "\n\"Count\" <- c(", sep="")
for(i in 1:(N-1)){
	datastring <- paste(datastring, counts[i], ", ", sep="")
}
datastring <- paste(datastring, counts[N], ")\n\n", sep="")

output <- file("data.txt", 'w')
cat(datastring, file=output,sep="")  
close(output)


scriptstring <- paste("model in <\"model.txt\">
data in <\"data.txt\">
compile
parameters in <\"inits.txt\">

initialize
update <", ini.runs, ">
monitor set <mean>
monitor set <b>
monitor set <prob>
update <", real.runs, ">
coda <*>
exit
\n", sep="")

output <- file("script.cmd", 'w')
cat(scriptstring, file=output,sep="")  
close(output)

totalupdates <- 3 * real.runs
crashed <- FALSE
converged <- NA
error <- 0
achieved <- NA

dataformean <- data
#dataformean[dataformean==0] <- NA
meanofnonzeros <- as.integer(mean(na.omit(dataformean)))

smallestmean <- as.integer(meanofnonzeros / 10)
largestmean <- as.integer(meanofnonzeros * 10)

if(is.na(smallestmean)==TRUE){
	smallestmean <- 1
}
if(is.na(largestmean)==TRUE){
	largestmean <- 10
}

if((smallestmean < 1)==TRUE){
	smallestmean <- 1
}
if((largestmean < 10)==TRUE){
	largestmean <- 10
}
if((smallestmean > 20)==TRUE){
	smallestmean <- 20
}
if((largestmean > 200)==TRUE){
	largestmean <- 200
}


#############  Write first init files and call JAGS to run simulation, then input results


probpos <- counts
zeros <- probpos==0
probpos[] <- 1
probpos[zeros] <- 0

initsone <- paste("\"mean\"  <-  ", smallestmean, "\n\n", sep="")

output <- file("inits.txt", 'w')
cat(initsone, file=output,sep="")
close(output)


print("Running the first simulation", quote=FALSE)
JAGScall("script", jags, quiet=TRUE)


print("Loading first set of coda files", quote=FALSE)
suppressWarnings(inputonesuccess <- try(input.data.one <- read.coda("jags.out", "jags.ind", quiet=TRUE), silent=TRUE))
if((class(inputonesuccess)=="try-error")==FALSE){
	if(length(input.data.one)!=totalupdates){
		crashed <- TRUE
		print(paste("First simulation crashed at ", as.numeric(length(input.data.one))/3, " iterations", sep=""), quote=FALSE)
	}
}else{
	setwd(save.directory)
	unlink(temp.directory, recursive = TRUE)
	print("Unable to load coda files", quote=FALSE)
	error <- 1
	crashed <- as.integer(crashed)
	results <- c(jobname, converged, crashed, error, achieved, errorpadding)
	return(results)
}



#############  Write second init files and call JAGS to run simulation, then input results


probpos <- counts
zeros <- probpos==0
probpos[] <- 1
probpos[zeros] <- 1


initsone <- paste("\"mean\"  <-  ", largestmean, "\n\n", sep="")

output <- file("inits.txt", 'w')
cat(initsone, file=output,sep="")
close(output)


print("Running the second simulation", quote=FALSE)
JAGScall("script", jags, quiet=TRUE)


print("Loading second set of coda files", quote=FALSE)
suppressWarnings(inputtwosuccess <- try(input.data.two <- read.coda("jags.out", "jags.ind", quiet=TRUE), silent=TRUE))
if((class(inputtwosuccess)=="try-error")==FALSE){
	if(length(input.data.two)!=totalupdates){
		crashed <- TRUE
		print(paste("Second simulation crashed at ", as.numeric(length(input.data.two))/3, " iterations", sep=""), quote=FALSE)
	}
}else{
	setwd(save.directory)
	unlink(temp.directory, recursive = TRUE)
	crashed <- as.integer(crashed)
	error <- 1
	results <- c(jobname, converged, crashed, error, achieved, errorpadding)
	print("Unable to load coda files", quote=FALSE)
	return(results)
}

if(crashed==TRUE){
	if(length(input.data.one) > length(input.data.two)){
		print("Chain one was shortened to match the length of chain two", quote=FALSE)
		new.data <- input.data.two
		new.data[,1] <- input.data.one[1:(length(input.data.two)/3),1]
		new.data[,2] <- input.data.one[1:(length(input.data.two)/3),2]
		new.data[,3] <- input.data.one[1:(length(input.data.two)/3),3]
		input.data.one <- as.mcmc(new.data)
	}
	if(length(input.data.two) > length(input.data.one)){
		print("Chain two was shortened to match the length of chain one", quote=FALSE)
		new.data <- input.data.one
		new.data[,1] <- input.data.two[1:(length(input.data.one)/3),1]
		new.data[,2] <- input.data.two[1:(length(input.data.one)/3),2]
		new.data[,3] <- input.data.two[1:(length(input.data.one)/3),3]
		input.data.two <- as.mcmc(new.data)
	}
}

achieved <- length(input.data.one) / 3

if(length(input.data.one)!=length(input.data.two)){
	setwd(save.directory)
	unlink(temp.directory, recursive = TRUE)
	print("Chains were returned as different lengths", quote=FALSE)
	crashed <- as.integer(crashed)
	error <- 1
	results <- c(jobname, converged, crashed, error, achieved, errorpadding)
	return(results)
}

if(updates > 999){
	if(length(input.data.one) < 3000){
		setwd(save.directory)
		unlink(temp.directory, recursive = TRUE)
		crashed <- as.integer(crashed)
		print("The model crashed before 1000 sampled iterations", quote=FALSE)
		results <- c(jobname, converged, crashed, error, achieved, errorpadding)
		return(results)
	}
}else{
	if(length(input.data.one) < (3 * updates)){
		setwd(save.directory)
		unlink(temp.directory, recursive = TRUE)
		message <- paste("The model crashed before ", updates, " sampled iterations", sep="")
		print(message, quote=FALSE)
		crashed <- as.integer(crashed)
		results <- c(jobname, converged, crashed, error, achieved, errorpadding)
		return(results)
	}
}

achieved <- length(input.data.one) / 3

print("Assessing convergence", quote=FALSE)

mcmclist <- try(input.data <- mcmc.list(input.data.one, input.data.two), silent = TRUE)

if((class(mcmclist)=="try-error")==TRUE){
	setwd(save.directory)
	unlink(temp.directory, recursive = TRUE)
	if(bugcheck==TRUE){
		assign("input.data.one", input.data.one, pos=".GlobalEnv")
		assign("input.data.two", input.data.two, pos=".GlobalEnv")
	}
	print("There was an error combining the chains", quote=FALSE)
	crashed <- as.integer(crashed)
	error <- 1
	results <- c(jobname, converged, crashed, error, achieved, errorpadding)
	return(results)
}

one <- as.mcmc(input.data[1])
two <- as.mcmc(input.data[2])

gelmans <- matrix(ncol=2, nrow=3)

parameternames <- c("mean", "b", "prob")

unconverged <- ""
n.unconv <- 0
converged <- matrix(NA, nrow=3)

for (i in 1:3){
	thing <- mcmc.list(as.mcmc(one[,i]), as.mcmc(two[,i]))
	gelmans[i,1] <- gelman.diag(thing)$psrf[1]
	gelmans[i,2] <- gelman.diag(thing)$psrf[2]
	if(is.na(gelmans[i,1])==FALSE){
		if(gelmans[i,1] > 1.05){
			unconverged <- paste(unconverged, parameternames[i], ", ", sep="")
			n.unconv <- n.unconv + 1
			converged[i] <- "No"
		}else{
			converged[i] <- "Yes"
		}
	}else{
		converged[i] <- "Yes"
	}
}

if(n.unconv > 0){
	setwd(save.directory)
	unlink(temp.directory, recursive = TRUE)
	print(paste("The following parameters failed to converge:  ", unconverged, sep=""), quote=FALSE)
	converged <- 0
	
}else{
	converged <- 1
}


###############  Analyse coda files

print("Calculating results", quote=FALSE)
print("", quote=FALSE)
mean <- c(as.matrix(input.data.one[,1]), as.matrix(input.data.two[,1]))
#overdis <- c(as.matrix(1 / input.data.one[,2]), as.matrix(1 / input.data.two[,2]))
#zi <- (1 - c(as.matrix(input.data.one[,3]),as.matrix(input.data.two[,3]))) * 100

#zi.an <- quantile(zi, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
#od.an <- quantile(overdis, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
mean.an <- quantile(mean, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)

od.an <- c(NA, NA, NA)
zi.an <- c(NA, NA, NA)

crashed <- as.integer(crashed)
results <- c(jobname, converged, crashed, error, achieved, mean.an[1], mean.an[2], mean.an[3], zi.an[1], zi.an[2], zi.an[3], od.an[1], od.an[2], od.an[3])

if(ml.equivalent==TRUE){
	suppressWarnings(ml.mean <- log(mean))
	#suppressWarnings(ml.overdis <- log(mean / overdis))
	#suppressWarnings(ml.zi <- log(zi / (1-zi)))

	ml.mean.an <- quantile(ml.mean, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
	ml.zi.an <- c(NA, NA, NA)
	ml.od.an <- c(NA, NA, NA)

	results <- c(results, ml.mean.an[1], ml.mean.an[2], ml.mean.an[3], ml.zi.an[1], ml.zi.an[2], ml.zi.an[3], ml.od.an[1], ml.od.an[2], ml.od.an[3])
}

setwd(save.directory)
unlink(temp.directory, recursive = TRUE)

print("DISCLAIMER:  *these results are intended for educational purposes only and should not be relied upon for real world applications*", quote=FALSE)
print ("", quote=FALSE)
print(paste("--- Completed dataset '", jobname, "' ---", sep=""), quote=FALSE)
print ("", quote=FALSE)

return(results)

}


##########################################################################################################################

#######  GP MODEL:

##########################################################################################################################


GP.model <- function (data, name, jags = "jags", burnin = 5000, updates = 10000, log.prior = TRUE, ml.equivalent = FALSE)
{

bugcheck <- FALSE

datana <- data

length(datana) <- 1

if(is.na(name)==TRUE){
	print("This job was not given a name", quote=FALSE)
	return(c("Error", "Name"))
}
if(is.na(datana)==TRUE){
	print("No data specified", quote=FALSE)
	return(c("Error", "Data"))
}


real.runs <- as.integer(updates)
ini.runs <- as.integer(burnin)
jobname <- name


if(!require(lattice)){
	stop("The required library 'lattice' is not installed")
}
if(!require(coda)){
	stop("The required library 'coda' is not installed")
}
if(!require(bayesmix)){
	stop("The required library 'bayesmix' is not installed")
}


test <- haveJAGS(jags)
if(test==FALSE){
	print(paste("Unable to call JAGS using '", jags, "'", sep=""))
	return(c("Error", "JAGS call"))
}

save.directory <- getwd()
on.exit(setwd(save.directory))

temp.directory <- new_unique("tempfiles")
if((temp.directory=="Directory not writable")==TRUE){
	print("Directory not writable", quote=FALSE)
	return(c("Error", "Write permissions"))
}

dir.create(temp.directory)
setwd(temp.directory)

if(ml.equivalent==FALSE){
	errorpadding <- c("NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA")
}else{
	errorpadding <- c("NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA")
}

print ("", quote=FALSE)
print(paste("--- Running GP model for dataset '", jobname, "' at ", updates, " iterations ---", sep=""), quote=FALSE)
print ("", quote=FALSE)


counts <- data

N <- length(counts)


##### Write model, data and script file

if(log.prior==TRUE){
modelstring <- paste("model {

for(row in 1 : N){
Count[row] ~ dpois(lambda[row]);
lambda[row] <- mean * gamma[row];

gamma[row] ~ dgamma(wb, wb);

}

prob <- 1;
wb <- mean * b;
b <- exp(logb);

# Priors
logb ~ dunif(-7,7);
mean ~ dunif(0.001,1000);
}")
}else{
modelstring <- paste("model {

for(row in 1 : N){
Count[row] ~ dpois(lambda[row]);
lambda[row] <- mean * gamma[row];

gamma[row] ~ dgamma(wb, wb);

}

prob <- 1;
wb <- mean * b;

# Priors
b ~ dunif(0.001,1000);
mean ~ dunif(0.001,1000);
}")
}
output <- file("model.txt", 'w')
cat(modelstring, file=output,sep="")  
close(output)

datastring <- paste("\"N\" <- ", N, "\n\"Count\" <- c(", sep="")
for(i in 1:(N-1)){
	datastring <- paste(datastring, counts[i], ", ", sep="")
}
datastring <- paste(datastring, counts[N], ")\n\n", sep="")

output <- file("data.txt", 'w')
cat(datastring, file=output,sep="")  
close(output)


scriptstring <- paste("model in <\"model.txt\">
data in <\"data.txt\">
compile
parameters in <\"inits.txt\">

initialize
update <", ini.runs, ">
monitor set <mean>
monitor set <b>
monitor set <prob>
update <", real.runs, ">
coda <*>
exit
\n", sep="")

output <- file("script.cmd", 'w')
cat(scriptstring, file=output,sep="")  
close(output)

totalupdates <- 3 * real.runs
crashed <- FALSE
converged <- NA
error <- 0
achieved <- NA

dataformean <- data
#dataformean[dataformean==0] <- NA
meanofnonzeros <- as.integer(mean(na.omit(dataformean)))

smallestmean <- as.integer(meanofnonzeros / 10)
largestmean <- as.integer(meanofnonzeros * 10)

if(is.na(smallestmean)==TRUE){
	smallestmean <- 1
}
if(is.na(largestmean)==TRUE){
	largestmean <- 10
}

if((smallestmean < 1)==TRUE){
	smallestmean <- 1
}
if((largestmean < 10)==TRUE){
	largestmean <- 10
}
if((smallestmean > 20)==TRUE){
	smallestmean <- 20
}
if((largestmean > 200)==TRUE){
	largestmean <- 200
}


#############  Write first init files and call JAGS to run simulation, then input results


probpos <- counts
zeros <- probpos==0
probpos[] <- 1
probpos[zeros] <- 0

if(log.prior==TRUE){
initsone <- paste("\"mean\"  <-  ", smallestmean, "\n\"logb\" <- -2.3\n\"gamma\"  <-  c(", sep="")
}else{
initsone <- paste("\"mean\"  <-  ", smallestmean, "\n\"b\" <- 0.1\n\"gamma\"  <-  c(", sep="")
}
for(i in 1:(N-1)){
	initsone <- paste(initsone, "1", ", ", sep="")
}
initsone <- paste(initsone, "1", ")\n\n", sep="")

output <- file("inits.txt", 'w')
cat(initsone, file=output,sep="")
close(output)


print("Running the first simulation", quote=FALSE)
JAGScall("script", jags, quiet=TRUE)


print("Loading first set of coda files", quote=FALSE)
suppressWarnings(inputonesuccess <- try(input.data.one <- read.coda("jags.out", "jags.ind", quiet=TRUE), silent=TRUE))
if((class(inputonesuccess)=="try-error")==FALSE){
	if(length(input.data.one)!=totalupdates){
		crashed <- TRUE
		print(paste("First simulation crashed at ", as.numeric(length(input.data.one))/3, " iterations", sep=""), quote=FALSE)
	}
}else{
	setwd(save.directory)
	unlink(temp.directory, recursive = TRUE)
	print("Unable to load coda files", quote=FALSE)
	error <- 1
	crashed <- as.integer(crashed)
	results <- c(jobname, converged, crashed, error, achieved, errorpadding)
	return(results)
}



#############  Write second init files and call JAGS to run simulation, then input results


probpos <- counts
zeros <- probpos==0
probpos[] <- 1
probpos[zeros] <- 1

if(log.prior==TRUE){
initsone <- paste("\"mean\"  <-  ", largestmean, "\n\"logb\" <- 2.3\n\"gamma\"  <-  c(", sep="")
}else{
initsone <- paste("\"mean\"  <-  ", largestmean, "\n\"b\" <- 10\n\"gamma\"  <-  c(", sep="")
}
for(i in 1:(N-1)){
	initsone <- paste(initsone, "1", ", ", sep="")
}
initsone <- paste(initsone, "1", ")\n\n", sep="")

output <- file("inits.txt", 'w')
cat(initsone, file=output,sep="")
close(output)


print("Running the second simulation", quote=FALSE)
JAGScall("script", jags, quiet=TRUE)


print("Loading second set of coda files", quote=FALSE)
suppressWarnings(inputtwosuccess <- try(input.data.two <- read.coda("jags.out", "jags.ind", quiet=TRUE), silent=TRUE))
if((class(inputtwosuccess)=="try-error")==FALSE){
	if(length(input.data.two)!=totalupdates){
		crashed <- TRUE
		print(paste("Second simulation crashed at ", as.numeric(length(input.data.two))/3, " iterations", sep=""), quote=FALSE)
	}
}else{
	setwd(save.directory)
	unlink(temp.directory, recursive = TRUE)
	crashed <- as.integer(crashed)
	error <- 1
	results <- c(jobname, converged, crashed, error, achieved, errorpadding)
	print("Unable to load coda files", quote=FALSE)
	return(results)
}

if(crashed==TRUE){
	if(length(input.data.one) > length(input.data.two)){
		print("Chain one was shortened to match the length of chain two", quote=FALSE)
		new.data <- input.data.two
		new.data[,1] <- input.data.one[1:(length(input.data.two)/3),1]
		new.data[,2] <- input.data.one[1:(length(input.data.two)/3),2]
		new.data[,3] <- input.data.one[1:(length(input.data.two)/3),3]
		input.data.one <- as.mcmc(new.data)
	}
	if(length(input.data.two) > length(input.data.one)){
		print("Chain two was shortened to match the length of chain one", quote=FALSE)
		new.data <- input.data.one
		new.data[,1] <- input.data.two[1:(length(input.data.one)/3),1]
		new.data[,2] <- input.data.two[1:(length(input.data.one)/3),2]
		new.data[,3] <- input.data.two[1:(length(input.data.one)/3),3]
		input.data.two <- as.mcmc(new.data)
	}
}

achieved <- length(input.data.one) / 3

if(length(input.data.one)!=length(input.data.two)){
	setwd(save.directory)
	unlink(temp.directory, recursive = TRUE)
	print("Chains were returned as different lengths", quote=FALSE)
	crashed <- as.integer(crashed)
	error <- 1
	results <- c(jobname, converged, crashed, error, achieved, errorpadding)
	return(results)
}

if(updates > 999){
	if(length(input.data.one) < 3000){
		setwd(save.directory)
		unlink(temp.directory, recursive = TRUE)
		crashed <- as.integer(crashed)
		print("The model crashed before 1000 sampled iterations", quote=FALSE)
		results <- c(jobname, converged, crashed, error, achieved, errorpadding)
		return(results)
	}
}else{
	if(length(input.data.one) < (3 * updates)){
		setwd(save.directory)
		unlink(temp.directory, recursive = TRUE)
		message <- paste("The model crashed before ", updates, " sampled iterations", sep="")
		print(message, quote=FALSE)
		crashed <- as.integer(crashed)
		results <- c(jobname, converged, crashed, error, achieved, errorpadding)
		return(results)
	}
}

achieved <- length(input.data.one) / 3

print("Assessing convergence", quote=FALSE)

mcmclist <- try(input.data <- mcmc.list(input.data.one, input.data.two), silent = TRUE)

if((class(mcmclist)=="try-error")==TRUE){
	setwd(save.directory)
	unlink(temp.directory, recursive = TRUE)
	if(bugcheck==TRUE){
		assign("input.data.one", input.data.one, pos=".GlobalEnv")
		assign("input.data.two", input.data.two, pos=".GlobalEnv")
	}
	print("There was an error combining the chains", quote=FALSE)
	crashed <- as.integer(crashed)
	error <- 1
	results <- c(jobname, converged, crashed, error, achieved, errorpadding)
	return(results)
}

one <- as.mcmc(input.data[1])
two <- as.mcmc(input.data[2])

gelmans <- matrix(ncol=2, nrow=3)

parameternames <- c("mean", "b", "prob")

unconverged <- ""
n.unconv <- 0
converged <- matrix(NA, nrow=3)

for (i in 1:3){
	thing <- mcmc.list(as.mcmc(one[,i]), as.mcmc(two[,i]))
	gelmans[i,1] <- gelman.diag(thing)$psrf[1]
	gelmans[i,2] <- gelman.diag(thing)$psrf[2]
	if(is.na(gelmans[i,1])==FALSE){
		if(gelmans[i,1] > 1.05){
			unconverged <- paste(unconverged, parameternames[i], ", ", sep="")
			n.unconv <- n.unconv + 1
			converged[i] <- "No"
		}else{
			converged[i] <- "Yes"
		}
	}else{
		converged[i] <- "Yes"
	}
}

if(n.unconv > 0){
	setwd(save.directory)
	unlink(temp.directory, recursive = TRUE)
	print(paste("The following parameters failed to converge:  ", unconverged, sep=""), quote=FALSE)
	converged <- 0
	
}else{
	converged <- 1
}


###############  Analyse coda files

print("Calculating results", quote=FALSE)
print("", quote=FALSE)
mean <- c(as.matrix(input.data.one[,1]), as.matrix(input.data.two[,1]))
overdis <- c(as.matrix(1 / input.data.one[,2]), as.matrix(1 / input.data.two[,2]))
#zi <- (1 - c(as.matrix(input.data.one[,3]),as.matrix(input.data.two[,3]))) * 100

#zi.an <- quantile(zi, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
od.an <- quantile(overdis, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
mean.an <- quantile(mean, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)

zi.an <- c(NA, NA, NA)

crashed <- as.integer(crashed)
results <- c(jobname, converged, crashed, error, achieved, mean.an[1], mean.an[2], mean.an[3], zi.an[1], zi.an[2], zi.an[3], od.an[1], od.an[2], od.an[3])

if(ml.equivalent==TRUE){
	suppressWarnings(ml.mean <- log(mean))
	suppressWarnings(ml.overdis <- log(mean / overdis))
	#suppressWarnings(ml.zi <- log(zi / (1-zi)))

	ml.mean.an <- quantile(ml.mean, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
	ml.zi.an <- c(NA, NA, NA)
	ml.od.an <- quantile(ml.overdis, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)

	results <- c(results, ml.mean.an[1], ml.mean.an[2], ml.mean.an[3], ml.zi.an[1], ml.zi.an[2], ml.zi.an[3], ml.od.an[1], ml.od.an[2], ml.od.an[3])
}

setwd(save.directory)
unlink(temp.directory, recursive = TRUE)

if(log.prior==TRUE){
if((od.an[1] < 0.002)==TRUE){
if((od.an[3] > 10)==TRUE){
test <- try(gp.largeod <- gp.largeod + 1, silent=TRUE)
if((class(test)=="try-error")==FALSE){
	assign("gp.largeod", gp.largeod, pos=".GlobalEnv")
}
print("*WARNING*  The 95% confidence interval for overdispersion is very large, indicating a lack of information about this parameter in the data.  The model should be re-run with a uniform prior distribution for overdispersion and the two sets of results compared to ensure that the results for all parameters are reliable", quote=FALSE)
print("", quote=FALSE)
}
}
}else{
if((od.an[1] < 0.002)==TRUE){
if((od.an[3] > 0.02)==TRUE){
test <- try(gp.largeod <- gp.largeod + 1, silent=TRUE)
if((class(test)=="try-error")==FALSE){
	assign("gp.largeod", gp.largeod, pos=".GlobalEnv")
}
print("*WARNING*  The 95% confidence interval for overdispersion is very large, indicating a lack of information about this parameter in the data.  The model should be re-run with a log uniform prior distribution for overdispersion and the two sets of results compared to ensure that the results for all parameters are reliable", quote=FALSE)
print("", quote=FALSE)
}
}
}
print("DISCLAIMER:  *these results are intended for educational purposes only and should not be relied upon for real world applications*", quote=FALSE)
print ("", quote=FALSE)

return(results)

}


##########################################################################################################################

#######  ZIP MODEL:

##########################################################################################################################


ZIP.model <- function (data, name, jags = "jags", burnin = 5000, updates = 10000, ml.equivalent = FALSE, adjust.mean = FALSE)
{

bugcheck <- FALSE

datana <- data

length(datana) <- 1

if(is.na(name)==TRUE){
	print("This job was not given a name", quote=FALSE)
	return(c("Error", "Name"))
}
if(is.na(datana)==TRUE){
	print("No data specified", quote=FALSE)
	return(c("Error", "Data"))
}


real.runs <- as.integer(updates)
ini.runs <- as.integer(burnin)
jobname <- name


if(!require(lattice)){
	stop("The required library 'lattice' is not installed")
}
if(!require(coda)){
	stop("The required library 'coda' is not installed")
}
if(!require(bayesmix)){
	stop("The required library 'bayesmix' is not installed")
}


test <- haveJAGS(jags)
if(test==FALSE){
	print(paste("Unable to call JAGS using '", jags, "'", sep=""))
	return(c("Error", "JAGS call"))
}

save.directory <- getwd()
on.exit(setwd(save.directory))

temp.directory <- new_unique("tempfiles")
if((temp.directory=="Directory not writable")==TRUE){
	print("Directory not writable", quote=FALSE)
	return(c("Error", "Write permissions"))
}

dir.create(temp.directory)
setwd(temp.directory)

if(ml.equivalent==FALSE){
	errorpadding <- c("NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA")
}else{
	errorpadding <- c("NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA")
}

print ("", quote=FALSE)
print(paste("--- Running ZIP model for dataset '", jobname, "' at ", updates, " iterations ---", sep=""), quote=FALSE)
print ("", quote=FALSE)


counts <- data

N <- length(counts)


##### Write model, data and script file


modelstring <- paste("model {

for(row in 1 : N){
Count[row] ~ dpois(lambda[row]);
lambda[row] <- probpos[row] * mean;
probpos[row] ~ dbern(prob);

}

b <- 1;

# Priors
prob ~ dunif(0,1);
mean ~ dunif(0.001,1000);
}")

output <- file("model.txt", 'w')
cat(modelstring, file=output,sep="")  
close(output)

datastring <- paste("\"N\" <- ", N, "\n\"Count\" <- c(", sep="")
for(i in 1:(N-1)){
	datastring <- paste(datastring, counts[i], ", ", sep="")
}
datastring <- paste(datastring, counts[N], ")\n\n", sep="")

output <- file("data.txt", 'w')
cat(datastring, file=output,sep="")  
close(output)


scriptstring <- paste("model in <\"model.txt\">
data in <\"data.txt\">
compile
parameters in <\"inits.txt\">

initialize
update <", ini.runs, ">
monitor set <mean>
monitor set <b>
monitor set <prob>
update <", real.runs, ">
coda <*>
exit
\n", sep="")

output <- file("script.cmd", 'w')
cat(scriptstring, file=output,sep="")  
close(output)

totalupdates <- 3 * real.runs
crashed <- FALSE
converged <- NA
error <- 0
achieved <- NA

dataformean <- data
dataformean[dataformean==0] <- NA
meanofnonzeros <- as.integer(mean(na.omit(dataformean)))

smallestmean <- as.integer(meanofnonzeros / 10)
largestmean <- as.integer(meanofnonzeros * 10)

if(is.na(smallestmean)==TRUE){
	smallestmean <- 1
}
if(is.na(largestmean)==TRUE){
	largestmean <- 10
}

if((smallestmean < 1)==TRUE){
	smallestmean <- 1
}
if((largestmean < 10)==TRUE){
	largestmean <- 10
}
if((smallestmean > 20)==TRUE){
	smallestmean <- 20
}
if((largestmean > 200)==TRUE){
	largestmean <- 200
}


#############  Write first init files and call JAGS to run simulation, then input results


probpos <- counts
zeros <- probpos==0
probpos[] <- 1
probpos[zeros] <- 0

initsone <- paste("\"prob\"  <-  0.05\n\"mean\"  <-  ", smallestmean, "\n\"probpos\"  <-  c(", sep="")
for(i in 1:(N-1)){
	initsone <- paste(initsone, probpos[i], ", ", sep="")
}
initsone <- paste(initsone, probpos[N], ")\n\n", sep="")

output <- file("inits.txt", 'w')
cat(initsone, file=output,sep="")
close(output)


print("Running the first simulation", quote=FALSE)
JAGScall("script", jags, quiet=TRUE)


print("Loading first set of coda files", quote=FALSE)
suppressWarnings(inputonesuccess <- try(input.data.one <- read.coda("jags.out", "jags.ind", quiet=TRUE), silent=TRUE))
if((class(inputonesuccess)=="try-error")==FALSE){
	if(length(input.data.one)!=totalupdates){
		crashed <- TRUE
		print(paste("First simulation crashed at ", as.numeric(length(input.data.one))/3, " iterations", sep=""), quote=FALSE)
	}
}else{
	setwd(save.directory)
	unlink(temp.directory, recursive = TRUE)
	print("Unable to load coda files", quote=FALSE)
	error <- 1
	crashed <- as.integer(crashed)
	results <- c(jobname, converged, crashed, error, achieved, errorpadding)
	return(results)
}



#############  Write second init files and call JAGS to run simulation, then input results


probpos <- counts
zeros <- probpos==0
probpos[] <- 1
probpos[zeros] <- 1


initsone <- paste("\"prob\"  <-  0.95\n\"mean\"  <-  ", largestmean, "\n\"probpos\"  <-  c(", sep="")
for(i in 1:(N-1)){
	initsone <- paste(initsone, probpos[i], ", ", sep="")
}
initsone <- paste(initsone, probpos[N], ")\n\n", sep="")

output <- file("inits.txt", 'w')
cat(initsone, file=output,sep="")
close(output)


print("Running the second simulation", quote=FALSE)
JAGScall("script", jags, quiet=TRUE)


print("Loading second set of coda files", quote=FALSE)
suppressWarnings(inputtwosuccess <- try(input.data.two <- read.coda("jags.out", "jags.ind", quiet=TRUE), silent=TRUE))
if((class(inputtwosuccess)=="try-error")==FALSE){
	if(length(input.data.two)!=totalupdates){
		crashed <- TRUE
		print(paste("Second simulation crashed at ", as.numeric(length(input.data.two))/3, " iterations", sep=""), quote=FALSE)
	}
}else{
	setwd(save.directory)
	unlink(temp.directory, recursive = TRUE)
	crashed <- as.integer(crashed)
	error <- 1
	results <- c(jobname, converged, crashed, error, achieved, errorpadding)
	print("Unable to load coda files", quote=FALSE)
	return(results)
}

if(crashed==TRUE){
	if(length(input.data.one) > length(input.data.two)){
		print("Chain one was shortened to match the length of chain two", quote=FALSE)
		new.data <- input.data.two
		new.data[,1] <- input.data.one[1:(length(input.data.two)/3),1]
		new.data[,2] <- input.data.one[1:(length(input.data.two)/3),2]
		new.data[,3] <- input.data.one[1:(length(input.data.two)/3),3]
		input.data.one <- as.mcmc(new.data)
	}
	if(length(input.data.two) > length(input.data.one)){
		print("Chain two was shortened to match the length of chain one", quote=FALSE)
		new.data <- input.data.one
		new.data[,1] <- input.data.two[1:(length(input.data.one)/3),1]
		new.data[,2] <- input.data.two[1:(length(input.data.one)/3),2]
		new.data[,3] <- input.data.two[1:(length(input.data.one)/3),3]
		input.data.two <- as.mcmc(new.data)
	}
}

achieved <- length(input.data.one) / 3

if(length(input.data.one)!=length(input.data.two)){
	setwd(save.directory)
	unlink(temp.directory, recursive = TRUE)
	print("Chains were returned as different lengths", quote=FALSE)
	crashed <- as.integer(crashed)
	error <- 1
	results <- c(jobname, converged, crashed, error, achieved, errorpadding)
	return(results)
}

if(updates > 999){
	if(length(input.data.one) < 3000){
		setwd(save.directory)
		unlink(temp.directory, recursive = TRUE)
		crashed <- as.integer(crashed)
		print("The model crashed before 1000 sampled iterations", quote=FALSE)
		results <- c(jobname, converged, crashed, error, achieved, errorpadding)
		return(results)
	}
}else{
	if(length(input.data.one) < (3 * updates)){
		setwd(save.directory)
		unlink(temp.directory, recursive = TRUE)
		message <- paste("The model crashed before ", updates, " sampled iterations", sep="")
		print(message, quote=FALSE)
		crashed <- as.integer(crashed)
		results <- c(jobname, converged, crashed, error, achieved, errorpadding)
		return(results)
	}
}

achieved <- length(input.data.one) / 3

print("Assessing convergence", quote=FALSE)

mcmclist <- try(input.data <- mcmc.list(input.data.one, input.data.two), silent = TRUE)

if((class(mcmclist)=="try-error")==TRUE){
	setwd(save.directory)
	unlink(temp.directory, recursive = TRUE)
	if(bugcheck==TRUE){
		assign("input.data.one", input.data.one, pos=".GlobalEnv")
		assign("input.data.two", input.data.two, pos=".GlobalEnv")
	}
	print("There was an error combining the chains", quote=FALSE)
	crashed <- as.integer(crashed)
	error <- 1
	results <- c(jobname, converged, crashed, error, achieved, errorpadding)
	return(results)
}

one <- as.mcmc(input.data[1])
two <- as.mcmc(input.data[2])

gelmans <- matrix(ncol=2, nrow=3)

parameternames <- c("mean", "b", "prob")

unconverged <- ""
n.unconv <- 0
converged <- matrix(NA, nrow=3)

for (i in 1:3){
	thing <- mcmc.list(as.mcmc(one[,i]), as.mcmc(two[,i]))
	gelmans[i,1] <- gelman.diag(thing)$psrf[1]
	gelmans[i,2] <- gelman.diag(thing)$psrf[2]
	if(is.na(gelmans[i,1])==FALSE){
		if(gelmans[i,1] > 1.05){
			unconverged <- paste(unconverged, parameternames[i], ", ", sep="")
			n.unconv <- n.unconv + 1
			converged[i] <- "No"
		}else{
			converged[i] <- "Yes"
		}
	}else{
		converged[i] <- "Yes"
	}
}

if(n.unconv > 0){
	setwd(save.directory)
	unlink(temp.directory, recursive = TRUE)
	print(paste("The following parameters failed to converge:  ", unconverged, sep=""), quote=FALSE)
	converged <- 0
	
}else{
	converged <- 1
}


###############  Analyse coda files

print("Calculating results", quote=FALSE)
print("", quote=FALSE)
mean <- c(as.matrix(input.data.one[,1]), as.matrix(input.data.two[,1]))
#overdis <- c(as.matrix(1 / input.data.one[,2]), as.matrix(1 / input.data.two[,2]))
prob <- (c(as.matrix(input.data.one[,3]),as.matrix(input.data.two[,3])))
zi <- 1 - prob

if(adjust.mean==TRUE){
	mean <- mean * prob
}

zi.an <- quantile((zi*100), probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
#od.an <- quantile(overdis, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
mean.an <- quantile(mean, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)

od.an <- c(NA, NA, NA)

crashed <- as.integer(crashed)
results <- c(jobname, converged, crashed, error, achieved, mean.an[1], mean.an[2], mean.an[3], zi.an[1], zi.an[2], zi.an[3], od.an[1], od.an[2], od.an[3])

if(ml.equivalent==TRUE){
	suppressWarnings(ml.mean <- log(mean))
	#suppressWarnings(ml.overdis <- log(mean / overdis))
	suppressWarnings(ml.zi <- log(zi / (1-zi)))

	ml.mean.an <- quantile(ml.mean, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
	ml.zi.an <- quantile(ml.zi, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
	ml.od.an <- c(NA, NA, NA)

	results <- c(results, ml.mean.an[1], ml.mean.an[2], ml.mean.an[3], ml.zi.an[1], ml.zi.an[2], ml.zi.an[3], ml.od.an[1], ml.od.an[2], ml.od.an[3])
}

setwd(save.directory)
unlink(temp.directory, recursive = TRUE)

print("DISCLAIMER:  *these results are intended for educational purposes only and should not be relied upon for real world applications*", quote=FALSE)
print ("", quote=FALSE)
print(paste("--- Completed dataset '", jobname, "' ---", sep=""), quote=FALSE)
print ("", quote=FALSE)

return(results)

}


##########################################################################################################################

#######  ZIGP MODEL:

##########################################################################################################################


ZIGP.model <- function (data, name, jags = "jags", burnin = 5000, updates = 10000, log.prior = TRUE, ml.equivalent = FALSE, adjust.mean = FALSE)
{

bugcheck <- FALSE

datana <- data

length(datana) <- 1

if(is.na(name)==TRUE){
	print("This job was not given a name", quote=FALSE)
	return(c("Error", "Name"))
}
if(is.na(datana)==TRUE){
	print("No data specified", quote=FALSE)
	return(c("Error", "Data"))
}


real.runs <- as.integer(updates)
ini.runs <- as.integer(burnin)
jobname <- name


if(!require(lattice)){
	stop("The required library 'lattice' is not installed")
}
if(!require(coda)){
	stop("The required library 'coda' is not installed")
}
if(!require(bayesmix)){
	stop("The required library 'bayesmix' is not installed")
}


test <- haveJAGS(jags)
if(test==FALSE){
	print(paste("Unable to call JAGS using '", jags, "'", sep=""))
	return(c("Error", "JAGS call"))
}

save.directory <- getwd()
on.exit(setwd(save.directory))

temp.directory <- new_unique("tempfiles")
if((temp.directory=="Directory not writable")==TRUE){
	print("Directory not writable", quote=FALSE)
	return(c("Error", "Write permissions"))
}

dir.create(temp.directory)
setwd(temp.directory)

if(ml.equivalent==FALSE){
	errorpadding <- c("NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA")
}else{
	errorpadding <- c("NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA")
}


print ("", quote=FALSE)
print(paste("--- Running ZIGP model for dataset '", jobname, "' at ", updates, " iterations ---", sep=""), quote=FALSE)
print ("", quote=FALSE)


counts <- data

N <- length(counts)


##### Write model, data and script file

if(log.prior==TRUE){
modelstring <- paste("model {

for(row in 1 : N){
Count[row] ~ dpois(lambda[row]);
lambda[row] <- probpos[row] * mean * gamma[row];
probpos[row] ~ dbern(prob);
gamma[row] ~ dgamma(wb, wb);

}

wb <- mean * b;
b <- exp(logb);

# Priors
logb ~ dunif(-7,7);
prob ~ dunif(0,1);
mean ~ dunif(0.001,1000);
}")
}else{
modelstring <- paste("model {

for(row in 1 : N){
Count[row] ~ dpois(lambda[row]);
lambda[row] <- probpos[row] * mean * gamma[row];
probpos[row] ~ dbern(prob);
gamma[row] ~ dgamma(wb, wb);

}

wb <- mean * b;

# Priors
b ~ dunif(0.001,1000);
prob ~ dunif(0,1);
mean ~ dunif(0.001,1000);
}")
}
output <- file("model.txt", 'w')
cat(modelstring, file=output,sep="")  
close(output)

datastring <- paste("\"N\" <- ", N, "\n\"Count\" <- c(", sep="")
for(i in 1:(N-1)){
	datastring <- paste(datastring, counts[i], ", ", sep="")
}
datastring <- paste(datastring, counts[N], ")\n\n", sep="")

output <- file("data.txt", 'w')
cat(datastring, file=output,sep="")  
close(output)


scriptstring <- paste("model in <\"model.txt\">
data in <\"data.txt\">
compile
parameters in <\"inits.txt\">

initialize
update <", ini.runs, ">
monitor set <mean>
monitor set <b>
monitor set <prob>
update <", real.runs, ">
coda <*>
exit
\n", sep="")

output <- file("script.cmd", 'w')
cat(scriptstring, file=output,sep="")  
close(output)

totalupdates <- 3 * real.runs
crashed <- FALSE
converged <- NA
error <- 0
achieved <- NA

dataformean <- data
dataformean[dataformean==0] <- NA
meanofnonzeros <- as.integer(mean(na.omit(dataformean)))

smallestmean <- as.integer(meanofnonzeros / 10)
largestmean <- as.integer(meanofnonzeros * 10)

if(is.na(smallestmean)==TRUE){
	smallestmean <- 1
}
if(is.na(largestmean)==TRUE){
	largestmean <- 10
}

if((smallestmean < 1)==TRUE){
	smallestmean <- 1
}
if((largestmean < 10)==TRUE){
	largestmean <- 10
}
if((smallestmean > 20)==TRUE){
	smallestmean <- 20
}
if((largestmean > 200)==TRUE){
	largestmean <- 200
}


#############  Write first init files and call JAGS to run simulation, then input results


probpos <- counts
zeros <- probpos==0
probpos[] <- 1
probpos[zeros] <- 0
if(log.prior==TRUE){
initsone <- paste("\"prob\"  <-  0.05\n\"mean\"  <-  ", smallestmean, "\n\"logb\" <- -2.3\n\"probpos\"  <-  c(", sep="")
}else{
initsone <- paste("\"prob\"  <-  0.05\n\"mean\"  <-  ", smallestmean, "\n\"b\" <- 0.1\n\"probpos\"  <-  c(", sep="")
}

for(i in 1:(N-1)){
	initsone <- paste(initsone, probpos[i], ", ", sep="")
}
initsone <- paste(initsone, probpos[N], ")\n\"gamma\"  <-  c(", sep="")
for(i in 1:(N-1)){
	initsone <- paste(initsone, "1", ", ", sep="")
}
initsone <- paste(initsone, "1", ")\n\n", sep="")

output <- file("inits.txt", 'w')
cat(initsone, file=output,sep="")
close(output)


print("Running the first simulation", quote=FALSE)
JAGScall("script", jags, quiet=TRUE)


print("Loading first set of coda files", quote=FALSE)
suppressWarnings(inputonesuccess <- try(input.data.one <- read.coda("jags.out", "jags.ind", quiet=TRUE), silent=TRUE))
if((class(inputonesuccess)=="try-error")==FALSE){
	if(length(input.data.one)!=totalupdates){
		crashed <- TRUE
		print(paste("First simulation crashed at ", as.numeric(length(input.data.one))/3, " iterations", sep=""), quote=FALSE)
	}
}else{
	setwd(save.directory)
	unlink(temp.directory, recursive = TRUE)
	print("Unable to load coda files", quote=FALSE)
	error <- 1
	crashed <- as.integer(crashed)
	results <- c(jobname, converged, crashed, error, achieved, errorpadding)
	return(results)
}



#############  Write second init files and call JAGS to run simulation, then input results


probpos <- counts
zeros <- probpos==0
probpos[] <- 1
probpos[zeros] <- 1

if(log.prior==TRUE){
initsone <- paste("\"prob\"  <-  0.95\n\"mean\"  <-  ", largestmean, "\n\"logb\" <- 2.3\n\"probpos\"  <-  c(", sep="")
}else{
initsone <- paste("\"prob\"  <-  0.95\n\"mean\"  <-  ", largestmean, "\n\"b\" <- 10\n\"probpos\"  <-  c(", sep="")
}

for(i in 1:(N-1)){
	initsone <- paste(initsone, probpos[i], ", ", sep="")
}
initsone <- paste(initsone, probpos[N], ")\n\"gamma\"  <-  c(", sep="")
for(i in 1:(N-1)){
	initsone <- paste(initsone, "1", ", ", sep="")
}
initsone <- paste(initsone, "1", ")\n\n", sep="")

output <- file("inits.txt", 'w')
cat(initsone, file=output,sep="")
close(output)


print("Running the second simulation", quote=FALSE)
JAGScall("script", jags, quiet=TRUE)


print("Loading second set of coda files", quote=FALSE)
suppressWarnings(inputtwosuccess <- try(input.data.two <- read.coda("jags.out", "jags.ind", quiet=TRUE), silent=TRUE))
if((class(inputtwosuccess)=="try-error")==FALSE){
	if(length(input.data.two)!=totalupdates){
		crashed <- TRUE
		print(paste("Second simulation crashed at ", as.numeric(length(input.data.two))/3, " iterations", sep=""), quote=FALSE)
	}
}else{
	setwd(save.directory)
	unlink(temp.directory, recursive = TRUE)
	crashed <- as.integer(crashed)
	error <- 1
	results <- c(jobname, converged, crashed, error, achieved, errorpadding)
	print("Unable to load coda files", quote=FALSE)
	return(results)
}

if(crashed==TRUE){
	if(length(input.data.one) > length(input.data.two)){
		print("Chain one was shortened to match the length of chain two", quote=FALSE)
		new.data <- input.data.two
		new.data[,1] <- input.data.one[1:(length(input.data.two)/3),1]
		new.data[,2] <- input.data.one[1:(length(input.data.two)/3),2]
		new.data[,3] <- input.data.one[1:(length(input.data.two)/3),3]
		input.data.one <- as.mcmc(new.data)
	}
	if(length(input.data.two) > length(input.data.one)){
		print("Chain two was shortened to match the length of chain one", quote=FALSE)
		new.data <- input.data.one
		new.data[,1] <- input.data.two[1:(length(input.data.one)/3),1]
		new.data[,2] <- input.data.two[1:(length(input.data.one)/3),2]
		new.data[,3] <- input.data.two[1:(length(input.data.one)/3),3]
		input.data.two <- as.mcmc(new.data)
	}
}

achieved <- length(input.data.one) / 3

if(length(input.data.one)!=length(input.data.two)){
	setwd(save.directory)
	unlink(temp.directory, recursive = TRUE)
	print("Chains were returned as different lengths", quote=FALSE)
	crashed <- as.integer(crashed)
	error <- 1
	results <- c(jobname, converged, crashed, error, achieved, errorpadding)
	return(results)
}

if(updates > 999){
	if(length(input.data.one) < 3000){
		setwd(save.directory)
		unlink(temp.directory, recursive = TRUE)
		crashed <- as.integer(crashed)
		print("The model crashed before 1000 sampled iterations", quote=FALSE)
		results <- c(jobname, converged, crashed, error, achieved, errorpadding)
		return(results)
	}
}else{
	if(length(input.data.one) < (3 * updates)){
		setwd(save.directory)
		unlink(temp.directory, recursive = TRUE)
		message <- paste("The model crashed before ", updates, " sampled iterations", sep="")
		print(message, quote=FALSE)
		crashed <- as.integer(crashed)
		results <- c(jobname, converged, crashed, error, achieved, errorpadding)
		return(results)
	}
}

achieved <- length(input.data.one) / 3

print("Assessing convergence", quote=FALSE)

mcmclist <- try(input.data <- mcmc.list(input.data.one, input.data.two), silent = TRUE)

if((class(mcmclist)=="try-error")==TRUE){
	setwd(save.directory)
	unlink(temp.directory, recursive = TRUE)
	if(bugcheck==TRUE){
		assign("input.data.one", input.data.one, pos=".GlobalEnv")
		assign("input.data.two", input.data.two, pos=".GlobalEnv")
	}
	print("There was an error combining the chains", quote=FALSE)
	crashed <- as.integer(crashed)
	error <- 1
	results <- c(jobname, converged, crashed, error, achieved, errorpadding)
	return(results)
}

one <- as.mcmc(input.data[1])
two <- as.mcmc(input.data[2])

gelmans <- matrix(ncol=2, nrow=3)

parameternames <- c("mean", "b", "prob")

unconverged <- ""
n.unconv <- 0
converged <- matrix(NA, nrow=3)

for (i in 1:3){
	thing <- mcmc.list(as.mcmc(one[,i]), as.mcmc(two[,i]))
	gelmans[i,1] <- gelman.diag(thing)$psrf[1]
	gelmans[i,2] <- gelman.diag(thing)$psrf[2]
	if(is.na(gelmans[i,1])==FALSE){
		if(gelmans[i,1] > 1.05){
			unconverged <- paste(unconverged, parameternames[i], ", ", sep="")
			n.unconv <- n.unconv + 1
			converged[i] <- "No"
		}else{
			converged[i] <- "Yes"
		}
	}else{
		converged[i] <- "Yes"
	}
}

if(n.unconv > 0){
	setwd(save.directory)
	unlink(temp.directory, recursive = TRUE)
	print(paste("The following parameters failed to converge:  ", unconverged, sep=""), quote=FALSE)
	converged <- 0
	
}else{
	converged <- 1
}


###############  Analyse coda files

print("Calculating results", quote=FALSE)
print("", quote=FALSE)

mean <- c(as.matrix(input.data.one[,1]), as.matrix(input.data.two[,1]))
overdis <- c(as.matrix(1 / input.data.one[,2]), as.matrix(1 / input.data.two[,2]))
prob <- (c(as.matrix(input.data.one[,3]),as.matrix(input.data.two[,3])))
zi <- 1 - prob

if(adjust.mean==TRUE){
	mean <- mean * prob
}

zi.an <- quantile((zi*100), probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
od.an <- quantile(overdis, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
mean.an <- quantile(mean, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)

crashed <- as.integer(crashed)
results <- c(jobname, converged, crashed, error, achieved, mean.an[1], mean.an[2], mean.an[3], zi.an[1], zi.an[2], zi.an[3], od.an[1], od.an[2], od.an[3])

if(ml.equivalent==TRUE){
	suppressWarnings(ml.mean <- log(mean))
	suppressWarnings(ml.overdis <- log(mean / overdis))
	suppressWarnings(ml.zi <- log(zi / (1-zi)))

	ml.mean.an <- quantile(ml.mean, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
	ml.zi.an <- quantile(ml.zi, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
	ml.od.an <- quantile(ml.overdis, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)

	results <- c(results, ml.mean.an[1], ml.mean.an[2], ml.mean.an[3], ml.zi.an[1], ml.zi.an[2], ml.zi.an[3], ml.od.an[1], ml.od.an[2], ml.od.an[3])
}


setwd(save.directory)
unlink(temp.directory, recursive = TRUE)

if(log.prior==TRUE){
if((od.an[1] < 0.002)==TRUE){
if((od.an[3] > 10)==TRUE){
test <- try(zigp.largeod <- zigp.largeod + 1, silent=TRUE)
if((class(test)=="try-error")==FALSE){
	assign("zigp.largeod", zigp.largeod, pos=".GlobalEnv")
}
print("*WARNING*  The 95% confidence interval for overdispersion is very large, indicating a lack of information about this parameter in the data.  The model should be re-run with a uniform prior distribution for overdispersion and the two sets of results compared to ensure that the results for all parameters are reliable", quote=FALSE)
print("", quote=FALSE)
}
}
}else{
if((od.an[1] < 0.002)==TRUE){
if((od.an[3] > 0.02)==TRUE){
test <- try(zigp.largeod <- zigp.largeod + 1, silent=TRUE)
if((class(test)=="try-error")==FALSE){
	assign("zigp.largeod", zigp.largeod, pos=".GlobalEnv")
}
print("*WARNING*  The 95% confidence interval for overdispersion is very large, indicating a lack of information about this parameter in the data.  The model should be re-run with a log uniform prior distribution for overdispersion and the two sets of results compared to ensure that the results for all parameters are reliable", quote=FALSE)
print("", quote=FALSE)
}
}
}
print("DISCLAIMER:  *these results are intended for educational purposes only and should not be relied upon for real world applications*", quote=FALSE)
print ("", quote=FALSE)
print(paste("--- Completed dataset '", jobname, "' ---", sep=""), quote=FALSE)
print ("", quote=FALSE)

return(results)

}

##########################################################################################################################

zigp.model <- ZIGP.model
zip.model <- ZIP.model
gp.model <- GP.model
p.model <- P.model
BAYESCOUNT <- bayescount