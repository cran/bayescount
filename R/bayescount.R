bayescount <- function (name = NA, data = NA, setnames = NA, div = 1, model = c("ZILP"), burnin = 5000, updates = c(10000,100000,500000), jags = "jags", rownames = FALSE, remove.zeros = TRUE, remove.missing = TRUE, test = TRUE, alt.prior = FALSE, write.file = TRUE, adjust.mean = FALSE, crash.retry = 1, silent.jags = FALSE, likelihood=FALSE)
{

datanames <- setnames
omit.zeros <- remove.zeros
runname <- name

updates <- sort(updates)

if(!require(lattice)){
	stop("The required library 'lattice' is not installed")
}
if(!require(coda)){
	stop("The required library 'coda' is not installed")
}

test.jags <- testjags(jags, silent=TRUE)
if(test.jags[[2]][1]==FALSE){
	cat("Unable to call JAGS using '", jags, "'\n", sep="")
	stop("Unable to call JAGS")
}


testwritable <- new_unique("test")
if(testwritable=="Directory not writable"){
	cat("\nThe working directory is not writable.  Please change the working directory\n\n")
	stop("Directory not writable")
}

if((crash.retry != as.integer(crash.retry)) | (crash.retry < 0)){
	cat("\nThe value of 'crash.retry' is not valid.  Please provide a non-negative integer\n\n")
	stop("Parameter invalid")
}

model <- switch(model, P="SP", ZIP="ZISP", model)

models <- c("SP", "ZISP", "GP", "ZIGP", "LP", "ZILP", "WP", "ZIWP", "IP")
modelsfull <- c("single Poisson", "zero-inflated single Poisson", "gamma Poisson", "zero-inflated gamma Poisson", "lognormal Poisson", "zero-inflated lognormal Poisson", "Weibull Poisson", "zero-inflated Weibull Poisson", "independant Poisson")

if(class(alt.prior)=="character"){
	stop("'alt.prior' must be either TRUE or FALSE for bayescount().  Use bayescount.single() for custom priors")
}

if(model[1]=="all") model <- models

if((length(model) == 0) | (length(model) > length(models)) | sum(is.na(model)) > 0){
	cat("Invalid model selection.  Please choose from the following models: ", sep="")
	cat(models, sep=", ")
	cat("\n")
	stop("Invalid model selection")
}

for(i in 1:length(model)){
	if(sum(model[i]==models) != 1){
		cat("Invalid model selection \'", model[i], "\'.  Please choose from the following models: ", sep="")
		cat(models, sep=", ")
		cat("\n")
		stop("Invalid model selection")
	}
}

duplicated <- ""
for(i in 1:length(models)){
	if(sum(models[i]==model) > 1){
		cat("Duplicated model selection \'", models[i], "\'.  The model will be run only once!\n", sep="")
		duplicated <- c(duplicated, models[i])
	}
}

if(length(duplicated) > 1){
	for(i in 2:length(duplicated)){
		model[model==duplicated[i]] <- ""
		model <- c(model, duplicated[i])
	}
	model <- model[model!=""]	
}

models.nos <- 1:length(models)

model.no <- numeric(length=length(model))
modelfull <- character(length=length(model))

for(i in 1:length(model)){
	model.no[i] <- models.nos[models[models.nos]==model[i]]
	modelfull[i] <- modelsfull[model.no[i]]
}

if(sum(model=="ZISP") > 0 | sum(model=="ZIGP") > 0 | sum(model=="ZILP") > 0 | sum(model=="ZIWP") > 0 | sum(model=="ZIIP") > 0){
	any.zero.inflation <- TRUE
}else{
	any.zero.inflation <- FALSE
}

data <- as.matrix(data)

cat("\n--- 'Bayescount': Analyse count data using a Bayesian distributional simulation model implemented in JAGS ---\n\n")
cat("*PLEASE NOTE:  THIS SOFTWARE IS INTENDED FOR EDUCATIONAL PURPOSES ONLY AND SHOULD NOT BE RELIED UPON FOR REAL WORLD APPLICATIONS*\n")
cat("*ANALYSING DATA USING MCMC SAMPLING CAN PRODUCE MISLEADING RESULTS IF USED INAPPROPRIATELY*\n\n")

if(is.na(runname)){
	runname <- ask(prompt = "Please enter a name for this analysis:  ", type="character")
}

datain <- data
datana <- data
datana <- as.vector(na.omit(datana))
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
		cat("ERROR:  The path you have entered does not appear to be valid\n")
	}else{
		if(valid.data==FALSE){
			cat("ERROR:  Invalid path / data\n") 
		}else{
			dataok=TRUE
		}
	}
	cat("\n")
}
while(dataok==FALSE){
	datain <- ask(prompt = "Please enter the path to a (comma delimited) CSV file containing the data (type 'exit' to quit):  ", type="character")
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
				cat("ERROR:  The path you have entered does not appear to be valid\n")
			}else{
				if(valid.data==FALSE){
					cat("ERROR:  The data you have entered is of length less than 2\n")
				}else{
					dataok=TRUE
				}
			}
		}else{
			cat("ERROR:  The path you entered does not appear to be valid; please retry (path may be absolute or relative)\n")
		}
		}else{
			cat("ERROR:  The path you entered does not appear to be valid; please retry (path may be absolute or relative)\n")
		}
	}else{
		cat("ERROR:  The path you entered does not appear to be valid; please retry (path may be absolute or relative)\n")
	}
	cat("\n")
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
	suppressWarnings(try({
		if(any(is.na(as.numeric(na.omit(data[1,]))))){
			names <- data[1,]
			data <- data[2:length(data[,1]),]
			datanames <- TRUE
		}else{
			if(class(dimnames(data)[[2]]) == "character"){
				names <- dimnames(data)[[2]]
				datanames <- TRUE
			}else{
				datanames <- FALSE
			}
		}
	}, silent=TRUE))
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

if(rownames==TRUE && remove.missing==TRUE){
	cat("Missing data cannot be removed if rownames are used.  Missing data will not be removed\n")
	remove.missing <- FALSE
}

tiers <- length(updates)

cat("Settings are as follows:")
cat(paste("\nName of analysis:                             ", runname, sep=""))
cat(paste("\nFirst row used for dataset names:             ", datanames, sep=""))
cat(paste("\nFirst column used for individual names:       ", rownames, sep=""))
cat(paste("\nNumber of datasets:                           ", as.numeric(length(names)), sep=""))
cat(paste("\nNumber of animals in column one               ", as.numeric(animalsindata1), sep=""))
cat(paste("\nDivide data by                                ", div, sep=""))
if(length(model)==1){
	cat(paste("\nModel to use:                                 ", modelfull[1], sep=""))
}else{
	cat(paste("\nModels to use:                                ", modelfull[1], " model", sep=""))
	for(i in 2:length(model)){
		cat(paste("\n                                              ", modelfull[i], " model", sep=""))
	}
}
cat(paste("\nOmit datasets with all zero counts:           ", omit.zeros, sep=""))
cat(paste("\nRemove missing data:                          ", remove.missing, sep=""))
cat(paste("\nBurn-in updates:                              ", burnin, sep=""))
for(i in 1:tiers){
	cat(paste("\nSample updates (tier ", i, "):                      ", updates[i], sep=""))
}
suppressWarnings(if(sum(model==c("GP", "ZIGP", "LP",	 "ZILP", "WP", "ZIWP"))>0){
	cat(paste("\nUse alternative prior for variance :          ", alt.prior, sep=""))
})
if(any.zero.inflation==TRUE){
	cat(paste("\nAdjust z-i means to mean of whole population: ", adjust.mean, sep=""))
}
cat(paste("\nCalculate the likelihood for each model:      ", likelihood, sep=""))
cat(paste("\nSuppress JAGS output to screen:               ", silent.jags, sep=""))
cat(paste("\nNumber of times to retry after crash:         ", as.numeric(crash.retry), sep=""))
cat(paste("\nWrite the results to file:                    ", write.file, sep=""))
cat(paste("\nSystem call to activate JAGS:                 ", jags, "\n\n", sep=""))
test.jags <- testjags(jags, silent=FALSE)

proceed <- !test
if(test==TRUE){

setdata <- as.integer(na.omit(data[,1]))
cat("\n\nTesting the model function\n")
output <- bayescount.single(burnin=1000, model = model[1], updates=1000, jags = jags, data = setdata, alt.prior = alt.prior, silent.jags = TRUE)

if(output[2]==1){
	proceed <- ask(prompt = "The test (first) dataset crashed.  Continue with other datasets?", type="logical")
}else{
	if(output[3]==0){
		cat("Test completed successfully\n")
		proceed <- TRUE
	}else{
		stop("Possible JAGS call or permissions error.  Please check that you have the latest version of JAGS (>0.99.0) installed")
	}
}

}

all.models <- model
all.modelfull <- modelfull
GP.largeod = ZIGP.largeod = WP.largeod = ZIWP.largeod = LP.largeod = ZILP.largeod <- 0

start.time <- Sys.time()

total.crashed <- 0
total.unconverged <- 0
total.error <- 0
total.largeod <- 0

cat("\n")

for(j in 1:length(all.models)){

	model <- all.models[j]
	modelfull <- all.modelfull[j]

	headers <- c("converged", "crashed", "error", "iterations")
	resultscol <- 4

	headers <- c(headers, "mean.l.95", "mean.median", "mean.u.95")
	resultscol <- resultscol + 3
	
	headers <- c(headers, "variance.l.95", "variance.median", "variance.u.95")
	resultscol <- resultscol + 3
	
	if(model=="ZIGP" | model=="ZILP" | model=="ZIWP" | model=="ZISP"){
		headers <- c(headers, "zi.l.95", "zi.median", "zi.u.95")
		resultscol <- resultscol + 3
	}
	
	if(model=="GP" | model=="ZIGP" | model=="WP" | model=="ZIWP"){
		headers <- c(headers, "scale.l.95", "scale.median", "scale.u.95")
		resultscol <- resultscol + 3
		headers <- c(headers, "shape.l.95", "shape.median", "shape.u.95")
		resultscol <- resultscol + 3
	}
	
	if(model=="LP" | model=="ZILP"){
		headers <- c(headers, "log.mean.l.95", "log.mean.median", "log.mean.u.95")
		resultscol <- resultscol + 3
		headers <- c(headers, "log.variance.l.95", "log.variance.median", "log.variance.u.95")
		resultscol <- resultscol + 3
	}
	
	headers <- c(headers, "multivariate.psrf")
	resultscol <- resultscol + 1
	
	if(likelihood==TRUE){
		headers <- c(headers, "likelihood.l.95", "likelihood.median", "likelihood.u.95")
		resultscol <- resultscol + 3
	}
	
	headers <- c(headers, "time.taken")
	resultscol <- resultscol + 1
	
	zeros <- character(length=resultscol-1)
	zeros[] <- "NA"


	results <- matrix(NA, ncol=resultscol, nrow=length(names), dimnames =list(names, headers))
	headers <- paste(paste(headers, collapse=","), "\n", sep="")
	crashed <- 0
	unconverged <- 0
	error <- 0
	largeod <- 0
	
	name <- new_unique(name=paste(runname, ".", model, sep=""), suffix=".csv", ask = test, prompt="A results file with this name already exists.  Overwrite?")
	outfile <- file(name, 'w')
	cat("dataset,", headers, file = outfile, sep = "")
	
	for(i in 1:length(names)){
		done <- FALSE
		if(remove.missing==TRUE){
			setdata <- as.integer(na.omit(data[,i]))
		}else{
			setdata <- as.integer(data[,i])
		}
		if(omit.zeros==TRUE){
			if((sum(na.omit(setdata)) < 1)==TRUE){
				output <- zeros
				cat("Dataset '", names[i], "' contained all zeros and was therefore omitted\n", sep="")
				done <- TRUE
			}
		}
		tries.left <- crash.retry + 1
		pre.time <- Sys.time()
		while(done==FALSE){
			cat("\nRunning the ", model, " model for dataset '", names[i], "'...\n", sep="")
			if(length(setdata)!=length(data[,i]) && remove.missing==TRUE){
				cat("\n*WARNING*  ", length(data[,i]) - length(setdata), " missing and/or non-numeric datapoints were removed from dataset '", names[i], "'\n", sep="")
			}
			if(sum(is.na(setdata)) > 0 && remove.missing==FALSE){
				cat("\n*WARNING*  There are ", sum(is.na(setdata)), " missing and/or non-numeric datapoints in dataset '", names[i], "'\n", sep="")
			}
			output <- bayescount.single(model=model, burnin=burnin, updates=updates, jags = jags, data = setdata, alt.prior = alt.prior, adjust.mean = adjust.mean, silent.jags=silent.jags, likelihood=likelihood)
			
			if(((output[2]=="1") | (output[3]=="1")) && (sum(na.omit(output[1] == "1")) == 0)){
				cat("The ", model, " model for dataset '", names[i], "' crashed", sep="")
				tries.left <- tries.left-1
				if(tries.left==0){
					cat("\n\n")
					done <- TRUE
				}else{
					cat(".  Re-trying:\n")
				}
			}else{
				done <- TRUE
			}
		}
		post.time <- Sys.time()
		
		if(output[2]=="1"){
			crashed <- crashed + 1
			total.crashed <- total.crashed + 1
		}
		if(output[3]=="1"){
			error <- error + 1
			total.error <- total.error + 1
		}
		if(sum(na.omit(output[1]=="0")) == 1 && (sum(na.omit(setdata)) > 0 || omit.zeros==FALSE)){
			cat("*WARNING*  The ", model, " model for dataset '", names[i], "' failed to achieve convergence\n", sep="")
			unconverged <- unconverged + 1
			total.unconverged <- total.unconverged + 1
		}
		
		if(sum(na.omit(output[1] == "1")) == 1){			
			largeod <- largeod + as.numeric(assess.variance(model=model, alt.prior=alt.prior, l.95=output["var.2.5%"], u.95=output["var.97.5%"]))
		}

		results[i,1:(resultscol-1)] <- as.numeric(output)
		results[i,resultscol] <- as.numeric(timestring(pre.time, post.time, units="secs", show.units=FALSE))
		
		time.taken <- timestring(pre.time, post.time, show.units=TRUE)
		total.time <- timestring(start.time, post.time, show.units=TRUE)
		number.models <- length(all.models) * length(names)
		current.number <- ((j-1) * length(names)) + i
		percent.complete <- current.number / number.models
		time.remaining <- timestring((as.integer(difftime(post.time, start.time, units="secs")) / percent.complete) - as.integer(difftime(post.time, start.time, units="secs")), show.units=TRUE)
		
		if(sum(na.omit(setdata)) > 0 || omit.zeros==FALSE){
			cat("Dataset '", names[i], "' completed in ", time.taken, " with the ", model, " model\n", sep="")
		}
		cat(i+(length(names)*(j-1)), " of ", number.models, " datasets (", round(percent.complete, digits=2)*100, "%) completed\n", sep="")
		cat(total.time, " elapsed, estimated time remaining: ", time.remaining, "\n", sep="") 
		if(j > 1){
			cat(unconverged, " failed convergence, ", crashed, " crashed, and ", error, " quit with an error with the ", model, " model\n", sep="")
			cat(total.unconverged, " failed convergence, ", total.crashed, " crashed, and ", total.error, " quit with an error for all models\n", sep="")
		}else{
			cat(unconverged, " failed convergence, ", crashed, " crashed, and ", error, " quit with an error\n", sep="")
		}
		cat("\n")
		cat(names[i], results[i,], file = outfile, sep = ",")
		cat("\n", file=outfile, sep="")
		gc()
	}
	close(outfile)
	if(write.file==FALSE){
		unlink(name)
	}
	assign(paste(runname, ".", model, ".results", sep=""), results, pos=".GlobalEnv")
	if((model=="GP") | (model=="ZIGP") | (model=="WP") | (model=="ZIWP")){
		assign(paste(model, ".largeod", sep=""), largeod)
	}
}

cat("All models completed.  Total time taken:  ", timestring(start.time, post.time), "\n\n", sep="")
cat("*PLEASE NOTE:  THIS SOFTWARE IS INTENDED FOR EDUCATIONAL PURPOSES ONLY AND SHOULD NOT BE RELIED UPON FOR REAL WORLD APPLICATIONS*\n")
cat("*ANALYSING DATA USING MCMC SAMPLING CAN PRODUCE MISLEADING RESULTS IF USED INAPPROPRIATELY*\n\n")
cat("--- End ---\n\n")


largeods <- c(GP.largeod, ZIGP.largeod, WP.largeod, ZIWP.largeod, LP.largeod, ZILP.largeod)
odmodels <- c("gamma Poisson", "zero-inflated gamma Poisson", "Weibull Poisson", "zero-inflated Weibull Poisson", "lognormal Poisson", "zero-inflated lognormal Poisson")

for(i in 1:length(largeods)){
	if((alt.prior==TRUE | is.character(alt.prior)) && largeods[i]>0){
		cat("*WARNING*  The 95% confidence interval for the variance was very large in a total of ", largeods[i], " datasets using the ", odmodels[i], " model, indicating a lack of information about this parameter in these datasets.  These models should be re-run with the standard prior distribution and the two sets of results compared to ensure that the results for all parameters are reliable\n")
	}
	if(alt.prior==FALSE && largeods[i]>0){
		cat("*WARNING*  The 95% confidence interval for the variance was very large in a total of ", largeods[i], " datasets using the ", odmodels[i], " model, indicating a lack of information about this parameter in these datasets.  These models should be re-run with an alternative prior distribution and the two sets of results compared to ensure that the results for all parameters are reliable\n")
	}
}

}