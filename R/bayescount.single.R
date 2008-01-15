bayescount.single <- function(data = stop("Data must be specified"), model="ZILP", burnin=5000, updates=c(10000), jags = findjags(), alt.prior = FALSE, adjust.mean = FALSE, silent.jags = FALSE, raw.output = FALSE, likelihood=FALSE){

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

updates <- sort(updates)

model <- toupper(model)

model <- switch(model, P="SP", ZIP="ZISP", model)

models <- c("SP", "ZISP", "GP", "ZIGP", "LP", "ZILP", "WP", "ZIWP", "IP")
modelsfull <- c("single Poisson", "zero-inflated single Poisson", "gamma Poisson", "zero-inflated gamma Poisson", "lognormal Poisson", "zero-inflated lognormal Poisson", "Weibull Poisson", "zero-inflated Weibull Poisson", "independant Poisson")

if((length(model) != 1) |  sum(is.na(model)) > 0 | sum(model==models)!=1){
	cat("Invalid model selection.  Please choose from ONE of the following models: ", sep="")
	cat(models, sep=", ")
	cat("\n")
	stop("Invalid model selection")
}

params.names <- switch(model, SP=c("mean"), ZISP=c("mean", "prob"), IP=c("mean", "variance"), LP=c("lmean", "lvariance"), ZILP=c("lmean", "lvariance", "prob"), GP=c("mean", "b"), ZIGP=c("mean", "b", "prob"), WP=c("a", "b"), ZIWP=c("a", "b", "prob"))

n.params <- switch(model, SP=1, ZISP=2, IP=2, LP=2, ZILP=3, GP=2, ZIGP=3, WP=2, ZIWP=3) 

total.updates <- n.params * updates

errorpad <- quantile(0, probs=c(0.025, 0.5, 0.975))
errorpad[] <- NA
errorpadding <- c(mean=errorpad[1], mean=errorpad[2], mean=errorpad[3], var=errorpad[1], var=errorpad[2], var=errorpad[3])

if(sum(strsplit(model, split="")[[1]][1:2] != c("Z", "I"))==0){
	errorpadding <- c(errorpadding, zi=errorpad[1], zi=errorpad[2], zi=errorpad[3])
}
if((model=="GP") | (model=="WP") | (model=="ZIGP") | (model=="ZIWP")){
	errorpadding <- c(errorpadding, scale=errorpad[1], scale=errorpad[2], scale=errorpad[3], shape=errorpad[1], shape=errorpad[2], shape=errorpad[3])
}
if((model=="LP") | (model=="ZILP")){
	errorpadding <- c(errorpadding, lmean=errorpad[1], lmean=errorpad[2], lmean=errorpad[3], lvar=errorpad[1], lvar=errorpad[2], lvar=errorpad[3])
}
errorpadding <- c(errorpadding, mpsrf=NA)
if(likelihood==TRUE){
	errorpadding <- c(errorpadding, likelihood=errorpad[1], likelihood=errorpad[2], likelihood=errorpad[3])
}

length(burnin) <- length(updates)
burnin[is.na(burnin)] <- 0
alreadydone <- 0

if(sum(strsplit(model, split="")[[1]][1:2] != c("Z", "I"))==0){
	zero.inflation <- TRUE
}

#if(model=="SP" | model=="ZISP"){
#	n.params.inc.lambda <- n.params #  FOR PROBPOS:  + (as.numeric(likelihood) * (as.numeric(zero.inflation) * length(data)))
#}else{
#	n.params.inc.lambda <- n.params #  NOT USING LAMBDA NOW + (as.numeric(likelihood) ) # FOR PROBPOS: * (length(data) + (as.numeric(zero.inflation) * length(data))))
#}

if(model=="IP" && likelihood==TRUE){
	#  The IP model is the only one that has something monitored for lambda
	n.params.inc.lambda <- n.params + length(data)
	params.names <- c(params.names, paste("Lambda", 1:length(data)))
}else{
	n.params.inc.lambda <- n.params
}

if(likelihood==TRUE){
	params.names <- c(params.names, "likelihood")
}	

oldtwo = oldone <- matrix(NA, ncol=n.params.inc.lambda, nrow=0)

strings <- run.model(model=model, data=data, alt.prior=alt.prior, jags=jags, silent.jags = silent.jags, call.jags=FALSE, monitor.lambda=likelihood) # Only using lambda for IP model - run.model corrects for this
modelstring <- strings[[1]]
datastring <- strings[[2]]
new.inits <- strings[[3]]
monitors <- strings[[4]]

running.string <- ""
cat("\n")

for(i in 1:length(updates)){
	
	cat("Running the ", model, " model at ", running.string, updates[i]-alreadydone, " iterations\n", sep="")
	
	success <- try(output <- run.jags(model=modelstring, burnin=burnin[i], updates=updates[i]-alreadydone, jags = jags, data = datastring, monitor=monitors, silent.jags = silent.jags, inits=new.inits, check.conv=FALSE), silent=FALSE)
	
	if(class(success)=="try-error"){
		cat("There was an error during the simulation\n")
		if(raw.output==TRUE){
			return(list(matrix(NA, ncol=n.params.inc.lambda+likelihood, nrow=1, dimnames=list(NULL, params.names)), matrix(NA, ncol=n.params.inc.lambda+likelihood, nrow=1, dimnames=list(NULL, params.names))))
		}else{
			return(c(conv=NA, crash=0, error=1, iterations=0, errorpadding))
		}
	}
	
	if(output[[1]][1]=="Unable to load coda files"){
		if(raw.output==TRUE){
			return(list(matrix(NA, ncol=n.params.inc.lambda+likelihood, nrow=1, dimnames=list(NULL, params.names)), matrix(NA, ncol=n.params.inc.lambda+likelihood, nrow=1, dimnames=list(NULL, params.names))))
		}else{
			return(c(conv=NA, crash=0, error=1, iterations=0, errorpadding))
		}
	}

	crashed <- 0
	error <- 0
	converged <- 0
	
	achieved <- length(output[[1]][,1])
	
	if(achieved != updates[i]-alreadydone){
		cat("The simulation crashed at ", as.numeric(achieved + alreadydone), " iterations", sep="")
		crashed <- 1
	}
	
	suppressWarnings(success <- try({	
	
	two = one <- matrix(NA, ncol=n.params.inc.lambda+likelihood, nrow=achieved+alreadydone, dimnames=list(NULL, params.names))
	
	for(j in 1:n.params.inc.lambda){
		
		one[,j] <- c(oldone[,j], output[[1]][,j])
		two[,j] <- c(oldtwo[,j], output[[2]][,j])
	}
			
	mcmc <- mcmc.list(as.mcmc(one[,1:n.params]), as.mcmc(two[,1:n.params]))
	new.inits <- c(output[[3]], output[[4]])

	convergence <- gelman.diag(mcmc)
	unconverged <- FALSE
	unconverged.string <- ""
	
	for(j in 1:n.params){
		param.conv <- convergence$psrf[j, 1]
		if(param.conv > 1.05){
			unconverged.string <- paste(unconverged.string, params.names[j], " (psrf=", round(param.conv, digits=2), "), ", sep="")
			unconverged <- TRUE
		}
	}
		
	if(unconverged==FALSE){
		achieved <- alreadydone + achieved
		cat("Convergence was achieved for all parameters at ", as.character(updates[i]), " iterations\n", sep="")
		converged <- 1
		break
	}
	}, silent=FALSE))
	
	if(class(success)=="try-error"){
		cat("There was an error when assessing convergence\n")
		return(c(NA, 0, 1, 0, errorpadding))
	}
	
	cat("The following parameters failed to converge at ", as.character(updates[i]), " iterations:  ", unconverged.string, "\n", sep="")	
	
	if(crashed==1){
		achieved <- alreadydone + achieved
		break
	}
	
	alreadydone <- alreadydone + achieved
	achieved <- alreadydone
	running.string <- "an additional "
	oldone <- one
	oldtwo <- two
}


if(raw.output==TRUE && likelihood==FALSE){
	if(unconverged==TRUE){
		cat("Returning UNCONVERGED results\n")
	}else{
		cat("Returning results\n")
	}
	cat("*PLEASE NOTE:  THIS SOFTWARE IS INTENDED FOR EDUCATIONAL PURPOSES ONLY AND SHOULD NOT BE RELIED UPON FOR REAL WORLD APPLICATIONS*\n\n")
	return(mcmc.list(as.mcmc(one), as.mcmc(two)))
}

if(updates[i] >= 1000){
	if(length(one[,1]) < (1000)){
		cat("The model crashed before 1000 sampled iterations\n")
		return(c(conv=!unconverged, crash=crashed, error=error, iterations=achieved, errorpadding))
	}
}else{
	if(length(one[,1]) < (updates[i])){
		cat("The model crashed before ", updates, " sampled iterations\n", sep="")
		return(c(conv=!unconverged, crash=crashed, error=error, iterations=achieved, errorpadding))
	}
}



###############  Analyse coda files

cat("Calculating results\n")

if(sum(strsplit(model, split="")[[1]][1:2] != c("Z", "I"))==0){
	zero.inflation <- TRUE
	model <- switch(model, ZISP="SP", ZILP="LP", ZIGP="GP", ZIWP="WP", ZILP="LP", ZIIP="IP")
}else{
	zero.inflation <- FALSE
}

shape = scale <- NA

lambda <- matrix(ncol=length(data), nrow=achieved*2)

if(model=="IP"){
	mean <- c(as.matrix(one[,1]), as.matrix(two[,1]))
	variance <- c(as.matrix(one[,2]), as.matrix(two[,2]))^2
	if(likelihood==TRUE){
		for(i in 1:length(data)){
			lambda[,i] <- c(as.matrix(one[,i+2+as.numeric(zero.inflation)]), as.matrix(two[,i+2+as.numeric(zero.inflation)]))
			#probpos[,i] <- c(as.matrix(one[,i+length(data)+2+as.numeric(zero.inflation)]), as.matrix(two[,i+length(data)+2+as.numeric(zero.inflation)]))
		}
	}
}

if(model=="SP"){
	mean <- c(as.matrix(one[,1]), as.matrix(two[,1]))
	variance <- 0
	#if(likelihood==TRUE){
	#	lambda <- data
	#}
}

if(model=="LP"){
	lmean <- c(as.matrix((one[,1])), as.matrix((two[,1])))
	lvariance <- c(as.matrix((one[,2])), as.matrix((two[,2])))
	
	normal.meanvar <- normal.params(lmean, sqrt(lvariance))
	mean <- normal.meanvar[,1]
	variance <- normal.meanvar[,2]^2
	#if(likelihood==TRUE){
		#for(i in 1:length(data)){
			#lambda[,i] <- c(as.matrix(one[,i+2+as.numeric(zero.inflation)]), as.matrix(two[,i+2+as.numeric(zero.inflation)]))
			#probpos[,i] <- c(as.matrix(one[,i+length(data)+2+as.numeric(zero.inflation)]), as.matrix(two[,i+length(data)+2+as.numeric(zero.inflation)]))
		#}
	#}
}

if(model=="GP"){
	mean <- c(as.matrix(one[,1]), as.matrix(two[,1]))
	scale <- c(as.matrix(1 / one[,2]), as.matrix(1 / two[,2]))
	shape <- mean / scale
	variance <- shape*scale^2
	#if(likelihood==TRUE){
	#	for(i in 1:length(data)){
	#		lambda[,i] <- c(as.matrix(one[,i+2+as.numeric(zero.inflation)]), as.matrix(two[,i+2+as.numeric(zero.inflation)])) * mean
			#probpos[,i] <- c(as.matrix(one[,i+length(data)+2+as.numeric(zero.inflation)]), as.matrix(two[,i+length(data)+2+as.numeric(zero.inflation)]))
	#	}
	#}
}

if(model=="WP"){
	#  b is different in JAGS than R, but translation done in model to make priors easier
	a <- c(as.matrix(one[,1]), as.matrix(two[,1]))
	b <- c(as.matrix(one[,2]), as.matrix(two[,2]))
	mean <- b/a * gamma(1/a)
	variance <- b^2/a * (2 * gamma(2/a) - (1 / a * gamma(1/a)^2))
	shape <- a
	scale <- b
	#if(likelihood==TRUE){
	#	for(i in 1:length(data)){
	#		lambda[,i] <- c(as.matrix(one[,i+2+as.numeric(zero.inflation)]), as.matrix(two[,i+2+as.numeric(zero.inflation)]))
			#probpos[,i] <- c(as.matrix(one[,i+length(data)+2+as.numeric(zero.inflation)]), as.matrix(two[,i+length(data)+2+as.numeric(zero.inflation)]))
	#	}
	#}
}

if(zero.inflation==TRUE & model=="SP"){
	prob <- c(as.matrix(one[,2]), as.matrix(two[,2]))
	zi <- (1 - prob) * 100
}
if(zero.inflation==TRUE & model!="SP"){
	prob <- c(as.matrix(one[,3]), as.matrix(two[,3]))
	zi <- (1- prob) * 100
}
		
if(adjust.mean==TRUE && zero.inflation==TRUE){
	mean <- mean * prob
}

try(zi.an <- quantile((zi), probs = c(0.025, 0.5, 0.975), na.rm = TRUE), silent=TRUE)
try(scale.an <- quantile(scale, probs = c(0.025, 0.5, 0.975), na.rm = TRUE), silent=TRUE)
try(shape.an <- quantile(shape, probs = c(0.025, 0.5, 0.975), na.rm = TRUE), silent=TRUE)
try(mean.an <- quantile(mean, probs = c(0.025, 0.5, 0.975), na.rm = TRUE), silent=TRUE)
try(variance.an <- quantile(variance, probs = c(0.025, 0.5, 0.975), na.rm = TRUE), silent=TRUE)
try(lmean.an <- quantile(lmean, probs = c(0.025, 0.5, 0.975), na.rm = TRUE), silent=TRUE)
try(lvariance.an <- quantile(lvariance, probs = c(0.025, 0.5, 0.975), na.rm = TRUE), silent=TRUE)

results <- c(conv=as.integer(converged), crash=as.integer(crashed), error=as.integer(error), iterations=as.integer(achieved), mean=mean.an[1], mean=mean.an[2], mean=mean.an[3], var=variance.an[1], var=variance.an[2], var=variance.an[3])

if(zero.inflation==TRUE){
	results <- c(results, zi=zi.an[1], zi=zi.an[2], zi=zi.an[3])
}

if((model=="GP") | (model=="WP")){
	results <- c(results, scale=scale.an[1], scale=scale.an[2], scale=scale.an[3], shape=shape.an[1], shape=shape.an[2], shape=shape.an[3])
}
if(model=="LP"){
	results <- c(results, lmean=lmean.an[1], lmean=lmean.an[2], lmean=lmean.an[3], lvar=lvariance.an[1], lvar=lvariance.an[2], lvar=lvariance.an[3])
}

if(n.params==1) mpsrf <- convergence$psrf[1,1] else mpsrf <- convergence$mpsrf

results <- c(results, mpsrf=mpsrf)

if(likelihood==TRUE){
	success <- try({
	if(zero.inflation==FALSE){
		l.zi <- NA
	}else{
		l.zi <- zi
	}
		
	if(model=="WP" | model=="GP"){
		likeli <- likelihood(model=model, data=data, shape=shape, scale=scale, zi=l.zi, silent=paste("noziwarn", silent.jags, sep=""), log=TRUE, raw.output=TRUE)  # Leave iterations as default
	}
	if(model=="LP"){
		newmeanvar <- normal.params(lmean, sqrt(lvariance))
		likemean <- newmeanvar[,1]
		likevar <- newmeanvar[,2]^2
		likeli <- likelihood(model=model, data=data, mean=likemean, variance=likevar, zi=l.zi, silent=paste("noziwarn", silent.jags, sep=""), log=TRUE, raw.output=TRUE)  # Leave iterations as default
		
	}
	if(model=="SP"){
		likeli <- likelihood(model=model, data=data, mean=mean, zi=l.zi, silent=paste("noziwarn", silent.jags, sep=""), log=TRUE, raw.output=TRUE)  # Leave iterations as default
	}
	if(model=="IP"){
		likeli <- likelihood(model=model, data=data, mean=lambda, zi=l.zi, silent=paste("noziwarn", silent.jags, sep=""), log=TRUE, raw.output=TRUE)
	}}, silent=FALSE)
	
	if(class(success)=="try-error"){
		cat("An error occured while computing the likelihood\n")
		likeli <- list(NA, NA)
	}
	#print(sum(is.na(likeli[[1]])))
	#assign('like', likeli, pos=.GlobalEnv)
	if(any(is.na(likeli[[1]]))){
		likeli.an <- quantile(NA, probs=c(0.025, 0.5, 0.975), na.rm=TRUE)
	}else{
		likeli.an <- quantile(likeli[[1]], probs=c(0.025, 0.5, 0.975))
	}
	results <- c(results, likelihood=likeli.an[1], likelihood=likeli.an[2], likelihood=likeli.an[3])
}

if(raw.output==TRUE && likelihood==TRUE){
	
	success <- try({
	
	one.backup <- one
	two.backup <- two
	
	sequence <- numeric(length=(length(one[,1]) + length(two[,1])))
	
	sequence[] <- NA
	
	sequence[likeli[[2]]] <- likeli[[1]]
		
	one[,length(one[1,])] <- sequence[1:length(one[,1])]
	two[,length(two[1,])] <- sequence[(length(one[,1])+1):(length(one[,1])+length(two[,1]))]
	
	}, silent=FALSE)
	
	if(class(success)=="try-error"){
		cat("An error occured while combining the chains with the calculated likelihoods.  The likelihoods will not be returned\n")
		one <- one.backup
		two <- two.backup
	}
	
	if(unconverged==TRUE){
		cat("Returning UNCONVERGED results\n")
	}else{
		cat("Returning results\n")
	}
	cat("*PLEASE NOTE:  THIS SOFTWARE IS INTENDED FOR EDUCATIONAL PURPOSES ONLY AND SHOULD NOT BE RELIED UPON FOR REAL WORLD APPLICATIONS*\n\n")
	
	return(mcmc.list(as.mcmc(one), as.mcmc(two)))
}

largeod <- assess.variance(model=model, alt.prior=alt.prior, l.95 = variance.an[1], u.95 = variance.an[3])

if((alt.prior==TRUE | is.character(alt.prior)) && largeod==TRUE){
	cat("*WARNING*  The 95% confidence interval for the variance is very large, indicating a lack of information about this parameter in the data.  The model should be re-run with the standard prior distribution and the two sets of results compared to ensure that the results for all parameters are reliable\n")
}
if(alt.prior==FALSE && largeod==TRUE){
	cat("*WARNING*  The 95% confidence interval for the variance is very large, indicating a lack of information about this parameter in the data.  The model should be re-run with an alternative prior distribution and the two sets of results compared to ensure that the results for all parameters are reliable\n")
}


cat("*PLEASE NOTE:  THIS SOFTWARE IS INTENDED FOR EDUCATIONAL PURPOSES ONLY AND SHOULD NOT BE RELIED UPON FOR REAL WORLD APPLICATIONS*\n")
cat("Finished running the model\n\n", sep="")

results[] <- as.numeric(results[])

return(results)

}	