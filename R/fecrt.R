fecrt <- function(name = NA, pre.data = NA, post.data = NA, data = list(pre=pre.data, post=post.data), animal.names = FALSE, efficacy=95, restrict.efficacy = TRUE, zero.inflation = FALSE, divide.data = 1, record.chains = FALSE, write.file = FALSE, bootstrap.iters=10000, plot.graph = TRUE, skip.mcmc = FALSE, individual.analysis=TRUE, ...){

runname <- name

if(!individual.analysis & zero.inflation){
	warning("Individual analysis must be used for the zero inflated model.  individual.analysis has been set to TRUE")
	individual.analysis <- TRUE
}
speed <- !individual.analysis

div <- divide.data  # Need both to prevent list data being divided twice

if(!require(runjags)){
	stop("The required library 'runjags' is not installed")
}
if(!require(lattice)){
	stop("The required library 'lattice' is not installed")
}
if(!require(coda)){
	stop("The required library 'coda' is not installed")
}

passthrough <- list(...)
if(is.null(passthrough$max.time)) passthrough$max.time <- "1hr"
if(is.null(passthrough$interactive)) passthrough$interactive <- FALSE
if(is.null(passthrough$plots)) passthrough$plots <- FALSE

arguments <- formals(autorun.jags)

for(i in 1:length(passthrough)){
	if(is.null(arguments[[names(passthrough)[i]]])){
		warning(paste("'", names(passthrough)[i], "' is not an argument to autorun.jags and was ignored"))
	}else{
		arguments[[names(passthrough)[i]]] <- passthrough[[i]]
	}
}

jags <- eval(arguments$jags)
silent.jags <- arguments$silent.jags

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

if(class(data)=="function") stop("The class of the object supplied for data is a function")

if(class(data)=="list"){
	if(all(is.na(data[[1]])) & all(is.na(data[[2]]))){
		data <- NA
	}else{
		if(all(is.na(data[[1]])) | all(is.na(data[[2]]))) stop("Both pre and post treatment counts must be supplied")
		pre.data <- data[[1]]
		post.data <- data[[2]]

		f <- function(data){
if(class(data)=="list"){
	
	if(suppressWarnings(class(data$totals) != "NULL" & class(data$repeats) != "NULL")){
		# data is a list of counts and repeats
		totals <- data$totals
		repeats <- data$repeats
		
		if(length(totals)!=length(repeats)) stop("The number of repeats does not match the length of totals")
		
		if(any(is.na(totals)) | any(is.na(repeats))) stop("Missing data is not allowed in repeats or totals")
		
		totals <- totals / div
		
		data <- matrix(NA, ncol=max(repeats), nrow=length(totals))
		
		for(i in 1:length(totals)){
				
			if(round(totals[i])!=totals[i]) stop("All totals must be an integer")
			if(round(repeats[i])!=repeats[i] | repeats[i] < 1) stop("All repeats must be an integer greater than 0")
			
			if(repeats[i] > 1){

				done <- FALSE
	
				while(!done){
					data[i,1:repeats[i]] <- rpois(repeats[i], totals[i]/repeats[i])
					if(sum(data[i,],na.rm=TRUE)==totals[i] & all(na.omit(data[i,]) >= 0)) done <- TRUE
		
				}
			}else{
				data[i,1] <- totals[i]
			}
		}
		
	}else{

		delist <- function(list){
		
		N <- length(list)
		R <- numeric(length=N)
	
		for(i in 1:N){
			R[i] <- length(na.omit(as.numeric(list[[i]])))
		
		}
	
		if(sum(R>0,na.rm=TRUE)==0) stop("No real values in the data")
	
		data <- matrix(NA, ncol=max(R), nrow=N)
	
		for(i in 1:N){
			data[i,(min(R[i], 1):R[i])] <- na.omit(as.numeric(list[[i]]))
		}
	
		data <- data[apply(!is.na(data), 1, sum)>0,]

		if(sum(R>0)==1) return(matrix(data, nrow=1)) else return(data)
		}
	
		data <- delist(data)	
	}
}
return(data)
}
		pre.data <- f(pre.data)
		post.data <- f(post.data)
		
		div <- 1
		
		data <- list(pre.counts=as.matrix(pre.data), post.counts=as.matrix(post.data))
	}
}

if(class(data)=="numeric" | class(data)=="integer") stop("Data cannot be provided as a single vector of data")

if(class(data)=="array"){
	if(dim(data)[3] == 1){
		data <- matrix(data, ncol=(dim(data)[2]))
	}else{
		if(dim(data)[3] !=2){
			stop("Data cannot be provided as an array of 3rd dimension length greater than 2")
		}else{
			data <- list(pre.counts=as.matrix(data[,,1]), post.counts=as.matrix(data[,,2]))
		}
	}
}
	
cat("\n--- FECRT: Analyse faecal egg cout reduction test data using a Bayesian distributional simulation model ---\n\n")
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

if(class(data)=="matrix" | class(data)=="list"){
	dataok <- TRUE
}else{
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

if(class(data)=="matrix") if(length(data[1,]) !=(suppressWarnings(sum(2, as.numeric(animal.names[1]), na.rm=TRUE)))) stop("If data is provided as a matrix, there must be either 2 columns if animal.names is FALSE, or 3 columns if animal.names is TRUE")



if(class(data)!="matrix" & animal.names[1]==TRUE){
	cat("Warning:  'animal.names' cannot be TRUE if the data is not provided as a matrix.  'animal.names' will be set to FALSE\n")
	animal.names <- FALSE
}

setnames <- animal.names

namesdone <- FALSE

if(setnames==TRUE){
	setnames <- data[,1]
	data <- matrix(data[2:length(data[1,]),], ncol=2)
	namesdone <- TRUE
}

if(class(data) == "matrix") data <- list(pre.counts=as.matrix(data[,1]), post.counts=as.matrix(data[,2]))

if(length(data$pre.counts[,1])!=length(data$post.counts[,1])) stop("There were an unequal number of animals in the pre and post treatment data")

if(class(data$pre.counts) == "NULL" | class(data$post.counts) == "NULL") stop("An error occured while transforming the data to an appropriate format")


names <- NA

if(namesdone==FALSE){
	if(length(setnames) != 1 | setnames[1] != FALSE){
		if(length(setnames) != length(data$pre.counts[,1])) stop("The length of the character vector animal.names does not match the lenth of the data provided")
		names <- setnames
	}else{
		names <- paste("Animal ", 1:length(data$pre.counts[,1]), sep="")
	}
}

N <- length(data$pre.counts)

if(N > 100 & !skip.mcmc){
	if(arguments$interactive){
		if(!speed){
			cat("There are a large number of animals in the data which will take a long time to analyse using MCMC, and may cause memory issues with individual analysis.  ")
			returned <- ask("Please choose from the following options:\n1 Continue with analysis (this may cause a crash)\n2 Set individual.analysis to FALSE\n3 Skip the MCMC analysis\n:", "numeric", c(1,3))
			cat("\n")
			if(returned==2){
				individual.analysis <- FALSE
				speed <- TRUE
			}
			if(returned==3) skip.mcmc <- TRUE
		}else{
			cat("There are a large number of animals in the data which will take a long time to analyse using MCMC.  ")
			returned <- ask("Please choose from the following options:\n1 Continue with analysis (this may take a while)\n2 Skip the MCMC analysis\n:", "numeric", c(1,2))
			cat("\n")
			if(returned==2) skip.mcmc <- TRUE
		}
	}else{
		if(!speed){
			individual.analysis <- FALSE
			speed <- TRUE
			warning("There are a large number of animals in the data which will take a long time to analyse using MCMC, and may cause memory issues with individual analysis.  individual analysis was therefore not performed.")
		}
	}
}

cat("Assessing the faecal egg count reduction for ", N, " animals.  This will take some time...\n", sep="")


pre <- data$pre.counts / div
post <- data$post.counts / div


#if(sum(pre+post,na.rm=TRUE)==0) stop("The data contains all 0 counts")

model <- paste("model{
	
	for(row in 1:N){
		for(repeat in 1:Pre.R[N]){
			Pre[row,repeat] ~ dpois(xpre.lambda[row])
		}
		for(repeat in 1:Post.R[N]){
			Post[row,repeat] ~ dpois(xpost.lambda[row])
		}
		
		xpre.lambda[row] <- ", if(zero.inflation) "probpos[row] * ", "pre.mean * pre.gamma[row]
		xpost.lambda[row] <- ", if(zero.inflation) "probpos[row] * ", "post.mean * post.gamma[row]
		
		ind.delta[row] <- xpost.lambda[row] / xpre.lambda[row]
		
		pre.gamma[row] ~ dgamma(pre.disp, pre.disp)
		post.gamma[row] ~ dgamma(post.disp, post.disp)
		
", if(zero.inflation) "probpos[row] ~ dbern(prob)","
	}
	
	pre.disp <- 1 / ia
	
	post.mean <- pre.mean * delta.mean
	post.disp <- pre.disp * delta.disp

	
	# Priors
	pre.mean ~ dunif(0.1, 1000)
", if(zero.inflation) "prob ~ dunif(0,1)","
	ia <- exp(logia)
	logia ~ ", {if(FALSE) "dunif(-9.21,5)#4.6" else "dunif(-4,4);"}, "
	
	ddispl <- ", {if(FALSE) exp(-9.21) else 0.02}, " / pre.disp;
	ddispu <- ", {if(FALSE) exp(5) else 100}, " / pre.disp;
	
	delta.mean ~ ", {if(restrict.efficacy) "dbeta(1,1)" else "dunif(0, 10)"}, "
	delta.disp ~ dlnorm(0, 0.1)T(ddispl,ddispu)
	
	}", sep="")


probposinit <- as.numeric(pmax(pre,post,na.rm=TRUE) > 0)
probinit <- max(sum(probposinit) / length(probposinit), 0.1)

inits1 <- dump.format(list(probpos=probposinit, pre.mean=max(mean(pre,na.rm=TRUE)*2, 1), delta.mean=0.01, logia=log(0.1), delta.disp=2, pre.gamma=replicate(length(pre), 1), post.gamma=replicate(length(pre), 1), prob=probinit))


inits2 <- dump.format(list(probpos=replicate(length(pre), 1), pre.mean=max(mean(pre,na.rm=TRUE)*0.5, 1), delta.mean=0.99, logia=log(10), delta.disp=0.5, pre.gamma=replicate(length(pre), 1), post.gamma=replicate(length(pre), 1), prob=1))

datastring <- dump.format(list(Pre=as.matrix(pre), Post=as.matrix(post), N=length(pre), Pre.R=replicate(length(pre),1), Post.R=replicate(length(post),1)))


monitor=c("pre.mean", "delta.mean", "delta.disp")
if(zero.inflation) monitor <- c(monitor, "prob", "probpos")
if(!speed) monitor <- c(monitor, "ind.delta")#, "pre.disp", "xpre.lambda", "xpost.lambda")


if(!skip.mcmc){
	
	arguments$model <- model
	arguments$inits <- c(inits1, inits2)
	arguments$n.chains <- length(arguments$inits)
	arguments$data <- datastring
	arguments$monitor <- monitor
	arguments$silent.jags <- list(silent.jags=arguments$silent.jags, killautocorr=TRUE)
	arguments$plots <- FALSE


	class(arguments) <- "list"
	
	results <- do.call(autorun.jags, arguments, quote=FALSE)
	
	#results <- autorun.jags(data=datastring, model=model, monitor=monitor, n.chains=2, inits=c(inits1, inits2), silent.jags = list(silent.jags=silent.jags, killautocorr=TRUE), plots = FALSE, thin.sample = TRUE, interactive=interactive, max.time=max.time, ...)



if(results[1]=="Error"){
	#print(results)
	cat("An unexpected error occured during the simulation.  Ensure that the values of divide.data, pre.data and post.data provided are correct and re-run the functon.\n\n--- Simulation aborted ---\n\n")
	stop("An error occured during the simulation")
}
if(any(names(results)=="pilot.mcmc")){
	warning("The chains did not achieve convergence during the simulation, you should interpret the Bayesian results with extreme caution")
	converged <- FALSE
	results$mcmc <- results$pilot.mcmc
	results$summary <- results$pilot.summary
}else{
	converged <- TRUE
}

}else{
	results <- list(mcmc="MCMC analysis not performed")
}

pre <- pre.data
post <- post.data
#return(results)

blankquant <- quantile(0, probs=c(0.025, 0.5, 0.975))
	
if(any(class(pre)==c("array", "matrix", "list")) | any(class(post)==c("array", "matrix", "list"))){
	warning("Bootstrap and WAAVP calculations are not available for repeated pre and/or post treatment egg counts")
	boot.reductions = bootquant = method.boot = waavpquant = method.waavp <- NA
}else{
	if(length(pre)!=length(post)) warning("Pre and post treatment data are of different lengths")
	
	cat("Calculating the Bootstrap and WAAVP method analysis...\n")
	# Bootstrap:
	
	boot.pre <- matrix(data=sample(pre, length(pre)*bootstrap.iters, replace=TRUE), ncol=bootstrap.iters, nrow=length(pre))
	boot.post <- matrix(data=sample(post, length(post)*bootstrap.iters, replace=TRUE), ncol=bootstrap.iters, nrow=length(post))
	
	boot.reductions <- (1 - apply(boot.post, 2, mean) / apply(boot.pre, 2, mean)) *100

	boot.reductions[boot.reductions==-Inf] <- NA
	bootprob <- sum(boot.reductions < efficacy, na.rm=TRUE) / sum(!is.na(boot.reductions)) * 100
	
	bootquant <- quantile(boot.reductions, probs=c(0.025, 0.5, 0.975), na.rm=TRUE)
	if(is.na(bootprob)){
		method.boot <- "Returned error"
	}else{
		method.boot <- if(bootprob >= 97.5) "Confirmed resistant" else if(bootprob >= 50) "Probable resistant" else if(bootprob >= 2.5) "Possible resistant" else "Confirmed susceptible"
	}

	# WAAVP:
	
	waavpquant <- blankquant
	
	pre.data <- pre
	post.data <- post

	N <- length(pre.data)
	arith.pre <- mean(pre.data)
	arith.post <- mean(post.data)
	var.pre <- var(pre.data)
	var.post <- var(post.data)

	red <- 100*(1-arith.post/arith.pre)
	var.red <- var.post/(N*arith.post^2) + var.pre/(N*arith.pre^2)
	try(if(var.post==0 & arith.post==0) var.red <- var.pre/(N*arith.pre^2))

	upper.ci <- 100 * (1-(arith.post/arith.pre * exp(-2.048*sqrt(var.red))))
	lower.ci <- 100 * (1-(arith.post/arith.pre * exp(2.048*sqrt(var.red))))

	waavpquant[1] <- lower.ci
	waavpquant[2] <- red
	waavpquant[3] <- upper.ci

	if(any(is.na(c(red,efficacy,lower.ci)))){
		method.waavp <- "Returned error"
	}else{
	method.waavp <- "susceptible"
	if(red < 95 | lower.ci < 90){
		method.waavp <- "Suspected resistant"
		if(red < 95 & lower.ci < 90){
			method.waavp <- "Confirmed resistant"
		}
	}
}
}

if(!skip.mcmc){

cat("Calculating the Bayesian method analysis...\n")

reduction <- (1-unlist(results$mcmc[,"delta.mean"]))*100
pre.mean <- (unlist(results$mcmc[,"pre.mean"]))
deltashape <- unlist(results$mcmc[,"delta.disp"])
if(zero.inflation) zi <- (1-unlist(results$mcmc[,"prob"]))*100 else zi <- NA

vars <- nvar(results$mcmc)

n <- (vars-(length(monitor)-(1+zero.inflation)))/(1+zero.inflation)

probpos = delta <- matrix(nrow=niter(results$mcmc)*2, ncol=n)

if(!speed){
for(i in 1:n){
	
	probpos[,i] <- if(zero.inflation) unlist(results$mcmc[,(i+4)]) else NA
	delta[,i] <- (1-unlist(results$mcmc[,(i+(zero.inflation*n)+3+zero.inflation)]))*100	
	
}
}

ind.prob.inf <- apply(delta, 2, function(x) return(sum(x < efficacy) / length(x))) * 100

ind.zi.probs <- (1-apply(probpos, 2, mean, na.rm=TRUE))*100

indredquant <- matrix(ncol=3, nrow=length(pre), dimnames=list(names, names(blankquant)))
indredquan <- if(!speed) HPDinterval(as.mcmc(delta)) else matrix(ncol=2,nrow=n)
indredmedians <- if(!speed) apply(delta, 2, median) else replicate(n, NA)

indredquant[,2] <- indredmedians
indredquant[,1] <- indredquan[,1]
indredquant[,3] <- indredquan[,2]

medred <- median(reduction)
ci <- HPDinterval(as.mcmc(reduction))
l95red <- ci[1]
u95red <- ci[2]
mcmcquant <- blankquant
mcmcquant[1] <- l95red
mcmcquant[3] <- u95red
mcmcquant[2] <- medred

meanquant <- blankquant
meanquant[2] <- median(pre.mean*divide.data)
ci <- HPDinterval(as.mcmc(pre.mean*divide.data))
meanquant[1] <- ci[1]
meanquant[3] <- ci[2]

ziquant <- blankquant
ziquant[2] <- if(zero.inflation) median(zi) else NA
ci <- if(zero.inflation) HPDinterval(as.mcmc(zi)) else c(NA,NA)
ziquant[1] <- ci[1]
ziquant[3] <- ci[2]

# post.disp = pre.disp * delta.disp
# 1/pod^2 = 1/prd^2 * dd
# 1/pod^2 = 1/prd^2 * dd
# prd^2 = pod^2 * dd

dshapequant <- blankquant
dshapequant[2] <- median(deltashape)
ci <- HPDinterval(as.mcmc(deltashape))
dshapequant[1] <- ci[1]
dshapequant[3] <- ci[2]


mcmcprob <- sum(reduction < efficacy) / length(reduction) * 100

method.mcmc <- if(mcmcprob >= 97.5) "Confirmed resistant" else if(mcmcprob >= 50) "Probable resistant" else if(mcmcprob >= 2.5) "Possible resistant" else "Confirmed susceptible"

#class <- (prob.res > 0.025) + (prob.res > 0.50) + (prob.res > 0.975) + 1
#method3 <- switch(class, "1"="susceptible", "2"="possible", "3"="probable", "4"="resistant")

}else{
	method.mcmc = mcmcquant = mcmcprob = ziquant = dshapequant = meanquant = indredquant = ind.prob.inf = ind.zi.probs = converged = efficacy = restrict.efficacy <- NA
}

cat("Finished calculations\n")

output <- list(results.mcmc=method.mcmc, quant.mcmc=mcmcquant, results.boot=method.boot, quant.boot=bootquant, results.waavp=method.waavp, quant.waavp = waavpquant, prob.mcmc=mcmcprob, prob.boot=bootprob, ziquant=ziquant, dshapequant=dshapequant, meanquant=meanquant, indredquant=indredquant, ind.prob.inf=ind.prob.inf, ind.zi.prob=ind.zi.probs, converged=converged, efficacy=efficacy, restrict.efficacy=restrict.efficacy, animal.names=names, name=name)
if(record.chains) output <- c(output, list(mcmc=results$mcmc))


if(write.file){
	cat("Writing results to file\n")
	filename <- new_unique(paste(name, ".results",sep=""), ".txt", ask=arguments$interactive)
	print.fecrt.results(output, filename=filename)
	filename <- new_unique(paste(name, ".graph",sep=""), ".pdf", ask=arguments$interactive)
	pdf(file=filename)
	if(!skip.mcmc) hist(reduction, col="red", breaks="fd", main="Posterior distribution for the true reduction", freq=FALSE, ylab="Probability", xlab="True FEC reduction (%)", xlim=c(0, 100))
	abline(v=efficacy)
	dev.off()
	if(record.chains) save(results, file=new_unique(paste("fecrt.", name, sep=""), ".Rsave", ask=FALSE))
	cat("Results file written successfully\n")
}else{
	if(plot.graph){
		if(!skip.mcmc) 	hist(reduction, col="red", breaks="fd", main="Posterior distribution for the true reduction", freq=FALSE, ylab="Probability", xlab="True FEC reduction (%)", xlim=c(0, 100))
		abline(v=efficacy)
	}
}


cat("\nAnalysis complete\n\n", sep="")
cat("*PLEASE NOTE:  THIS SOFTWARE IS INTENDED FOR EDUCATIONAL PURPOSES ONLY AND SHOULD NOT BE RELIED UPON FOR REAL WORLD APPLICATIONS*\n")
cat("*ANALYSING DATA USING MCMC SAMPLING CAN PRODUCE MISLEADING RESULTS IF USED INAPPROPRIATELY*\n\n")
cat("--- End ---\n\n")

class(output) <- "fecrt.results"
return(output)


}

FECRT <- fecrt
