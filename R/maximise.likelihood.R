maximise.likelihood <- function(data=stop("Data must be specified"), model=stop("Please specify a distribution"), mean=NA, variance=NA, zi=NA, shape=NA, scale=NA, silent=FALSE){
	
	model <- toupper(model)
	testdata <- data
	
	models <- c("P", "ZIP", "G", "ZIG", "L", "ZIL", "W", "ZIW", "GP", "ZIGP", "LP", "ZILP", "WP", "ZIWP")
		
	if((length(model) != 1) |  sum(is.na(model)) > 0 | sum(model==models)!=1){
		if(silent==FALSE){
			cat("Invalid model selection.  Please choose from ONE of the following distributions: ", sep="")
			cat(models, sep=", ")
			cat("\n")
		}
		stop("Invalid model selection")
	}
	
	if(silent==FALSE) cat("\n")
	
	if(sum(is.na(data)) > 0){
		if(silent==FALSE) cat("WARNING:  Missing data was removed before calculating the likelihood\n\n")
		data <- na.omit(data)
	}
	
	guitest <- testJAGS(silent=TRUE)
	
	#if(guitest$os=='windows' | (guitest$R.GUI == "AQUA" & guitest$R.package.type == "mac.binary")) eol <- "" else eol <- "                  \n"
	if(guitest$R.GUI == "AQUA" & guitest$R.package.type == "mac.binary"){
		clearline <- ""
		eol <- "\n"
	}else{
		clearline <- "\r                                                         \r"
		eol <- "\t\t"
	}
	
	if(silent==FALSE){
		cat("Maximising the likelihood for the '", model, "' model.  This may take some time\n")
	}
	
	if(model=="P"){
		if(is.na(mean)) mean <- mean(data)
		
		if(silent==TRUE){
			f <- function(mean) return(likelihood(model=model, data=data, mean=mean, silent=TRUE))
		}else{
			cat("\n\tmean\n", eol, "\t", signif(mean+(mean/100000), 4), eol, sep="")
			f <- function(mean){
				cat(clearline, "\t", signif(mean+(mean/100000), 4), eol, sep="")
				return(likelihood(model=model, data=data, mean=mean, silent=TRUE))
			}
		}
		maxim <- optimise(f, mean, lower=0, upper=max(data)*10, maximum=TRUE)
		results <- c(mean=maxim$maximum)
	}

	if(model=="ZIP"){
		if(is.na(mean)) mean <- mean(data[data!=0])
		if(is.na(zi)) zi <- sum(data==0)/length(data)*100
		
		if(silent==TRUE){
			f <- function(pars=c(NA, NA)){
				if(pars[2] < 0 | pars[2] > 100 | pars[1] < 0) return(-Inf)
				return(likelihood(model=model, data=data, mean=pars[1], zi=pars[2], silent=TRUE))
			}
		}else{
			cat("\n\tmean\t\tzi\n", eol, "\t", signif(mean+(mean/100000), 4), "\t\t", signif(zi+(zi/100000), 4), eol, sep="")
			f <- function(pars=c(NA, NA)){
				if(pars[2] < 0 | pars[2] > 100 | pars[1] < 0) return(-Inf)
				cat(clearline, "\t", signif(pars[1]+(pars[1]/100000), 4), "\t\t", signif(pars[2]+(pars[2]/100000), 4), eol, sep="")
				return(likelihood(model=model, data=data, mean=pars[1], zi=pars[2], silent=TRUE))
			}
		}
		maxim <- optim(c(mean, zi), f, control=list(fnscale=-1))
		results <- c(mean=maxim$par[1], zi=maxim$par[2])
	}
	
	if(model=="G" | model=="GP"){
		if(is.na(scale)) scale <- var(data)/mean(data)
		if(is.na(shape)) shape <- mean(data) / scale
		
		if(silent==TRUE){
			f <- function(pars=c(NA, NA)){
				if(pars[2] <= 0 | pars[1] <= 0) return(-Inf)
				return(likelihood(model=model, data=data, shape=pars[1], scale=pars[2], silent=TRUE))
			}
		}else{
			cat("\n\tshape\t\tscale\n", eol, "\t", signif(shape, 4), "\t\t", signif(scale, 4), eol, sep="")
			f <- function(pars=c(NA, NA)){
				if(pars[2] <= 0 | pars[1] <= 0) return(-Inf)
				cat(clearline, "\t", signif(pars[1], 4), "\t\t", signif(pars[2], 4), eol, sep="")
				return(likelihood(model=model, data=data, shape=pars[1], scale=pars[2], silent=TRUE))
			}
		}
		maxim <- optim(c(shape, scale), f, control=list(fnscale=-1))
		results <- c(shape=maxim$par[1], scale=maxim$par[2])
	}
	
	if(model=="ZIG" | model=="ZIGP"){
		if(is.na(scale)) scale <- var(data[data!=0])/mean(data[data!=0])
		if(is.na(shape)) shape <- mean(data[data!=0]) / scale
		if(is.na(zi)) zi <- sum(data==0)/length(data)*100
		
		if(silent==TRUE){
			f <- function(pars=c(NA, NA, NA)){
				if(pars[2] <= 0 | pars[1] <= 0 | pars[3] < 0 | pars[3] > 100) return(-Inf)
				return(likelihood(model=model, data=data, shape=pars[1], scale=pars[2], zi=pars[3], silent=TRUE))
			}
		}else{
			cat("\n\tshape\t\tscale\t\tzi\n", eol, "\t", signif(shape, 4), "\t\t", signif(scale, 4), "\t\t", signif(zi, 4), eol, sep="")
			f <- function(pars=c(NA, NA, NA)){
				if(pars[2] <= 0 | pars[1] <= 0 | pars[3] < 0 | pars[3] > 100) return(-Inf)
				cat(clearline, "\t", signif(pars[1], 4), "\t\t", signif(pars[2], 4), "\t\t", signif(pars[3], 4), eol, sep="")
				return(likelihood(model=model, data=data, shape=pars[1], scale=pars[2], zi=pars[3], silent=TRUE))
			}
		}
		maxim <- optim(c(shape, scale, zi), f, control=list(fnscale=-1))
		results <- c(shape=maxim$par[1], scale=maxim$par[2], zi=maxim$par[3])
	}
	
	if(model=="W" | model=="WP"){
		if(is.na(scale)) scale <- 4
		if(is.na(shape)) shape <- mean(data)
		
		if(silent==TRUE){
			f <- function(pars=c(NA, NA)){
				if(pars[2] <= 0 | pars[1] <= 0) return(-Inf)
				return(likelihood(model=model, data=data, shape=pars[1], scale=pars[2], silent=TRUE))
			}
		}else{
			cat("\n\tshape\t\tscale\n", eol, "\t", signif(shape, 4), "\t\t", signif(scale, 4), eol, sep="")
			f <- function(pars=c(NA, NA)){
				if(pars[2] <= 0 | pars[1] <= 0) return(-Inf)
				cat(clearline, "\t", signif(pars[1], 4), "\t\t", signif(pars[2], 4), eol, sep="")
				return(likelihood(model=model, data=data, shape=pars[1], scale=pars[2], silent=TRUE))
			}
		}
		maxim <- optim(c(shape, scale), f, control=list(fnscale=-1))
		results <- c(shape=maxim$par[1], scale=maxim$par[2])
	}
	
	if(model=="ZIW" | model=="ZIWP"){
		if(is.na(scale)) scale <- 4
		if(is.na(shape)) shape <- mean(data[data!=0])
		if(is.na(zi)) zi <- sum(data==0)/length(data)*100
		
		if(silent==TRUE){
			f <- function(pars=c(NA, NA, NA)){
				if(pars[2] <= 0 | pars[1] <= 0 | pars[3] < 0 | pars[3] > 100) return(-Inf)
				return(likelihood(model=model, data=data, shape=pars[1], scale=pars[2], zi=pars[3], silent=TRUE))
			}
		}else{
			cat("\n\tshape\t\tscale\t\tzi\n", eol, "\t", signif(shape, 4), "\t\t", signif(scale, 4), "\t\t", signif(zi, 4), eol, sep="")
			f <- function(pars=c(NA, NA, NA)){
				if(pars[2] <= 0 | pars[1] <= 0 | pars[3] < 0 | pars[3] > 100) return(-Inf)
				cat(clearline, "\t", signif(pars[1], 4), "\t\t", signif(pars[2], 4), "\t\t", signif(pars[3], 4), eol, sep="")
				return(likelihood(model=model, data=data, shape=pars[1], scale=pars[2], zi=pars[3], silent=TRUE))
			}
		}
		maxim <- optim(c(shape, scale, zi), f, control=list(fnscale=-1))
		results <- c(shape=maxim$par[1], scale=maxim$par[2], zi=maxim$par[3])
	}
	
	if(model=="L" | model=="LP"){
		if(is.na(mean)) mean <- mean(data)
		if(is.na(variance)) variance <- var(data)
		
		if(silent==TRUE){
			f <- function(pars=c(NA, NA)){
				if(pars[2] < 0 | pars[1] < 0) return(-Inf)
				return(likelihood(model=model, data=data, mean=pars[1], variance=pars[2], silent=TRUE))
			}
		}else{
			cat("\n\tmean\t\tvariance\n", eol, "\t", signif(mean, 4), "\t\t", signif(variance, 4), eol, sep="")
			f <- function(pars=c(NA, NA)){
				if(pars[2] <= 0 | pars[1] < 0) return(-Inf)
				cat(clearline, "\t", signif(pars[1], 4), "\t\t", signif(pars[2], 4), eol, sep="")
				return(likelihood(model=model, data=data, mean=pars[1], variance=pars[2], silent=TRUE))
			}
		}
		maxim <- optim(c(mean, variance), f, control=list(fnscale=-1))
		results <- c(mean=maxim$par[1], variance=maxim$par[2])
	}
	
	if(model=="ZIL" | model=="ZILP"){
		if(is.na(mean)) mean <- mean(data[data!=0])
		if(is.na(variance)) variance <- var(data[data!=0])
		if(is.na(zi)) zi <- sum(data==0)/length(data)*100
		
		if(silent==TRUE){
			f <- function(pars=c(NA, NA, NA)){
				if(pars[2] < 0 | pars[1] < 0 | pars[3] < 0 | pars[3] > 100) return(-Inf)
				return(likelihood(model=model, data=data, mean=pars[1], variance=pars[2], zi=pars[3], silent=TRUE))
			}
		}else{
			cat("\n\tmean\t\tvariance\t\tzi\n", eol, "\t", signif(mean, 4), "\t\t", signif(variance, 4), "\t\t", signif(zi, 4), eol, sep="")
			f <- function(pars=c(NA, NA, NA)){
				if(pars[2] < 0 | pars[1] < 0 | pars[3] < 0 | pars[3] > 100) return(-Inf)
				cat(clearline, "\t", signif(pars[1], 4), "\t\t", signif(pars[2], 4), "\t\t", signif(pars[3], 4), eol, sep="")
				return(likelihood(model=model, data=data, mean=pars[1], variance=pars[2], zi=pars[3], silent=TRUE))
			}
		}
		maxim <- optim(c(mean, variance, zi), f, control=list(fnscale=-1))
		results <- c(mean=maxim$par[1], variance=maxim$par[2], zi=maxim$par[3])
	}
	
	if(silent==FALSE){
		#if(guitest$R.GUI == "AQUA" & guitest$R.package.type == "mac.binary") cat("\n")
		cat("\n\nFinished maximising the likelihood\n\n")
	}
	return(results)
}