print.fecrt.results <- function(x=list(), filename=FALSE, ...){
		
	if(filename!=FALSE) file=filename else file=""
	
	N <- length(x$animal.names)
	
	digits <- 1
	
	if(filename==FALSE) cat("\n", file=file,append=TRUE)
	
	cat(paste("Results of the feacal egg count reduction test analysis for '", x$name, "' with ", N, " animals:\n\n", sep=""), file=file,append=TRUE)
	
	if(all(is.na(x$results.mcmc))) cat("The Bayesian MCMC analysis method was not used\n\n", file=file,append=TRUE) else cat(paste("The Bayesian MCMC method concluded that these data have a ", round(x$prob.mcmc, digits), "% probability of coming from a resistant herd, which is defined as '", x$results.mcmc, "'.  The 95% credible interval and median result for the mean efficacy according to the Bayesian method is:\nl95: ", round(x$quant.mcmc[1], digits), "%    median: ", round(x$quant.mcmc[2], digits), "%    u95: ", round(x$quant.mcmc[3], digits), "%\n\n", sep=""), file=file,append=TRUE)
	
	if(identical(NA, x$results.boot)){
		cat("Results from bootstrapping and the WAAVP method were not available\n\n", file=file,append=TRUE)
		
	}else{
		cat(paste("The bootstrap method concluded that these data have a ", round(x$prob.boot,digits), "% probability of coming from a resistant herd, which is defined as '", x$results.boot, "'.  The 95% confidence interval and median result for the mean efficacy according to the bootstrap method is:\nl95: ", round(x$quant.boot[1], digits), "%    median: ", round(x$quant.boot[2],digits), "%    u95: ", round(x$quant.boot[3],digits), "%\n\n", sep=""), file=file,append=TRUE)
		cat(paste("The WAAVP method defined the data as '", x$results.boot, "'.\nThe 95% confidence interval and median result for the mean efficacy according to the WAAVP method is:\nl95: ", round(x$quant.waavp[1],digits), "%    median: ", round(x$quant.waavp[2],digits), "%    u95: ", round(x$quant.waavp[3],digits), "%\n\n", sep=""), file=file,append=TRUE)
		
	}
	
	if(all(!is.na(x$results.mcmc))){
	cat("\nThe following additional results were obtained using the Bayesian MCMC method:\n\n", file=file,append=TRUE)
	cat(paste("The median and 95% credible intervals for the true mean egg count before treatment:\nl95: ", round(x$meanquant[1],digits+1), "    median: ", round(x$meanquant[2],digits+1), "    u95: ", round(x$meanquant[3],digits+1), "\n\n", sep=""), file=file,append=TRUE)
	cat(paste("The median and 95% credible intervals for the change in shape parameter (<1 denotes an increase in variability and >1 a decrease in variability):\nl95: ", round(x$dshapequant[1],digits+1), "    median: ", round(x$dshapequant[2],digits+1), "    u95: ", round(x$dshapequant[3],digits+1), "\n\n", sep=""), file=file,append=TRUE)
	if(all(is.na(x$ziquant))) cat("The zero-inflated model was not used\n\n", file=file, append=TRUE) else cat(paste("The median and 95% credible intervals for the zero-inflation:\nl95: ", round(x$ziquant[1],digits+1), "%    median: ", round(x$ziquant[2],digits+1), "%    u95: ", round(x$ziquant[3],digits+1), "%\n\n", sep=""), file=file,append=TRUE)
	if(!all(is.na(x$indredquant))){
		cat(paste("Individual animal probabilities of 'efficacy < ", x$efficacy, "%', along with median and 95% credible intervals for mean egg count reduction:\n\n", sep=""), file=file,append=TRUE)
	for(i in 1:N){
		cat(paste(x$animal.names[i], ":  prob reduction < ", x$efficacy, "%: ", round(x$ind.prob.inf[i],digits), "%    l95: ", round(x$indredquant[i,1],digits), "%    median:", round(x$indredquant[i,2],digits), "%    u95: ", round(x$indredquant[i,3],digits), "%\n", sep=""), file=file,append=TRUE)
	}
	cat("\n", file=file, append=TRUE)
	}else{
		cat("Individual animal efficacies were not assessed\n\n", file=file, append=TRUE)
	}
	
	if(!x$converged) cat("*WARNING* The chains did not achieve convergence during the simulation, you should interpret the Bayesian MCMC results with extreme caution\n\n")
	cat(paste("Note:  The prior for efficacy with the Bayesian MCMC method was ", if(!x$restrict.efficacy) "not ", "restricted to values of greater than 0%\n\n", sep=""), file=file,append=TRUE)
	
}

	cat("*THIS SOFTWARE IS INTENDED FOR EDUCATIONAL PURPOSES ONLY AND SHOULD NOT BE RELIED UPON FOR REAL WORLD APPLICATIONS*\n", file=file,append=TRUE)
	
    cat("\n", file=file,append=TRUE)
	
}

assess.variance <- function(model){
	
	alt.prior <- model$alt.prior
	l.95 <- model$l.95
	u.95 <- model$u.95
	model <- toupper(model$model)

	largeod <- FALSE

	if(model=="GP" | model=="ZIGP"){				
		if(alt.prior==TRUE | is.character(alt.prior)){
			if((l.95 < 0.002) && (u.95 > 0.1)){
				largeod <- TRUE
			}
		}else{
			if((l.95 < 0.01) && (u.95 > 100)){
				largeod <- TRUE
			}
		}
	}
	
	if(model=="WP" | model=="ZIWP"){				
		if(alt.prior==TRUE | is.character(alt.prior)){
			if((l.95 < 0.0001) && (u.95 > 0.01)){
				largeod <- TRUE
			}
		}else{
			if((l.95 < 0.001) && (u.95 > 10)){
				largeod <- TRUE
			}
		}
	}

	if(model=="LP" | model=="ZILP"){				
		if(alt.prior==TRUE | is.character(alt.prior)){
			if((l.95 < 0.01) && (u.95 > 1000)){
				largeod <- TRUE
			}
		}else{
			if((l.95 < 0.001) && (u.95 > 10)){
				largeod <- TRUE
			}
		}
	}

	return(largeod)
}