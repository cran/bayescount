normal.params <- function(log.mean=stop("No log mean specified"), log.sd=stop("No log standard deviation specified")){

lsd <- log.sd
lmu <- log.mean

mu <- exp(lmu + ((lsd^2) / 2))
sd <- sqrt(exp((2*lmu)+(lsd^2)) * (exp(lsd^2) - 1))

if(length(log.mean)>1 | length(log.sd) >1){
	results <- matrix(NA, ncol=2, nrow=length(log.mean))
	results[,1] <- mu
	results[,2] <- sd
}else{
	results <- c(mu, sd)
}
return(results)
}

lnormal.params <- function(mean=stop("No mean specified"), sd=stop("No standard deviation specified")){

tau <- sd / mean
lsd <- sqrt(log(tau^2 +1))
lmu <- log(mean) - ((lsd^2) / 2)

if(length(mean)>1 | length(sd) >1){
	results <- matrix(NA, ncol=2, nrow=length(mean))
	results[,1] <- lmu
	results[,2] <- lsd
}else{
	results <- c(lmu, lsd)
}
return(results)
}