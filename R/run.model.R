run.model <- function(data=stop("No data supplied"), model=stop("No model specified"), burnin = 5000, updates = 10000, call.jags = TRUE, alt.prior=FALSE, jags="jags", silent.jags = FALSE, check.conv=TRUE, monitor.lambda=FALSE){


###  Currently only the IP model requires gamma monitored for the likelihood bit - others are integrated

if(monitor.lambda==TRUE && model!="IP"){
	monitor.lambda <- FALSE
}

counts <- data
N <- length(counts)

models <- c("SP", "ZISP", "GP", "ZIGP", "LP", "ZILP", "WP", "ZIWP", "IP")
modelsfull <- c("single Poisson", "zero-inflated single Poisson", "gamma Poisson", "zero-inflated gamma Poisson", "lognormal Poisson", "zero-inflated lognormal Poisson", "Weibull Poisson", "zero-inflated Weibull Poisson", "independant Poisson")

if(sum(model==models) != 1){
	cat("Invalid model selection.  Please choose from ONE of the following models: ", sep="")
	cat(models, sep=", ")
	cat("\n")
	stop("Invalid model selection")
}

chains <- 2

inits <- matrix("", ncol=5, nrow=2)
		
dataformean <- data
dataformean[dataformean==0] <- NA
meanofnonzeros <- as.integer(mean(na.omit(dataformean)))

initstring <- character(length=chains)

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

if(model=="SP" | model=="ZISP"){
	initstring[1] <- dump.format("mean", smallestmean)
	initstring[2] <- dump.format("mean", largestmean)
	
}

if(model=="LP" | model=="ZILP"){
	
	if(class(alt.prior)=="character"){
		priorstring <- paste("sd ~ ", alt.prior, ";\n", sep="")
	}else{
		if(alt.prior==FALSE){
			priorstring <- "sd ~ dexp(0.8);\n"
		}else{
			priorstring <- "sd ~ dunif(0,10);\n"
		}	
	}
	
	lgammas1 = lgammas2 <- numeric(length=length(data))
	lgammas1[] <- 1
	lgammas2[] <- 1
	
	initstring[1] <- dump.format(c("pmean", "lgamma", "sd"), list(smallestmean, lgammas1, 0.1))
	initstring[2] <- dump.format(c("pmean", "lgamma", "sd"), list(largestmean, lgammas2, 5))			
}

if(model=="IP"){

	gammas1 <- data
	gammas2 <- data
	gammas1[gammas1 == 0] <- 1
	gammas2[gammas2 == 0] <- 1

	initstring[1] <- dump.format("gamma", gammas1)
	initstring[2] <- dump.format("gamma", gammas2)
}

if(model=="GP" | model=="ZIGP"){

	gammas1 = gammas2 <- numeric(length=length(data))
	gammas1[] <- 1
	gammas2[] <- 1
	
	initstring[1] <- dump.format(c("mean", "gamma", "sd"), list(smallestmean, gammas1, 0.1))
	initstring[2] <- dump.format(c("mean", "gamma", "sd"), list(largestmean, gammas2, 10))			
	
	if(class(alt.prior)=="character"){
		priorstring <- paste("b ~ ", alt.prior, ";\n", sep="")
		initstring[1] <- paste(initstring[1], dump.format("b", 0.1), sep="")
		initstring[2] <- paste(initstring[2], dump.format("b", 10), sep="")
	}else{
		if(alt.prior==FALSE){
			priorstring <- "b <- exp(logb);\nlogb ~ dunif(-7,7);\n"
			initstring[1] <- paste(initstring[1], dump.format("logb", -2.3), sep="")
			initstring[2] <- paste(initstring[2], dump.format("logb", 2.3), sep="")
		}else{
			priorstring <- "b ~ dunif(0.001,1000);\n"
			initstring[1] <- paste(initstring[1], dump.format("b", 0.1), sep="")
			initstring[2] <- paste(initstring[2], dump.format("b", 10), sep="")
		}
	}
	
}

if(model=="WP" | model=="ZIWP"){

	if(class(alt.prior)=="character"){
		priorstring <- paste("a ~ ", alt.prior, ";\nb ~ ", alt.prior, ";\n", sep="")
	}else{
		if(alt.prior==FALSE){
			priorstring <- "a ~ dgamma(1,0.01);\nb ~ dgamma(1,0.01);\n"
		}else{
			priorstring <- "a ~ dunif(0.001,1000);\nb ~ dunif(0.001,1000);\n"
		}
	}
	
	gammas1 = gammas2 <- numeric(length=length(data))
	gammas1[] <- 1
	gammas2[] <- 100
	
	initstring[1] <- dump.format(c("a", "gamma", "b"), list(0.01, gammas1, 0.1))
	initstring[2] <- dump.format(c("a", "gamma", "b"), list(100, gammas2, 100))
	
}

if(model=="ZISP" | model=="ZILP" | model=="ZIGP" | model=="ZIWP"){

	probpos <- data
	zeros <- probpos==0
	probpos[] <- 1
	probpos[zeros] <- 0
	
	initstring[1] <- paste(initstring[1], dump.format(c("prob", "probpos"), list(0.05, probpos)))
	probpos[zeros] <- 1
	initstring[2] <- paste(initstring[2], dump.format(c("prob", "probpos"), list(0.95, probpos)))
	
}


##### make model string

if(model=="SP"){
modelstring <- paste("model {
for(row in 1 : N){
Count[row] ~ dpois(mean);
}

# Priors
mean ~ dunif(0.001,1000);
}")
monitors <- "mean"
}
if(model=="ZISP"){
modelstring <- paste("model {
for(row in 1 : N){
Count[row] ~ dpois(lambda[row]);
lambda[row] <- probpos[row] * mean;
probpos[row] ~ dbern(prob);
}

# Priors
mean ~ dunif(0.001,1000);
prob ~ dunif(0,1);
}")
monitors <- c("mean", "prob")
}

if(model=="IP"){
modelstring <- paste("model {
for(row in 1 : N){
Count[row] ~ dpois(gamma[row]);
# Priors
gamma[row] ~ dunif(0.001,1000);
}

mean <- mean(gamma[])
sd <- sd(gamma[])

}")
monitors <- c("mean", "sd")
if(monitor.lambda==TRUE){
	monitors <- c(monitors, "gamma")
}
}


if(model=="LP"){
modelstring <- paste("model {
for(row in 1 : N){
Count[row] ~ dpois(gamma[row]);
gamma[row] ~ dlnorm(lmu, lprec);

}

lprec <- 1/log((sd/pmean)^2 + 1);
lvar <- 1/lprec;

# Priors
lmu <- log(pmean);
pmean ~ dunif(0.001,1000);
", priorstring, "}", sep="")
monitors <- c("lmu", "lvar")
if(monitor.lambda==TRUE){
	monitors <- c(monitors, "gamma")
}
}

if(model=="ZILP"){
modelstring <- paste("model {
for(row in 1 : N){
Count[row] ~ dpois(lambda[row]);
lambda[row] <- probpos[row] * gamma[row];
gamma[row] ~ dlnorm(lmu, lprec);
probpos[row] ~ dbern(prob);
}

lprec <- 1/log((sd/pmean)^2 + 1);
lvar <- 1/lprec;

# Priors
lmu <- log(pmean);
pmean ~ dunif(0.001,1000);
prob ~ dunif(0,1);
", priorstring, "}", sep="")
monitors <- c("lmu", "lvar", "prob")
if(monitor.lambda==TRUE){
	monitors <- c(monitors, "gamma")#, "probpos")
}
}


if(model=="GP"){
modelstring <- paste("model {

for(row in 1 : N){
Count[row] ~ dpois(lambda[row]);
lambda[row] <- mean * gamma[row];
gamma[row] ~ dgamma(wb, wb);

}

wb <- mean * b;

# Priors
mean ~ dunif(0.001,1000);
", priorstring, "}", sep="")
monitors <- c("mean", "b")
if(monitor.lambda==TRUE){
	monitors <- c(monitors, "gamma")
}
}

if(model=="ZIGP"){
modelstring <- paste("model {

for(row in 1 : N){
Count[row] ~ dpois(lambda[row]);
lambda[row] <- probpos[row] * mean * gamma[row];
probpos[row] ~ dbern(prob);
gamma[row] ~ dgamma(wb, wb);

}

wb <- mean * b;

# Priors
prob ~ dunif(0,1);
mean ~ dunif(0.001,1000);
", priorstring, "}", sep="")
monitors <- c("mean", "b", "prob")
if(monitor.lambda==TRUE){
	monitors <- c(monitors, "gamma")#, "probpos")
}
}

if(model=="WP"){
modelstring <- paste("model {

for(row in 1 : N){
Count[row] ~ dpois(gamma[row]);
gamma[row] ~ dweib(a,nb);
}

nb <- exp(-(log(b) * a));

# Priors
", priorstring, "}", sep="")
monitors <- c("a", "b")
if(monitor.lambda==TRUE){
	monitors <- c(monitors, "gamma")
}
}

if(model=="ZIWP"){
modelstring <- paste("model {

for(row in 1 : N){
Count[row] ~ dpois(lambda[row]);
lambda[row] <- probpos[row] * gamma[row];
probpos[row] ~ dbern(prob);
gamma[row] ~ dweib(a, nb);

}

nb <- exp(-(log(b) * a));

# Priors
prob ~ dunif(0,1);
", priorstring, "}", sep="")
monitors <- c("a", "b", "prob")
if(monitor.lambda==TRUE){
	monitors <- c(monitors, "gamma")#, "probpos")
}
}

datastring <- dump.format(c("N", "Count"), list(N, counts))

if(call.jags==FALSE){
	return(list(modelstring, datastring, initstring, monitors))
}else{
	return(run.jags(data=datastring, model=modelstring, inits=initstring, monitor=monitors, burnin=burnin, updates=10000, jags=jags, silent.jags=silent.jags, check.conv=TRUE))
}

}
