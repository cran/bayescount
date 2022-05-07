library('bayescount')

data <- rpois(100, 10)
data[1:15] <- 0
likelihood('ZIP', data, mean=10, zi=15)
# now calculate the likelihood for the same data using an MCMC object
# to provide the values for mean and zero-inflation
## Not run:
values <- fec.analysis(data, model='ZISP', raw.output=TRUE)$mcmc
means <- c(values[,'mean'][[1]], values[,'mean'][[2]])
zis <- (1-c(values[,'prob'][[1]], values[,'prob'][[2]]))*100
# The function outputs the prevalence of disease when raw.ouput is
# TRUE, so zero-inflation must be calculated from this
likes <- likelihood('ZIP', data, mean=means, zi=zis,
raw.output=TRUE)$likelihood
hist(likes, breaks='fd', col='red')