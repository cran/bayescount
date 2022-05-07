library('bayescount')
data <- rpois(100, rlnorm(3, 0.2))
model <- run.model(model="LP", data=data, call.jags=FALSE)
library('runjags')
results <- extend.jags(model, burnin=5000, sample=10000)

