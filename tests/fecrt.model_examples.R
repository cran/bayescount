library('bayescount')

data <- data.frame(Count=rpois(80,rep(c(10,10,2,2), 20)),
Subject=rep(1:20, each=4), Time=rep(rep(1:2,each=2),40),
Sample=1:2, Control=rep(c(0,1), each=40))
# Compile the model - a paired model is required because
# there are replicate samples within an individual:
model <- fecrt.model(data, paired.model=TRUE)
# Update the model - requires runjags:
library('runjags')
results <- extend.jags(model, burnin=5000)
