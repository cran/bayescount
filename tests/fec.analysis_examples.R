library('bayescount')

results <- fec.analysis(data=c(0,5,3,7,0,4,3,8,0),
                        model="ZILP", silent.jags=TRUE)
