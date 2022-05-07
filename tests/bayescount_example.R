library("bayescount")

data <- array(dim=c(20,10,2))
means1 <- rgamma(20, 10, 1)
means2 <- rgamma(20, 5, 1)
for(i in 1:20){
data[i,,1] <- rpois(10, means1[i])
data[i,,2] <- rpois(10, means2[i])
}
# Missing data is permissible but means the likelihood cannot be
# calculated - a warning will be printed:
data[sample(1:(20*10*2), 10)] <- NA
try(unlink("analysis.ZILP.csv"), silent=TRUE)
# Run the analysis:
bayescount(name="analysis", data=data, model = "ZILP",
setnames=c("Simulated group A", "Simulated group B"), likelihood=TRUE)
unlink("analysis.ZILP.csv")
