dump.format <- function(variable, value){		# usage eg jags.string(c("f", "g"), list(c(1:10), 5))
	if(length(variable)==1){
		value <- list(value)
	}
	
	output.string <- ""
	for(i in 1:length(variable)){
		if(length(value[[i]])==1 && length(dim(value[[i]]))==0){
			value.string <- as.character(value[[i]])
		}else{
			dims <- dim(value[[i]])
			if(length(dims) > 1){
				value.string <- "structure(c("			
			}else{
				value.string <- "c("
			}
			value.string <- paste(value.string, paste(value[[i]], collapse=", "), ")", sep="")
			if(length(dims) > 1){
				value.string <- paste(value.string, ", .Dim = c(", paste(dims, collapse=", "), "))", sep="")
			}
		}
		output.string <- paste(output.string, "\"", variable[[i]], "\" <- ", value.string, "\n", sep="")
	}
	
	
	return(output.string)
}
