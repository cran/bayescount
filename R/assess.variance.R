assess.variance <- function(model, alt.prior, l.95, u.95){

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