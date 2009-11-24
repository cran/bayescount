#include <math.h>
#include <R.h>
#include <Rmath.h>
#include <stdio.h>

int arel(int position[], int dimlengths[], int dims){
	
	/* Function to find the element number of a multi-dimensional array (probably obtained from R).  Needs to know the dimlengths and the number of dimensions */
	
	int element = position[0];
	int i;
	
	int mult[dims];
	
	mult[0] = 1;
	
	for(i=1; i<dims; i++){
		mult[i]=mult[i-1]*dimlengths[i-1];
		element = element + (mult[i]*(position[i]-1));
	}
	
	return(element-1);
	
}

/*
double sumdarray(double value[], int N){
	double sum=0;
	int i;
	for(i=0; i<N; i++){
		sum=sum+value[i];
	}
	return(sum);
}

int sumiarray(int value[], int N){
	int sum=0;
	int i;
	for(i=0; i<N; i++){
		sum=sum+value[i];
	}
	return(sum);
}
*/

double myround(double number, int decimals){
	number = number * pow(10,decimals);
	return ((double)((number >= 0) ? (int)(number + 0.5) : (int)(number - 0.5))) / pow(10,decimals);
}


/*
void poweranalysispopulation(double *meanepg, double *gfaeces, double *sensitivity, int *animals, double *coeffvarind, double *coeffvargroup, double *accuracy, int *miniterations, int *maxiterations, int *skip, int *precision, int *nin, int *ntotal){

double shapegp; 
double shapein;
	
shapegp = 1/(coeffvargroup[0]*coeffvargroup[0]);
shapein = 1/(coeffvarind[0]*coeffvarind[0]);

//printf("%f, %f", shapegp, coeffvargroup[0]);

nin[0] = 0;
ntotal[0] = 0;

double indmeans[animals[0]];
double samplemeans[animals[0]];

int sumcount;
int set,a, skipset;

double lci, uci, meancount;

int decimals = precision[0];

double upper = meanepg[0]*(1. + accuracy[0]);
double lower = meanepg[0]*(1. - accuracy[0]);

int newminiters;

//printf("%f, %f, \n", shapegp, (meanepg[0] / shapegp));

newminiters = (int)ceil((double)miniterations[0] / (double)skip[0]) * skip[0];

for(set=0; set<newminiters; set++){
	
	sumcount = 0;
	
	for(a=0; a<animals[0]; a++){
		GetRNGstate();
		indmeans[a] = rgamma(shapegp, (meanepg[0] / shapegp));
		PutRNGstate(); 
		GetRNGstate();
		samplemeans[a] = rgamma(shapein*gfaeces[0], indmeans[a]/(shapein*gfaeces[0]));
		PutRNGstate(); 
		GetRNGstate();
		sumcount = sumcount + rpois(samplemeans[a]);
		PutRNGstate(); 
	}
	
	ntotal[0]++;
	meancount = ((double)sumcount)/animals[0];
	nin[0] = nin[0] + ((meancount <= upper) && (meancount >= lower));
	//printf("%f, %f, %f\n", lower, meancount, upper);
	
	//printf("%f\n", sumdarray(indmeans, animals[0])/animals[0]);
}

//printf("\n\n");

lci = qbeta(0.025, (ntotal[0]-nin[0])+1, nin[0]+1, 1, 0);//myround(qbeta(0.025, (ntotal[0]-nin[0])+1, nin[0]+1, 0, 0), decimals);
uci = qbeta(0.975, (ntotal[0]-nin[0])+1, nin[0]+1, 1, 0);//myround(qbeta(0.975, (ntotal[0]-nin[0])+1, nin[0]+1, 0, 0), decimals);

printf("%f, %f\n", lci, uci);
printf("%f, %f\n", myround(lci, decimals), myround(uci,decimals));


for(set=miniterations[0]; set<(maxiterations[0]-skip[0]); set=set+skip[0]){
//for(set=0; set<maxiterations[0]; set++){
		
	sumcount = 0;
	
	for(skipset=0; skipset<skip[0]; skipset++){
		for(a=0; a<animals[0]; a++){
			GetRNGstate();
			indmeans[a] = rgamma(shapegp, (meanepg[0] / shapegp));
			PutRNGstate(); 
			GetRNGstate();
			samplemeans[a] = rgamma(shapein*gfaeces[0], indmeans[a]/(shapein*gfaeces[0]));
			PutRNGstate(); 
			GetRNGstate();
			sumcount = sumcount + rpois(samplemeans[a]);
			PutRNGstate(); 
		}
		ntotal[0]++;
		meancount = ((double)sumcount)/animals[0];
		nin[0] = nin[0] + ((meancount <= upper) && (meancount >= lower));
		//printf("%f, %f, %f\n", lower, meancount, upper);
	}
	
	lci = qbeta(0.025, (ntotal[0]-nin[0])+1, nin[0]+1, 1, 0);//myround(qbeta(0.025, (ntotal[0]-nin[0])+1, nin[0]+1, 0, 0), decimals);
	uci = qbeta(0.975, (ntotal[0]-nin[0])+1, nin[0]+1, 1, 0);//myround(qbeta(0.975, (ntotal[0]-nin[0])+1, nin[0]+1, 0, 0), decimals);
	
	//printf("%f, %f\n", lci, uci);
	if(myround(lci, decimals) == myround(uci,decimals)){
		break;
	}
	
	//qbeta(double, double, double, int, int);
	//qbeta(c(0.025, 0.5, 0.975), iterations-ci, ci, lower.tail=FALSE))
	//printf("%f\n", sumdarray(indmeans, animals[0])/animals[0]);
}


printf("%f, %f\n", lci, uci);
printf("%f, %f\n", myround(lci, decimals), myround(uci,decimals));

}

*/

void poweranalysispopulation(double *meanepg, double *gfaeces, double *sensitivity, int *replicates, int *animals, double *coeffvarrep, double *coeffvarind, double *coeffvargroup, double *lowerl, double *upperl, int *maxiterations, int *precision, double *lcil, double *ucil, int *print, int *nin, int *ntotal){

double shapegp;
double shapeindt;

shapeindt = 1.0 / (pow(coeffvarind[0],2) + pow(coeffvarrep[0]/sqrt(gfaeces[0]),2) + pow(coeffvarind[0],2)*pow(coeffvarrep[0]/sqrt(gfaeces[0]),2));
shapegp = 1/(coeffvargroup[0]*coeffvargroup[0]);

//Rprintf("CVS:  %f, %f\n", shapeindt, shapegp);

double indmeans;//[animals[0]];
double replicatemeans;//[animals[0]];

register unsigned int set, a, skipset;
double lci, uci, meancount, sumcount;

register unsigned int decimals = precision[0];

double lower = lowerl[0];
double upper = upperl[0];

if(print[0]){
	Rprintf("< Determining power >\n   l95       u95    iteration\n");
}

GetRNGstate();

for(set=maxiterations[0]; set--; ){
//Same as for(set=0; set<maxiterations[0]; set++){ but quicker - INDEXING IS BACKWARDS	
	sumcount = 0.;
	
	for(a=animals[0]; a--; ){
	//Same as for(a=0; a<animals[0]; a++){ but quicker - INDEXING IS BACKWARDS
		indmeans = rgamma(shapegp, (meanepg[0] / shapegp));
		replicatemeans = rgamma(shapeindt*replicates[0], indmeans/(shapeindt*replicates[0]));
		sumcount = sumcount + (((double)rpois(replicatemeans*replicates[0]*sensitivity[0]))*(1/sensitivity[0]));
	}
	ntotal[0]++;
	meancount = sumcount/(animals[0]*replicates[0]);
	
	nin[0] = nin[0] + ((meancount <= upper) && (meancount >= lower));
	
	lci = qbeta(lcil[0], nin[0]+1, (ntotal[0]-nin[0])+1, 1, 0);
	uci = qbeta(ucil[0], nin[0]+1, (ntotal[0]-nin[0])+1, 1, 0);
	
	if(print[0]){
		Rprintf("%f, %f, %i\r", lci, uci, ntotal[0]);
	}
	if(myround(lci, decimals) == myround(uci,decimals)){
		break;
	}
}

PutRNGstate(); 

if(print[0]){
	if(myround(lci, decimals) == myround(uci,decimals)){
		Rprintf("%f, %f, %i\n", lci, uci, ntotal[0]);
		Rprintf("%f, %f, (rounded)\n", myround(lci, decimals), myround(uci,decimals));
		Rprintf("< Power determined >\n");
	}else{
		Rprintf("%f, %f, %i\n", lci, uci, ntotal[0]);
		Rprintf("%f, %f, (rounded)\n", myround(lci, decimals), myround(uci,decimals));
		Rprintf("< Power not determined >\n");
	}
}

}


void poweranalysispopulationfixed(double *meanepg, double *gfaeces, double *sensitivity, int *replicates, int *animals, double *coeffvarrep, double *coeffvarind, double *coeffvargroup, int *maxiterations, int *print, double *meancounts){

double shapegp;
double shapeindt;

shapeindt = 1.0 / (pow(coeffvarind[0],2) + pow(coeffvarrep[0]/sqrt(gfaeces[0]),2) + pow(coeffvarind[0],2)*pow(coeffvarrep[0]/sqrt(gfaeces[0]),2));
shapegp = 1/(coeffvargroup[0]*coeffvargroup[0]);

//Rprintf("CVS:  %f, %f\n", shapeindt, shapegp);

double indmeans;//[animals[0]];
double replicatemeans;//[animals[0]];

register unsigned int set, a, skipset;
double lci, uci, meancount, sumcount;

if(print[0]){
	Rprintf("< Running simulation >\n0%% complete\n");
}

GetRNGstate();

for(set=0; set<maxiterations[0]; set++){
//Same as for(set=0; set<maxiterations[0]; set++){ but quicker - INDEXING IS BACKWARDS	
	sumcount = 0.;
	
	for(a=animals[0]; a--; ){
	//Same as for(a=0; a<animals[0]; a++){ but quicker - INDEXING IS BACKWARDS
		indmeans = rgamma(shapegp, (meanepg[0] / shapegp));
		replicatemeans = rgamma(shapeindt*replicates[0], indmeans/(shapeindt*replicates[0]));
		sumcount = sumcount + (((double)rpois(replicatemeans*replicates[0]*sensitivity[0]))*(1/sensitivity[0]));
	}

	meancount = sumcount/(animals[0]*replicates[0]);
	
	meancounts[set] = meancount;
	
	if(print[0]){
		Rprintf("%f%% complete\r", set/maxiterations[0]);
	}
}

PutRNGstate(); 

if(print[0]){
	Rprintf("< Finished >\n");
}

}





void poweranalysissample(double *meanepg, double *gfaeces, double *sensitivity, int *replicates, int *animals, double *coeffvarrep, double *coeffvarind, double *coeffvargroup, double *lowerl, double *upperl, int *maxiterations, int *precision, double *lcil, double *ucil, int *print, int *nin, int *ntotal){

double shapegp;
double shapeindt;

shapeindt = 1.0 / (pow(coeffvarind[0],2) + pow(coeffvarrep[0]/sqrt(gfaeces[0]),2) + pow(coeffvarind[0],2)*pow(coeffvarrep[0]/sqrt(gfaeces[0]),2));
shapegp = 1/(coeffvargroup[0]*coeffvargroup[0]);

//Rprintf("CVS:  %f, %f\n", shapeindt, shapegp);

double indmeans[animals[0]];
double replicatemeans;//[animals[0]];

register unsigned int set, a, skipset, done;
double lci, uci, meancount, sumcount, samplesum, adjust;

register unsigned int decimals = precision[0];

double lower = lowerl[0];
double upper = upperl[0];

if(print[0]){
	Rprintf("< Determining power >\n   l95       u95    iteration\n");
}

GetRNGstate();

for(set=maxiterations[0]; set--; ){
//Same as for(set=0; set<maxiterations[0]; set++){ but quicker - INDEXING IS BACKWARDS	
	
	for(;;){
		samplesum = 0;
		done = 1;
	
		for(a=animals[0]; a--; ){
		//Same as for(a=0; a<animals[0]; a++){ but quicker - INDEXING IS BACKWARDS
			indmeans[a] = rgamma(shapegp, (meanepg[0] / shapegp));
			samplesum = samplesum+indmeans[a];
		}
	
		adjust = meanepg[0] - (samplesum / animals[0]);
	
		for(a=animals[0]; a--; ){
			indmeans[a] = indmeans[a] + adjust;
			if(indmeans[a] < 0){
				done = 0;
				break;
			}
		}
		
		if(done==1){
			break;
		}
	}
	
	sumcount = 0.;
		
	for(a=animals[0]; a--; ){
		replicatemeans = rgamma(shapeindt*replicates[0], indmeans[a]/(shapeindt*replicates[0]));
		sumcount = sumcount + (((double)rpois(replicatemeans*replicates[0]*sensitivity[0]))*(1/sensitivity[0]));
	}
	
	ntotal[0]++;
	meancount = sumcount/(animals[0]*replicates[0]);
	
	nin[0] = nin[0] + ((meancount <= upper) && (meancount >= lower));
	
	lci = qbeta(lcil[0], nin[0]+1, (ntotal[0]-nin[0])+1, 1, 0);
	uci = qbeta(ucil[0], nin[0]+1, (ntotal[0]-nin[0])+1, 1, 0);
	
	if(print[0]){
		Rprintf("%f, %f, %i\r", lci, uci, ntotal[0]);
	}
	if(myround(lci, decimals) == myround(uci,decimals)){
		break;
	}
}

PutRNGstate();

if(print[0]){
	if(myround(lci, decimals) == myround(uci,decimals)){
		Rprintf("%f, %f, %i\n", lci, uci, ntotal[0]);
		Rprintf("%f, %f, (rounded)\n", myround(lci, decimals), myround(uci,decimals));
		Rprintf("< Power determined >\n");
	}else{
		Rprintf("%f, %f, %i\n", lci, uci, ntotal[0]);
		Rprintf("%f, %f, (rounded)\n", myround(lci, decimals), myround(uci,decimals));
		Rprintf("< Power not determined >\n");
	}
}

}


void poweranalysissamplefixed(double *meanepg, double *gfaeces, double *sensitivity, int *replicates, int *animals, double *coeffvarrep, double *coeffvarind, double *coeffvargroup, int *maxiterations, int *print, double *meancounts){

double shapegp;
double shapeindt;

shapeindt = 1.0 / (pow(coeffvarind[0],2) + pow(coeffvarrep[0]/sqrt(gfaeces[0]),2) + pow(coeffvarind[0],2)*pow(coeffvarrep[0]/sqrt(gfaeces[0]),2));
shapegp = 1/(coeffvargroup[0]*coeffvargroup[0]);

//Rprintf("CVS:  %f, %f\n", shapeindt, shapegp);

double indmeans[animals[0]];
double replicatemeans;//[animals[0]];

register unsigned int set, a, skipset;
double lci, uci, meancount, sumcount, done, samplesum, adjust;

if(print[0]){
	Rprintf("< Running simulation >\n0%% complete\n");
}

GetRNGstate();

for(set=maxiterations[0]; set--; ){
//Same as for(set=0; set<maxiterations[0]; set++){ but quicker - INDEXING IS BACKWARDS	
	
	for(;;){
		samplesum = 0;
		done = 1;
	
		for(a=animals[0]; a--; ){
		//Same as for(a=0; a<animals[0]; a++){ but quicker - INDEXING IS BACKWARDS
			indmeans[a] = rgamma(shapegp, (meanepg[0] / shapegp));
			samplesum = samplesum+indmeans[a];
		}
	
		adjust = meanepg[0] - (samplesum / animals[0]);
	
		for(a=animals[0]; a--; ){
			indmeans[a] = indmeans[a] + adjust;
			if(indmeans[a] < 0){
				done = 0;
				break;
			}
		}
		
		if(done==1){
			break;
		}
	}
	
	sumcount = 0.;
		
	for(a=animals[0]; a--; ){
		replicatemeans = rgamma(shapeindt*replicates[0], indmeans[a]/(shapeindt*replicates[0]));
		sumcount = sumcount + (((double)rpois(replicatemeans*replicates[0]*sensitivity[0]))*(1/sensitivity[0]));
	}
	
	meancount = sumcount/(animals[0]*replicates[0]);
	meancounts[set] = meancount;
		
	if(print[0]){
		Rprintf("%f%% complete\r", set/maxiterations[0]);
	}
}

PutRNGstate();

if(print[0]){
	Rprintf("< Finished >\n");
}

}





//out <- .C("fecrtpowerpopulation", as.numeric(meanepg), as.numeric(delta), as.numeric(g.faeces), as.numeric(sensitivity), as.integer(replicates), as.integer(animals), as.numeric(pre.coeffvarrep), as.numeric(pre.coeffvarind), as.numeric(pre.coeffvargroup), as.numeric(post.coeffvarrep), as.numeric(post.coeffvarind), as.numeric(post.coeffvargroup), as.numeric(lowerl), as.numeric(upperl), as.integer(maxiterations), as.integer(precision), as.numeric(lci), as.numeric(uci), as.integer(feedback), as.integer(0), as.integer(0), PACKAGE="bayescount")

void fecrtpowerpopulation(double *meanepg, double *reduction, double *gfaeces, double *sensitivity, int *replicates, int *animals, double *precoeffvarrep, double *precoeffvarind, double *precoeffvargroup, double *postcoeffvarrep, double *postcoeffvarind, double *postcoeffvargroup, double *lowerl, double *upperl, int *maxiterations, int *precision, double *lcil, double *ucil, int *print, int *nin, int *ntotal){

double preshapegp;
double preshapeindt;
double postshapegp;
double postshapeindt;

preshapeindt = 1.0 / (pow(precoeffvarind[0],2) + pow(precoeffvarrep[0]/sqrt(gfaeces[0]),2) + pow(precoeffvarind[0],2)*pow(precoeffvarrep[0]/sqrt(gfaeces[0]),2));
preshapegp = 1/(precoeffvargroup[0]*precoeffvargroup[0]);

postshapeindt = 1.0 / (pow(postcoeffvarind[0],2) + pow(postcoeffvarrep[0]/sqrt(gfaeces[0]),2) + pow(postcoeffvarind[0],2)*pow(postcoeffvarrep[0]/sqrt(gfaeces[0]),2));
postshapegp = 1/(postcoeffvargroup[0]*postcoeffvargroup[0]);

double indmeans;//[animals[0]];
double replicatemeans;//[animals[0]];

register unsigned int set, a, skipset;
double lci, uci, meanred, presumcount, postsumcount;

register unsigned int decimals = precision[0];

double lower = lowerl[0];
double upper = upperl[0];

if(print[0]){
	Rprintf("< Determining power >\n   l95       u95    iteration\n");
}

GetRNGstate();

for(set=maxiterations[0]; set--; ){
//Same as for(set=0; set<maxiterations[0]; set++){ but quicker - INDEXING IS BACKWARDS	
	presumcount = 0.;
	postsumcount = 0.;
	
	for(a=animals[0]; a--; ){
	//Same as for(a=0; a<animals[0]; a++){ but quicker - INDEXING IS BACKWARDS
		indmeans = rgamma(preshapegp, (meanepg[0] / preshapegp));
		replicatemeans = rgamma(preshapeindt*replicates[0], indmeans/(preshapeindt*replicates[0]));
		presumcount = presumcount + (((double)rpois(replicatemeans*replicates[0]*sensitivity[0]))*(1/sensitivity[0]));
		
		indmeans = rgamma(postshapegp, ((meanepg[0]*reduction[0]) / postshapegp));
		replicatemeans = rgamma(postshapeindt*replicates[0], indmeans/(postshapeindt*replicates[0]));
		postsumcount = postsumcount + (((double)rpois(replicatemeans*replicates[0]*sensitivity[0]))*(1/sensitivity[0]));
	}
	ntotal[0]++;
	
	//Will give NaN if pre=post=0 or Inf if pre=0 & post>0.  The reductions for these should be 0, so delta =1:
	if(presumcount==0){
		meanred=1;
	}else{
		meanred = postsumcount / presumcount;
	}
	
	nin[0] = nin[0] + ((meanred <= upper) && (meanred >= lower));
	
	lci = qbeta(lcil[0], nin[0]+1, (ntotal[0]-nin[0])+1, 1, 0);
	uci = qbeta(ucil[0], nin[0]+1, (ntotal[0]-nin[0])+1, 1, 0);
	
	if(print[0]){
		Rprintf("%f, %f, %i\r", lci, uci, ntotal[0]);
	}
	if(myround(lci, decimals) == myround(uci,decimals)){
		break;
	}
}

PutRNGstate(); 

if(print[0]){
	if(myround(lci, decimals) == myround(uci,decimals)){
		Rprintf("%f, %f, %i\n", lci, uci, ntotal[0]);
		Rprintf("%f, %f, (rounded)\n", myround(lci, decimals), myround(uci,decimals));
		Rprintf("< Power determined >\n");
	}else{
		Rprintf("%f, %f, %i\n", lci, uci, ntotal[0]);
		Rprintf("%f, %f, (rounded)\n", myround(lci, decimals), myround(uci,decimals));
		Rprintf("< Power not determined >\n");
	}
}

}



//out <- .C("fecrtpowersample", as.numeric(meanepg), as.numeric(delta), as.numeric(g.faeces), as.numeric(sensitivity), as.integer(replicates), as.integer(animals), as.numeric(pre.coeffvarrep), as.numeric(pre.coeffvarind), as.numeric(pre.coeffvargroup), as.numeric(post.coeffvarrep), as.numeric(post.coeffvarind), as.numeric(post.coeffvargroup), as.numeric(lowerl), as.numeric(upperl), as.integer(maxiterations), as.integer(precision), as.numeric(lci), as.numeric(uci), as.integer(feedback), as.integer(0), as.integer(0), PACKAGE="bayescount")

void fecrtpowersample(double *meanepg, double *reduction, double *gfaeces, double *sensitivity, int *replicates, int *animals, double *precoeffvarrep, double *precoeffvarind, double *precoeffvargroup, double *postcoeffvarrep, double *postcoeffvarind, double *postcoeffvargroup, double *lowerl, double *upperl, int *maxiterations, int *precision, double *lcil, double *ucil, int *print, int *nin, int *ntotal){

double preshapegp;
double preshapeindt;
double postshapegp;
double postshapeindt;

preshapeindt = 1.0 / (pow(precoeffvarind[0],2) + pow(precoeffvarrep[0]/sqrt(gfaeces[0]),2) + pow(precoeffvarind[0],2)*pow(precoeffvarrep[0]/sqrt(gfaeces[0]),2));
preshapegp = 1/(precoeffvargroup[0]*precoeffvargroup[0]);

postshapeindt = 1.0 / (pow(postcoeffvarind[0],2) + pow(postcoeffvarrep[0]/sqrt(gfaeces[0]),2) + pow(postcoeffvarind[0],2)*pow(postcoeffvarrep[0]/sqrt(gfaeces[0]),2));
postshapegp = 1/(postcoeffvargroup[0]*postcoeffvargroup[0]);

double preindmeans[animals[0]];
double postindmeans[animals[0]];
double replicatemeans;//[animals[0]];

register unsigned int set, a, skipset, done;
double lci, uci, meanred, presumcount, postsumcount, samplesum, adjust;

register unsigned int decimals = precision[0];

double lower = lowerl[0];
double upper = upperl[0];

if(print[0]){
	Rprintf("< Determining power >\n   l95       u95    iteration\n");
}

GetRNGstate();

for(set=maxiterations[0]; set--; ){
//Same as for(set=0; set<maxiterations[0]; set++){ but quicker - INDEXING IS BACKWARDS	
	
	for(;;){
		samplesum = 0;
		done = 1;
	
		for(a=animals[0]; a--; ){
		//Same as for(a=0; a<animals[0]; a++){ but quicker - INDEXING IS BACKWARDS
			preindmeans[a] = rgamma(preshapegp, (meanepg[0] / preshapegp));
			samplesum = samplesum+preindmeans[a];
		}
	
		adjust = meanepg[0] - (samplesum / animals[0]);
	
		for(a=animals[0]; a--; ){
			preindmeans[a] = preindmeans[a] + adjust;
			if(preindmeans[a] < 0){
				done = 0;
				break;
			}
		}
		
		if(done==1){
			break;
		}
	}
	
	for(;;){
		samplesum = 0;
		done = 1;
	
		for(a=animals[0]; a--; ){
		//Same as for(a=0; a<animals[0]; a++){ but quicker - INDEXING IS BACKWARDS
			postindmeans[a] = rgamma(postshapegp, ((meanepg[0]*reduction[0]) / postshapegp));
			samplesum = samplesum+postindmeans[a];
		}
	
		adjust = (meanepg[0]*reduction[0]) - (samplesum / animals[0]);
	
		for(a=animals[0]; a--; ){
			postindmeans[a] = postindmeans[a] + adjust;
			if(postindmeans[a] < 0){
				done = 0;
				break;
			}
		}
		
		if(done==1){
			break;
		}
	}
	
	presumcount = 0.;
	postsumcount = 0.;
	
	for(a=animals[0]; a--; ){
	//Same as for(a=0; a<animals[0]; a++){ but quicker - INDEXING IS BACKWARDS
		replicatemeans = rgamma(preshapeindt*replicates[0], preindmeans[a]/(preshapeindt*replicates[0]));
		presumcount = presumcount + (((double)rpois(replicatemeans*replicates[0]*sensitivity[0]))*(1/sensitivity[0]));

		replicatemeans = rgamma(postshapeindt*replicates[0], postindmeans[a]/(postshapeindt*replicates[0]));
		postsumcount = postsumcount + (((double)rpois(replicatemeans*replicates[0]*sensitivity[0]))*(1/sensitivity[0]));

	}
	
	ntotal[0]++;
	
	//Will give NaN if pre=post=0 or Inf if pre=0 & post>0.  The reductions for these should be 0, so delta =1:
	if(presumcount==0){
		meanred=1;
	}else{
		meanred = postsumcount / presumcount;
	}
	
	nin[0] = nin[0] + ((meanred <= upper) && (meanred >= lower));
	
	lci = qbeta(lcil[0], nin[0]+1, (ntotal[0]-nin[0])+1, 1, 0);
	uci = qbeta(ucil[0], nin[0]+1, (ntotal[0]-nin[0])+1, 1, 0);
	
	if(print[0]){
		Rprintf("%f, %f, %i\r", lci, uci, ntotal[0]);
	}
	if(myround(lci, decimals) == myround(uci,decimals)){
		break;
	}
}

PutRNGstate(); 

if(print[0]){
	if(myround(lci, decimals) == myround(uci,decimals)){
		Rprintf("%f, %f, %i\n", lci, uci, ntotal[0]);
		Rprintf("%f, %f, (rounded)\n", myround(lci, decimals), myround(uci,decimals));
		Rprintf("< Power determined >\n");
	}else{
		Rprintf("%f, %f, %i\n", lci, uci, ntotal[0]);
		Rprintf("%f, %f, (rounded)\n", myround(lci, decimals), myround(uci,decimals));
		Rprintf("< Power not determined >\n");
	}
}

}


//out <- .C("fecrtpowerpopulationfixed", as.numeric(meanepg), as.numeric(delta), as.numeric(g.faeces), as.numeric(sensitivity), as.integer(replicates), as.integer(animals), as.numeric(pre.coeffvarrep), as.numeric(pre.coeffvarind), as.numeric(pre.coeffvargroup), as.numeric(post.coeffvarrep), as.numeric(post.coeffvarind), as.numeric(post.coeffvargroup), as.integer(maxiterations), as.integer(feedback), numeric(maxiterations), PACKAGE="bayescount")

void fecrtpowerpopulationfixed(double *meanepg, double *reduction, double *gfaeces, double *sensitivity, int *replicates, int *animals, double *precoeffvarrep, double *precoeffvarind, double *precoeffvargroup, double *postcoeffvarrep, double *postcoeffvarind, double *postcoeffvargroup, int *maxiterations, int *print, double *meanreds){

double preshapegp;
double preshapeindt;
double postshapegp;
double postshapeindt;

preshapeindt = 1.0 / (pow(precoeffvarind[0],2) + pow(precoeffvarrep[0]/sqrt(gfaeces[0]),2) + pow(precoeffvarind[0],2)*pow(precoeffvarrep[0]/sqrt(gfaeces[0]),2));
preshapegp = 1/(precoeffvargroup[0]*precoeffvargroup[0]);

postshapeindt = 1.0 / (pow(postcoeffvarind[0],2) + pow(postcoeffvarrep[0]/sqrt(gfaeces[0]),2) + pow(postcoeffvarind[0],2)*pow(postcoeffvarrep[0]/sqrt(gfaeces[0]),2));
postshapegp = 1/(postcoeffvargroup[0]*postcoeffvargroup[0]);

double indmeans;//[animals[0]];
double replicatemeans;//[animals[0]];

register unsigned int set, a, skipset;
double lci, uci, meanred, presumcount, postsumcount;

if(print[0]){
	Rprintf("< Running simulation >\n0%% complete\n");
}

GetRNGstate();

for(set=0; set<maxiterations[0]; set++){
//Same as for(set=0; set<maxiterations[0]; set++){ but quicker - INDEXING IS BACKWARDS	
	presumcount = 0.;
	postsumcount = 0.;
	
	for(a=animals[0]; a--; ){
	//Same as for(a=0; a<animals[0]; a++){ but quicker - INDEXING IS BACKWARDS
		indmeans = rgamma(preshapegp, (meanepg[0] / preshapegp));
		replicatemeans = rgamma(preshapeindt*replicates[0], indmeans/(preshapeindt*replicates[0]));
		presumcount = presumcount + (((double)rpois(replicatemeans*replicates[0]*sensitivity[0]))*(1/sensitivity[0]));
		
		indmeans = rgamma(postshapegp, ((meanepg[0]*reduction[0]) / postshapegp));
		replicatemeans = rgamma(postshapeindt*replicates[0], indmeans/(postshapeindt*replicates[0]));
		postsumcount = postsumcount + (((double)rpois(replicatemeans*replicates[0]*sensitivity[0]))*(1/sensitivity[0]));
	}

	//Will give NaN if pre=post=0 or Inf if pre=0 & post>0.  The reductions for these should be 0, so delta =1:
	if(presumcount==0){
		meanred=1;
	}else{
		meanred = postsumcount / presumcount;
	}
	
	meanreds[set] = meanred;
	
	if(print[0]){
		Rprintf("%f%% complete\r", set/maxiterations[0]);
	}
}

PutRNGstate(); 

if(print[0]){
	Rprintf("< Finished >\n");
}

}



//out <- .C("fecrtpowersamplefixed", as.numeric(meanepg), as.numeric(delta), as.numeric(g.faeces), as.numeric(sensitivity), as.integer(replicates), as.integer(animals), as.numeric(pre.coeffvarrep), as.numeric(pre.coeffvarind), as.numeric(pre.coeffvargroup), as.numeric(post.coeffvarrep), as.numeric(post.coeffvarind), as.numeric(post.coeffvargroup), as.integer(maxiterations), as.integer(feedback), numeric(maxiterations), PACKAGE="bayescount")


void fecrtpowersamplefixed(double *meanepg, double *reduction, double *gfaeces, double *sensitivity, int *replicates, int *animals, double *precoeffvarrep, double *precoeffvarind, double *precoeffvargroup, double *postcoeffvarrep, double *postcoeffvarind, double *postcoeffvargroup, int *maxiterations, int *print, double *meanreds){

double preshapegp;
double preshapeindt;
double postshapegp;
double postshapeindt;

preshapeindt = 1.0 / (pow(precoeffvarind[0],2) + pow(precoeffvarrep[0]/sqrt(gfaeces[0]),2) + pow(precoeffvarind[0],2)*pow(precoeffvarrep[0]/sqrt(gfaeces[0]),2));
preshapegp = 1/(precoeffvargroup[0]*precoeffvargroup[0]);

postshapeindt = 1.0 / (pow(postcoeffvarind[0],2) + pow(postcoeffvarrep[0]/sqrt(gfaeces[0]),2) + pow(postcoeffvarind[0],2)*pow(postcoeffvarrep[0]/sqrt(gfaeces[0]),2));
postshapegp = 1/(postcoeffvargroup[0]*postcoeffvargroup[0]);

double preindmeans[animals[0]];
double postindmeans[animals[0]];
double replicatemeans;//[animals[0]];

register unsigned int set, a, skipset, done;
double lci, uci, meanred, presumcount, postsumcount, samplesum, adjust;

if(print[0]){
	Rprintf("< Running simulation >\n0%% complete\n");
}

GetRNGstate();

for(set=0; set<maxiterations[0]; set++){
//Same as for(set=0; set<maxiterations[0]; set++){ but quicker - INDEXING IS BACKWARDS	
	
	for(;;){
		samplesum = 0;
		done = 1;
	
		for(a=animals[0]; a--; ){
		//Same as for(a=0; a<animals[0]; a++){ but quicker - INDEXING IS BACKWARDS
			preindmeans[a] = rgamma(preshapegp, (meanepg[0] / preshapegp));
			samplesum = samplesum+preindmeans[a];
		}
	
		adjust = meanepg[0] - (samplesum / animals[0]);
	
		for(a=animals[0]; a--; ){
			preindmeans[a] = preindmeans[a] + adjust;
			if(preindmeans[a] < 0){
				done = 0;
				break;
			}
		}
		
		if(done==1){
			break;
		}
	}
	
	for(;;){
		samplesum = 0;
		done = 1;
	
		for(a=animals[0]; a--; ){
		//Same as for(a=0; a<animals[0]; a++){ but quicker - INDEXING IS BACKWARDS
			postindmeans[a] = rgamma(postshapegp, ((meanepg[0]*reduction[0]) / postshapegp));
			samplesum = samplesum+postindmeans[a];
		}
	
		adjust = (meanepg[0]*reduction[0]) - (samplesum / animals[0]);
	
		for(a=animals[0]; a--; ){
			postindmeans[a] = postindmeans[a] + adjust;
			if(postindmeans[a] < 0){
				done = 0;
				break;
			}
		}
		
		if(done==1){
			break;
		}
	}
	
	presumcount = 0.;
	postsumcount = 0.;
	
	for(a=animals[0]; a--; ){
	//Same as for(a=0; a<animals[0]; a++){ but quicker - INDEXING IS BACKWARDS
		replicatemeans = rgamma(preshapeindt*replicates[0], preindmeans[a]/(preshapeindt*replicates[0]));
		presumcount = presumcount + (((double)rpois(replicatemeans*replicates[0]*sensitivity[0]))*(1/sensitivity[0]));

		replicatemeans = rgamma(postshapeindt*replicates[0], postindmeans[a]/(postshapeindt*replicates[0]));
		postsumcount = postsumcount + (((double)rpois(replicatemeans*replicates[0]*sensitivity[0]))*(1/sensitivity[0]));

	}
	
	//Will give NaN if pre=post=0 or Inf if pre=0 & post>0.  The reductions for these should be 0, so delta =1:
	if(presumcount==0){
		meanred=1;
	}else{
		meanred = postsumcount / presumcount;
	}
	meanreds[set] = meanred;

	if(print[0]){
		Rprintf("%f%% complete\r", set/maxiterations[0]);
	}
}

PutRNGstate(); 

if(print[0]){
	Rprintf("< Finished >\n");
}

}


/* 
dyn.load("power.so")
start <- Sys.time()
out <- .C("poweranalysispopulation", as.numeric(10), as.numeric(1), as.numeric(1), as.integer(100), as.numeric(0.1), as.numeric(0.1), as.numeric(0.1), as.integer(100), as.integer(100000), as.integer(10000), as.integer(3), as.integer(0), as.integer(0))
timestring(start, Sys.time())


power <- function(meanepg=10, g.faeces=3, sensitivity=1/25, animals=10, coeffvarind=0.5, coeffvargroup=0.5, accuracy=0.1, maxiterations=1000000, precision=2, confidence = 0.99, feedback=TRUE, true.sample=FALSE){
dyn.load("power.so")

if(confidence >= 1) stop("Confidence must be < 1")
conf <- (1-confidence)/2
lci <- 0+conf
uci <- 1-conf

start <- Sys.time()
if(true.sample){

}else{
out <- .C("poweranalysispopulation", as.numeric(meanepg), as.numeric(g.faeces), as.numeric(sensitivity), as.integer(animals), as.numeric(coeffvarind), as.numeric(coeffvargroup), as.numeric(accuracy), as.integer(maxiterations), as.integer(precision), as.numeric(lci), as.numeric(uci), as.integer(feedback), as.integer(0), as.integer(0))
time <- timestring(start, Sys.time())
}

lo <- length(out)

return(list(within=out[[lo-1]], without=out[[lo]]-out[[lo-1]], total=out[[lo]], roundedci = round(qbeta(c(0.025, 0.975), out[[lo-1]]+1, (out[[lo]]-out[[lo-1]])+1), precision), ci=qbeta(c(0.025, 0.975), out[[lo-1]]+1, (out[[lo]]-out[[lo-1]])+1), time.taken=time))
}

f <- function(mean){
	return

*/


/*
data.matx <- numeric(length = iterations)
for (set in 1:iterations)  {

individual.means <- rgamma(animals, shape.gp, scale=mean.epg/shape.gp)

sample.means <- replicate(animals, mean(rgamma(g.faeces, shape.in, scale=individual.means/shape.in)))
	
obs.counts <- rpois(animals, sample.means*sensitivity)*(1/sensitivity)
obs.mean <- mean(obs.counts, na.rm=TRUE)

data.matx[set] <- obs.mean
}
ci <- (sum(na.omit(data.matx) > mean.epg*(1+accuracy)) + sum(na.omit(data.matx) < mean.epg*(1-accuracy)))
suppressWarnings(probs <- qbeta(c(0.025, 0.5, 0.975), iterations-ci, ci, lower.tail=FALSE))
probs[ci==0 | ci==iterations] <- NA
return(probs)
}

\end{verbatim}


The R code representing this model with mean.epg equal to the true sample mean is as follows.  In this case, the individual means are corrected so that the true sample mean reflects exactly the population mean input:

\begin{verbatim}
power.analysis.true.sample <- function(mean.epg=NA, g.faeces=NA, sensitivity=NA, animals=NA, coeff.var.ind=NA, coeff.var.group=NA, accuracy=NA, iterations=NA){

shape.gp <- 1/coeff.var.group^2
shape.in <- 1/coeff.var.ind^2
	
data.matx <- numeric(length = iterations)
for (set in 1:iterations)  {

repeat{
	individual.means <- rgamma(animals, shape.gp, scale=mean.epg/shape.gp)
	individual.means <- individual.means+mean.epg-mean(individual.means)
	if(all(individual.means > 0)) break
}

!!!!!!!!!!!better to use true.mean[set] <- mean(individual.means) and use true.mean rather than mean.epg below

sample.means <- rgamma(animals, shape.in*g.faeces, scale=individual.means/(shape.in*g.faeces))

obs.counts <- rpois(animals, sample.means*sensitivity)*(1/sensitivity)
obs.mean <- mean(obs.counts, na.rm=TRUE)

data.matx[set] <- obs.mean
}

ci <- (sum(na.omit(data.matx) > mean.epg*(1+accuracy)) + sum(na.omit(data.matx) < mean.epg*(1-accuracy)))
#return(ci)
suppressWarnings(probs <- qbeta(c(0.025, 0.5, 0.975), iterations-ci, ci, lower.tail=FALSE))
probs[ci==0 | ci==iterations] <- NA
return(probs)
}

*/
