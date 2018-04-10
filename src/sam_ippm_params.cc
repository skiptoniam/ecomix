#include"sam_ippm.h"

sam_ippm_params::sam_ippm_params(){};
sam_ippm_params::~sam_ippm_params(){};

void sam_ippm_params::setVals( const sam_ippm_data &dat, SEXP &Ralpha, SEXP &Rbeta, SEXP &Rpi){
//	double *tmpD;

	Alpha = REAL( Ralpha);
	Beta = REAL( Rbeta);
	Pi = REAL( Rpi);

	nalpha = dat.nS;
	nbeta = dat.nG*dat.nP;
	npi = (dat.nG-1);

	nTot = nalpha + nbeta + npi; 
}

void sam_ippm_params::getArray(double *parArr, const sam_ippm_data &dat){
	int kount=0;
	for( int i=0; i<dat.nS; i++){
		parArr[kount] = Alpha[i];
		kount++;
	}
	for( int i=0; i<((dat.nG*dat.nP)); i++){
		parArr[kount] = Beta[i];
		kount++;
	}
	for( int i=0; i<((dat.nG-1)); i++){ 
		parArr[kount] = Pi[i];
		//Rprintf("%d\n", i);
		kount++;
	}

}

void sam_ippm_params::update( double *parArr, const sam_ippm_data &dat){
	int kount=0;
	for( int i=0; i<dat.nS; i++){
		Alpha[i] = parArr[kount];
		kount++;
	}
	for( int i=0; i<((dat.nG*dat.nP)); i++){
		Beta[i] = parArr[kount];
		kount++;
	}
	for( int i=0; i<((dat.nG-1)); i++){
		Pi[i] = parArr[kount];
		kount++;
	}
}

void sam_ippm_params::printParms( const sam_ippm_data &dat){
	
	Rprintf( "ALPHA:\n");
	for( int i=0; i<dat.nS; i++)
		Rprintf( "%3.2f\t", Alpha[i]);
		Rprintf( "\n");
		Rprintf( "BETA:\n");
	for( int g=0; g<(dat.nG); g++){
		for( int i=0; i<dat.nP; i++)
			Rprintf( "%3.2f\t", Beta[MATREF2D(g,i,(dat.nG))]);
			Rprintf( "\n");
	}
	Rprintf( "PI:\n");
	for( int g=0; g<(dat.nG-1); g++){
		Rprintf( "%3.2f\t", Pi[g]);
		Rprintf( "\n");
	}
		
		
}
