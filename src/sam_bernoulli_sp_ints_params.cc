#include"sam_bernoulli_sp_ints.h"

sam_bernoulli_sp_ints_params::sam_bernoulli_sp_ints_params(){};
sam_bernoulli_sp_ints_params::~sam_bernoulli_sp_ints_params(){};

void sam_bernoulli_sp_ints_params::setVals( const sam_bernoulli_sp_ints_data &dat, SEXP &Ralpha, SEXP &Rbeta, SEXP &Reta){
//	double *tmpD;

	Alpha = REAL( Ralpha);
	Beta = REAL( Rbeta);
	Eta = REAL( Reta);

	nalpha = dat.nS;
	nbeta = dat.nG*dat.nP;
	npi = (dat.nG-1);

	nTot = nalpha + nbeta + npi; 
}

void sam_bernoulli_sp_ints_params::getArray(double *parArr, const sam_bernoulli_sp_ints_data &dat){
	int kount=0;
	for( int i=0; i<((dat.nS)); i++){
		parArr[kount] = Alpha[i];
		kount++;
	}
	for( int i=0; i<((dat.nG*dat.nP)); i++){
		parArr[kount] = Beta[i];
		kount++;
	}
	for( int i=0; i<((dat.nG-1)); i++){ 
		parArr[kount] = Eta[i];
		kount++;
	}

}

void sam_bernoulli_sp_ints_params::update( double *parArr, const sam_bernoulli_sp_ints_data &dat){
	int kount=0;
	for( int i=0; i<((dat.nS)); i++){
		Alpha[i] = parArr[kount];
		kount++;
	}
	for( int i=0; i<((dat.nG*dat.nP)); i++){
		Beta[i] = parArr[kount];
		kount++;
	}
	for( int i=0; i<((dat.nG-1)); i++){
		Eta[i] = parArr[kount];
		kount++;
	}
}

void sam_bernoulli_sp_ints_params::printParms( const sam_bernoulli_sp_ints_data &dat){
	
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
	Rprintf( "PI (transformed Pi):\n");
	for( int g=0; g<(dat.nG-1); g++){
		Rprintf( "%3.2f\t", Eta[g]);
		Rprintf( "\n");
	}
		
		
}
