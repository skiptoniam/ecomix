#include"sam_cpp2.h"

sam_params::sam_params(){};
sam_params::~sam_params(){};

void sam_params::setVals(const sam_data &dat, SEXP &Ralpha, SEXP &Rbeta,
						SEXP &Reta, SEXP &Rtheta, &Rgamma){
//	double *tmpD;

	Alpha = REAL( Ralpha);
	Beta = REAL( Rbeta);
	Gamma = REAL( Rgamma);
	Eta = REAL( Reta);
	Theta = REAL( Rtheta);

	nalpha = dat.nS;
	nbeta = dat.nG*dat.nPX;
	ngamma = dat.nS*dat.nPW;
	neta = (dat.nG-1);
	if(dat.isDispersion())
		ntheta = dat.nS;
	else
		ntheta = 0;

	nTot = nalpha + nbeta + ngamma + neta + ntheta; 
}

void sam_params::getArray(double *parArr, const sam_data &dat){
	int kount=0;
	for( int i=0; i<((dat.nS)); i++){
		parArr[kount] = Alpha[i];
		kount++;
	}
	for( int i=0; i<((dat.nG*dat.nPX)); i++){
		parArr[kount] = Beta[i];
		kount++;
	}
	
	for( int i=0; i<((dat.nS*dat.nPW)); i++){
		parArr[kount] = Gamma[i];
		kount++;
	}
	for( int i=0; i<((dat.nG-1)); i++){ 
		parArr[kount] = Eta[i];
		kount++;
	}
	if( dat.isDispersion())
		for( int i=0; i<dat.nS; i++){
				parArr[kount] = Theta[i];
				kount++;
		}

}

void sam_params::update( double *parArr, const sam_data &dat){
	int kount=0;
	for( int i=0; i<((dat.nS)); i++){
		Alpha[i] = parArr[kount];
		kount++;
	}
	for( int i=0; i<((dat.nG*dat.nPX)); i++){
		Beta[i] = parArr[kount];
		kount++;
	}
	for( int i=0; i<((dat.nS*dat.nPW)); i++){
		Gamma[i] = parArr[kount];
		kount++;
	}	
	for( int i=0; i<((dat.nG-1)); i++){
		Eta[i] = parArr[kount];
		kount++;
	}
	if( dat.isDispersion() & dat.doOptiDisp())
		for( int i=0; i<dat.nS; i++){
			Theta[i] = parArr[kount];
			kount++;
		}
		
}

void sam_params::printParms( const sam_data &dat){
	
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
	Rprintf( "GAMMA:\n");
	for( int g=0; g<(dat.nS); g++){
		for( int i=0; i<dat.nPW; i++)
			Rprintf( "%3.2f\t", Gamma[MATREF2D(s,i,(dat.nS))]);
			Rprintf( "\n");
	}
	Rprintf( "PI (transformed Pi):\n");
	for( int g=0; g<(dat.nG-1); g++){
		Rprintf( "%3.2f\t", Eta[g]);
		Rprintf( "\n");
	}
	if( dat.isDispersion()==true){
		Rprintf("DISPERSION:\n");
		for( int i=0; i<dat.nS; i++)
		Rprintf( "%3.2f\t", Theta[i]);
		Rprintf( "\n");
	}
		
}
