#include"sam_cpp.h"

sam_params::sam_params(){};
sam_params::~sam_params(){};

void sam_params::setVals( const sam_data &dat, SEXP &Ralpha, SEXP &Rbeta, SEXP &Reta, SEXP &Rdisp){
//	double *tmpD;

	Alpha = REAL( Ralpha);
	Beta = REAL( Rbeta);
	Eta = REAL( Reta);
	Disp = REAL( Rdisp);

	nalpha = dat.nS;
	nbeta = dat.nG*dat.nP;
	neta = (dat.nG-1);
	if(dat.isDispersion())
		ndisp = dat.nS;
	else
		ndisp = 0;

	nTot = nalpha + nbeta + neta + ndisp; 
}

void sam_params::getArray(double *parArr, const sam_data &dat){
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
	if( dat.isDispersion())
		for( int i=0; i<dat.nS; i++){
				parArr[kount] = Disp[i];
				kount++;
		}

}

void sam_params::update( double *parArr, const sam_data &dat){
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
	if( dat.isDispersion() & dat.doOptiDisp())
		for( int i=0; i<dat.nS; i++){
			Disp[i] = parArr[kount];
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
	Rprintf( "PI (transformed Pi):\n");
	for( int g=0; g<(dat.nG-1); g++){
		Rprintf( "%3.2f\t", Eta[g]);
		Rprintf( "\n");
	}
	if( dat.isDispersion()==true){
		Rprintf("DISPERSION:\n");
		for( int i=0; i<dat.nS; i++)
		Rprintf( "%3.2f\t", Disp[i]);
		Rprintf( "\n");
	}
		
		
}
