#include"sam.h"

spmix_params::spmix_params(){};
spmix_params::~spmix_params(){};

void spmix_params::setVals( const spmix_data &dat, SEXP &Ralpha, SEXP &Rbeta, SEXP &Rtau, SEXP &Rdisps){ //, SEXP &Rpowers, SEXP &Rconc, SEXP &Rsd, SEXP &RsdGamma, SEXP &RdispLocat, SEXP &RdispScale)
{
//	double *tmpD;

	Alpha = REAL( Ralpha);
	Beta = REAL( Rbeta);
	Tau = REAL( Rtau);
	Disp = REAL( Rdisps);
	//Power = REAL( Rpowers);
	//conc = *(REAL( Rconc));
	//sd = *(REAL( Rsd));
	//sdGamma = *(REAL( RsdGamma));
	//dispLocat = *(REAL( RdispLocat));
	//dispScale = *(REAL( RdispScale));

	nalpha = dat.nS;
	nbeta = (dat.nG-1)*dat.np;
	ntau = (dat.nG-1)*dat.nS;
	if( dat.isDispersion())
		ndisp = dat.nS;
	else
		ndisp = 0;
	nTot = nalpha + ntau + nbeta + ngamma + ndisp;
}

void spmix_params::getArray(double *parArr, const spmix_data &dat) const
{
	int kount=0;
	for( int i=0; i<dat.nS; i++){
		parArr[kount] = Alpha[i];
		kount++;
	}
	for( int i=0; i<((dat.nG-1)*dat.np); i++){
		parArr[kount] = Beta[i];
		kount++;
	}
	for( int i=0; i<((dat.nG-1)*dat.nS); i++){
		parArr[kount] = Tau[i];
		kount++;
	}
	if( dat.isDispersion())
		for( int i=0; i<dat.nS; i++){
			parArr[kount] = Disp[i];
			kount++;
		}
}

void spmix_params::update( double *parArr, const spmix_data &dat)
{
	int kount=0;
	for( int i=0; i<dat.nS; i++){
		Alpha[i] = parArr[kount];
		kount++;
	}
	for( int i=0; i<((dat.nG-1)*dat.np); i++){
		Beta[i] = parArr[kount];
		kount++;
	}
	for( int i=0; i<((dat.nG-1)*dat.nS); i++){
		Tau[i] = parArr[kount];
		kount++;
	}
	if( dat.isDispersion() & dat.doOptiDisp())
		for( int i=0; i<dat.nS; i++){
			Disp[i] = parArr[kount];
			kount++;
		}
}

void spmix_params::getAllTaus( vector<double> &newTau, const spmix_data &dat) const
{
	double su;

	newTau.assign(dat.nG*dat.nS, dat.NAnum);
	//calculate sum-to-zero taus
	for( int s=0; s<dat.nS; s++){
		su = 0.0;
		for( int g=0; g<(dat.nG-1); g++){
			newTau.at( MATREF2D(g,s,dat.nG)) = Tau[MATREF2D(g,s,(dat.nG-1))];
			su += Tau[MATREF2D(g,s,(dat.nG-1))];
		}
		newTau.at( MATREF2D((dat.nG-1),s,dat.nG)) = -su;
	}

}

void spmix_params::printParms( const spmix_data &dat){
	Rprintf( "ALPHA:\n");
	for( int i=0; i<dat.nS; i++)
		Rprintf( "%3.2f\t", Alpha[i]);
	Rprintf( "\n");
	Rprintf( "BETA:\n");
	for( int g=0; g<(dat.nG-1); g++){
		for( int i=0; i<dat.np; i++)
			Rprintf( "%3.2f\t", Beta[MATREF2D(g,i,(dat.nG-1))]);
		Rprintf( "\n");
	}
	Rprintf( "TAU:\n");
	for( int g=0; g<(dat.nG-1); g++){
		for( int i=0; i<dat.nS; i++)
			Rprintf( "%3.2f\t", Tau[MATREF2D(g,i,(dat.nG-1))]);
		Rprintf( "\n");
	}
	if( dat.isDispersion()==true){
		Rprintf("DISPERSION:\n");
		for( int i=0; i<dat.nS; i++)
			Rprintf( "%3.2f\t", Disp[i]);
		Rprintf( "\n");
	}
		
		
}
