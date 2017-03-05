#include"rcp4.h"

myParms::myParms(){};
myParms::~myParms(){};

void myParms::setVals( const myData &dat, SEXP &Ralpha, SEXP &Rbeta, SEXP &Rtau, SEXP &Rgamma, SEXP &Rdisps, SEXP &Rpowers, SEXP &Rconc, SEXP &Rsd, SEXP &RsdGamma, SEXP &RdispLocat, SEXP &RdispScale)
{
//	double *tmpD;

	Alpha = REAL( Ralpha);
	Tau = REAL( Rtau);
	Beta = REAL( Rbeta);
	Gamma = REAL( Rgamma);
	Disp = REAL( Rdisps);
	Power = REAL( Rpowers);
	conc = *(REAL( Rconc));
	sd = *(REAL( Rsd));
	sdGamma = *(REAL( RsdGamma));
	dispLocat = *(REAL( RdispLocat));
	dispScale = *(REAL( RdispScale));

	nalpha = dat.nS;
	ntau = (dat.nG-1)*dat.nS;
	nbeta = (dat.nG-1)*dat.np;
	ngamma = dat.nS*dat.npw;
	if( dat.isDispersion())
		ndisp = dat.nS;
	else
		ndisp = 0;
	nTot = nalpha + ntau + nbeta + ngamma + ndisp;

}

void myParms::getArray(double *parArr, const myData &dat) const
{
	int kount=0;
	for( int i=0; i<dat.nS; i++){
		parArr[kount] = Alpha[i];
		kount++;
	}
	for( int i=0; i<((dat.nG-1)*dat.nS); i++){
		parArr[kount] = Tau[i];
		kount++;
	}
	for( int i=0; i<((dat.nG-1)*dat.np); i++){
		parArr[kount] = Beta[i];
		kount++;
	}
	for( int i=0; i<(dat.nS*dat.npw); i++){
		parArr[kount] = Gamma[i];
		kount++;
	}
	if( dat.isDispersion())
		for( int i=0; i<dat.nS; i++){
			parArr[kount] = Disp[i];
			kount++;
		}
}

void myParms::update( double *parArr, const myData &dat)
{
	int kount=0;
	for( int i=0; i<dat.nS; i++){
		Alpha[i] = parArr[kount];
		kount++;
	}
	for( int i=0; i<((dat.nG-1)*dat.nS); i++){
		Tau[i] = parArr[kount];
		kount++;
	}
	for( int i=0; i<((dat.nG-1)*dat.np); i++){
		Beta[i] = parArr[kount];
		kount++;
	}
	for( int i=0; i<(dat.nS*dat.npw); i++){
		Gamma[i] = parArr[kount];
		kount++;
	}
	if( dat.isDispersion() & dat.doOptiDisp())
		for( int i=0; i<dat.nS; i++){
			Disp[i] = parArr[kount];
			kount++;
		}
}

void myParms::getAllTaus( vector<double> &newTau, const myData &dat) const
{
	double su;

	newTau.assign(dat.nG*dat.nS, dat.NAnum);
	//calculate sum-to-zero taus
	for( int s=0; s<dat.nS; s++){
		su = 0.0;
		for( int g=0; g<(dat.nG-1); g++){
			newTau.at( MATREF(g,s,dat.nG)) = Tau[MATREF(g,s,(dat.nG-1))];
			su += Tau[MATREF(g,s,(dat.nG-1))];
		}
		newTau.at( MATREF((dat.nG-1),s,dat.nG)) = -su;
	}

}

void myParms::printParms( const myData &dat){
	Rprintf( "ALPHA:\n");
	for( int i=0; i<dat.nS; i++)
		Rprintf( "%3.2f\t", Alpha[i]);
	Rprintf( "\n");
	Rprintf( "TAU:\n");
	for( int g=0; g<(dat.nG-1); g++){
		for( int i=0; i<dat.nS; i++)
			Rprintf( "%3.2f\t", Tau[MATREF(g,i,(dat.nG-1))]);
		Rprintf( "\n");
	}
	Rprintf( "BETA:\n");
	for( int g=0; g<(dat.nG-1); g++){
		for( int i=0; i<dat.np; i++)
			Rprintf( "%3.2f\t", Beta[MATREF(g,i,(dat.nG-1))]);
		Rprintf( "\n");
	}
	if( dat.npw > 0){
		Rprintf( "GAMMA:\n");
		for( int g=0; g<dat.nS; g++){
			for( int i=0; i<dat.npw; i++)
				Rprintf( "%3.2f\t", Gamma[MATREF(g,i,dat.nS)]);
			Rprintf( "\n");
		}
	}
	if( dat.isDispersion()==true){
		Rprintf("DISPERSION:\n");
		for( int i=0; i<dat.nS; i++)
			Rprintf( "%3.2f\t", Disp[i]);
		Rprintf( "\n");
	}
		
		
}
