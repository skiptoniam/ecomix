#include"rcp4.h"

myData::myData(){};
myData::~myData(){};

void myData::setVals( const SEXP &Ry, const SEXP &RX, const SEXP &RW, const SEXP &Roffset, const SEXP &RS, const SEXP &RG, const SEXP &Rp, const SEXP &Rpw, const SEXP &RnObs, const SEXP &Rdisty, const SEXP &RdoOptiDisp, const SEXP &Rwts){

//	double *tmpD;

	nS = *(INTEGER( RS));
	nG = *(INTEGER( RG));
	np = *(INTEGER( Rp));
	npw = *(INTEGER( Rpw));
	nObs = *(INTEGER( RnObs));
	disty = *(INTEGER( Rdisty));
	optiDisp = *(INTEGER( RdoOptiDisp));
	NAnum = -999999;

/*	tmpD = REAL( Ry); y.assign( tmpD, tmpD + LENGTH( Ry));
	tmpD = REAL( RX); X.assign( tmpD, tmpD + LENGTH( RX));*/
	y = REAL( Ry);
	X = REAL( RX);
	W = REAL( RW);
	offset = REAL( Roffset);
	wts = REAL( Rwts);
}

bool myData::doOptiDisp() const{
	if( optiDisp == 1)
		return( TRUE);
	return( FALSE);
}

bool myData::isDispersion() const{
	if( (disty == 3) | (disty == 4) | (disty == 5))
		return( true);
	return( false);
}

void myData::printVals( int printX=0, int printW = 0, int printy = 0){
	if( printy == 1)
		for( int i=0; i<nObs; i++){
			for( int j=0; j<nS; j++)
				Rprintf( "%3.2f\t", y[MATREF(i,j,nObs)]);
			Rprintf( "\n");
		}
	
	if( printX == 1)
		for( int i=0; i<nObs; i++){
			for( int j=0; j<np; j++)
				Rprintf( "%3.2f\t", X[MATREF(i,j,nObs)]);
			Rprintf( "\n");
		}

	if( printW == 1)
		for( int i=0; i<nObs; i++){
			for( int j=0; j<npw; j++)
				Rprintf( "%3.2f\t", W[MATREF(i,j,nObs)]);
			Rprintf( "\n");
		}
}
