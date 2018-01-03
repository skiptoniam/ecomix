#include"sam.h"

void spmix_data::setVals( const SEXP &Ry, const SEXP &RX, const SEXP &Roffset, const SEXP &Rwts, const SEXP &Ry_not_na, const SEXP &RS, const SEXP &RG, const SEXP &Rp, const SEXP &RnObs, const SEXP &Rdisty, const SEXP &RdoOptiDisp){

//	double *tmpD;

	nS = *(INTEGER( RS));
	nG = *(INTEGER( RG));
	np = *(INTEGER( Rp));
	nObs = *(INTEGER( RnObs));
	disty = *(INTEGER( Rdisty));
	optiDisp = *(INTEGER( RdoOptiDisp));
	NAnum = -999999;

/*	tmpD = REAL( Ry); y.assign( tmpD, tmpD + LENGTH( Ry));
	tmpD = REAL( RX); X.assign( tmpD, tmpD + LENGTH( RX));*/
	y = REAL( Ry);
	X = REAL( RX);
	offset = REAL( Roffset);
	wts = REAL( Rwts);
	y_is_not_na = INTEGER(Ry_not_na));		
}

bool spmix_data::doOptiDisp() const{
	if( optiDisp == 1)
		return( TRUE);
	return( FALSE);
}

bool spmix_data::isDispersion() const{
	if( (disty == 2) | (disty == 3))
		return( true);
	return( false);
}

void spmix_data::printVals( int printX=0, int printy = 0){
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
}
