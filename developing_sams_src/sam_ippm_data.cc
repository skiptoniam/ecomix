#include"sam_ippm.h"

void sam_ippm_data::setVals( const SEXP &Ry, const SEXP &RX, const SEXP &Roffset, const SEXP &Rwts, const SEXP &Ry_not_na, const SEXP &RS, const SEXP &RG, const SEXP &Rp, const SEXP &RnObs){//, const SEXP &Rdisty, const SEXP &RdoOptiDisp){

	nS = *(INTEGER( RS));
	nG = *(INTEGER( RG));
	np = *(INTEGER( Rp));
	nObs = *(INTEGER( RnObs));
	NAnum = -999999;
	y = REAL( Ry);
	X = REAL( RX);
	offset = REAL( Roffset);
	wts = REAL( Rwts);
	y_not_na = INTEGER(Ry_not_na);	

}

void sam_ippm_data::printVals( int printX=0, int printy = 0){
	if( printy == 1)
		for( int i=0; i<nObs; i++){
			for( int j=0; j<nS; j++)
				Rprintf( "%3.2f\t", y[MATREF2D(i,j,nObs)]);
				Rprintf( "\n");
		}
	
	if( printX == 1)
		for( int i=0; i<nObs; i++){
			for( int j=0; j<np; j++)
				Rprintf( "%3.2f\t", X[MATREF2D(i,j,nObs)]);
				Rprintf( "\n");
		}
}
