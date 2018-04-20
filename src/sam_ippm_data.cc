#include"sam_ippm.h"

sam_ippm_data::sam_ippm_data(){};
sam_ippm_data::~sam_ippm_data(){};

void sam_ippm_data::setVals( SEXP &Ry, SEXP &RX, SEXP &Roffset, SEXP &Rwts, SEXP &Ry_not_na,
                             SEXP &RS, SEXP &RG, SEXP &Rp, SEXP &RnObs){

	nS = *(INTEGER( RS));
	nG = *(INTEGER( RG));
	nP = *(INTEGER( Rp));
	nObs = *(INTEGER( RnObs));
	NAnum = -99999;
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
			for( int j=0; j<nP; j++)
				Rprintf( "%3.2f\t", X[MATREF2D(i,j,nObs)]);
				Rprintf( "\n");
		}
}
