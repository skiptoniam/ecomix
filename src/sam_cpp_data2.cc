#include"sam_cpp2.h"

sam_data::sam_data(){};
sam_data::~sam_data(){};


//slowly adding the bias formula design matrix. 

void sam_data::setVals( SEXP &Ry, SEXP &RX, SEXP &RW, //SEXP &RU, 
						SEXP &Roffset, SEXP &Rspp_wts, 
						SEXP &Rsite_spp_wts, SEXP &Ry_not_na, SEXP &Rbinsize,
						SEXP &RS, SEXP &RG, SEXP &Rpx, SEXP &Rpw, //SEXP &Rpu,
						SEXP &RnObs, SEXP &Rdisty, SEXP &RoptiDisp, SEXP &RoptiPart){

	nS = *(INTEGER( RS));
	nG = *(INTEGER( RG));
	nPX = *(INTEGER( Rpx));
	nPW = *(INTEGER( Rpw));
	//nPU = *(INTEGER( Rpu));
	nObs = *(INTEGER( RnObs));
	disty = *(INTEGER( Rdisty));
	optiDisp = *(INTEGER( RoptiDisp));
	optiPart = *(INTEGER( RoptiPart));
	NAnum = -999999;

	y = REAL( Ry);
	X = REAL( RX);
	W = REAL( RW);
	//U = REAL( RU);
	offset = REAL( Roffset);
	spp_wts = REAL( Rspp_wts); // this is for the bayesian bootstrap
	site_spp_wts = REAL( Rsite_spp_wts); //this is for the ippm
	binsize = REAL(Rbinsize);
	y_not_na = INTEGER(Ry_not_na);


}

// bool sam_data::doOptiPart() const{
// 	if( optiPart == 1)
// 		return( TRUE);
// 	return( FALSE);
// }

bool sam_data::isDispersion() const{ // currently just check if negative binomial.
	if((disty == 4) | (disty == 5) | (disty == 6) )
		return( true);
	return( false);
}

bool sam_data::doOptiDisp() const{
	if( optiDisp == 1)
		return( TRUE);
	return( FALSE);
}

// bool sam_data::isPartial() const{ // currently just check if negative binomial.
// 	if(optiPart==1)
// 		return( true);
// 	return( false);
// }


// void sam_data::printVals( int printX=0, int printW=0, int printy = 0){
// 	if( printy == 1)
// 		for( int i=0; i<nObs; i++){
// 			for( int j=0; j<nS; j++)
// 				Rprintf( "%3.2f\t", y[MATREF2D(i,j,nObs)]);
// 				Rprintf( "\n");
// 		}
//
// 	if( printX == 1)
// 		for( int i=0; i<nObs; i++){
// 			for( int j=0; j<nPX; j++)
// 				Rprintf( "%3.2f\t", X[MATREF2D(i,j,nObs)]);
// 				Rprintf( "\n");
// 		}
//
// 	if( printW == 1)
// 		for( int i=0; i<nObs; i++){
// 			for( int j=0; j<nPW; j++)
// 				Rprintf( "%3.2f\t", W[MATREF2D(i,j,nObs)]);
// 				Rprintf( "\n");
// 		}
// }
