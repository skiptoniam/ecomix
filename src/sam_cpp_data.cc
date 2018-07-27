#include"sam_cpp.h"

sam_data::sam_data(){};
sam_data::~sam_data(){};

void sam_data::setVals( SEXP &Ry, SEXP &RX, SEXP &Roffset, SEXP &Rspp_wts, SEXP &Rsite_spp_wts, SEXP &Ry_not_na,
                             SEXP &RS, SEXP &RG, SEXP &Rp, SEXP &RnObs, SEXP &Rdisty, SEXP &RoptiDisp){

	nS = *(INTEGER( RS));
	nG = *(INTEGER( RG));
	nP = *(INTEGER( Rp));
	nObs = *(INTEGER( RnObs));
	disty = *(INTEGER( Rdisty));
	optiDisp = *(INTEGER( RoptiDisp));
	NAnum = -99999;
	
	y = REAL( Ry);
	X = REAL( RX);
	offset = REAL( Roffset);
	spp_wts = REAL( Rspp_wts); // this is for the bayesian bootstrap
	site_spp_wts = REAL( Rsite_spp_wts); //this is for the ippm
	y_not_na = INTEGER(Ry_not_na);
	

}

bool sam_data::doOptiDisp() const{
	if( optiDisp == 1)
		return( TRUE);
	return( FALSE);
}

//"bernoulli" = 1,"poisson" = 2,"negative_binomial" = 3,"tweedie" = 4,'ippm' = 6
//"bernoulli_sp" = 1, "poisson" = 2, "ippm" = 3, "negative_binomial"=4, //Ignore these for now "tweedie","gaussian"

bool sam_data::isDispersion() const{ // currently just check if negative binomial.
	if((disty == 4) | (disty == 5) | (disty == 6) )
		return( true);
	return( false);
}



void sam_data::printVals( int printX=0, int printy = 0){
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
