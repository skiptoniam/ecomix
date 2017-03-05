#ifndef rcpPred4_hh
#define rcpPred4_hh

#include<R.h>
#include<Rmath.h>
#include<Rinternals.h>
#include<R_ext/Applic.h>
#include<vector>

#undef length
#include <iostream>

using namespace std;
using std::vector;         // use vector as abbreviation for std::vector
using std::cout;

////////////////////////////////////////////////////////
/////////////	Function Definitions	////////////////////
////////////////////////////////////////////////////////
extern "C" SEXP RCP_predict_C( SEXP Ry, SEXP RX, SEXP RW, SEXP Roffset, SEXP Rwts,
				SEXP RS, SEXP RG, SEXP Rpx, SEXP Rpw, SEXP RnObs, SEXP Rdisty,
				SEXP Ralpha, SEXP Rtau, SEXP Rbeta, SEXP Rgamma, SEXP Rdisps, SEXP Rpowers,
				SEXP Rconc, SEXP Rsd, SEXP RsdGamma, SEXP RdispLocat, SEXP RdispScale,
				SEXP RalphaBoot, SEXP RtauBoot, SEXP RbetaBoot,
				SEXP Rnboot,
				SEXP RptPreds, SEXP RbootPreds, SEXP RoptiDisp);	
//void calcMargFits( double *ptPreds, int bootCount, const vector<double> &allMus, const vector< vector<double> > &allPis, const myData &dat);


#endif
