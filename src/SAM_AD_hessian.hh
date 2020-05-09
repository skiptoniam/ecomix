#ifndef SAM_AD_hessian_hh
#define SAM_AD_hessian_hh

#include <R_ext/Applic.h>
#include <vector>
#include <iostream>
#include <cppad/cppad.hpp> 
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>

#undef length

using CppAD::AD;           // use AD as abbreviation for CppAD::AD
using CppAD::Value;
using std::vector;         // use vector as abbreviation for std::vector

////////////////////////////////////////////////////////
/////////////	Inline Function Definitions	//////////////
////////////////////////////////////////////////////////

#define MATREF( i, j, nrow) i + j*nrow	//indices start at 0

////////////////////////////////////////////////////////
/////////////	Function Definitions	////////////////////
////////////////////////////////////////////////////////

AD<double> lgammaAD( AD<double> xx);
AD<double> ldnegbin(double yi, AD<double> mui, AD<double> phii);
void invMultLogit( vector< AD<double> > &pi2, const vector< AD<double> > alpha2, int nG2);
AD<double> ldbeta( AD<double> theta1, double *thetaRange2, double penParm2);
AD<double> calcLogl( vector< AD<double> > allPars1, double *y1, double *offy1, double *XnoMix1, double *XMix1, int nS1, int nG1, int nn1, int npNoMix1, int npMix1, double *thetaRange1, double penParm1);
extern "C" SEXP calcDerHess( SEXP Ry, SEXP Roffy, SEXP RallPars, SEXP RXnoMix, SEXP RXMix, SEXP RS, SEXP RG, SEXP Rn, SEXP RnpNoMix, SEXP RnpMix, SEXP RthetaRange, SEXP RpenParm);

#endif
