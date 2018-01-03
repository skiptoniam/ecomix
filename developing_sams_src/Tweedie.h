//Header for Tweedie C++ functions

#ifndef Tweedie_hh
#define Tweedie_hh

#include<vector>
//#include<algorithm>
//#include<time.h>

//#include<iostream>


#include<Rinternals.h>
#include<Rmath.h>
#include<R.h>
//#include<R_ext/Utils.h>



using std::vector;         // use vector as abbreviation for std::vector

//////////////////////////////////////////
//Macros
//////////////////////////////////////////

#define MATREF( i, j, nrow) i + j*nrow	//indices start at 0

//////////////////////////////////////////
//Headers
//////////////////////////////////////////

extern "C" SEXP dTweedie( SEXP, SEXP, SEXP, SEXP, SEXP);
extern "C" SEXP dTweedieDeriv( SEXP, SEXP, SEXP, SEXP);

double dTweedie( double, double, double, double, int);
void dTGLM( vector<double> &outDerivs, vector<double> &X, const double y, const double offs, const vector<double> &beta, const double phi, const double p);

double logWfun( double, double, double, double, double);
double logWderivApprox( double, double, double);
double findW( double, double, double, double, double, double, double);
double findjMax( double y3, double muN3, double muZ3, double alpha3, double beta3, double z13, double z23, double &logWmax3);
void findlogWjs( double, double, double, double, double, double, double, double, double &, double &, double, vector<double> &);

void findWDeriv( double y2, double muN2, double muZ2, double alpha2, double beta2, double z12, double z22, vector<double> &jmax2, vector<double> &jlims2, vector<double> &derivsW);
void findjMaxDerivs( double y3, double muN3, double muZ3, double alpha3, double beta3, double z13, double z23, vector<double> &logMaxs3, vector<double> &jmax3);
void ddjOFlogdWjdLambda( const vector<double> &jnr4, double z14, double alpha4, vector<double> &deriv4);
/*void ddjOFdWjdAlpha( const vector<double> &jnr5, double y5, double muN5, double muZ5, double z15, double alpha5, double beta5, vector<double> &deriv5);*/
void findLogWjsForDeriv( double y4, double muN4, double muZ4, double alpha4, double beta4, double z14, double z24, const vector<double> &jmax4, vector<double> &jlims4, const vector<double> &logMaxs4, vector<double> &logWjs4, vector<double> &logdlambda4, vector<double> &logdmuZ4, vector<double> &logdalpha4, vector<double> &signalpha4);
inline bool checkTol( const vector<double> &maxes, const double &currLogW, const double &currdlambda, const double &currdmuZ, const double &currdalpha, const double &eps1, const double &expeps1);
void findEachDeriv( const double y8, const double muN8, const double muZ8, const double alpha8, const double beta8, const double z18, const double z228, const vector<double> &logWjs8, const vector<double> &logdlambda8, const vector<double> &logdmuZ8, const vector<double> &logdalpha8, const vector<double> &signalpha8, const vector<double> &logMaxs8, vector<double> &derivsW8);

//headers for mixture derivatives

extern "C" SEXP mixUtilTauBetaAlpha( SEXP dtau, SEXP dbeta, SEXP dlalpha, SEXP maxLocs, SEXP mults, SEXP Rn, SEXP RS, SEXP RG, SEXP Rptau, SEXP Rpbeta, SEXP Rzero);
extern "C" SEXP mixUtilDelta( SEXP ddelta, SEXP maxLocs, SEXP mults, SEXP Rn, SEXP RS, SEXP RG);
extern "C" SEXP mixUtilKappaConstPi( SEXP dpi, SEXP maxLocs, SEXP mults, SEXP Rn, SEXP RS, SEXP RG);
extern "C" SEXP mixUtilKappa( SEXP dpi, SEXP maxLocs, SEXP mults, SEXP RX, SEXP Rn, SEXP RS, SEXP RG, SEXP Rnvar);

///////////////////////////////////////////
//Classes
///////////////////////////////////////////

#endif
