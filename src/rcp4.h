#ifndef rcp4_hh
#define rcp4_hh

#include<R.h>
#include<Rmath.h>
#include<Rinternals.h>
#include<R_ext/Applic.h>
#include<vector>
#include"Tweedie.h"

#undef length
#include <iostream>

using namespace std;
using std::vector;         // use vector as abbreviation for std::vector
using std::cout;

////////////////////////////////////////////////////////
/////////////	Inline Function Definitions	//////////////
////////////////////////////////////////////////////////

#define MATREF( i, j, nrow) i + j*nrow	//indices start at 0
#define MATREF3D( i, j, k, nrow, ncol) i + j*nrow + k*nrow*ncol	//indices start at 0 and nrow (ncol) is the number of rows (cols) of a single face of 3D array.

////////////////////////////////////////////////////////
////////////	Class Definitions	////////////////////////
////////////////////////////////////////////////////////

class myData
{
	public:
		myData();
		~myData();
		void setVals( const SEXP &Ry, const SEXP &RX, const SEXP &RW, const SEXP &Roffset, const SEXP &RS, const SEXP &RG, const SEXP &Rp, const SEXP &Rpw, const SEXP &RnObs, const SEXP &Rdisty, const SEXP &RoptiDisp, const SEXP &RlogWts);
		bool isDispersion() const;		
		bool doOptiDisp() const;
		void printVals( int printX, int printW, int printy);

		int np, 		//the number of parameters in each of the (G-1) habitat lps
			npw,		//the number of parameters in the spp level model
			nG,			//the number of habitats
			nS, 		//the number of species
			nObs,		//the number of observations
			disty,		//the distrbution code
			optiDisp,	//should the dispersions be optimised
			NAnum;	//a common number to insert for NAs

		double 	*X, //the design matrix in vector form (nObs x np)
				*W,	//the design matrix for the spp model, in vector form (nObs x npw)
				*y,	//the outcome matrix, in vector form (nObs x nS)
				*offset, //the offset vector (length nObs)
				*wts;  //the wts for the logl (typically all zero and of length nObs).
};

class myParms
{
	public:
		myParms();
		~myParms();
		void setVals( const myData &dat, SEXP &Ralpha, SEXP &Rbeta, SEXP &Rtau, SEXP &Rgamma, SEXP &Rdisps, SEXP &Rpowers, SEXP &Rconc, SEXP &Rsd, SEXP &RsdGamma, SEXP &RdispLocat, SEXP &RdispScale);
		void getArray(double *parArr, const myData &dat) const;
		void update( double *parArr, const myData &dat);
		void getAllTaus( vector<double> &newTau, const myData &dat) const;
		void printParms( const myData &dat);

//		void updatePars( double *newPar);

		double 	*Alpha, //the species' prevalences (Gx1)==(1xG)
						*Tau, 	//the species' habitat adjustments ((G-1)xS) -- only free params!
						*Beta,	//the habitats' free covariate params ((G-1)xp)
						*Gamma,	//the spp free covariate params (S x pw)
						*Disp,	//the spp dispersion parameters (Sx1) == (1xS)
						*Power,	//the spp power parameters (for Tweedie) (Sx1)==(1xS)
						conc,		//the concentration for the dirichlet penalty for pi
						sd,			//the sd for the normal penalty for tau
						sdGamma,	//the sd for the normal penalty on gammas
						dispLocat,	//the location parameter for the pnealty on the dispersion parameters	(default is 0.0)
						dispScale;	//the scale parameter for the pnealty on the dispersion parameters	(default is 10)
		int nalpha, ntau, nbeta, ngamma, ndisp, nTot;
//		int *myMask;
};

class myDerivs
{
	public:
		myDerivs();
		~myDerivs();
		void setVals( const myData &dat, SEXP &RderivsAlpha, SEXP &RderivsTau, SEXP &RderivsBeta, SEXP &RderivsGamma, SEXP &RderivsDisp, SEXP &RgetScores, SEXP &Rscores);
		void zeroDerivs( const myData &dat);
		void updateDerivs( const myData &dat, const vector<double> &alphaDerivsI, const vector<double> &tauDerivsI, const vector<double> &betaDerivsI, const vector<double> &gammaDerivsI, const vector<double> &dispDerivsI, const int &i);
		void update( double *grArr, const myData &dat);
		void getArray( double *grArr, const myData &dat);

		int getScoreFlag;	//Should the scores be calculated for empirical information
		double 	*Alpha, //the derivatives of logl w.r.t. alpha
				*Tau, 	//the derivatives of logl w.r.t. tau
				*Beta,	//the derivatives of logl w.r.t. beta
				*Gamma,	//the derivatives of logl w.r.t. gamma
				*Disp,	//the derivative of logl w.r.t. dispersions
				*Scores;	//the score contribution for each site (for empirical information)
};

class myOptContr
{
	public:
		myOptContr();
		~myOptContr();
		void setVals( const SEXP &Rmaxit, const SEXP &Rtrace, const SEXP &RnReport, const SEXP &Rabstol, const SEXP &Rreltol, SEXP &Rconv);

		int maxitQN, traceQN, nReport, fnKount, grKount, ifail, *conv;
		double abstol, reltol, denomEps;
};

class myFits
{
	public: 
		myFits();
		~myFits();
		void initialise( const int &nObs, const int &nG, const int &nS, const int &NAnum);
		void zero(const int &NAnum);
		
		vector< vector<double> > allPis;	//2D array for the fitted pis
		vector<double> allMus; 	//3D array for the fitted mus (note that indexing must be done with MATREF3D)
		vector< vector<double> > allLogDens;	//2D array for the logls, conditional on RCP type
		vector<double> allLogls; //Vector for marginal logls

};

class allClasses
{
	public:
	allClasses();
	~allClasses();

	myData data;
	myParms parms;
	myDerivs derivs;
	myOptContr contr;
	myFits fits;
};

////////////////////////////////////////////////////////
/////////////	Function Definitions	////////////////////
////////////////////////////////////////////////////////
extern "C" SEXP RCP_C( SEXP Ry, SEXP RX, SEXP RW, SEXP Roffset, SEXP Rwts,
						 SEXP RS, SEXP RG, SEXP Rp, SEXP Rpw, SEXP RnObs, SEXP Rdisty,
						 SEXP Ralpha, SEXP Rtau, SEXP Rbeta, SEXP Rgamma, SEXP Rdisps, SEXP Rpowers, 
						 SEXP Rconc, SEXP Rsd, SEXP RsdGamma, SEXP RdispLocat, SEXP RdispScale,
						 SEXP RderivsAlpha, SEXP RderivsTau, SEXP RderivsBeta, SEXP RderivsGamma, SEXP RderivsDisps, SEXP Rscores,
						 SEXP Rpis, SEXP Rmus, SEXP RlogDens, SEXP Rlogli,
						 SEXP Rmaxit, SEXP Rtrace, SEXP RnReport, SEXP Rabstol, SEXP Rreltol, SEXP Rconv,
						 SEXP Roptimise, SEXP RloglOnly, SEXP RderivsOnly, SEXP RoptiDisp, SEXP RgetScores);


double invLogit( double x);
double mixLogl( const myData &dat, const myParms &parms, myFits &fits);
double calcDispPen( const myData &dat, const myParms &parms);
double calcTauPen( const myData &dat, const myParms &parms);
double calcGammaPen( const myData &dat, const myParms &parms);
double calcPiPen( const vector<double> &logPis, const myData &dat, const myParms &parms);
void calcLogPis( vector<double> &logPis, vector<double> &pis, const myData &dat, const myParms &parms, 		int i);
void calcLogCondDens( vector<double> &condDens, const vector<double> &fits, const myData &dat, const myParms &parms, int i);
void calcMuFits( vector<double> &fits, const myData &dat, const myParms &parms);
double logBernoulli( const double &y, double const &mu);
double logPoisson( const double &y, const double &mu);
double logNegBin( const double &y, const double &mu, const double &od);
double logTweedie( const double &y, const double &mu, const double &phi, const double &p);
double logNormal( const double &y, const double &mu, const double &sig);
double calcMixSum( const vector<double> &logPis, const vector<double> &condDens, double &wi, vector<double> &wij, int &maxg);
void loglDerivs( const myData &dat, const myParms &parms, myDerivs &derivs, myFits &fits);
void weightDerivs( vector<double> &alphaDerivI, vector<double> &tauDerivsI, vector<double> &gammaDerivsI, vector<double> &BetaDerivsI, vector<double> &dispDerivsI, const myData &dat, const int &i);
void calcDispDeriv( vector<double> &dispDerivsI, const vector<double> &fits, const myData &dat, const myParms &parms, const double &wi, const vector<double> &wij, const int &m, const int &i);
void calcDispPenDeriv( vector<double> &dispDerivsI, const myData &dat, const myParms &parms);
double logNegBinDispDer( double y, double fit, double dispParm);
double logTweedieDispDer( double y, double fit, double dispParm , double p);
double logNormalDispDer( double y, double fit, double dispParm);
void calcDerivMu( vector<double> &muDerivs, const vector<double> &fits, const myData &dat, const myParms &parms, const double wi, const vector<double> &wij, const int &m, const int &i);
void calcDerivEtaMu( vector<double> &etaDerivsI, const myData &dat, const vector<double> &muDerivsI, const vector<double> &fits, const int &i);
double logBernDer( double y, double mu);
double logPoissonDer( const double &y, const double &mu);
double logNegBinLocatDer( const double &y, const double &mu, const double &od);
double logNormalLocatDer( const double &y, const double &mu, const double &sig);
void calcAlphaDeriv( vector<double> &alphaDerivsI, const vector<double> &etaDerivs, const myData &dat);
void calcTauDeriv( vector<double> &tauDerivsI, const vector<double> &etaDerivs, const myData &dat, const myParms &parms);
void calcGammaDeriv( vector<double> &gammaDerivsI, const vector<double> &etaDerivsI, const myData &dat, const myParms &parms, const int &i);
void calcTauPenDeriv( vector<double> &tauDerivsI, const myData &dat, const myParms &parms);
void calcGammaPenDeriv( vector<double> &gammaDerivsI, const myData &dat, const myParms &parms);
void calcPiDeriv( vector<double> &piDerivs, const myData &dat, const myParms &parms, const vector<double> &pis, const double wi, const vector<double> &wij, int m);
void calcBetaDeriv( vector<double> &betaDerivsI, const vector<double> &piDerivsI, const vector<double> &pis, const myData &dat, int i);
double ALLoptimise( allClasses &all);
bool converged( double *oldP, double *newP, const myOptContr &contr, int nTot);
double optimise_function(int n, double *par, void *ex);
void gradient_function(int n, double *par, double *gr, void *ex);



#endif
