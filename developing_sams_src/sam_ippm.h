#ifndef sam_ippm_hh
#define sam_ippm_hh

#include<R.h>
#include<Rinternals.h>
#include<Rmath.h>
#include<R_ext/Applic.h>
#undef length

#include <vector>
#include<algorithm>
#include <iostream>


//set up matref
#define MATREF2D(i,j,nx) i+nx*j
#define MATREF3D(i,j,k,nx,ny)  i + nx*j + k*(nx*ny)
using std::vector;         // use vector as abbreviation for std::vector

/* Setting up new data structure for species mix. 
 * I'm trying to make it someway consistent with RCP so it'll be easier to maintain throught time. Fingers crossed :P 
 */
 
class sam_ippm_data {
	public:
		sam_ippm_data();
		~sam_ippm_data();
		void setVals( const SEXP &Ry, const SEXP &RX, const SEXP &Roffset, const SEXP &Rwts, const SEXP &Ry_not_na,
					  const SEXP &RS, const SEXP &RG, const SEXP &Rnp, const SEXP &Rnobs);
		//bool isDispersion() const;
		//bool doOptiDisp() const;
		void printVals( int printX, int printy);

		int np, 		//the number of parameters in each of the (G-1) habitat lps, same as lpar
			nG,			//the number of habitats
			nS, 		//the number of species
			nobs,		//the number of observations
			NAnum;	//a common number to insert for NAs

		double 	*X, //the design matrix in vector form (nObs x np)
				*y,	//the outcome matrix, in vector form (nObs x nS)
				*offset, //the offset vector (length nObs)
				*wts;  //the wts for the logl (typically all zero and of length nObs).
		int 	*y_not_na; //a matrix which keeps track of NAs in ippm data. If non-ippm model all == 1.		
		
};

class sam_ippm_derivs
{
	public:
		sam_ippm_derivs();
		~sam_ippm_derivs();
		void setVals( const sam_ippm_data &dat, SEXP &RderivsAlpha, SEXP &RderivsBeta, SEXP &RderivsTau, SEXP &RgetScores, SEXP &Rscores);
		void zeroDerivs( const sam_ippm_data &dat);
		void updateDerivs( const sam_ippm_data &dat, const vector<double> &alphaDerivsI, const vector<double> &betaDerivsI, const vector<double> &tauDerivsI, const int &i);
		void update( double *grArr, const sam_ippm_data &dat);
		void getArray( double *grArr, const sam_ippm_data &dat);

		int getScoreFlag;	//Should the scores be calculated for empirical information
		double 	*Alpha, //the derivatives of logl w.r.t. alpha
				*Tau, 	//the derivatives of logl w.r.t. tau
				*Beta,	//the derivatives of logl w.r.t. beta
				*Scores;//the score contribution for each site (for empirical information)
};

// control functions for species mix.
class sam_ippm_opt_contr
{
	public:
		sam_ippm_opt_contr();
		~sam_ippm_opt_contr();
		void setVals( const SEXP &Rmaxit, const SEXP &Rtrace, const SEXP &RnReport, const SEXP &Rabstol, const SEXP &Rreltol, SEXP &Rconv);

		int maxitQN, traceQN, nReport, fnKount, grKount, ifail, *conv;
		double abstol, reltol, denomEps;
};

class sam_ippm_fits
{
	public:
		sam_ippm_fits();
		~sam_ippm_fits();
		void initialise( const int &nObs, const int &nG, const int &nS, const int &NAnum);
		void zero(const int &NAnum);

		vector< vector<double> > allPis;	//2D array for the fitted pis
		vector<double> allMus; 	//3D array for the fitted mus (note that indexing must be done with MATREF3D)
		vector< vector<double> > allLogDens;	//2D array for the logls, conditional on SAM type
		vector<double> allLogls; //Vector for marginal logls

};

// a wrapper around all data classes. 
class sam_ippm_all_classes
{
	public:
	sam_ippm_all_classes();
	~sam_ippm_all_classes();

	sam_ippm_data data;
	sam_ippm_params parms;
	sam_ippm_derivs derivs;
	sam_ippm_opt_contr contr;
	sam_ippm_fits fits;
};

////////////////////////////////////////////////////////
/////////////	Function Definitions	////////////////
////////////////////////////////////////////////////////

extern "C" SEXP species_mix_ippm_cpp(SEXP Ry, SEXP RX, SEXP Roffset, SEXP Rwts, SEXP Ry_not_na,
									 SEXP RnS, SEXP RnG, SEXP RnP, SEXP RnObs,
									 SEXP Ralpha, SEXP Rbeta, SEXP Rtau, 
									 SEXP RderivsAlpha, SEXP RderivsBeta, SEXP RderivsTau, SEXP Rscores,
									 SEXP Rpis, SEXP Rmus, SEXP RlogDens, SEXP Rloglike,
									 SEXP Rmaxit, SEXP Rtrace, SEXP RnReport, SEXP Rabstol, SEXP Rreltol, SEXP Rconv,
									 SEXP Roptimise, SEXP RloglOnly, SEXP RderivsOnly, SEXP RoptiDisp, SEXP RgetScores);
									 
// functions for optimisation.									 
double sam_optimise(sam_ippm_all_classes &all);

// functions for calculating ippm likelihood.
double optimise_function_sam_ippm(int n, double *par, void *ex)
double sam_ippm_mix_loglike( const sam_ippm_data &dat, const sam_ippm_params &parms, sam_ippm_fits &fits)
void calc_mu_fits( vector<double> &fits, const sam_ipp_data &dat, const sam_ippm_params &parms)
void calc_log_cond_dens( vector<double> &condDens, const vector<double> &fits, const sam_ippm_data &dat, const sam_ippm_params &parms, int i)
double log_ippm(const double &y, const double &mu, const double &wts)

//// functions for calculating ippm derivatives.
//void gradient_function_sam_ippm(int n, double *par, double *gr, void *ex)
//double log_ippm_derivative( const double &y, const double &mu, const double &wts)

#endif
