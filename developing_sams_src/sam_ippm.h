#include<R.h>
#include<Rinternals.h>
#include<Rmath.h>
#include<R_ext/Applic.h>
#undef length

#include <vector>
#include<algorithm>
#include <iostream>
#include"Tweedie.h"

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
			nobs,		//the number of observations ly
			//disty,		//the distrbution code
			//optiDisp,	//should the dispersions be optimised
			NAnum;	//a common number to insert for NAs

		double 	*X, //the design matrix in vector form (nObs x np)
				*y,	//the outcome matrix, in vector form (nObs x nS)
				*offset, //the offset vector (length nObs)
				*wts;  //the wts for the logl (typically all zero and of length nObs).
		int 	*y_not_na; //a matrix which keeps track of NAs in ippm data. If non-ippm model all == 1.		
		
		//// this will all be setup in params and derivs.
		
		//double logl;
		//vector<double> lp;//linear predictor.
		//vector<double> sum_f_species; // log( f(yi,Bi)) S*G long
		//vector<double> deriv_f_B;  //  d(log f(yi,Bi)) / d( Bi) G*Xc*S long
        //vector<double> deriv_f_alphaS; //d(log(yi,alphai,Bi)/d(alphai) S long
        //vector<double> parpi; // vector containing calculated pi's G long
        //vector<double> species_l_contrib; // vector of each species likelihood contribution // S long
        //vector<double> species_group_l_contrib; // matrix holding likelihood contribution for each combination of species & group
        //vector<double> log_y_factorial; //calculation of factorial for nbinom
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
				//*Disp,	//the derivative of logl w.r.t. dispersions
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

extern "C"  SEXP species_mix_bernoulli_cpp(SEXP R_pars, SEXP R_y, SEXP R_X, SEXP R_ID,SEXP R_tau, SEXP R_gradient, SEXP R_offset);//, SEXP R_model_type);
extern "C"  SEXP species_mix_bernoulli_gradient_cpp(SEXP R_pars, SEXP R_y, SEXP R_X, SEXP R_ID,SEXP R_tau, SEXP R_gradient, SEXP R_offset);//, SEXP R_model_type);

double optimise_function_sam(int n, double *pars, void *ex);
void gradient_function_sam(int n, double *pars, double *gr, void *ex );

vector<double> calc_logl(const vector<double> &pars, Optimise_data &data);
double like_function(vector <double> &estpi, vector < double > &coef, const double *y, const double *X, int Xr, int Xc, int start, int end, double *tau, int s, vector<double> &sum_f_species, vector<double> &deriv_f_B);
double link_function(double p, int link_type);
double inv_link_function(double lpre, int link_type);

void additive_logistic(vector< double> &x,int inv);

// negative binom 
extern "C"  SEXP species_mix_negative_binomial_cpp(SEXP R_pars, SEXP R_X, SEXP R_y, SEXP R_w, SEXP R_offset, SEXP R_gradient, SEXP R_fitted_values);
extern "C"  SEXP species_mix_negative_binomial_gradient_cpp(SEXP R_pars, SEXP R_X, SEXP R_y, SEXP R_w, SEXP R_offset, SEXP R_gradient);
double optimise_nbinom(int n, double *pars, void *ex);
void gradient_nbinom(int n, double *pars, double *gr, void *ex );
double negative_binomial_logl(vector<double> &pars, Optimise_data_nbinom &data );

// nbinom mix
double optimise_mixnbinom_function(int n, double *pars, void *ex);
void gradient_mixnbinom_function(int n, double *pars, double *gr, void *ex );
vector <double> calc_mixnbinom_logl(const vector<double> &pars, Optimise_data &data);
double like_mixnbinom_function(vector< double > &estpi, vector < double > &coef, vector < double > &sp_int, vector < double > &sp_dispersion, const double *y, const double *X, int Xr, int Xc, int start, int end, double *tau, int s, vector<double> &sum_f_species, vector<double> &deriv_f_B, vector<double> &deriv_f_alphaS, vector<double> &deriv_f_dispersionS, vector<double> &log_y_factorial, const double *offset);

//ippm functions
extern "C"  SEXP species_mix_ippm_cpp(SEXP R_pars, SEXP R_X, SEXP R_y, SEXP R_w, SEXP R_offset, SEXP R_y_is_not_na, SEXP R_gradient, SEXP R_fitted_values);
extern "C"  SEXP species_mix_ippm_gradient_cpp(SEXP R_pars, SEXP R_X, SEXP R_y, SEXP R_w, SEXP R_offset, SEXP R_y_is_not_na, SEXP R_gradient);
double optimise_ippm(int n, double *pars, void *ex);
void gradient_ippm(int n, double *pars, double *gr, void *ex );
double ippm_logl(vector<double> &pars, Optimise_data_ippm &data );

// ippm mix functions
double optimise_mix_ippm_function(int n, double *pars, void *ex);
void gradient_mix_ippm_function(int n, double *pars, double *gr, void *ex );
vector <double> calc_mix_ippm_logl(const vector<double> &pars, Optimise_data_ippm &data);
double like_mix_ippm_function(vector< double > &estpi, vector < double > &coef, vector < double > &sp_int, const double *y, const double *X, int Xr, int Xc, int start, int end, double *tau, int s, vector<double> &sum_f_species, vector<double> &deriv_f_B, vector<double> &deriv_f_alphaS, vector<double> &log_y_factorial, const double *offset);
double log_ippm_derivative(double y, double mu, double wts);


// tweedie mix
double optimise_tweedie_function(int n, double *pars, void *ex);
void gradient_tweedie_function(int n, double *pars, double *gr, void *ex );
vector <double> calc_tweedie_logl(const vector<double> &pars, Optimise_data &data);
double like_tweedie_function(vector< double > &estpi, vector < double > &coef, vector < double > &sp_int, vector < double > &sp_dispersion, const double *y, const double *X, int Xr, int Xc, int start, int end, double *tau, int s, vector<double> &sum_f_species, vector<double> &deriv_f_B, vector<double> &deriv_f_alphaS, vector<double> &deriv_f_dispersionS, vector<double> &log_y_factorial, const double *offset);
