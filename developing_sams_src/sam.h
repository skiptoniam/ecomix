#include<R.h>
#include<Rinternals.h>
#include<Rmath.h>
#include<R_ext/Applic.h>
#undef length

#include <vector>
#include<algorithm>
#include <iostream>
#include"Tweedie.h"

/* setup MATREF
 * 
 * 
 */

#define MATREF2D(i,j,nx) i+nx*j
#define MATREF3D(i,j,k,nx,ny)  i + nx*j + k*(nx*ny)
using std::vector;         // use vector as abbreviation for std::vector

//class Optimise_data{
//public:
  //Optimise_data();
  //~Optimise_data();
  
  //double *y, *X, *tau, *offset;
  //int *ID;
  //int S, G,  Xr, Xc, lpar, ly;
  //vector<int> StartEndPos;
  
  //void SetVars(double *ty, double *tX, int *tID, int tS, int tG, int tXr, int tXc, int tlpar, int tlobs, double *ttau, double *toffset);
  
  //// holders for derivitive calculation values

  //double logl;
  //vector<double> sum_f_species; // log( f(yi,Bi)) S*G long
  //vector<double> deriv_f_B;  //  d(log f(yi,Bi)) / d( Bi) G*Xc*S long
  //vector<double> deriv_f_alphaS; //d(log(yi,alphai,Bi)/d(alphai) S long
  //vector<double> deriv_f_dispersionS; // d(log(yi,Bi,thetai)/d(thetai) S long
  //vector<double> parpi; // vector containing calculated pi's G long
  //vector<double> species_l_contrib; // vector of each species likelihood contribution // S long
  //vector<double> species_group_l_contrib; // matrix holding likelihood contribution for each combination of species & group
  //vector<double> log_y_factorial; //calculation of factorial for nbinom
//};

//class Optimise_data_nbinom{
//public:
  //Optimise_data_nbinom();
  //~Optimise_data_nbinom();
  
  //double *y, *X,*w, offset;
  //int  Xr, Xc;
  //vector<double> lp;
  //vector<double> log_y_factorial;
  
  //void SetVars(double *ty, double *tX, double *tw, double toffset, int tXr, int tXc);

 
//};

//class Optimise_data_ippm{
//public:
  //Optimise_data_ippm();
  //~Optimise_data_ippm();
  
  ////double *y, *X, *w, offset;
  //int *y_is_not_na, *ID, Xr, Xc, S, G, lpar, ly;
  //vector<double> lp;
  //double *y, *X, *w, *tau, *offset;
  
  //void set_vars_ippm(double *ty, double *tX, double *tw, double toffset, int *y_is_not_na, int *tID, int tS, int tG, int tXr, int tXc, int tlpar, int tlobs, double *ttau);
    
  //double logl;
  //vector<double> sum_f_species; // log( f(yi,Bi)) S*G long
  //vector<double> deriv_f_B;  //  d(log f(yi,Bi)) / d( Bi) G*Xc*S long
  //vector<double> deriv_f_alphaS; //d(log(yi,alphai,Bi)/d(alphai) S long
  //vector<double> parpi; // vector containing calculated pi's G long
  //vector<double> species_l_contrib; // vector of each species likelihood contribution // S long
  //vector<double> species_group_l_contrib; // matrix holding likelihood contribution for each combination of species & group
 
//};

/* Setting up new data structure for species mix. 
 * I'm trying to make it someway consistent with RCP so it'll be easier to maintain throught time. Fingers crossed :P 
 */
 
class spmix_data {
	public:
		spmix_data();
		~spmix_data();
		void setVals( const SEXP &Ry, const SEXP &RX, const SEXP &Roffset, const SEXP &Rwts, const SEXP &Ry_not_na,
					  const SEXP &RS, const SEXP &RG, const SEXP &Rnp, const SEXP &Rnobs, const SEXP &Rdisty, const SEXP &RoptiDisp);
		bool isDispersion() const;
		bool doOptiDisp() const;
		void printVals( int printX, int printy);

		int np, 		//the number of parameters in each of the (G-1) habitat lps, same as lpar
			nG,			//the number of habitats
			nS, 		//the number of species
			nobs,		//the number of observations ly
			disty,		//the distrbution code
			optiDisp,	//should the dispersions be optimised
			//NAnum;	//a common number to insert for NAs

		double 	*X, //the design matrix in vector form (nObs x np)
				*y,	//the outcome matrix, in vector form (nObs x nS)
				*offset, //the offset vector (length nObs)
				*wts;  //the wts for the logl (typically all zero and of length nObs).
		int 	*y_not_na; //a matrix which keeps track of NAs in ippm data. If non-ippm model all == 1.		
		
		double logl;
		vector<double> lp;//linear predictor.
		vector<double> sum_f_species; // log( f(yi,Bi)) S*G long
		vector<double> deriv_f_B;  //  d(log f(yi,Bi)) / d( Bi) G*Xc*S long
        vector<double> deriv_f_alphaS; //d(log(yi,alphai,Bi)/d(alphai) S long
        vector<double> parpi; // vector containing calculated pi's G long
        vector<double> species_l_contrib; // vector of each species likelihood contribution // S long
        vector<double> species_group_l_contrib; // matrix holding likelihood contribution for each combination of species & group
        vector<double> log_y_factorial; //calculation of factorial for nbinom
};

class spmix_derivs
{
	public:
		spmix_derivs();
		~spmix_derivs();
		void setVals( const myData &dat, SEXP &RderivsAlpha, SEXP &RderivsBeta, SEXP &RderivsTau, SEXP &RderivsDisp, SEXP &RgetScores, SEXP &Rscores);
		void zeroDerivs( const myData &dat);
		void updateDerivs( const myData &dat, const vector<double> &alphaDerivsI, const vector<double> &betaDerivsI, const vector<double> &tauDerivsI, const vector<double> &dispDerivsI, const int &i);
		void update( double *grArr, const myData &dat);
		void getArray( double *grArr, const myData &dat);

		int getScoreFlag;	//Should the scores be calculated for empirical information
		double 	*Alpha, //the derivatives of logl w.r.t. alpha
				*Tau, 	//the derivatives of logl w.r.t. tau
				*Beta,	//the derivatives of logl w.r.t. beta
				*Disp,	//the derivative of logl w.r.t. dispersions
				*Scores;//the score contribution for each site (for empirical information)
};

// control functions for species mix.
class spmix_opt_contr
{
	public:
		spmix_opt_contr();
		~spmix_opt_contr();
		void setVals( const SEXP &Rmaxit, const SEXP &Rtrace, const SEXP &RnReport, const SEXP &Rabstol, const SEXP &Rreltol, SEXP &Rconv);

		int maxitQN, traceQN, nReport, fnKount, grKount, ifail, *conv;
		double abstol, reltol, denomEps;
};

class spmix_fits
{
	public:
		spmix_fits();
		~spmix_fits();
		void initialise( const int &nObs, const int &nG, const int &nS, const int &NAnum);
		void zero(const int &NAnum);

		vector< vector<double> > allPis;	//2D array for the fitted pis
		vector<double> allMus; 	//3D array for the fitted mus (note that indexing must be done with MATREF3D)
		vector< vector<double> > allLogDens;	//2D array for the logls, conditional on RCP type
		vector<double> allLogls; //Vector for marginal logls

};

// a wrapper around all data classes. 
class spmix_all_classes
{
	public:
	spmix_all_classes();
	~spmix_all_classes();

	spmix_data data;
	//myParms parms;
	//spmix_derivs derivs;
	spmix_opt_contr contr;
	//spmix_fits fits;
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
