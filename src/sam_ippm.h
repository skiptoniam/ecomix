#ifndef sam_ippm_hh
#define sam_ippm_hh

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

//set up matref
#define MATREF2D(i,j,nx) i+nx*j
#define MATREF3D(i,j,k,nx,ny)  i + nx*j + k*(nx*ny)

// settting up the data class which should hold the data.
class sam_ippm_data {
	public:
		sam_ippm_data();
		~sam_ippm_data();
		void setVals( SEXP &Ry, SEXP &RX, SEXP &Roffset, SEXP &Rwts, SEXP &Ry_not_na, SEXP &RS, SEXP &RG, SEXP &Rp, SEXP &RnObs);
		//bool isDispersion() const;
		//bool doOptiDisp() const;
		void printVals( int printX, int printy);

		int nP, 	//the number of parameters in each of the (G-1) habitat lps, same as lpar
			nG,			//the number of habitats
			nS, 		//the number of species
			nObs,		//the number of observations
			NAnum;	//a common number to insert for NAs

		double 	*X, //the design matrix in vector form (nObs x nP)
				*y,	//the outcome matrix, in vector form (nObs x nS)
				*offset, //the offset vector (length nObs)
				*wts;  //the wts for the logl (typically all zero and of length nObs).
		int 	*y_not_na; //a matrix which keeps track of NAs in ippm data. If non-ippm model all == 1.

};

// Doh' I forgot to add this!
class sam_ippm_params {
	public:
		sam_ippm_params();
		~sam_ippm_params();
		void setVals( const sam_ippm_data &dat, SEXP &Ralpha, SEXP &Rbeta, SEXP &Rpi);
		void getArray(double *parArr, const sam_ippm_data &dat);
		void update( double *parArr, const sam_ippm_data &dat);
		void printParms( const sam_ippm_data &dat);


		double 	*Alpha, //the species' prevalences (Gx1)==(1xG)
				*Beta,	//the habitats' free covariate params ((G-1)xp)
				*Pi;	//the pis - mmmmm pies.
		int nalpha, nbeta, npi, nTot;
};



// setting up the derivatives which should setup, hold and update derivatives.
class sam_ippm_derivs{
	public:
		sam_ippm_derivs();
		~sam_ippm_derivs();
		void setVals( const sam_ippm_data &dat, SEXP &RderivsAlpha, SEXP &RderivsBeta, SEXP &RderivsPi, SEXP &RgetScores, SEXP &Rscores);
		void zeroDerivs( const sam_ippm_data &dat);
		void updateDerivs( const sam_ippm_data &dat, const vector<double> &alphaDerivs, const vector<double> &betaDerivs, const vector<double> &piDerivs);
		void update( double *grArr, const sam_ippm_data &dat);
		void getArray( double *grArr, const sam_ippm_data &dat);

        int getScoreFlag;	//Should the scores be calculated for empirical information
		double 	*Alpha, //the derivatives of logl w.r.t. alpha
				*Beta,	//the derivatives of logl w.r.t. beta
				*Pi, 	//the derivatives of logl w.r.t. pi
				*Scores;//the score contribution for each site (for empirical information)
	
};

// control functions for species mix.
class sam_ippm_opt_contr {
	public:
		sam_ippm_opt_contr();
		~sam_ippm_opt_contr();
		void setVals( const SEXP &Rmaxit, const SEXP &Rtrace, const SEXP &RnReport, const SEXP &Rabstol, const SEXP &Rreltol, SEXP &Rconv);

		int maxitQN, traceQN, nReport, fnKount, grKount, ifail, *conv;
		double abstol, reltol, denomEps;
};

class sam_ippm_fits {
	public:
		sam_ippm_fits();
		~sam_ippm_fits();
		void initialise( const int &nObs, const int &nG, const int &nS, const int &nP, const int &NAnum);
		void zero(const int &NAnum);

		vector<double> pis_loglike;	
		//vector<double> pis_grads;
		//vector<double> grad_pis_trans;	
		//vector<double> grad_pis_non_trans;
		vector<double> allMus; 	//3D array for the fitted mus (note that indexing must be done with MATREF3D)
		vector<double> log_like_species_group_contrib; // 2D array of loglikes species and groups.
		vector<double> log_like_species_contrib; //vector of species logls.
		vector<double> dlogdalpha;
		vector<double> dlogdbeta;
		vector<double> dlogdpi;


};

// a wrapper around all data classes.
class sam_ippm_all_classes
{
	public:
	sam_ippm_all_classes();
	~sam_ippm_all_classes();

	sam_ippm_data data;
	sam_ippm_params params;
	sam_ippm_derivs derivs;
	sam_ippm_opt_contr contr;
	sam_ippm_fits fits;
};

////////////////////////////////////////////////////////
/////////////	Function Definitions	////////////////
////////////////////////////////////////////////////////

extern "C" SEXP species_mix_ippm_cpp(SEXP Ry, SEXP RX, SEXP Roffset, SEXP Rwts, SEXP Ry_not_na,
									 SEXP RnS, SEXP RnG, SEXP Rp, SEXP RnObs,
									 SEXP Ralpha, SEXP Rbeta, SEXP Rpi,
									 SEXP RderivsAlpha, SEXP RderivsBeta, SEXP RderivsPi, SEXP RgetScores, SEXP Rscores,
									 SEXP Rpis, SEXP Rmus, SEXP logliS, SEXP logliSG,
									 SEXP Rmaxit, SEXP Rtrace, SEXP RnReport, SEXP Rabstol, SEXP Rreltol, SEXP Rconv,
									 SEXP Roptimise, SEXP RloglOnly, SEXP RderivsOnly);


// functions for optimisation.
double sam_ippm_optimise(sam_ippm_all_classes &all);
bool converged_ippm( double *oldP, double *newP, const sam_ippm_opt_contr &contr, int nTot);

// functions for calculating ippm likelihood.
double optimise_function_ippm(int n, double *par, void *ex);
double sam_ippm_mix_loglike(sam_ippm_data &dat, sam_ippm_params &params, sam_ippm_fits &fits);
void calc_mu_fits( vector<double> &fits, const sam_ippm_data &dat, const sam_ippm_params &params);
double calc_ippm_loglike_per_species(sam_ippm_data &dat, sam_ippm_params &params, sam_ippm_fits &fits, int i);
double log_ippm(const double &y, const double &mu, const double &wts);
void additive_logistic_ippm(vector< double > &x,int inv, int G);

// functions for calculating the gradient.
void gradient_function_ippm(int n, double *par, double *gr, void *ex);
void sam_ippm_mix_gradient_function( const sam_ippm_data &dat, const sam_ippm_params &params, sam_ippm_derivs &derivs, sam_ippm_fits &fits);
double log_poisson_deriv( const double &y, const double &mu, const double &wts);
void calc_dlog_dalpha(const sam_ippm_data &dat, sam_ippm_fits &fits);
void calc_dlog_dbeta(const sam_ippm_data &dat, sam_ippm_fits &fits);
void calc_alpha_deriv( vector<double> &alphaDerivs, const sam_ippm_data &dat, sam_ippm_fits &fits);
void calc_beta_deriv( vector<double> &betaDerivs, const sam_ippm_data &dat, sam_ippm_fits &fits);
void calc_pi_deriv( vector<double> &piDerivs, const sam_ippm_data &dat, sam_ippm_derivs &derivs, sam_ippm_fits &fits);


#endif
