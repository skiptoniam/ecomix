#ifndef sam_bernoulli_sp_ints_hh
#define sam_bernoulli_sp_ints_hh

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
class sam_bernoulli_sp_ints_data {
	public:
		sam_bernoulli_sp_ints_data();
		~sam_bernoulli_sp_ints_data();
		void setVals( SEXP &Ry, SEXP &RX, SEXP &Roffset, SEXP &Rwts, SEXP &RS, SEXP &RG, SEXP &Rp, SEXP &RnObs);
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
				*wts; //a vector of weights (length nObs) - useful for BB.
};

// Doh' I forgot to add this!
class sam_bernoulli_sp_ints_params {
	public:
		sam_bernoulli_sp_ints_params();
		~sam_bernoulli_sp_ints_params();
		void setVals( const sam_bernoulli_sp_ints_data &dat, SEXP &Ralpha, SEXP &Rbeta, SEXP &Rpi);
		void getArray(double *parArr, const sam_bernoulli_sp_ints_data &dat);
		void update(double *parArr, const sam_bernoulli_sp_ints_data &dat);
		void printParms( const sam_bernoulli_sp_ints_data &dat);


		double 	*Alpha, //the species' prevalences 
				*Beta,	//the habitats' free covariate params (G*xp)
				*Eta;	//the pis - mmmmm pies.
		int nalpha, nbeta, npi, nTot;
};



// setting up the derivatives which should setup, hold and update derivatives.
class sam_bernoulli_sp_ints_derivs{
	public:
		sam_bernoulli_sp_ints_derivs();
		~sam_bernoulli_sp_ints_derivs();
		void setVals( const sam_bernoulli_sp_ints_data &dat, SEXP &RderivsAlpha, SEXP &RderivsBeta, SEXP &RderivsEta, SEXP &RgetScores, SEXP &Rscores);
		void zeroDerivs( const sam_bernoulli_sp_ints_data &dat);
		void updateDerivs( const sam_bernoulli_sp_ints_data &dat, const vector<double> &alphaDerivs, const vector<double> &betaDerivs, const vector<double> &etaDerivs);
		void update( double *grArr, const sam_bernoulli_sp_ints_data &dat);
		void getArray( double *grArr, const sam_bernoulli_sp_ints_data &dat);

        int getScoreFlag;	//Should the scores be calculated for empirical information
		double 	*Alpha, //the derivatives of logl w.r.t. alpha
				*Beta,	//the derivatives of logl w.r.t. beta
				*Eta, 	//the derivatives of logl w.r.t. eta (transformed pi)
				*Scores;//the score contribution for each site (for empirical information)
	
};

// control functions for species mix.
class sam_bernoulli_sp_ints_opt_contr {
	public:
		sam_bernoulli_sp_ints_opt_contr();
		~sam_bernoulli_sp_ints_opt_contr();
		void setVals( const SEXP &Rmaxit, const SEXP &Rtrace, const SEXP &RnReport, const SEXP &Rabstol, const SEXP &Rreltol, SEXP &Rconv, const SEXP &Rprintparam);

		int maxitQN, traceQN, nReport, fnKount, grKount, ifail, *conv, printparams;
		double abstol, reltol, denomEps;
};

class sam_bernoulli_sp_ints_fits {
	public:
		sam_bernoulli_sp_ints_fits();
		~sam_bernoulli_sp_ints_fits();
		void initialise( const int &nObs, const int &nG, const int &nS, const int &nP, const int &NAnum);
		void zero(const int &NAnum);

		//vector<double> par_pis;
		//vector<double> par_etas;	
		vector<double> allMus; 	//3D array for the fitted mus (note that indexing must be done with MATREF3D)
		vector<double> log_like_species_group_contrib; // 2D array of loglikes species and groups.
		vector<double> log_like_species_contrib; //vector of species logls.
		vector<double> dlogdalpha;
		vector<double> dlogdbeta;
		vector<double> dlogdpi;


};

// a wrapper around all data classes.
class sam_bernoulli_sp_ints_all_classes
{
	public:
	sam_bernoulli_sp_ints_all_classes();
	~sam_bernoulli_sp_ints_all_classes();

	sam_bernoulli_sp_ints_data data;
	sam_bernoulli_sp_ints_params params;
	sam_bernoulli_sp_ints_derivs derivs;
	sam_bernoulli_sp_ints_opt_contr contr;
	sam_bernoulli_sp_ints_fits fits;
};

////////////////////////////////////////////////////////
/////////////	Function Definitions	////////////////
////////////////////////////////////////////////////////

extern "C" SEXP species_mix_bernoulli_cpp(SEXP Ry, SEXP RX, SEXP Roffset, SEXP Rwts,
									 SEXP RnS, SEXP RnG, SEXP Rp, SEXP RnObs,
									 SEXP Ralpha, SEXP Rbeta, SEXP Reta,
									 SEXP RderivsAlpha, SEXP RderivsBeta, SEXP RderivsEta, SEXP RgetScores, SEXP Rscores,
									 SEXP Rpis, SEXP Rmus, SEXP logliS, SEXP logliSG,
									 SEXP Rmaxit, SEXP Rtrace, SEXP RnReport, SEXP Rabstol, SEXP Rreltol, SEXP Rconv, SEXP Rprintparams,
									 SEXP Roptimise, SEXP RloglOnly, SEXP RderivsOnly);


// functions for optimisation.
double sam_bernoulli_sp_ints_optimise(sam_bernoulli_sp_ints_all_classes &all);
bool converged_bernoulli( double *oldP, double *newP, const sam_bernoulli_sp_ints_opt_contr &contr, int nTot);

// functions for calculating bernoulli likelihood.
double optimise_function_bernoulli_sp_ints(int n, double *par, void *ex);
double sam_bernoulli_sp_ints_mix_loglike(const sam_bernoulli_sp_ints_data &dat, const sam_bernoulli_sp_ints_params &params, sam_bernoulli_sp_ints_fits &fits);
void calc_mu_fits(vector<double> &fits, const sam_bernoulli_sp_ints_params &params, const sam_bernoulli_sp_ints_data &dat);
void calc_bernoulli_loglike_SG(vector<double> &loglSG, vector<double> &fits, const sam_bernoulli_sp_ints_data &dat);
double calc_bernoulli_loglike_S(vector<double> &fits, vector<double> const &pis, const sam_bernoulli_sp_ints_data &dat, int s);
double log_bernoulli(const double &y, const double &mu);
double dmu_deta_bernoulli(const double &mu);
double invLogit_bern(const double eta);
void additive_logistic_bernoulli_sp_ints(vector< double > &x,int inv, int G);

// functions for calculating the gradient.
void gradient_function_bernoulli_sp_ints(int n, double *par, double *gr, void *ex);
void sam_bernoulli_sp_ints_mix_gradient_function(const sam_bernoulli_sp_ints_data &dat, const sam_bernoulli_sp_ints_params &params, sam_bernoulli_sp_ints_derivs &derivs, sam_bernoulli_sp_ints_fits &fits);
double log_bernoulli_deriv(const double &y, const double &mu);
void calc_dlog_dalpha(vector<double> &dlda, vector<double> const &mus, const sam_bernoulli_sp_ints_data &dat);
void calc_dlog_dbeta(vector<double> &dldb, vector<double> const &mus, const sam_bernoulli_sp_ints_data &dat);
void calc_dlog_dpi(vector<double> &dldpi, vector<double> const &llS, vector<double> const &llSG, const sam_bernoulli_sp_ints_data &dat);
void calc_alpha_deriv( vector<double> &alphaDerivs, vector<double> const &dlogdalpha, vector<double> const &llSG, vector<double> const &llS, vector<double> const &pis, const sam_bernoulli_sp_ints_data &dat);
void calc_beta_deriv( vector<double> &betaDerivs, vector<double> const &dlogdbeta, vector<double> const &llSG, vector<double> const &llS, vector<double> const &pis, const sam_bernoulli_sp_ints_data &dat);
void calc_eta_deriv( vector<double> &etaDerivs, vector<double> const &dlogdpi, vector<double> const eta, const sam_bernoulli_sp_ints_data &dat);


#endif
