#ifndef sam_cpp_hh
#define sam_cpp_hh

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
class sam_data {
	public:
		sam_data();
		~sam_data();
		void setVals( SEXP &Ry, SEXP &RX, SEXP &Roffset, SEXP &Rspp_wts, SEXP &Rsite_spp_wts, SEXP &Ry_not_na,
		 SEXP &RS, SEXP &RG, SEXP &Rp, SEXP &RnObs, SEXP &Rdisty, SEXP &RoptiDisp);
		bool isDispersion() const;
		bool doOptiDisp() const;
		void printVals( int printX, int printy);

		int nP,       //the number of parameters in each of the (G-1) habitat lps, same as lpar
			//nPW,    //the number of parameters in the species level moodel *** yet to be implemented ***
			nG,	      //the number of habitats
			nS,       //the number of species
			nObs,     //the number of observations
			disty,    //the distribution code
			optiDisp, //should the dispersion parameter be optimised (this is for negative binomial).
			NAnum;    //a common number to insert for NAs

		double 	*X, //the design matrix in vector form (nObs x nP)
				*y,	//the outcome matrix, in vector form (nObs x nS)
				//*W, //the design matix in vector for for the species model (nObs x npw) ** yet to be implemented. 
				*offset, //the offset vector (length nObs)
				*spp_wts,  //the wts for the logl dependent on the species - used for the bayesian bootstrap (typically all zero and of length nObs).
				*site_spp_wts; // the wts for the ippm.
		int 	*y_not_na; //a matrix which keeps track of NAs in ippm data. If non-ippm model all == 1.

};

// Doh' I forgot to add this!
class sam_params {
	public:
		sam_params();
		~sam_params();
		void setVals( const sam_data &dat, SEXP &Ralpha, SEXP &Rbeta, SEXP &Reta, SEXP &Rdisp);
		void getArray(double *parArr, const sam_data &dat);
		void update(double *parArr, const sam_data &dat);
		void printParms( const sam_data &dat);


		double 	*Alpha, //the species' prevalences 
				*Beta,	//the habitats' free covariate params (G*xp)
				*Eta,	//the pis - mmmmm pies.
				*Disp;  //the dispersion parameter for negative binomial model nspp long.
		int nalpha, nbeta, neta, ndisp, nTot;
};



// setting up the derivatives which should setup, hold and update derivatives.
class sam_derivs{
	public:
		sam_derivs();
		~sam_derivs();
		void setVals( const sam_data &dat, SEXP &RderivsAlpha, SEXP &RderivsBeta, SEXP &RderivsEta, SEXP &RderivsDisp, SEXP &RgetScores, SEXP &Rscores);
		void zeroDerivs( const sam_data &dat);
		void updateDerivs( const sam_data &dat, const vector<double> &alphaDerivs, const vector<double> &betaDerivs, const vector<double> &etaDerivs, const vector<double> &dispDerivs);
		void update( double *grArr, const sam_data &dat);
		void getArray( double *grArr, const sam_data &dat);

        int getScoreFlag;	//Should the scores be calculated for empirical information
		double 	*Alpha, //the derivatives of logl w.r.t. alpha
				*Beta,	//the derivatives of logl w.r.t. beta
				*Eta, 	//the derivatives of logl w.r.t. eta (transformed pi)
				*Disp,  //the derivatives of logl w.r.t. dispersion parameters.
				*Scores;//the score contribution for each site (for empirical information)
	
};

// control functions for species mix.
class sam_opt_contr {
	public:
		sam_opt_contr();
		~sam_opt_contr();
		void setVals( const SEXP &Rmaxit, const SEXP &Rtrace, const SEXP &RnReport, const SEXP &Rabstol, const SEXP &Rreltol, SEXP &Rconv, const SEXP &Rprintparam);

		int maxitQN, traceQN, nReport, fnKount, grKount, ifail, *conv, printparams;
		double abstol, reltol, denomEps;
};

class sam_fits {
	public:
		sam_fits();
		~sam_fits();
		void initialise( const int &nObs, const int &nG, const int &nS, const int &nP, const int &NAnum);
		void zero(const int &NAnum);

		//vector<double> par_pis;
		//vector<double> par_etas;	
		vector<double> allMus; 	//3D array for the fitted mus (note that indexing must be done with MATREF3D)
		vector<double> log_like_species_group_contrib; // 2D array of loglikes species and groups.
		vector<double> log_like_species_contrib; //vector of species logls.
		vector<double> all_derivs_mu;
		vector<double> dlogdalpha;
		vector<double> dlogdbeta;
		vector<double> dlogdpi;
		vector<double> dlogddispersion;


};

// a wrapper around all data classes.
class sam_cpp_all_classes
{
	public:
	sam_cpp_all_classes();
	~sam_cpp_all_classes();

	sam_data data;
	sam_params params;
	sam_derivs derivs;
	sam_opt_contr contr;
	sam_fits fits;
};

////////////////////////////////////////////////////////
/////////////	Function Definitions	////////////////
////////////////////////////////////////////////////////

extern "C" SEXP species_mix_cpp(SEXP Ry, SEXP RX, SEXP Roffset, SEXP Rspp_wts, SEXP Rsite_spp_wts, SEXP Ry_not_na,
									 SEXP RnS, SEXP RnG, SEXP Rp, SEXP RnObs, SEXP Rdisty,
									 SEXP Ralpha, SEXP Rbeta, SEXP Reta, SEXP Rdisp,
									 SEXP RderivsAlpha, SEXP RderivsBeta, SEXP RderivsEta, SEXP RderivsDisp, SEXP RgetScores, SEXP Rscores,
									 SEXP Rpis, SEXP Rmus, SEXP logliS, SEXP logliSG,
									 SEXP Rmaxit, SEXP Rtrace, SEXP RnReport, SEXP Rabstol, SEXP Rreltol, SEXP Rconv, SEXP Rprintparams,
									 SEXP Roptimise, SEXP RloglOnly, SEXP RderivsOnly, SEXP RoptiDisp);


// functions for optimisation.
double SAM_optimise(sam_cpp_all_classes &all);
bool converged_sam( double *oldP, double *newP, const sam_opt_contr &contr, int nTot);

// functions for calculating ippm likelihood.
double optimise_function_sam(int n, double *par, void *ex);
double sam_cpp_mix_loglike(const sam_data &dat, const sam_params &params, sam_fits &fits);
void calc_mu_fits(vector<double> &fits, const sam_params &params, const sam_data &dat);
void calc_sam_loglike_SG(vector<double> &loglSG, vector<double> &fits, const sam_data &dat, const sam_params &params);
double calc_sam_loglike_S(vector<double> &fits, vector<double> const &pis, const sam_data &dat, int s);
void additive_logistic_sam(vector< double > &x,int inv, int G);
double inverse_logit(const double eta);

// functions for calculating the gradient.
void gradient_function_sam(int n, double *par, double *gr, void *ex);
void sam_cpp_mix_gradient(const sam_data &dat, const sam_params &params, sam_derivs &derivs, sam_fits &fits);
double log_poisson_deriv( const double &y, const double &mu, const double &wts);
void calc_mu_deriv( vector<double> &mu_derivs, const vector<double> &fits, const sam_data &dat, const sam_params &params);
void calc_eta_mu_deriv( vector<double> &etaDerivs, const sam_data &dat, const vector<double> &muDerivs, const vector<double> &fits);
void calc_dlog_dalpha(vector<double> &dlda, vector<double> const &mu_eta_derivs, const sam_data &dat);
void calc_dlog_dbeta(vector<double> &dldb, vector<double> const &mu_eta_derivs, const sam_data &dat);
void calc_dlog_ddispersionS(vector<double> &dldd, vector<double> const &mus, const sam_data &dat, const sam_params &params);
void calc_dlog_dpi(vector<double> &dldpi, vector<double> const &llS, vector<double> const &llSG, const sam_data &dat);
void calc_alpha_deriv( vector<double> &alphaDerivs, vector<double> const &dlogdalpha, vector<double> const &llSG, vector<double> const &llS, vector<double> const &pis, const sam_data &dat);
void calc_beta_deriv( vector<double> &betaDerivs, vector<double> const &dlogdbeta, vector<double> const &llSG, vector<double> const &llS, vector<double> const &pis, const sam_data &dat);
void calc_dispersion_deriv( vector<double> &dispDerivs, vector<double> const &dlogddispersionS, vector<double> const &llSG, vector<double> const &llS, vector<double> const &pis, const sam_data &dat);
void calc_eta_deriv( vector<double> &etaDerivs, vector<double> const &dlogdpi, vector<double> const eta, const sam_data &dat);

//distribution functions
double log_bernoulli_sam( const double &y, const double &mu);
double log_bernoulli_deriv_sam(const double &y, const double &mu);
double log_poisson_sam( const double &y, const double &mu);
double log_poisson_deriv_sam( const double &y, const double &mu);
double log_ippm_sam(const double &y, const double &mu, const double &st_sp_wts);
double log_ippm_deriv_sam( const double &y, const double &mu, const double &st_sp_wts);
double log_negative_binomial_sam( const double &y, const double &mu, const double &od);
double log_negative_binomial_deriv_mu_sam( const double &y, const double &mu, const double &od);
double log_negative_binomial_deriv_disp_sam(const double &y, const double &mu, const double &od);
double log_normal_sam( const double &y, const double &mu, const double &sig);
double log_normal_deriv_mu_sam( const double &y, const double &mu, const double &sig);
double log_normal_deriv_disp_sam( const double &y, const double &mu, const double &sig);
#endif
