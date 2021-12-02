#ifndef sam_cpp2_hh
#define sam_cpp2_hh

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

//set up matref
#define MATREF2D(i,j,nx) i+nx*j
#define MATREF3D(i,j,k,nx,ny)  i + nx*j + k*(nx*ny)

// settting up the data class which should hold the data.
class sam_data {
	public:
		sam_data();
		~sam_data();
		void setVals( SEXP &Ry, SEXP &RX, SEXP &RW, SEXP &RU, SEXP &Roffset, SEXP &Rspp_wts, SEXP &Rsite_spp_wts, SEXP &Rbinsize, //SEXP &Ry_not_na, 
		 SEXP &RS, SEXP &RG, SEXP &Rpx, SEXP &Rpw, SEXP &Rpu, SEXP &RnObs, SEXP &Rdisty, SEXP &Rlinky,
		 SEXP &RoptiDisp, SEXP &RoptiPart, SEXP &RoptiAll, SEXP &RdoPenalties);
		//bool isDispersion() const;
		//bool doOptiDisp() const;
		//bool isPartial() const;
		//bool doOptiPart() const;
		//void printVals( int printX, int printW, int printy);

		int nPX,      //the number of parameters in each of the (G-1) habitat lps, same as lpar
			nPW,      //the number of parameters in the species level model
			nPU,      //the number of parameters in the bias model
			nG,	      //the number of habitats
			nS,       //the number of species
			nObs,     //the number of observations
			disty,    //the distribution code
			linky,    //which link function to use # 0 = logit, 1 = cloglog.
			optiDisp, //should the dispersion parameter be optimised (this is for negative binomial, tweedie, gaussian).
			optiPart, //should the partial sam parameters be estimated.
			optiAll, //should the all/bias sam parameters be estimated.
			doPenalties, // should we estimate penalties?
			NAnum;    //a common number to insert for NAs

		double 	*y,	//the outcome matrix, in vector form (nObs x nS)
				*X, //the design matrix in vector form (nObs x nP)
				*W, //the design matrix in vector form for the species model (nObs x npw).
				*U, //the design matrix in vector form for the bias model (nObs x npu).
				*offset, //the offset vector (length nObs)
				*spp_wts,  //the wts for the logl dependent on the species - used for the bayesian bootstrap (typically all zero and of length nObs).
				*site_spp_wts, // the wts for the ippm.
				*binsize; //Size for a binomial.
		//int 	*y_not_na; //a matrix which keeps track of NAs in ippm data. If non-ippm model all == 1. Removed 17/11/2021

};

// Set the parameters for the model.
class sam_params {
	public:
		sam_params();
		~sam_params();
		void setVals( const sam_data &dat, SEXP &Ralpha, SEXP &Rbeta, SEXP &Reta, SEXP &Rgamma, SEXP &Rdelta, SEXP &Rtheta, SEXP &Rpowers,
		SEXP &RalphaPen, SEXP &RbetaPen, SEXP &RpiPen, SEXP &RgammaPen, SEXP &RdeltaPen, SEXP &RthetaLocatPen, SEXP &RthetaScalePen);
		void getArray(double *parArr, const sam_data &dat);
		void update(double *parArr, const sam_data &dat);
		void printParms( const sam_data &dat);


		double 	*Alpha, //the species' prevalences
				*Beta,	//the archetype' free covariate params (G*xp)
				*Eta,	//the pis - mmmmm pies.
				*Gamma, //species x npw parameters form partial SAMs
				*Delta, //bias
				*Theta, //species specific dispersion parameter for negative binomial and guassian model nspp long.
				*Power, // powers for tweedie
				AlphaPen, // penalites for alpha, beta, gamma, delta, dispersion.
				BetaPen,
				PiPen,
				GammaPen,
				DeltaPen,
				ThetaLocatPen,
				ThetaScalePen;
		int nalpha, nbeta, ngamma, neta, ntheta, ndelta, nTot;
};



// setting up the derivatives which should setup, hold and update derivatives.
class sam_derivs{
	public:
		sam_derivs();
		~sam_derivs();
		void setVals( const sam_data &dat, SEXP &RderivsAlpha,  SEXP &RderivsBeta, SEXP &RderivsEta, SEXP &RderivsGamma, SEXP &RderivsDelta, SEXP &RderivsTheta, SEXP &RgetScores, SEXP &Rscores);
		void zeroDerivs( const sam_data &dat);
		void updateDerivs( const sam_data &dat, const vector<double> &alphaDerivs, const vector<double> &betaDerivs, const vector<double> &etaDerivs, const vector<double> &gammaDerivs, const vector<double> &deltaDerivs, const vector<double> &thetaDerivs);
		void update( double *grArr, const sam_data &dat);
		void getArray( double *grArr, const sam_data &dat);

        int getScoreFlag;	//Should the scores be calculated for empirical information
		double 	*Alpha, //the derivatives of logl w.r.t. alpha
				*Beta,	//the derivatives of logl w.r.t. beta
				*Eta, 	//the derivatives of logl w.r.t. eta (transformed pi)
				*Gamma, //the derivatives of logl w.r.t. species specific parameters.
				*Delta, //the derivatives of logl w.r.t. bias delta.
				*Theta, //the derivatives of logl w.r.t. dispersion parameters.
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
		void initialise( const int &nObs, const int &nG, const int &nS, const int &nPX, const int &nPW, const int &nPU, const int &NAnum);
		void zero(const int &NAnum);

		vector<double> allMus; 	//3D array for the fitted mus (note that indexing must be done with MATREF3D)
		vector<double> log_like_species_group_contrib; // 2D array of loglikes species and groups.
		vector<double> log_like_species_contrib; //vector of species logls.
		vector<double> all_derivs_mu;
		vector<double> dlogdalpha;
		vector<double> dlogdbeta;
		vector<double> dlogdgamma;
		vector<double> dlogddelta;
		vector<double> dlogdpi;
		vector<double> dlogdtheta;

};

// a wrapper around all data classes.
class sam_cpp_all_classes {
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

extern "C" SEXP species_mix_cpp(SEXP Ry, SEXP RX, SEXP RW,  SEXP RU, SEXP Roffset, SEXP Rspp_wts, SEXP Rsite_spp_wts, SEXP Rbinsize, //SEXP Ry_not_na,
								SEXP RnS, SEXP RnG, SEXP Rpx, SEXP Rpw, SEXP Rpu, SEXP RnObs, SEXP Rdisty,
								SEXP RoptiDisp, SEXP RoptiPart, SEXP RoptiAll, SEXP RdoPenalties, SEXP Rlinky,
								SEXP Ralpha, SEXP Rbeta, SEXP Reta, SEXP Rgamma, SEXP Rdelta, SEXP Rtheta, SEXP Rpowers,
								SEXP RalphaPen, SEXP RbetaPen, SEXP RpiPen, SEXP RgammaPen, SEXP RdeltaPen, SEXP RthetaLocatPen, SEXP RthetaScalePen,
								SEXP RderivsAlpha, SEXP RderivsBeta,  SEXP RderivsEta, SEXP RderivsGamma, SEXP RderivsDelta, SEXP RderivsTheta, SEXP RgetScores, SEXP Rscores,
								SEXP Rpis, SEXP Rmus, SEXP logliS, SEXP logliSG,
								SEXP Rmaxit, SEXP Rtrace, SEXP RnReport, SEXP Rabstol, SEXP Rreltol, SEXP Rconv, SEXP Rprintparams,
								SEXP Roptimise, SEXP RloglOnly, SEXP RderivsOnly);


// functions for optimisation.
double SAM_optimise(sam_cpp_all_classes &all);
bool converged_sam( double *oldP, double *newP, const sam_opt_contr &contr, int nTot);

// functions for calculating likelihoods.
double optimise_function_sam(int n, double *par, void *ex);
double sam_cpp_mix_loglike(const sam_data &dat, const sam_params &params, sam_fits &fits);
void calc_mu_fits(vector<double> &fits, const sam_params &params, const sam_data &dat);
void calc_sam_loglike_SG(vector<double> &loglSG, vector<double> &fits, const sam_data &dat, const sam_params &params);
double calc_sam_loglike_S(vector<double> &fits, vector<double> const &pis, const sam_data &dat, int s);
void additive_logistic_sam(vector< double > &x,int inv, int G);
double inverse_logit(const double eta);
double inverse_cloglog(const double eta);
double inverse_log(const double eta);

// functions for calculating the gradient.
void gradient_function_sam(int n, double *par, double *gr, void *ex);
void sam_cpp_mix_gradient(const sam_data &dat, const sam_params &params, sam_derivs &derivs, sam_fits &fits);
void calc_mu_deriv( vector<double> &mu_derivs, const vector<double> &fits, const sam_data &dat, const sam_params &params);
void calc_eta_mu_deriv( vector<double> &etaDerivs, const sam_data &dat, const vector<double> &muDerivs, const vector<double> &fits);
void calc_dlog_dalpha(vector<double> &dlda, vector<double> const &mu_eta_derivs, const sam_data &dat);
void calc_dlog_dbeta(vector<double> &dldb, vector<double> const &mu_eta_derivs, const sam_data &dat);
void calc_dlog_ddelta(vector<double> &dldd, vector<double> const &mu_eta_derivs, const sam_data &dat);
void calc_dlog_dgamma(vector<double> &dldg, vector<double> const &mu_eta_derivs, const sam_data &dat);
void calc_dlog_dtheta(vector<double> &dldt, vector<double> const &mus, const sam_data &dat, const sam_params &params);
void calc_dlog_dpi(vector<double> &dldpi, vector<double> const &llS, vector<double> const &llSG, const sam_data &dat, const sam_params &params);
void calc_alpha_deriv( vector<double> &alphaDerivs, vector<double> const &dlogdalpha, vector<double> const &llSG, vector<double> const &llS, vector<double> const &pis, const sam_data &dat);
double calc_alpha_pen( const sam_data &dat, const sam_params &params);
void calc_alpha_pen_deriv( vector<double> &alphaDerivs, const sam_data &dat, const sam_params &params);
void calc_beta_deriv( vector<double> &betaDerivs, vector<double> const &dlogdbeta, vector<double> const &llSG, vector<double> const &llS, vector<double> const &pis, const sam_data &dat);
double calc_beta_pen(  const sam_data &dat, const sam_params &params);
void calc_beta_pen_deriv( vector<double> &betaDerivs, const sam_data &dat, const sam_params &params);
void calc_gamma_deriv( vector<double> &gammaDerivs, vector<double> const &dlogdgamma, vector<double> const &llSG, vector<double> const &llS, vector<double> const &pis, const sam_data &dat);
double calc_gamma_pen(  const sam_data &dat, const sam_params &params);
void calc_gamma_pen_deriv( vector<double> &gammaDerivs, const sam_data &dat, const sam_params &params);
void calc_delta_deriv( vector<double> &deltaDerivs, vector<double> const &dlogddelta, vector<double> const &llSG, vector<double> const &llS, vector<double> const &pis, const sam_data &dat);
double calc_delta_pen(  const sam_data &dat, const sam_params &params);
void calc_delta_pen_deriv( vector<double> &deltaDerivs, const sam_data &dat, const sam_params &params);
void calc_theta_deriv( vector<double> &thetaDerivs, vector<double> const &dlogdtheta, vector<double> const &llSG, vector<double> const &llS, vector<double> const &pis, const sam_data &dat);
double calc_theta_pen( const sam_data &dat, const sam_params &params);
void calc_theta_pen_deriv(vector<double> &thetaDerivs, const sam_data &dat, const sam_params &params);
void calc_eta_deriv( vector<double> &etaDerivs, vector<double> const &dlogdpi, vector<double> const eta, const sam_data &dat);
double calc_pi_pen(const sam_data &dat, const sam_params &params);

//distribution functions
double log_bernoulli_sam( const double &y, const double &mu);
double log_bernoulli_deriv_sam(const double &y, const double &mu);
double log_poisson_sam( const double &y, const double &mu);
double log_poisson_deriv_sam( const double &y, const double &mu);
double log_negative_binomial_sam( const double &y, const double &mu, const double &od);
double log_negative_binomial_deriv_mu_sam( const double &y, const double &mu, const double &od);
double log_negative_binomial_deriv_theta_sam(const double &y, const double &mu, const double &od);
double log_tweedie_sam( const double &y, const double &mu, const double &phi, const double &p);
double log_tweedie_deriv_sam( double y, double fit, double dispParm , double p);
double log_normal_sam( const double &y, const double &mu, const double &sig);
double log_normal_deriv_mu_sam( const double &y, const double &mu, const double &sig);
double log_normal_deriv_theta_sam( const double &y, const double &mu, const double &sig);
double log_binomial_sam( const double &y, const double &mu, const double &n);
double log_binomial_deriv_sam(const double &y, const double &mu, const double &n);

#endif
