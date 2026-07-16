#ifndef sam_cpp_pred_v1_hh
#define sam_cpp_pred_v1_hh

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
#define MATREF4D(i,j,k,l,nx,ny,nz)  i + nx*j + k*(nx*ny) + l*(nx*ny*nz)

// Classes
// setting up the data class which should hold the data.
class sam_pred_data {
public:
  sam_pred_data();
  ~sam_pred_data();
  void setVals( SEXP &RX, SEXP &RW, SEXP &Roffset,
                SEXP &RS, SEXP &RG, SEXP &Rpx, SEXP &Rpw, SEXP &RnObs,
                SEXP &Rdisty, SEXP &Rlinky,  SEXP &Rtype, SEXP &Rpredtype, SEXP &Rnboot);

  int nPX,  //the number of parameters in the archetype model
  nPW,      //the number of parameters in the species level model
  nG,	      //the number of habitats
  nS,       //the number of species
  nObs,     //the number of observations
  disty,    //the distribution code
  linky,    //which link function to use # 0 = logit, 1 = cloglog.
  type,     //which prediction version? 0=response, 1=link.
  predtype, //which prediction level? 0=archetype, 1=species.
  nboot;    //Number of bootstraps if being used.


  double 	*X, //the design matrix in vector form (nObs x nP)
  *W, //the design matrix in vector form for the species model (nObs x npw).
  *offset; //the offset vector (length nObs)

};

// setting up the data class which should hold the data.
class sam_pred_params {
public:
  sam_pred_params();
  ~sam_pred_params();
void setParams(const sam_pred_data &dat,
               SEXP &Ralpha,
               SEXP &Rbeta,
               SEXP &Rgamma,
               SEXP &Rtau);

  double 	*Alpha, //the species' prevalences
  *Beta,	//the archetype' free covariate params (G*xp)
  *Gamma, //species x npw parameters form partial SAMs
  *Tau; // Derived taus from pis and loglike (G*S)


  int nalpha, nbeta, ngamma, ntau, nTot;


};

// setting up the data class which should hold the data.
class sam_pred_bootparams {
public:
  sam_pred_bootparams();
  ~sam_pred_bootparams();

  void setParams(const sam_pred_data &dat,
                     SEXP &Rbootalpha,
                     SEXP &Rbootbeta,
                     SEXP &Rbootgamma,
                     SEXP &Rboottau);

  double 	*bootAlpha, //the species' prevalences
  *bootBeta,	//the archetype' free covariate params (G*xp)
  *bootGamma, //species x npw parameters form partial SAMs
  *bootTau; // Derived taus from pis and loglike, one set PER bootstrap draw (nboot*G*S)

  int nbootalpha, nbootbeta, nbootgamma, nboottau, nbootTot;

};

class sam_pred_classes{
public:
  sam_pred_classes();
  ~sam_pred_classes();

  sam_pred_data data;
  sam_pred_params params;
  sam_pred_bootparams bootparams;

};

////////////////////////////////////////////////////////
/////////////	Function Definitions	////////////////////
////////////////////////////////////////////////////////

extern "C" SEXP sam_cpp_pred(SEXP RX, SEXP RW, SEXP Roffset,
                             SEXP RG, SEXP RS, SEXP RnObs, SEXP Rpx, SEXP Rpw,
                             SEXP Rdisty, SEXP Rlinky, SEXP Rtype, SEXP Rpredtype,
                             SEXP Ralpha, SEXP Rbeta, SEXP Rgamma, SEXP Rtau,
                             SEXP Rbootalpha, SEXP Rbootbeta, SEXP Rbootgamma, SEXP Rboottau,
                             SEXP Rnboot);

// functions to help with prediction
void pt_predict_fun(const sam_pred_data &dat, const sam_pred_params &params, vector<double> &preds);
void boot_predict_fun(const sam_pred_data &dat, const sam_pred_bootparams &bootparams, vector<double> &bootpreds);
void pt_predict_species_fun(const sam_pred_data &dat, const sam_pred_params &params, vector<double> &preds);
void boot_predict_species_fun(const sam_pred_data &dat, const sam_pred_bootparams &bootparams, vector<double> &bootpreds);

#endif
