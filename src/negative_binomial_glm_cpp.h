#include<R.h>
#include<Rinternals.h>
#include<Rmath.h>
#include<R_ext/Applic.h>
#undef length

#include <vector>
#include<algorithm>
#include <iostream>

#define MAT_RF(i,j,nx) i+nx*j
#define MAT_3D(i,j,k,nx,ny)  i + nx*j + k*(nx*ny)
using std::vector;         // use vector as abbreviation for std::vector

class Optimise_data_nbinom{
public:
  Optimise_data_nbinom();
  ~Optimise_data_nbinom();

  double *y, *X, *w, offset;
  int  Xr, Xc;
  vector<double> lp;
  vector<double> log_y_factorial;

  void SetVars(double *ty, double *tX, double *tw, double toffset, int tXr, int tXc);


};

// Negative binomial external calls
extern "C"  SEXP Neg_Bin(SEXP R_pars, SEXP R_X, SEXP R_y, SEXP R_w, SEXP R_offset, SEXP R_gradient, SEXP R_fitted_values);
extern "C"  SEXP Neg_Bin_Gradient(SEXP R_pars, SEXP R_X, SEXP R_y, SEXP R_w, SEXP R_offset, SEXP R_gradient);
double optimise_nbinom(int n, double *pars, void *ex);
void gradient_nbinom(int n, double *pars, double *gr, void *ex );
double NBlogl(vector<double> &pars, Optimise_data_nbinom &data );
