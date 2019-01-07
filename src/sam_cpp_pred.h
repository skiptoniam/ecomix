#ifndef sam_cpp_pred_hh
#define sam_cpp_pred_hh

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

////////////////////////////////////////////////////////
/////////////	Function Definitions	////////////////////
////////////////////////////////////////////////////////
extern "C" SEXP SAM_predict_C(SEXP Ry, SEXP RX, SEXP Roffset, SEXP Rspp_wts,
						   SEXP Rsite_spp_wts, SEXP Ry_not_na, SEXP Rtaus,
				           SEXP RnS, SEXP RnG, SEXP Rp, SEXP RnObs, SEXP Rdisty,
				           SEXP Ralpha, SEXP Rbeta, SEXP Reta, SEXP Rdisp,
				           SEXP RalphaBoot, SEXP RbetaBoot, //SEXP RetaBoot,// SEXP RdispBoot,
				           SEXP Rnboot,
				           SEXP Rspp_pt_preds, SEXP Rgrp_pt_preds,
				           SEXP Rspp_boot_preds, SEXP Rgrp_boot_preds,
				           SEXP RoptiDisp);

#endif
