#include"sam_bernoulli_sp_ints.h"

sam_bernoulli_sp_ints_opt_contr::sam_bernoulli_sp_ints_opt_contr(){};
sam_bernoulli_sp_ints_opt_contr::~sam_bernoulli_sp_ints_opt_contr(){};

void sam_bernoulli_sp_ints_opt_contr::setVals( const SEXP &Rmaxit, const SEXP &Rtrace, const SEXP &RnReport, 
											   const SEXP &Rabstol, const SEXP &Rreltol, SEXP &Rconv, const SEXP &Rprintparams)
{

	maxitQN = *(INTEGER( Rmaxit));
	traceQN = *(INTEGER(Rtrace));
	nReport = *(INTEGER(RnReport));
	abstol = *(REAL(Rabstol));
	reltol = *(REAL(Rreltol));
	conv = INTEGER( Rconv);
	printparams = *INTEGER(Rprintparams);

	denomEps = 1e-5;
}
