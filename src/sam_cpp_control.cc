#include"sam_cpp.h"

sam_opt_contr::sam_opt_contr(){};
sam_opt_contr::~sam_opt_contr(){};

void sam_opt_contr::setVals( const SEXP &Rmaxit, const SEXP &Rtrace, const SEXP &RnReport, const SEXP &Rabstol, 
const SEXP &Rreltol, SEXP &Rconv,const SEXP &Rprintparams)
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
