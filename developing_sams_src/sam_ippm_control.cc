#include"sam_ippm.h"

sam_ippm_opt_contr::sam_ippm_opt_contr(){};
sam_ippm_opt_contr::~sam_ippm_opt_contr(){};

void sam_ippm_opt_contr::setVals( const SEXP &Rmaxit, const SEXP &Rtrace, const SEXP &RnReport, const SEXP &Rabstol, const SEXP &Rreltol, SEXP &Rconv)
{

	maxitQN = *(INTEGER( Rmaxit));
	traceQN = *(INTEGER(Rtrace));
	nReport = *(INTEGER(RnReport));
	abstol = *(REAL(Rabstol));
	reltol = *(REAL(Rreltol));
	conv = INTEGER( Rconv);

	denomEps = 1e-5;
}
