#include"rcp4.h"

myOptContr::myOptContr(){};
myOptContr::~myOptContr(){};

void myOptContr::setVals( const SEXP &Rmaxit, const SEXP &Rtrace, const SEXP &RnReport, const SEXP &Rabstol, const SEXP &Rreltol, SEXP &Rconv)
{

	maxitQN = *(INTEGER( Rmaxit));
	traceQN = *(INTEGER(Rtrace));
	nReport = *(INTEGER(RnReport));
	abstol = *(REAL(Rabstol));
	reltol = *(REAL(Rreltol));
	conv = INTEGER( Rconv);

	denomEps = 1e-5;
}
