#include"sam_cpp2.h"

sam_params::sam_params(){};
sam_params::~sam_params(){};

void sam_params::setVals(const sam_data &dat, SEXP &Ralpha, SEXP &Rbeta, SEXP &Reta, SEXP &Rgamma, SEXP &Rtheta, SEXP &Rpowers, SEXP &RalphaPen, SEXP &RbetaPen,
						 SEXP &RpiPen, SEXP &RgammaPen, SEXP &RthetaLocatPen, SEXP &RthetaScalePen){
//	double *tmpD;

	Alpha = REAL( Ralpha);
	Beta = REAL( Rbeta);
	Eta = REAL( Reta);
	Gamma = REAL( Rgamma);
	Theta = REAL( Rtheta);
	Power = REAL( Rpowers);
	AlphaPen = *(REAL( RalphaPen));
	BetaPen = *(REAL( RbetaPen));
	PiPen = *(REAL( RpiPen)); //penalties on the pis, rather than the etas. Should represent a symetric dirichlet.
	GammaPen = *(REAL( RgammaPen));
	ThetaLocatPen = *(REAL( RthetaLocatPen));
	ThetaScalePen = *(REAL( RthetaScalePen));
	nalpha = dat.nS;
	nbeta = dat.nG*dat.nPX;
	neta = (dat.nG-1);
	ngamma = dat.nS*dat.nPW;

	if(dat.optiDisp>0)
		ntheta = dat.nS;
	else
		ntheta = 0;

	nTot = nalpha + nbeta + neta + ngamma + ntheta;
}

void sam_params::getArray(double *parArr, const sam_data &dat){
	int kount=0;
	for( int i=0; i<((dat.nS)); i++){
		parArr[kount] = Alpha[i];
		kount++;
	}
	for( int i=0; i<((dat.nG*dat.nPX)); i++){
		parArr[kount] = Beta[i];
		kount++;
	}
	for( int i=0; i<((dat.nG-1)); i++){
		parArr[kount] = Eta[i];
		kount++;
	}
	if( dat.optiPart>0){
	for( int i=0; i<((dat.nS*dat.nPW)); i++){
			parArr[kount] = Gamma[i];
			kount++;
		}
	} else {
	  kount = kount + dat.nS;
	}

	if( dat.optiDisp>0){
		for( int i=0; i<dat.nS; i++){
				parArr[kount] = Theta[i];
				kount++;
		}
	}
}

void sam_params::update( double *parArr, const sam_data &dat){

	const double etaBound = 30.0;
	int kount=0;
	for( int i=0; i<((dat.nS)); i++){
		Alpha[i] = parArr[kount];
		kount++;
	}
	for( int i=0; i<((dat.nG*dat.nPX)); i++){
		Beta[i] = parArr[kount];
		kount++;
	}
	for( int i=0; i<((dat.nG-1)); i++){
		Eta[i] = parArr[kount];
		if( Eta[i] > etaBound) Eta[i] = etaBound;
		if( Eta[i] < -etaBound) Eta[i] = -etaBound;
		kount++;
	}
	if( dat.optiPart>0){
		for( int i=0; i<((dat.nS*dat.nPW)); i++){
			Gamma[i] = parArr[kount];
			kount++;
		}
	} else {

	kount = kount + dat.nS;

	}
	if( dat.optiDisp>0){
		for( int i=0; i<dat.nS; i++){
			Theta[i] = parArr[kount];
			kount++;
		}
	}
}

