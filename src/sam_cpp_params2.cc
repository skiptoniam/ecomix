#include"sam_cpp2.h"

sam_params::sam_params(){};
sam_params::~sam_params(){};

void sam_params::setVals(const sam_data &dat, SEXP &Ralpha, SEXP &Rbeta, SEXP &Reta, SEXP &Rgamma, SEXP &Rdelta, SEXP &Rtheta, SEXP &Rpowers, SEXP &RalphaPen, SEXP &RbetaPen,
						 SEXP &RpiPen, SEXP &RgammaPen, SEXP &RdeltaPen, SEXP &RthetaLocatPen, SEXP &RthetaScalePen){
//	double *tmpD;

	Alpha = REAL( Ralpha);
	Beta = REAL( Rbeta);
	Eta = REAL( Reta);
	Gamma = REAL( Rgamma);
	Delta = REAL( Rdelta);
	Theta = REAL( Rtheta);
	Power = REAL( Rpowers);
	AlphaPen = *(REAL( RalphaPen));
	BetaPen = *(REAL( RbetaPen));
	PiPen = *(REAL( RpiPen)); //penalties on the pis, rather than the etas. Should represent a symetric dirichlet. 
	GammaPen = *(REAL( RgammaPen));
	DeltaPen = *(REAL( RdeltaPen));
	ThetaLocatPen = *(REAL( RthetaLocatPen));
	ThetaScalePen = *(REAL( RthetaScalePen));
	nalpha = dat.nS;
	nbeta = dat.nG*dat.nPX;
	neta = (dat.nG-1);
	ngamma = dat.nS*dat.nPW;
	ndelta = dat.nPU;

	if(dat.isDispersion())
		ntheta = dat.nS;
	else
		ntheta = 0;

	nTot = nalpha + nbeta + ngamma + neta + ntheta + ndelta;
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
	
	if( dat.optiAll>0){
	for( int i=0; i<((dat.nPU)); i++){
			parArr[kount] = Delta[i];
			kount++;
		}
	} else {
	   kount = kount + dat.nPU;	
	}
	if( dat.isDispersion() ){
		for( int i=0; i<dat.nS; i++){
				parArr[kount] = Theta[i];
				kount++;
		}
	}
}

void sam_params::update( double *parArr, const sam_data &dat){
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
	if( dat.optiAll>0){
		for( int i=0; i<(dat.nPU); i++){
			Delta[i] = parArr[kount];
			kount++;
		}
	} else {
	
	kount = kount + dat.nPU;
	 	
	}	
	if( dat.isDispersion() & dat.doOptiDisp()){
		for( int i=0; i<dat.nS; i++){
			Theta[i] = parArr[kount];
			kount++;
		}
	} 
}

//void sam_params::printParms( const sam_data &dat){

	//Rprintf( "ALPHA:\n");
	//for( int i=0; i<dat.nS; i++){
		//Rprintf( "%3.2f\t", Alpha[i]);
	//}
	//Rprintf( "\n");
	//Rprintf( "BETA:\n");
	//for( int g=0; g<(dat.nG); g++){
		//for( int i=0; i<dat.nPX; i++)
			//Rprintf( "%3.2f\t", Beta[MATREF2D(g,i,(dat.nG))]);
			//Rprintf( "\n");
	//}
	//Rprintf( "PI (transformed Pi):\n");
	//for( int g=0; g<(dat.nG-1); g++){
		//Rprintf( "%3.2f\t", Eta[g]);
		//Rprintf( "\n");
	//}
	//if( dat.optiPart>0){
	//Rprintf( "GAMMA:\n");
		//for( int s=0; s<(dat.nS); s++){
			//for( int i=0; i<dat.nPW; i++)
				//Rprintf( "%3.2f\t", Gamma[MATREF2D(s,i,(dat.nS))]);
				//Rprintf( "\n");
		//}
	//}
	//if( dat.optiAll>0){
	//Rprintf( "DELTA:\n");
		//for( int i=0; i<dat.nPW; i++)
				//Rprintf( "%3.2f\t", Delta[i]);
				//Rprintf( "\n");
	//}
	//if( dat.isDispersion()==true){
		//Rprintf("DISPERSION:\n");
		//for( int i=0; i<dat.nS; i++)
		//Rprintf( "%3.2f\t", Theta[i]);
		//Rprintf( "\n");
	//}

//}
