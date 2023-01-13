#include"sam_cpp2.h"

sam_derivs::sam_derivs(){};
sam_derivs::~sam_derivs(){};

void sam_derivs::setVals( const sam_data &dat, SEXP &RderivsAlpha, SEXP &RderivsBeta, SEXP &RderivsEta, SEXP &RderivsGamma, SEXP &RderivsDelta, SEXP &RderivsTheta, SEXP &RgetScores, SEXP &Rscores){
	AlphaDeriv = REAL(RderivsAlpha);
	BetaDeriv = REAL(RderivsBeta);
	EtaDeriv = REAL(RderivsEta);
	GammaDeriv = REAL(RderivsGamma);
	DeltaDeriv = REAL(RderivsDelta);
	ThetaDeriv = REAL(RderivsTheta);
	getScoreFlag = *INTEGER(RgetScores);
	Scores = REAL(Rscores);
}

void sam_derivs::zeroDerivs( const sam_data &dat){

	for( int i=0; i<(dat.nS); i++)
		AlphaDeriv[i] = 0.0;
	for( int i=0; i<((dat.nG*dat.nPX)); i++)
		BetaDeriv[i] = 0.0;
	for( int i=0; i<(dat.nG-1); i++)
		EtaDeriv[i] = 0.0;
	if( dat.optiPart>0)
		for( int i=0; i<((dat.nS*dat.nPW)); i++)
			GammaDeriv[i] = 0.0;
	if( dat.optiAll>0)
		for( int i=0; i<(dat.nPU); i++)
			DeltaDeriv[i] = 0.0;
	if(dat.optiDisp>0)
		for( int i=0; i<(dat.nS); i++)
			ThetaDeriv[i] = 0.0;
}

void sam_derivs::updateDerivs( const sam_data &dat, const vector<double> &alphaDerivs, const vector<double> &betaDerivs,
							   const vector<double> &etaDerivs, const vector<double> &gammaDerivs, const vector<double> &deltaDerivs,
							   const vector<double> &thetaDerivs)
{
	for(int s=0; s<(dat.nS); s++){
			AlphaDeriv[s] = alphaDerivs.at(s);
			//Rprintf( " %f", Alpha[s],"\n");
			}
	for(int g=0; g<(dat.nG); g++){
		for( int p=0; p<(dat.nPX); p++){
			BetaDeriv[MATREF2D(g,p,(dat.nG))] = betaDerivs.at(MATREF2D(g,p,(dat.nG)));
			//Rprintf( " %f", Beta[MATREF2D(g,p,(dat.nG))],"\n");
			}
		}
	for(int g=0; g<(dat.nG-1); g++){
			EtaDeriv[g] = etaDerivs.at(g);
			//Rprintf( " %f", Eta[g],"\n");
			}
	if( dat.optiPart>0){
		for(int s=0; s<(dat.nS); s++){
			for( int p=0; p<(dat.nPW); p++){
				GammaDeriv[MATREF2D(s,p,(dat.nS))] = gammaDerivs.at(MATREF2D(s,p,(dat.nS)));
			}
		}
	}
	if( dat.optiAll>0){
		for( int p=0; p<(dat.nPU); p++){
			DeltaDeriv[p] = gammaDerivs.at(p);
		}
	}
	if(dat.optiDisp>0){
		for( int s=0; s<dat.nS; s++){
			ThetaDeriv[s] += thetaDerivs.at(s);
			//Rprintf( " %f", Theta[s],"\n");
		}
	}
	////Updating the score contributions for empirical information, if requested.
	if( getScoreFlag != 1)
		return;
	int k = 0;
		for( int s=0; s<dat.nS; s++){
			Scores[k++] = alphaDerivs.at(s);
		}
		for( int p=0; p<dat.nPX; p++){
			for( int g=0; g<(dat.nG); g++){
				Scores[k++] = betaDerivs.at(MATREF2D(g,p,(dat.nG)));
			}
		}
		for( int g=0; g<(dat.nG-1); g++){
				Scores[k++] = etaDerivs.at(g);
		}
		if( dat.optiPart>0){		
			for( int p=0; p<dat.nPW; p++){
				for( int s=0; s<(dat.nS); s++){
					Scores[k++] = gammaDerivs.at(MATREF2D(s, p,(dat.nS)));
				}
			}
		}
		if( dat.optiPart>0){
			for( int p=0; p<dat.nPU; p++){
				Scores[k++] = deltaDerivs.at(p);
			}
		}
		if(dat.optiDisp>0){
			for( int s=0; s<dat.nS; s++){
				Scores[k++] = thetaDerivs.at(s);
			}
		}
}

void sam_derivs::update( double *grArr, const sam_data &dat){
	int kount=0;
	for( int i=0; i<dat.nS; i++){
		AlphaDeriv[i] = grArr[kount];
		kount++;
	}
	for( int i=0; i<(dat.nG*dat.nPX); i++){
		BetaDeriv[i] = grArr[kount];
		kount++;
	}
	for( int i=0; i<(dat.nG-1); i++){
		EtaDeriv[i] = grArr[kount];
		kount++;
	}
	if( dat.optiPart>0){
		for( int i=0; i<(dat.nS*dat.nPW); i++){
			GammaDeriv[i] = grArr[kount];
			kount++;
		}
	} else {
	
	kount = kount + dat.nS;
	 	
	}	
		
	if( dat.optiAll>0){
		for( int i=0; i<(dat.nPU); i++){
			DeltaDeriv[i] = grArr[kount];
			kount++;
		}
	} else {
		
	kount = kount + dat.nPU;	
		
	}	
	if(dat.optiDisp>0){
		for( int s=0; s<dat.nS; s++){
			ThetaDeriv[s] = grArr[kount];
			kount++;
		}
	} else {
	
	kount = kount + dat.nS;
	 	
	}	
}

void sam_derivs::getArray( double *grArr, const sam_data &dat){
	int kount=0;
	for( int i=0; i<(dat.nS); i++){
		grArr[kount] = AlphaDeriv[i];
		kount++;
	}
	for( int i=0; i<(dat.nG*dat.nPX); i++){
		grArr[kount] = BetaDeriv[i];
		kount++;
	}
	for( int i=0; i<(dat.nG-1); i++){
		grArr[kount] = EtaDeriv[i];
		kount++;
	}
	if( dat.optiPart>0){
	for( int i=0; i<(dat.nS*dat.nPW); i++){
		grArr[kount] = GammaDeriv[i];
		kount++;
	}
	} else {
	
	kount = kount + dat.nS;
	 	
	}	
	if( dat.optiAll>0){
	for( int i=0; i<(dat.nPU); i++){
		grArr[kount] = DeltaDeriv[i];
		kount++;
	}
	} else {
		
	kount = kount + dat.nPU;	
		
	}	
	if(dat.optiDisp>0){
		for( int s=0; s<dat.nS; s++){
			grArr[kount] = ThetaDeriv[s];
			kount++;
		}
	} else {
	
	kount = kount + dat.nS;
	 	
	}	
}
