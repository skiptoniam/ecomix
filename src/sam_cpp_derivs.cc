#include"sam_cpp.h"

sam_derivs::sam_derivs(){};
sam_derivs::~sam_derivs(){};

void sam_derivs::setVals( const sam_data &dat, SEXP &RderivsAlpha, SEXP &RderivsBeta, SEXP &RderivsEta, SEXP &RderivsDisp, SEXP &RgetScores, SEXP &Rscores){
	Alpha = REAL(RderivsAlpha);
	Beta = REAL(RderivsBeta);
	Eta = REAL(RderivsEta);
	Disp = REAL(RderivsDisp);
	getScoreFlag = *INTEGER(RgetScores);
	Scores = REAL(Rscores);
}

void sam_derivs::zeroDerivs( const sam_data &dat){
	
	for( int i=0; i<(dat.nS); i++)
		Alpha[i] = 0.0;
	for( int i=0; i<((dat.nG*dat.nP)); i++)
		Beta[i] = 0.0;
	for( int i=0; i<(dat.nG-1); i++)
		Eta[i] = 0.0;
	if(dat.isDispersion() & dat.doOptiDisp())
		for( int i=0; i<(dat.nS); i++)
			Disp[i] = 0.0;		
}

void sam_derivs::updateDerivs( const sam_data &dat, const vector<double> &alphaDerivs, const vector<double> &betaDerivs, const vector<double> &etaDerivs, const vector<double> &dispDerivs)
{
	for(int s=0; s<(dat.nS); s++){
			Alpha[s] = alphaDerivs.at(s);
			//Rprintf( " %f", dfdAlpha[s],"\n");
			}
	for(int g=0; g<(dat.nG); g++){
		for( int p=0; p<(dat.nP); p++){
			Beta[MATREF2D(g,p,(dat.nG))] = betaDerivs.at(MATREF2D(g,p,(dat.nG)));
			//Rprintf( " %f", dfdBeta[MATREF2D(g,p,(dat.nG))],"\n");
			}
		}
    for(int g=0; g<(dat.nG-1); g++){
			Eta[g] = etaDerivs.at(g); 	
			//Rprintf( " %f", dfdEta[g],"\n");
			}
	if( dat.isDispersion() & dat.doOptiDisp())
		for( int s=0; s<dat.nS; s++)
			Disp[s] += dispDerivs.at(s);						
	
	////Updating the score contributions for empirical information, if requested.
	if( getScoreFlag != 1)
		return;
	int k = 0;
		for( int s=0; s<dat.nS; s++)
			Scores[k++] = alphaDerivs.at(s);
		for( int p=0; p<dat.nP; p++)
			for( int g=0; g<(dat.nG); g++)
				Scores[k++] = betaDerivs.at(MATREF2D(g,p,(dat.nG)));
		for( int g=0; g<(dat.nG-1); g++)
				Scores[k++] = etaDerivs.at(g);
		if( dat.isDispersion())
			for( int s=0; s<dat.nS; s++)
				Scores[k++] = dispDerivs.at(s);		

}

void sam_derivs::update( double *grArr, const sam_data &dat){
	int kount=0;
	for( int i=0; i<dat.nS; i++){
		Alpha[i] = grArr[kount];
		kount++;
	}
	for( int i=0; i<(dat.nG*dat.nP); i++){
		Beta[i] = grArr[kount];
		kount++;
	}
	for( int i=0; i<(dat.nG-1); i++){
		Eta[i] = grArr[kount];
		kount++;
	}
	if( dat.isDispersion())
		for( int s=0; s<dat.nS; s++){
			Disp[s] = grArr[kount];
			kount++;
		}

}

void sam_derivs::getArray( double *grArr, const sam_data &dat){
	int kount=0;
	for( int i=0; i<(dat.nS); i++){
		grArr[kount] = Alpha[i];
		kount++;
	}
	for( int i=0; i<(dat.nG*dat.nP); i++){
		grArr[kount] = Beta[i];
		kount++;
	}
	for( int i=0; i<(dat.nG-1); i++){
		grArr[kount] = Eta[i];
		kount++;
	}
	if( dat.isDispersion())
		for( int s=0; s<dat.nS; s++){
			grArr[kount] = Disp[s];
			kount++;
		}

}
