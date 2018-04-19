#include"sam_ippm.h"

sam_ippm_derivs::sam_ippm_derivs(){};
sam_ippm_derivs::~sam_ippm_derivs(){};

void sam_ippm_derivs::setVals( const sam_ippm_data &dat, SEXP &RderivsAlpha, SEXP &RderivsBeta, SEXP &RderivsEta, SEXP &RgetScores, SEXP &Rscores){
	Alpha = REAL(RderivsAlpha);
	Beta = REAL(RderivsBeta);
	Eta = REAL(RderivsEta);
	getScoreFlag = *INTEGER(RgetScores);
	Scores = REAL(Rscores);
}

void sam_ippm_derivs::zeroDerivs( const sam_ippm_data &dat){
	
	for( int i=0; i<(dat.nS); i++)
		Alpha[i] = 0.0;
	for( int i=0; i<((dat.nG*dat.nP)); i++)
		Beta[i] = 0.0;
	for( int i=0; i<(dat.nG-1); i++)
		Eta[i] = 0.0;	
}

void sam_ippm_derivs::updateDerivs( const sam_ippm_data &dat, const vector<double> &alphaDerivs, const vector<double> &betaDerivs, const vector<double> &etaDerivs)
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

}

void sam_ippm_derivs::update( double *grArr, const sam_ippm_data &dat)
{
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

}

void sam_ippm_derivs::getArray( double *grArr, const sam_ippm_data &dat)
{
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

}
