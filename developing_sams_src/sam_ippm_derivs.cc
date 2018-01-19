#include"sam_ippm.h"

sam_ippm_derivs::sam_ippm_derivs(){};
sam_ippm_derivs::~sam_ippm_derivs(){};

void sam_ippm_derivs::setVals( const sam_ippm_data &dat, SEXP &RderivsAlpha, SEXP &RderivsBeta, SEXP &RderivsPi){
	Alpha = REAL(RderivsAlpha);
	Beta = REAL(RderivsBeta);
	Pi = REAL(RderivsPi);
	//getScoreFlag = *INTEGER(RgetScores);
	//Scores = REAL(Rscores);
}

void sam_ippm_derivs::zeroDerivs( const sam_ippm_data &dat){
	
	for( int i=0; i<(dat.nS-1); i++)
		Alpha[i] = 0.0;
	for( int i=0; i<((dat.nG*dat.nP)-1); i++)
		Beta[i] = 0.0;
	for( int i=0; i<(dat.nG-1); i++)
		Pi[i] = 0.0;	
}

void sam_ippm_derivs::updateDerivs( const sam_ippm_data &dat, const vector<double> &alphaDerivs, const vector<double> &betaDerivs, const vector<double> &piDerivs)
{
	for(int s=0; s<(dat.nS-1); s++)
			Alpha[s] += alphaDerivs.at(s);
	for(int g=0; g<(dat.nG-1); g++)
		for( int p=0; p<(dat.nP-1); p++)
			Beta[MATREF2D(g,p,(dat.nG-1))] += betaDerivs.at(MATREF2D(g,p,(dat.nG-1)));
    for(int g=0; g<(dat.nG-1); g++)
			Pi[g] += piDerivs.at(g); 					
	
	////Updating the score contributions for empirical information, if requested.
	//if( getScoreFlag != 1)
		//return;
	//int k = 0;
	//if( i != -1){
		//for( int s=0; s<dat.nS; s++)
			//Scores[MATREF2D(i,k++,dat.nObs)] = alphaDerivs.at(s);
		//for( int p=0; p<dat.nP; p++)
			//for( int g=0; g<(dat.nG-1); g++)
				//Scores[MATREF2D(i,k++,dat.nObs)] = betaDerivs.at(MATREF2D(g,p,(dat.nG-1)));
		//for( int g=0; g<(dat.nG-1); g++)
				//Scores[MATREF2D(i,k++,dat.nObs)] = piDerivs.at(g);
	//}
	//else{
		//for( int i=0; i<dat.nObs; i++){
			//k = 0;
			//for( int s=0; s<dat.nS; s++)
				//Scores[MATREF2D(i,k++,dat.nObs)] += dat.wts[i]*alphaDerivsI.at(s) / dat.nObs;
			//for( int p=0; p<dat.nP; p++)
				//for( int g=0; g<(dat.nG-1); g++)
					//Scores[MATREF2D(i,k++,dat.nObs)] += dat.wts[i]*betaDerivsI.at(MATREF2D(g,p,(dat.nG-1))) / dat.nObs;
			//for( int s=0; s<dat.nS; s++)
				//for( int g=0; g<(dat.nG-1); g++)
					//Scores[MATREF2D(i,k++,dat.nObs)] += dat.wts[i]*tauDerivsI.at(MATREF2D(g,s,(dat.nG-1))) / dat.nObs;		
		//}
	//}
}

void sam_ippm_derivs::update( double *grArr, const sam_ippm_data &dat)
{
	int kount=0;
	for( int i=0; i<dat.nS; i++){
		Alpha[i] = grArr[kount];
		kount++;
	}
	for( int i=0; i<((dat.nG*dat.nP)-1); i++){
		Beta[i] = grArr[kount];
		kount++;
	}
	for( int i=0; i<(dat.nG-1); i++){
		Pi[i] = grArr[kount];
		kount++;
	}

}

void sam_ippm_derivs::getArray( double *grArr, const sam_ippm_data &dat)
{
	int kount=0;
	for( int i=0; i<(dat.nS-1); i++){
		grArr[kount] = Alpha[i];
		kount++;
	}
	for( int i=0; i<((dat.nG*dat.nP)-1); i++){
		grArr[kount] = Beta[i];
		kount++;
	}
	for( int i=0; i<(dat.nG-1); i++){
		grArr[kount] = Pi[i];
		kount++;
	}

}
