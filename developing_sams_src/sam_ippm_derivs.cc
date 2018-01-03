#include"sam.h"

sam_ippm_derivs::sam_ippm_derivs(){};
sam_ippm_derivs::~sam_ippm_derivs(){};

void sam_ippm_derivs::setVals( const sam_ippm_data &dat, SEXP &RderivsAlpha, SEXP &RderivsBeta, SEXP &RderivsTau, SEXP &RgetScores, SEXP &Rscores)
{
	Alpha = REAL( RderivsAlpha);
	Beta = REAL( RderivsBeta);
	Tau = REAL( RderivsTau);
		//Disp = REAL( RderivsDisp);
	getScoreFlag = *INTEGER(RgetScores);
	Scores = REAL( Rscores);
}

void sam_ippm_derivs::zeroDerivs( const sam_ippm_data &dat)
{
	for( int i=0; i<dat.nS; i++)
		Alpha[i] = 0.0;
	for( int i=0; i<((dat.nG-1)*dat.np); i++)
		Beta[i] = 0.0;
	for( int i=0; i<((dat.nG-1)*dat.nS); i++)
		Tau[i] = 0.0;	
	//if( dat.isDispersion())
		//for( int i=0; i<dat.nS; i++)
			//Disp[i] = 0.0;
}

void sam_ippm_derivs::updateDerivs( const sam_ippm_data &dat, const vector<double> &alphaDerivsI, const vector<double> &tauDerivsI, const vector<double> &betaDerivsI, const vector<double> &gammaDerivsI, const vector<double> &dispDerivsI, const int &i)
{
	for( int s=0; s<dat.nS; s++)
		Alpha[s] += alphaDerivsI.at(s);
	for( int g=0; g<(dat.nG-1); g++)
		for( int s=0; s<dat.nS; s++)
			Tau[MATREF2D(g,s,(dat.nG-1))] += tauDerivsI.at(MATREF2D(g,s,(dat.nG-1)));
	for( int g=0; g<(dat.nG-1); g++)
		for( int p=0; p<dat.np; p++)
			Beta[MATREF2D(g,p,(dat.nG-1))] += betaDerivsI.at(MATREF2D(g,p,(dat.nG-1)));		
	//if( dat.isDispersion() & dat.doOptiDisp())
		//for( int s=0; s<dat.nS; s++)
			//Disp[s] += dispDerivsI.at(s);
	
	//Updating the score contributions for empirical information, if requested.
	if( getScoreFlag != 1)
		return;
	int k = 0;
	if( i != -1){
		for( int s=0; s<dat.nS; s++)
			Scores[MATREF2D(i,k++,dat.nObs)] = alphaDerivsI.at(s);
		for( int p=0; p<dat.np; p++)
			for( int g=0; g<(dat.nG-1); g++)
				Scores[MATREF2D(i,k++,dat.nObs)] = betaDerivsI.at(MATREF2D(g,p,(dat.nG-1)));
		for( int s=0; s<dat.nS; s++)
			for( int g=0; g<(dat.nG-1); g++)
				Scores[MATREF2D(i,k++,dat.nObs)] = tauDerivsI.at(MATREF2D(g,s,(dat.nG-1)));
		//if( dat.isDispersion())
			//for( int s=0; s<dat.nS; s++)
				//Scores[MATREF2D(i,k++, dat.nObs)] = dispDerivsI.at(s);
	}
	else{
		for( int i=0; i<dat.nObs; i++){
			k = 0;
			for( int s=0; s<dat.nS; s++)
				Scores[MATREF2D(i,k++,dat.nObs)] += dat.wts[i]*alphaDerivsI.at(s) / dat.nObs;
			for( int p=0; p<dat.np; p++)
				for( int g=0; g<(dat.nG-1); g++)
					Scores[MATREF2D(i,k++,dat.nObs)] += dat.wts[i]*betaDerivsI.at(MATREF2D(g,p,(dat.nG-1))) / dat.nObs;
			for( int s=0; s<dat.nS; s++)
				for( int g=0; g<(dat.nG-1); g++)
					Scores[MATREF2D(i,k++,dat.nObs)] += dat.wts[i]*tauDerivsI.at(MATREF2D(g,s,(dat.nG-1))) / dat.nObs;		
			//if( dat.isDispersion())
				//for( int s=0; s<dat.nS; s++)
					//Scores[MATREF2D(i,k++,dat.nObs)] += dat.wts[i]*dispDerivsI.at(s) / dat.nObs;
		}
	}
}

void sam_ippm_derivs::update( double *grArr, const sam_ippm_data &dat)
{
	int kount=0;
	for( int i=0; i<dat.nS; i++){
		Alpha[i] = grArr[kount];
		kount++;
	}
	for( int i=0; i<((dat.nG-1)*dat.np); i++){
		Beta[i] = grArr[kount];
		kount++;
	}
	for( int i=0; i<((dat.nG-1)*dat.nS); i++){
		Tau[i] = grArr[kount];
		kount++;
	}
	//if( dat.isDispersion())
		//for( int s=0; s<dat.nS; s++){
			//Disp[s] = grArr[kount];
			//kount++;
		//}
}

void sam_ippm_derivs::getArray( double *grArr, const sam_ippm_data &dat)
{
	int kount=0;
	for( int i=0; i<dat.nS; i++){
		grArr[kount] = Alpha[i];
		kount++;
	}
	for( int i=0; i<((dat.nG-1)*dat.np); i++){
		grArr[kount] = Beta[i];
		kount++;
	}
	for( int i=0; i<((dat.nG-1)*dat.nS); i++){
		grArr[kount] = Tau[i];
		kount++;
	}
	//if( dat.isDispersion())
		//for( int s=0; s<dat.nS; s++){
			//grArr[kount] = Disp[s];
			//kount++;
		//}
}
