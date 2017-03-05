#include"rcp4.h"

myDerivs::myDerivs(){};
myDerivs::~myDerivs(){};

void myDerivs::setVals( const myData &dat, SEXP &RderivsAlpha, SEXP &RderivsTau, SEXP &RderivsBeta, SEXP &RderivsGamma, SEXP &RderivsDisp, SEXP &RgetScores, SEXP &Rscores)
{
	Alpha = REAL( RderivsAlpha);
	Tau = REAL( RderivsTau);
	Beta = REAL( RderivsBeta);
	Gamma = REAL( RderivsGamma);
	Disp = REAL( RderivsDisp);
	getScoreFlag = *INTEGER(RgetScores);
	Scores = REAL( Rscores);
}

void myDerivs::zeroDerivs( const myData &dat)
{
	for( int i=0; i<dat.nS; i++)
		Alpha[i] = 0.0;
	for( int i=0; i<((dat.nG-1)*dat.nS); i++)
		Tau[i] = 0.0;
	for( int i=0; i<((dat.nG-1)*dat.np); i++)
		Beta[i] = 0.0;
	for( int i=0; i<(dat.nS*dat.npw); i++)
		Gamma[i] = 0.0;
	if( dat.isDispersion())
		for( int i=0; i<dat.nS; i++)
			Disp[i] = 0.0;
}

void myDerivs::updateDerivs( const myData &dat, const vector<double> &alphaDerivsI, const vector<double> &tauDerivsI, const vector<double> &betaDerivsI, const vector<double> &gammaDerivsI, const vector<double> &dispDerivsI, const int &i)
{
	for( int s=0; s<dat.nS; s++)
		Alpha[s] += alphaDerivsI.at(s);
	for( int g=0; g<(dat.nG-1); g++)
		for( int s=0; s<dat.nS; s++)
			Tau[MATREF(g,s,(dat.nG-1))] += tauDerivsI.at(MATREF(g,s,(dat.nG-1)));
	for( int g=0; g<(dat.nG-1); g++)
		for( int p=0; p<dat.np; p++)
			Beta[MATREF(g,p,(dat.nG-1))] += betaDerivsI.at(MATREF(g,p,(dat.nG-1)));
	for( int s=0; s<dat.nS; s++)
		for( int p=0; p<dat.npw; p++)
			Gamma[MATREF(s,p,dat.nS)] += gammaDerivsI.at(MATREF(s,p,dat.nS));
	if( dat.isDispersion() & dat.doOptiDisp())
		for( int s=0; s<dat.nS; s++)
			Disp[s] += dispDerivsI.at(s);
	
	//Updating the score contributions for empirical information, if requested.
	if( getScoreFlag != 1)
		return;
	int k = 0;
	if( i != -1){
		for( int s=0; s<dat.nS; s++)
			Scores[MATREF(i,k++,dat.nObs)] = alphaDerivsI.at(s);
		for( int s=0; s<dat.nS; s++)
			for( int g=0; g<(dat.nG-1); g++)
				Scores[MATREF(i,k++,dat.nObs)] = tauDerivsI.at(MATREF(g,s,(dat.nG-1)));
		for( int p=0; p<dat.np; p++)
			for( int g=0; g<(dat.nG-1); g++)
				Scores[MATREF(i,k++,dat.nObs)] = betaDerivsI.at(MATREF(g,p,(dat.nG-1)));
		for( int p=0; p<dat.npw; p++)
			for( int s=0; s<dat.nS; s++)
				Scores[MATREF(i,k++, dat.nObs)] = gammaDerivsI.at(MATREF(s,p,dat.nS));
		if( dat.isDispersion())
			for( int s=0; s<dat.nS; s++)
				Scores[MATREF(i,k++, dat.nObs)] = dispDerivsI.at(s);
	}
	else{
		for( int i=0; i<dat.nObs; i++){
			k = 0;
			for( int s=0; s<dat.nS; s++)
				Scores[MATREF(i,k++,dat.nObs)] += dat.wts[i]*alphaDerivsI.at(s) / dat.nObs;
			for( int s=0; s<dat.nS; s++)
				for( int g=0; g<(dat.nG-1); g++)
					Scores[MATREF(i,k++,dat.nObs)] += dat.wts[i]*tauDerivsI.at(MATREF(g,s,(dat.nG-1))) / dat.nObs;
			for( int p=0; p<dat.np; p++)
				for( int g=0; g<(dat.nG-1); g++)
					Scores[MATREF(i,k++,dat.nObs)] += dat.wts[i]*betaDerivsI.at(MATREF(g,p,(dat.nG-1))) / dat.nObs;
			for( int p=0; p<dat.npw; p++)
				for( int s=0; s<dat.nS; s++)
					Scores[MATREF(i,k++,dat.nObs)] += dat.wts[i]*gammaDerivsI.at(MATREF(s,p,dat.nS)) / dat.nObs;
			if( dat.isDispersion())
				for( int s=0; s<dat.nS; s++)
					Scores[MATREF(i,k++,dat.nObs)] += dat.wts[i]*dispDerivsI.at(s) / dat.nObs;
		}
	}
}

void myDerivs::update( double *grArr, const myData &dat)
{
	int kount=0;
	for( int i=0; i<dat.nS; i++){
		Alpha[i] = grArr[kount];
		kount++;
	}
	for( int i=0; i<((dat.nG-1)*dat.nS); i++){
		Tau[i] = grArr[kount];
		kount++;
	}
	for( int i=0; i<((dat.nG-1)*dat.np); i++){
		Beta[i] = grArr[kount];
		kount++;
	}
	for( int i=0; i<(dat.nS*dat.npw); i++){
		Gamma[i] = grArr[kount];
		kount++;
	}
	if( dat.isDispersion())
		for( int s=0; s<dat.nS; s++){
			Disp[s] = grArr[kount];
			kount++;
		}
}

void myDerivs::getArray( double *grArr, const myData &dat)
{
	int kount=0;
	for( int i=0; i<dat.nS; i++){
		grArr[kount] = Alpha[i];
		kount++;
	}
	for( int i=0; i<((dat.nG-1)*dat.nS); i++){
		grArr[kount] = Tau[i];
		kount++;
	}
	for( int i=0; i<((dat.nG-1)*dat.np); i++){
		grArr[kount] = Beta[i];
		kount++;
	}
	for( int i=0; i<(dat.nS*dat.npw); i++){
		grArr[kount] = Gamma[i];
		kount++;
	}
	if( dat.isDispersion())
		for( int s=0; s<dat.nS; s++){
			grArr[kount] = Disp[s];
			kount++;
		}
}
