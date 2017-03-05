#include"rcp4.h"


extern "C" { SEXP RCP_C( SEXP Ry, SEXP RX, SEXP RW, SEXP Roffset, SEXP Rwts,
				     SEXP RS, SEXP RG, SEXP Rp, SEXP Rpw, SEXP RnObs, SEXP Rdisty,
					 SEXP Ralpha, SEXP Rtau, SEXP Rbeta, SEXP Rgamma, SEXP Rdisps, SEXP Rpowers, 
					 SEXP Rconc, SEXP Rsd, SEXP RsdGamma, SEXP RdispLocat, SEXP RdispScale,
					 SEXP RderivsAlpha, SEXP RderivsTau, SEXP RderivsBeta, SEXP RderivsGamma, SEXP RderivsDisps, SEXP Rscores,
					 SEXP Rpis, SEXP Rmus, SEXP RlogDens, SEXP Rlogli,
					 SEXP Rmaxit, SEXP Rtrace, SEXP RnReport, SEXP Rabstol, SEXP Rreltol, SEXP Rconv,
					 SEXP Roptimise, SEXP RloglOnly, SEXP RderivsOnly, SEXP RoptiDisp, SEXP RgetScores)
{
	allClasses all;

	//initialise the data structures -- they are mostly just pointers to REAL()s...
	all.data.setVals( Ry, RX, RW, Roffset, RS, RG, Rp, Rpw, RnObs, Rdisty, RoptiDisp, Rwts);	//read in the data
	all.parms.setVals( all.data, Ralpha, Rbeta, Rtau, Rgamma, Rdisps, Rpowers, Rconc, Rsd, RsdGamma, RdispLocat, RdispScale);	//read in the parameters
	all.derivs.setVals( all.data, RderivsAlpha, RderivsTau, RderivsBeta, RderivsGamma, RderivsDisps, RgetScores, Rscores);
	all.contr.setVals( Rmaxit, Rtrace, RnReport, Rabstol, Rreltol, Rconv);
	all.fits.initialise( all.data.nObs, all.data.nG, all.data.nS, all.data.NAnum);

	double logl = -999999;

	//doing the optimising
	if( *INTEGER(Roptimise) == 1)
		logl = ALLoptimise( all);

	//re-running to get pis and mus
	if( *INTEGER(RloglOnly) == 1)
		logl = mixLogl( all.data, all.parms, all.fits);
	//and derivatives (inlcuding scores, for empirical info, if requested)
	if( *INTEGER(RderivsOnly) == 1)
		loglDerivs( all.data, all.parms, all.derivs, all.fits);

//	for( int i=0; i<(all.data.nObs*all.parms.nTot); i++)
//		Rprintf( " %f", (REAL(Rscores))[i]);

	//bundling up things to return
	//first the fitted pis
	double *tmpPi = REAL( Rpis);
	for( int i=0; i<all.data.nObs; i++)
		for( int g=0; g<all.data.nG; g++)
			tmpPi[MATREF(i,g,all.data.nObs)] = all.fits.allPis.at(i).at(g);
	//the fitted expectations
	double *tmpMus = REAL( Rmus);
	for( size_t i=0; i<all.fits.allMus.size(); i++)
		tmpMus[i] = all.fits.allMus.at(i);
	//the log conditional densities
	double *tmpDens = REAL( RlogDens);
	for( int g=0; g<all.data.nG;g++)
		for( int i=0; i<all.data.nObs; i++)
				tmpDens[MATREF(i,g,all.data.nObs)] = all.fits.allLogDens.at(i).at(g);
	//the logl contributions
	double *tmplogls = REAL( Rlogli);
	for( int i=0; i<all.data.nObs; i++)
		tmplogls[i] = all.fits.allLogls.at(i);
	//Convergence code
	int *tmpconv = INTEGER( Rconv);
	tmpconv[0] = all.contr.ifail;
	//the logl
	SEXP Rres;	//R object to return -- it is the logl!
	Rres = PROTECT( allocVector(REALSXP,1));
    REAL( Rres)[0] = logl;
	UNPROTECT(1);
	return( Rres);
}
}

double invLogit( double x)
{
	double tmp;
	tmp = exp( x);
	tmp = tmp / (1+tmp);
	return( tmp);
}

double mixLogl( const myData &dat, const myParms &parms, myFits &fits)
{
	vector<double> logPis(dat.nG, dat.NAnum);//, pis( dat.nG, dat.NAnum);
	double res=0.0, pen=0.0, penTau=0.0, penGamma = 0.0, penDisp = 0.0;
	double wi, logli, peni;	//for coding purposes really -- note that this is NOT related to the spp design matrix W
	vector<double> wij( dat.nG, dat.NAnum);	//for coding purposes -- note that this is NOT related to the spp design matrix W
	int m;	//for coding purposes

	fits.zero( dat.NAnum);
	//calculate fitted values (constant over i)
	calcMuFits( fits.allMus, dat, parms);
	for( int i=0; i<dat.nObs; i++){
		calcLogPis( logPis, fits.allPis.at(i), dat, parms, i);
		calcLogCondDens( fits.allLogDens.at(i), fits.allMus, dat, parms, i);
		logli = calcMixSum( logPis, fits.allLogDens.at(i), wi, wij, m);
		logli *= dat.wts[i];
		res += logli;
		peni = calcPiPen( logPis, dat, parms);
		peni *= dat.wts[i]*peni;
		pen += peni;
		fits.allLogls.at(i) = logli + peni;
	}
	penTau = calcTauPen( dat, parms);
	penGamma = calcGammaPen( dat, parms);
	
	res += pen;
	res += penTau;
	res += penGamma;
	if( dat.isDispersion()){
		penDisp = calcDispPen( dat, parms);
		res += penDisp;
	}

	return( res);
}

double calcDispPen( const myData &dat, const myParms &parms)
{
	double pen = 0.0, penContr = 0.0;//, sig;
	
	for( int s=0; s<dat.nS; s++){
		penContr = - (parms.Disp[s]-parms.dispLocat) * (parms.Disp[s]-parms.dispLocat) / (2*parms.dispScale*parms.dispScale);	//dispersions are log-normally distributed (params are normally distributed)
//		Rprintf( "Species: %i, param: %f, penalty: %f, Cumulative %f \n", s, parms.Disp[s], penContr, pen);
		pen += penContr;
	}
	return( pen);
}

double calcTauPen( const myData &dat, const myParms &parms)
{
	double penTau = 0.0;
	vector<double> newTau( dat.nG*dat.nS, dat.NAnum);

	parms.getAllTaus( newTau, dat);

	for( int g=0; g<dat.nG; g++)
		for( int s=0; s<dat.nS; s++)
			penTau += -newTau.at(MATREF(g,s,dat.nG))*newTau.at(MATREF(g,s,dat.nG)) / (2*parms.sd*parms.sd);
	return( penTau);
}

double calcGammaPen( const myData &dat, const myParms &parms)
{
	double penGamma = 0.0;
	
	for( int s=0; s<dat.nS; s++)
		for( int p=0; p<dat.npw; p++)
			penGamma += -parms.Gamma[MATREF(s,p,dat.nS)]*parms.Gamma[MATREF(s,p,dat.nS)] / (2*parms.sdGamma*parms.sdGamma);
	
	return( penGamma);
}

double calcPiPen( const vector<double> &logPis, const myData &dat, const myParms &parms)
{
	double pen=0.0;

	for( int g=0; g<dat.nG; g++)
		pen += logPis.at(g);
	pen *= parms.conc;

	return( pen);
}

void calcLogPis( vector<double> &logPis, vector<double> &pis, const myData &dat, const myParms &parms, int i)
{
	vector<double> lp((dat.nG-1),0.0);
	double sumlp=0.0, sumpi=0.0;

	lp.assign( (dat.nG-1), 0.0);
	for( int k=0; k<(dat.nG-1); k++){
		for( int p=0; p<dat.np; p++)
			lp.at(k) += parms.Beta[MATREF(k,p,(dat.nG-1))] * dat.X[MATREF(i,p,dat.nObs)];
		lp.at(k) = exp( lp.at(k));
		sumlp += lp.at(k);
	}
	for( int k=0; k<(dat.nG-1); k++){
		pis.at(k) = lp.at(k) / ( 1+sumlp);
		sumpi += pis.at(k);
	}
	pis.at(dat.nG-1) = 1-sumpi;
	for( int k=0; k<dat.nG; k++)
		logPis.at(k) = log( pis.at(k));
		
	for( int k=0; k<dat.nG; k++){
		if( logPis.at(k)>=0)
			logPis.at(k) = -DBL_MIN;	//Smallest (absolute) non-zero number on your machine
		if( !R_FINITE(logPis.at(k)))
			logPis.at(k) = -DBL_MAX;	//Smallest (most negative) number on your machine
	}

}

void calcLogCondDens( vector<double> &condDens, const vector<double> &fits, const myData &dat, const myParms &parms, int i)
{
	vector<double> condDensSG( dat.nG*dat.nS, dat.NAnum);

	//calcualte the G*S log conditional densities
	for( int g=0; g<dat.nG; g++){
		for( int s=0; s<dat.nS; s++)
			switch( dat.disty){
				case 1:
					condDensSG.at(MATREF(g,s,dat.nG)) = logBernoulli( dat.y[MATREF(i,s,dat.nObs)], fits.at(MATREF3D(i,s,g,dat.nObs,dat.nS)));
					break;
				case 2:
					condDensSG.at(MATREF(g,s,dat.nG)) = logPoisson( dat.y[MATREF(i,s,dat.nObs)], fits.at(MATREF3D(i,s,g,dat.nObs,dat.nS)));
					break;
				case 3:
					condDensSG.at(MATREF(g,s,dat.nG)) = logNegBin( dat.y[MATREF(i,s,dat.nObs)], fits.at(MATREF3D(i,s,g,dat.nObs,dat.nS)), parms.Disp[s]);
					break;
				case 4:
					condDensSG.at(MATREF(g,s,dat.nG)) = logTweedie( dat.y[MATREF(i,s,dat.nObs)], fits.at(MATREF3D(i,s,g,dat.nObs,dat.nS)), parms.Disp[s], parms.Power[s]);
//					Rprintf("%f\t",condDensSG.at(MATREF(g,s,dat.nG)));
					break;
				case 5:
					condDensSG.at(MATREF(g,s,dat.nG)) = logNormal( dat.y[MATREF(i,s,dat.nObs)], fits.at(MATREF3D(i,s,g,dat.nObs,dat.nS)), parms.Disp[s]);
					break;
			}
//		Rprintf("\n");
	}
	//calculate the G log conditional densities (under independence)
	for( int g=0; g<dat.nG; g++){
		condDens.at(g) = 0.0;
		for( int s=0; s<dat.nS; s++)
			condDens.at(g) += condDensSG.at(MATREF(g,s,dat.nG));
	}
}

void calcMuFits( vector<double> &fits, const myData &dat, const myParms &parms)
{
	//fits is a G*S matrix of the fitted values if dat.npw==0 and a G*S*nObs array if npw>0
	vector<double> newTau( dat.nG*dat.nS, dat.NAnum);
	vector<double> lps(dat.nG*dat.nS, dat.NAnum);	//the nG x nS intercepts
	double lp=0.0;	//the lin pred for the gth group, sth species and ith site

	//calculate sum-to-zero taus
	parms.getAllTaus( newTau, dat);
	//calcualte the G*S*n fits
	for( int g=0; g<dat.nG; g++)
		for( int s=0; s<dat.nS; s++){
			lps.at(MATREF(g,s,dat.nG)) = parms.Alpha[s] + newTau[MATREF(g,s,dat.nG)];
			for( int i=0; i<dat.nObs; i++){
				lp = lps.at(MATREF(g,s,dat.nG)) + dat.offset[i];
				for( int p=0; p<dat.npw; p++)
					lp += dat.W[MATREF(i,p,dat.nObs)] * parms.Gamma[MATREF(s,p,dat.nS)];
				switch( dat.disty){
					case 1: 
						fits.at(MATREF3D(i,s,g,dat.nObs,dat.nS)) = invLogit( lp);
						break;
					case 2:
					case 3:
					case 4:
						fits.at( MATREF3D(i,s,g,dat.nObs,dat.nS)) = exp( lp);
						break;
					case 5:
						fits.at( MATREF3D(i,s,g,dat.nObs,dat.nS)) = lp;
						break;
				}
			}
		}
}

double logBernoulli( const double &y, const double &mu)
{
	double tmp;
	if( y==1){
		tmp = log( mu);
		return( tmp);
	}
	tmp = log( 1-mu);
	return( tmp);
}

double logPoisson( const double &y, const double &mu)
{
	double tmp;
	tmp = y * log( mu);
	tmp -= lgammafn(y+1);
	tmp -= mu;
	return( tmp);	
}

double logNegBin( const double &y, const double &mu, const double &od)
{
	double tmp, theta;
	theta = 1/exp( od);
/*	tmp1 = lgamma( (theta+y));
	tmp1 -= lgamma( theta);
	tmp1 += y*log( mu);
	tmp1 += theta*log( theta);
	tmp1 -= ( theta+y)*log(mu+theta);
	tmp1 -= lgamma(y+1);// this shouldn't be needed, but it is...  I can't understand why.  Perhaps there is something funny going on with errors cancelling?  Can't be, I think.*/
	
	tmp = dnbinom_mu(y, theta, mu, 1);
	//Equivalent to before but possibly more stable -- use this!
	
	return( tmp);	
}

double logTweedie( const double &y, const double &mu, const double &phi, const double &p)
{
	double lambda, alpha, tau, muZ, tmp, phi1;
	phi1 = exp( phi);
	lambda = R_pow( mu, (2-p)) / ( phi1*(2-p));
	alpha = ( 2-p) / ( p-1);
	tau = phi1*(p-1)*R_pow(mu,(p-1));
	muZ = alpha * tau;
	
	tmp = dTweedie( y, lambda, muZ, alpha, 1);
	return( tmp);
}

double logNormal( const double &y, const double &mu, const double &sig)
{
	double tmp, sig1;
	sig1 = exp( sig);
	
	tmp = -log( sig1);
	tmp -= (y-mu)*(y-mu)/(2*sig1*sig1);
	return( tmp);
}

double calcMixSum( const vector<double> &logPis, const vector<double> &condDens, double &wi, vector<double> &wij, int &maxg)
{
//function used in both mixLogl and loglDerivs
	vector<double> summands(logPis.size(), 0.0);
	double max, res=0.0;

	max = logPis.at(0) + condDens.at(0);
	maxg = 0;
	for( size_t i=0; i<logPis.size(); i++){
		summands.at(i) = logPis.at(i) + condDens.at(i);
		if( summands.at(i) > max){
			max = summands.at(i);
			maxg = i;
		}
	}
	wi = 0.0;
	for( size_t g=0; g<summands.size(); g++){
		wij.at(g) = exp( summands.at(g) - max);
		wi += wij.at(g);
	}
	res = log( wi);
	res += max;

	return( res);
}

void loglDerivs( const myData &dat, const myParms &parms, myDerivs &derivs, myFits &fits)
{
	vector<double> logPis(dat.nG, dat.NAnum), pis( dat.nG, dat.NAnum);
	vector<double> logCondDens( dat.nG, dat.NAnum);
	double wi, tmp1;
	vector<double> wij( dat.nG, dat.NAnum);
	int m;	//location of the maximum group contribution
	
	vector<double> muDerivsI( dat.nG*dat.nS, dat.NAnum);
	vector<double> etaDerivsI( dat.nG*dat.nS, dat.NAnum);
	vector<double> alphaDerivsI( dat.nS, dat.NAnum);
	vector<double> tauDerivsI( (dat.nG-1)*dat.nS, dat.NAnum);
	vector<double> piDerivsI( dat.nG, dat.NAnum);
	vector<double> betaDerivsI( (dat.nG-1)*dat.np, dat.NAnum);
	vector<double> gammaDerivsI( dat.nS*dat.npw, dat.NAnum);
	vector<double> dispDerivsI( dat.nS, dat.NAnum);

	vector<double> tmpPiDerivs( dat.nG*dat.nS, 0.0);

	calcMuFits( fits.allMus, dat, parms);
	derivs.zeroDerivs( dat);
	for( int i=0; i<dat.nObs; i++){
		//calculating the w_i and the {w_ig}
		calcLogPis( logPis, pis, dat, parms, i);
		calcLogCondDens( logCondDens, fits.allMus, dat, parms, i);
		tmp1 = calcMixSum( logPis, logCondDens, wi, wij, m);
		//calc deriv w.r.t. mu and the eta (all of them)
		calcDerivMu( muDerivsI, fits.allMus, dat, parms, wi, wij, m, i);
		calcDerivEtaMu( etaDerivsI, dat, muDerivsI, fits.allMus, i);
		//calc deriv w.r.t. alpha and then tau and then gamma
		calcAlphaDeriv( alphaDerivsI, etaDerivsI, dat);
		calcTauDeriv( tauDerivsI, etaDerivsI, dat, parms);
		calcGammaDeriv( gammaDerivsI, etaDerivsI, dat, parms, i);
		//calc deriv w.r.t. beta
		calcPiDeriv( piDerivsI, dat, parms, pis, wi, wij, m);
		calcBetaDeriv( betaDerivsI, piDerivsI, pis, dat, i);
		//calc deriv w.r.t. dispersions
		calcDispDeriv( dispDerivsI, fits.allMus, dat, parms, wi, wij, m, i);	//if dat.disty specifies no dispersion then a vector of zeros is returned (and not used)
		//put on weights
		weightDerivs( alphaDerivsI, tauDerivsI, gammaDerivsI, betaDerivsI, dispDerivsI, dat, i);
		//update derivatives
		derivs.updateDerivs( dat, alphaDerivsI, tauDerivsI, betaDerivsI, gammaDerivsI, dispDerivsI, i);	//if dat.disty specifies no dispersion then no place for disp derivs (and are not updated, of course)
	}
	calcTauPenDeriv( tauDerivsI, dat, parms);
	calcGammaPenDeriv( gammaDerivsI, dat, parms);
	calcDispPenDeriv( dispDerivsI, dat, parms);

	alphaDerivsI.assign(alphaDerivsI.size(), 0.0);
	betaDerivsI.assign(betaDerivsI.size(), 0.0);

	derivs.updateDerivs( dat, alphaDerivsI, tauDerivsI, betaDerivsI, gammaDerivsI, dispDerivsI, -1);
	(void)tmp1;	//tmp1 is not used again, this is just a little trick to avoid a compile warning.
}

void weightDerivs( vector<double> &alphaDerivsI, vector<double> &tauDerivsI, vector<double> &gammaDerivsI, vector<double> &betaDerivsI, vector<double> &dispDerivsI, const myData &dat, const int &i)
{
	for( size_t s=0; s<alphaDerivsI.size(); s++)
		alphaDerivsI.at(s) *= dat.wts[i];
	for( size_t s=0; s<tauDerivsI.size(); s++)
		tauDerivsI.at(s) *= dat.wts[i];
	for( size_t s=0; s<gammaDerivsI.size(); s++)
		gammaDerivsI.at(s) *= dat.wts[i];
	for( size_t s=0; s<betaDerivsI.size(); s++)
		betaDerivsI.at(s) *= dat.wts[i];
	for( size_t s=0; s<dispDerivsI.size(); s++)
		dispDerivsI.at(s) *= dat.wts[i];		
}

void calcDispDeriv( vector<double> &dispDerivsI, const vector<double> &fits, const myData &dat, const myParms &parms, const double &wi, const vector<double> &wij, const int &m, const int &i)
{
	vector<double> tmpDerivs( dat.nS*dat.nG, 0.0);
	double summand;

	if( !dat.isDispersion())
		return;	//nothing to do here, move along please
			
	for( int s=0; s<dat.nS; s++)
		for( int g=0; g<dat.nG; g++){
			switch( dat.disty){
				case 3:
					tmpDerivs.at(MATREF(g,s,dat.nG)) = logNegBinDispDer( dat.y[MATREF(i,s,dat.nObs)], fits.at( MATREF3D(i,s,g,dat.nObs, dat.nS)), parms.Disp[s]);
					break;
				case 4:
					tmpDerivs.at(MATREF(g,s,dat.nG)) = logTweedieDispDer( dat.y[MATREF(i,s,dat.nObs)], fits.at( MATREF3D(i,s,g,dat.nObs, dat.nS)), parms.Disp[s], parms.Power[s]);
					break;
				case 5:
					tmpDerivs.at(MATREF(g,s,dat.nG)) = logNormalDispDer( dat.y[MATREF(i,s,dat.nObs)], fits.at( MATREF3D(i,s,g,dat.nObs, dat.nS)), parms.Disp[s]);
					break;
			}
		}
	
	dispDerivsI.assign(dispDerivsI.size(), 0.0);
	for( int s=0; s<dat.nS; s++){
		summand = 0.0;
		for( int g=0; g<dat.nG; g++)
			summand += wij.at(g) * (tmpDerivs.at(MATREF(g,s,dat.nG)) - tmpDerivs.at(MATREF(m,s,dat.nG)));
		dispDerivsI.at(s) = tmpDerivs.at(MATREF(m,s,dat.nG)) + summand / wi;
	}
}

void calcDispPenDeriv( vector<double> &dispDerivsI, const myData &dat, const myParms &parms)
{
	dispDerivsI.assign(dispDerivsI.size(), 0.0);
	if( dat.isDispersion())
		for( int s=0; s<dat.nS; s++)
			dispDerivsI.at(s) = -(parms.Disp[s]-parms.dispLocat)/(parms.dispScale*parms.dispScale);
}

double logNegBinDispDer( double y, double fit, double dispParm)
{
	double theta, sig, res=0.0;
	
	sig = exp( dispParm);
	theta = 1 / sig;
	
	res = digamma( theta+y);
	res -= digamma( theta);
	res += 1 + log( theta) - log( fit+theta) - (theta+y)/(theta+fit);
	res /= -sig*sig;	//for the change of variable sig --> r
	res *= sig;	//for the change of variable dispParm --> sig

	return( res);
}

double logTweedieDispDer( double y, double fit, double dispParm , double p)
{
	double phi, tmp;
	
	phi = exp( dispParm);
	
	tmp = dTweediePhi( y, fit, phi, p);
	tmp *= phi;
	
	return( tmp);	
}

double logNormalDispDer( double y, double fit, double dispParm)
{
	double sig, res=0.0;
	
	sig = exp( dispParm);
	res = (y-fit)*(y-fit) / (sig*sig*sig);
	res -= 1/sig;
	res *= sig; //for the change of variable dispParm --> sig
	
	return( res);
}

void calcDerivMu( vector<double> &muDerivs, const vector<double> &fits, const myData &dat, const myParms &parms, const double wi, const vector<double> &wij, const int &m, const int &i)
{
	//derivatives of conditional density w.r.t. its mean
	//muDerivs is a GxS matrix of first derivatives
	vector<double> tmpDerivs( dat.nG*dat.nS, 0.0);

	for( int g=0; g<dat.nG; g++)
		for( int s=0; s<dat.nS; s++)
			switch( dat.disty){
				case 1: 
					tmpDerivs.at(MATREF(g,s,dat.nG)) = logBernDer( dat.y[MATREF(i,s,dat.nObs)], fits.at(MATREF3D(i,s,g,dat.nObs,dat.nS)));
					break;
				case 2:
					tmpDerivs.at(MATREF(g,s,dat.nG)) = logPoissonDer( dat.y[MATREF(i,s,dat.nObs)], fits.at(MATREF3D(i,s,g,dat.nObs,dat.nS)));
					break;
				case 3:
					tmpDerivs.at(MATREF(g,s,dat.nG)) = logNegBinLocatDer( dat.y[MATREF(i,s,dat.nObs)], fits.at(MATREF3D(i,s,g,dat.nObs,dat.nS)), parms.Disp[s]);
					break;
				case 4:
					tmpDerivs.at(MATREF(g,s,dat.nG)) = 	dTweedieMu( dat.y[MATREF(i,s,dat.nObs)], fits.at(MATREF3D(i,s,g,dat.nObs,dat.nS)), exp( parms.Disp[s]), parms.Power[s]);
					break;
				case 5:
					tmpDerivs.at(MATREF(g,s,dat.nG)) = logNormalLocatDer( dat.y[MATREF(i,s,dat.nObs)], fits.at(MATREF3D(i,s,g,dat.nObs,dat.nS)), parms.Disp[s]);
					break;					
			}
		

	for( int s=0; s<dat.nS; s++){
		muDerivs.at(MATREF(m,s,dat.nG)) = 0.0;
		for( int g=0; g<dat.nG; g++){
			if( m!=g){
				muDerivs.at(MATREF(g,s,dat.nG)) = tmpDerivs.at(MATREF(g,s,dat.nG)) * wij.at(g) / wi;
				muDerivs.at(MATREF(m,s,dat.nG)) -= tmpDerivs.at(MATREF(m,s,dat.nG)) * wij.at(g) / wi;
			}
			else
				muDerivs.at(MATREF(g,s,dat.nG)) += tmpDerivs.at(MATREF(g,s,dat.nG));
		}
	}
}

void calcDerivEtaMu( vector<double> &etaDerivsI, const myData &dat, const vector<double> &muDerivsI, const vector<double> &fits, const int &i)
{
	//derivatives of conditional densities w.r.t. its linear predictor(s)
	for( int g=0; g<dat.nG; g++)
		for( int s=0; s<dat.nS; s++)
			switch( dat.disty){
				case 1: 
					etaDerivsI.at(MATREF(g,s,dat.nG)) = fits.at(MATREF3D(i,s,g,dat.nObs,dat.nS)) * (1-fits.at(MATREF3D(i,s,g,dat.nObs,dat.nS))) * muDerivsI.at(MATREF(g,s,dat.nG));	//logit link
					break;
				case 2:
				case 3:
				case 4:
					etaDerivsI.at(MATREF(g,s,dat.nG)) = fits.at(MATREF3D(i,s,g,dat.nObs,dat.nS)) * muDerivsI.at(MATREF(g,s,dat.nG));	//log link
					break;
				case 5:
					etaDerivsI.at(MATREF(g,s,dat.nG)) = muDerivsI.at(MATREF(g,s,dat.nG));	//identity link
					break;
			}
}

double logBernDer( double y, double mu)
{
	double tmp, negOne = -1.0;
	if( y==1){
		tmp = 1/mu;
		return( tmp);
	}
	if( y==0){
		tmp = -1/(1-mu);
		return( tmp);
	}
	return( log( negOne));	//to give an error
}

double logPoissonDer( const double &y, const double &mu)
{
	double tmp;
	tmp = y/mu;
	tmp -= 1;
	
	return( tmp);
}

double logNegBinLocatDer( const double &y, const double &mu, const double &od)
{
	double tmp, theta;
	theta = 1/exp( od);
	tmp = -(theta+y)/(theta+mu);
	tmp += y/mu;
	return( tmp);
}

double logNormalLocatDer( const double &y, const double &mu, const double &sig)
{
	double tmp, sig1;
	
	sig1 = exp( sig);	
	tmp = (y-mu) / (sig1*sig1);
	return( tmp);
}

void calcAlphaDeriv( vector<double> &alphaDerivsI, const vector<double> &etaDerivs, const myData &dat)
{
	alphaDerivsI.assign(alphaDerivsI.size(), 0.0);
	for( int s=0; s<dat.nS; s++)
		for( int g=0; g<dat.nG; g++)
			alphaDerivsI.at(s) += etaDerivs.at(MATREF(g,s,dat.nG));
}

void calcTauDeriv( vector<double> &tauDerivsI, const vector<double> &etaDerivs, const myData &dat, const myParms &parms)
{
	vector<double> newTau( dat.nG*dat.nS, dat.NAnum);

	tauDerivsI.assign(tauDerivsI.size(), 0.0);
	for( int s=0; s<dat.nS; s++){
		for( int g=0; g<(dat.nG-1); g++){
			tauDerivsI.at(MATREF(g,s,(dat.nG-1))) = etaDerivs.at(MATREF(g,s,dat.nG));
			tauDerivsI.at(MATREF(g,s,(dat.nG-1))) -= etaDerivs.at(MATREF((dat.nG-1),s,dat.nG));
		}
	}
}

void calcGammaDeriv( vector<double> &gammaDerivsI, const vector<double> &etaDerivs, const myData &dat, const myParms &parms, const int &i)
{
	//This function could be a future source of problems.
	//Reading that comment, I am not sure why.  Perhaps it is to do with non-identifiability? 10-may-2015.
	gammaDerivsI.assign(gammaDerivsI.size(), 0.0);
	for( int s=0; s<dat.nS; s++)
		for( int p=0; p<dat.npw;p++)
			for( int g=0; g<dat.nG; g++)
				gammaDerivsI.at(MATREF(s,p,dat.nS)) += etaDerivs.at(MATREF(g,s,dat.nG)) * dat.W[MATREF(i,p,dat.nObs)];
}

void calcTauPenDeriv( vector<double> &tauDerivsI, const myData &dat, const myParms &parms)
{
	vector<double> newTau( dat.nG*dat.nS, dat.NAnum);

	tauDerivsI.assign(tauDerivsI.size(), 0.0);
	parms.getAllTaus( newTau, dat);
	for( int s=0; s<dat.nS; s++)
		for( int g=0; g<(dat.nG-1); g++){
			tauDerivsI.at(MATREF(g,s,(dat.nG-1))) += - ( newTau.at( MATREF(g,s,dat.nG)) - newTau.at( MATREF((dat.nG-1),s,dat.nG))) / (parms.sd*parms.sd);
	}
}

void calcGammaPenDeriv( vector<double> &gammaDerivsI, const myData &dat, const myParms &parms)
{
	gammaDerivsI.assign(gammaDerivsI.size(), 0.0);
	for( int s=0; s<dat.nS; s++)
		for( int p=0; p<dat.npw; p++)
			gammaDerivsI.at( MATREF(s,p,dat.nS)) += - parms.Gamma[MATREF(s,p,dat.nS)] / (parms.sdGamma*parms.sdGamma);
}

void calcPiDeriv( vector<double> &piDerivsI, const myData &dat, const myParms &parms, const vector<double> &pis, const double wi, const vector<double> &wig, int m)
{
	vector<double> wigDerivs(dat.nG, 0.0);

	for( int g=0; g<dat.nG; g++){
		if( g!=m)
			piDerivsI.at(g) = wig.at(g) / (wi*pis.at(g));
	}
	piDerivsI.at(m) = 1 / pis.at(m);
	for( int g=0; g<dat.nG; g++)
		if( g!=m)
			piDerivsI.at(m) -= wig.at(g) / (wi*pis.at(m));

	//The penalty derivative
	for( int g=0; g<dat.nG; g++)
		piDerivsI.at(g) += parms.conc / pis.at(g);
}

void calcBetaDeriv( vector<double> &betaDerivsI, const vector<double> &piDerivsI, const vector<double> &pis, const myData &dat, int i)
{
	vector<double> dpideta( (dat.nG*(dat.nG-1)), 0.0);
	vector<double> dldeta( (dat.nG-1), 0.0);

	//derivs of pi w.r.t. eta (all (G-1) of 'em)
	for( int g=0; g<(dat.nG-1); g++){
		dpideta.at(MATREF(g,g,dat.nG)) += pis.at(g);
		for( int h=0; h<(dat.nG-1); h++)
			dpideta.at(MATREF(g,h,dat.nG)) += -pis.at(g) * pis.at(h);
	}
	for( int g=0; g<(dat.nG-1); g++){
		dpideta.at(MATREF((dat.nG-1),g,dat.nG)) = 0.0;
		for( int h=0; h<(dat.nG-1); h++)
			dpideta.at(MATREF((dat.nG-1),g,dat.nG)) -= dpideta.at(MATREF(h,g,dat.nG));
	}

	//deriv is dldpi X dpideta X detadbeta_h, for lp number h (of course)
	//logl_i w.r.t eta first: a 1x(G-1) vector
	for( int h=0; h<(dat.nG-1); h++)
		for( int g=0; g<dat.nG; g++)
			dldeta.at(h) += piDerivsI.at(g)*dpideta.at(MATREF(g,h,dat.nG));

	//now for each of the beta_h vectors
	betaDerivsI.assign( betaDerivsI.size(), 0.0);
	for( int h=0; h<(dat.nG-1); h++)
		for( int p=0; p<dat.np; p++)
			betaDerivsI.at(MATREF(h,p,(dat.nG-1))) += dldeta.at(h)*dat.X[MATREF(i,p,dat.nObs)];

}

double ALLoptimise( allClasses &all)
{
	double *vmminGrad, *vmminParms, *oldParms;	//arrays to pass to vmmin
	double *vmminParmsIn;	//del
	vmminParms = (double *) R_alloc(all.parms.nTot,sizeof(double));
	vmminParmsIn = (double *) R_alloc(all.parms.nTot,sizeof(double));	//Del
	oldParms = (double *) R_alloc(all.parms.nTot,sizeof(double));
	vmminGrad = (double *) R_alloc(all.parms.nTot,sizeof(double));
	int *myMask;
	vector<int> vecMask(all.parms.nTot, 1);
	double vmminLogl[1];

//	Rprintf( "Quasi-Newton iterations\n");
	all.parms.getArray( oldParms, all.data);
	myMask = &vecMask[0];
	//optimise
	all.parms.getArray( vmminParms, all.data);
	all.parms.getArray( vmminParmsIn, all.data);	//del
	vmmin( all.parms.nTot, vmminParms, vmminLogl, optimise_function, gradient_function, all.contr.maxitQN, all.contr.traceQN, myMask, all.contr.abstol, all.contr.reltol, all.contr.nReport, &all, &all.contr.fnKount, &all.contr.grKount, &all.contr.ifail);
//	nmmin( all.parms.nTot, vmminParmsIn, vmminParms, vmminLogl, optimise_function, &all.contr.ifail, all.contr.abstol, all.contr.reltol, &all, 1.0, 0.5, 2.0, all.contr.traceQN, &all.contr.fnKount, all.contr.maxitQN);	
	//update parameters
	all.parms.update( vmminParms, all.data);
	gradient_function(all.parms.nTot, vmminParms, vmminGrad, &all);
	all.derivs.update( vmminGrad, all.data);

	return( vmminLogl[0]);
}

bool converged( double *oldP, double *newP, const myOptContr &contr, int nTot)
{
	double tmp, eps=1e-5;

	for( int i=0; i<nTot; i++){
		tmp = fabs(newP[i]-oldP[i]);
		tmp /= fabs(oldP[i])+eps;
		if( tmp > contr.reltol)
			return( false);
	}
	return( true);

}

double optimise_function(int n, double *par, void *ex)
{
	allClasses *all = (allClasses *) ex;
	double logl;

	all->parms.update( par, all->data);
	logl = mixLogl( all->data, all->parms, all->fits);
	
	return( (0.0-logl));
}

void gradient_function(int n, double *par, double *gr, void *ex)
{
	allClasses *all = (allClasses *) ex;

	loglDerivs( all->data, all->parms, all->derivs, all->fits);

	all->derivs.getArray(gr, all->data);

	for( int i=0; i<n; i++)
		gr[i] = 0-gr[i];
}
