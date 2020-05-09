#include"SAM_AD_hessian.hh"
AD<double> lgammaAD( AD<double> xx){
	//log gamma function according to numerical recipes
	int j;
	AD<double> x,y,tmp,ser;
	static const double cof[6]= {76.18009172947146,-86.50532032941677,24.01409824083091,-1.231739572450155,0.1208650973866179e-2,-0.5395239384953e-5};

	y = x = xx;
	tmp = x+5+0.5;
	tmp -= (x + 0.5) * log( tmp);
	ser=1.000000000190015;
	for( j=0;j<6;j++){
		y = y+1;
		ser += cof[j]/y;
	}
	return -tmp+log( 2.5066282746310005*ser/x);
}

AD<double> ldnegbin(double yi, AD<double> mui, AD<double> phii){
	AD<double> tmp=0;

	tmp = lgammaAD(phii+yi) - lgammaAD(phii) - lgamma(yi+1) + yi*log(mui) + phii*log(phii) - (phii+yi)*log(mui+phii);
//	tmp = lgammaAD(phii+yi) - lgammaAD(phii) - lgamma(yi+1) + yi*log( mui/(mui+phii)) + phii*log( phii/(mui+phii));
	return( tmp);
}

// AD<double> ldnegbin(double yi, AD<double> mui, AD<double> phii){
//   AD<double> tmp=0, theta;
//   theta = 1/exp(phii);
//   tmp = dnbinom_mu(yi, theta, mui, 1);
//   return( tmp);
// }


void invMultLogit( vector< AD<double> > &pi2, const vector< AD<double> > alpha2, int nG2){
	vector< AD<double> > tmp(nG2, 0.0);
	AD<double> sum=0;

	for( size_t gg=0; gg<(nG2-1); gg++){
		tmp.at( gg) = exp( alpha2.at(gg));
		sum += tmp.at(gg);
	}
	tmp.at( nG2-1) = 1.0;
	sum += tmp.at( nG2-1);

	for( size_t gg=0; gg<nG2; gg++)
		pi2.push_back(tmp.at( gg) / sum);
}

AD<double> ldbeta( AD<double> theta1, double *thetaRange2, double penParm2){
	AD<double> thetaStar, tmp;

	thetaStar = ( theta1 - thetaRange2[0]) / (thetaRange2[1]-thetaRange2[0]);
	tmp = lgammafn( 2*penParm2) - 2*lgammafn(penParm2) + (penParm2-1)*(log(thetaStar) + log(1-thetaStar));

	return( tmp);
}

AD<double> calcLogl( vector< AD<double> > allPars1, double *y1, double *offy1,
 double *XnoMix1, double *XMix1, int nS1, int nG1, int nn1, int npNoMix1,
  int npMix1, double *thetaRange1, double penParm1){

	vector< AD<double> > sppEta, grpEta, pi, sppPars, grpPars, alpha, phi, sppLogl(nG1,0.0), summand(nG1,0.0), summandAlt(nG1, 0.0), wt(nG1, 0.0);
	AD<double> logl=0.0, tmp, tmp1, mu, avSummand, sppContr=0.0, penalty=0.0;
	double tau=100;//magic number
	int kount=0;

	double tmp_nb, diff, maxDiff=0.0;
	AD<double> tmp_myNB;

	for( size_t gg=0; gg<(nG1-1); gg++){
		alpha.push_back(allPars1.at(kount));
		kount++;
	}
	invMultLogit(pi, alpha, nG1);
	for( size_t pp=0; pp<npNoMix1; pp++){
		for( size_t ss=0; ss<nS1; ss++){
			sppPars.push_back(allPars1.at(kount));
			kount++;
		}
	}
	for( size_t ss=0; ss<nS1; ss++){
		phi.push_back(allPars1.at(kount));
		kount++;
	}
	for( size_t pp=0; pp<npMix1; pp++){
		for( size_t gg=0; gg<nG1; gg++){
			grpPars.push_back(allPars1.at(kount));
			kount++;
		}
	}
	logl = 0.0;
	//calc eta for each spp into a nn1 by nS1 matrix
	for( size_t ss=0; ss<nS1; ss++){
		for( size_t ii=0; ii<nn1; ii++){
			tmp = 0.0;
			for( size_t jj=0; jj<npNoMix1; jj++)
				tmp += XnoMix1[MATREF(ii,jj,nn1)] * sppPars.at( MATREF(ss,jj,nS1));
			sppEta.push_back( tmp);
		}
	}
	//calc eta for each grp into an nn1 by nG1 matrix
	for( size_t gg=0; gg<nG1; gg++){
		for( size_t ii=0; ii<nn1; ii++){
			tmp = 0.0;
			for( size_t jj=0; jj<npMix1; jj++)
				tmp += XMix1[MATREF(ii,jj,nn1)] * grpPars.at( MATREF(gg,jj,nG1));
			grpEta.push_back( tmp);
		}
	}
	logl = 0.0;
	//for each spp and then for each grp...: calc eta, mu and logl contribution (on log scale, summands are scaled by average summand)
	for( size_t ss=0; ss<nS1; ss++){
		avSummand = 0.0;
		for( size_t gg=0; gg<nG1; gg++){
			sppLogl.at(gg) = 0.0;
			for( size_t ii=0; ii<nn1; ii++){
				mu = exp( sppEta.at( MATREF(ii,ss,nn1)) + grpEta.at( MATREF( ii,gg,nn1)) + offy1[ii]);
				sppLogl.at(gg) += ldnegbin( y1[MATREF(ii,ss,nn1)], mu, phi.at(ss));
			}
			summand.at(gg) = log( pi.at(gg)) + sppLogl.at(gg);
		}
		tmp = 0.0;
		tau = 100;	//magic number
		for( size_t gg=0; gg<nG1; gg++){
			summandAlt.at(gg) = exp( summand.at(gg) / tau);
			tmp += summandAlt.at(gg);
		}
		for( size_t gg=0; gg<nG1; gg++)
			wt.at(gg) = summandAlt.at(gg) / tmp;
		avSummand = 0.0;
		for( size_t gg=0; gg<nG1; gg++)
			avSummand += wt.at(gg) * summand.at(gg);
		tmp = 0.0;
		for( size_t gg=0; gg<nG1; gg++)
			tmp += exp( summand.at(gg) - avSummand);
		sppContr = avSummand + log( tmp);
		logl += sppContr;
		penalty += ldbeta( phi.at(ss), thetaRange1, penParm1);
	}

	return( logl + penalty);

}


extern "C" {
	SEXP calcDerHess( SEXP Ry, SEXP Roffy, SEXP RallPars, SEXP RXnoMix, SEXP RXMix, SEXP RS, SEXP RG, SEXP Rn, SEXP RnpNoMix, SEXP RnpMix, SEXP RthetaRange, SEXP RpenParm){
	//calculates the logl, scores and Hessian for the partialSAM.nb model.
	//Ry contains the vectorised n X S data
	// RallPars contains all the parameters
	//RXnoMix contains the (vectorised) design matrix for the spp-specific part
	//RXMix contains the (vectorised) design matrix for the group-specific part
	//RS is the number of spp
	//RG is the number of groups
	//Rn is the number of obs
	//RnpNoMix is the number of parameters in each spp's model
	//RnpMix is the number of parameters in each group's model
	//RthetaRange is the allowable range for the dispersion parameter
	//RpenParm is the penalty parameter for the lubricating dispersion penalty.

	int nS, nG, nn, npNoMix, npMix;
	double *y, *XnoMix, *XMix, *offy, *thetaRange, penParm;
	vector< AD<double> > allPars, logl(1,0);


	nS = *(INTEGER( RS));
	nG = *(INTEGER( RG));
	nn = *(INTEGER( Rn));
	npNoMix = *(INTEGER( RnpNoMix));
	npMix = *(INTEGER( RnpMix));
	y = REAL( Ry);
	XnoMix = REAL( RXnoMix);
	XMix = REAL( RXMix);
	offy = REAL( Roffy);
	thetaRange = REAL( RthetaRange);
	penParm = *(REAL(RpenParm));
	for( size_t ii=0; ii < LENGTH(RallPars); ii++)
		allPars.push_back( REAL(RallPars)[ii]);

	CppAD::ADFun<double> F;
	CppAD::Independent( allPars);	//set the parameters to be the AD-parameters
	logl.at(0) = calcLogl( allPars, y, offy, XnoMix, XMix, nS, nG, nn, npNoMix, npMix, thetaRange, penParm);
	F.Dependent( logl);
	F.optimize();


	SEXP Rres;//R object to return -- it is the logl
	Rres = allocVector(REALSXP,1 + LENGTH(RallPars) * ( 1 + LENGTH(RallPars)));	//giving it some memory
	double *R_pointer = REAL(Rres);	//so we can look at something that isn't so scary
	R_pointer[0] = Value( logl.at(0));	//assigning the value we calculated

	//scores
	vector<double> pt, grads, Hess;
	for( size_t ii=0; ii<LENGTH(RallPars); ii++)
		pt.push_back(REAL(RallPars)[ii]);
	grads = F.RevOne( pt, 0);
	Hess = F.Hessian( pt, 0);

	size_t kount=1;	//already have logl in first position
	for( size_t ii=0; ii<LENGTH( RallPars); ii++){
		R_pointer[kount] = grads.at( ii);
		kount++;
	}
	kount = LENGTH( RallPars)+1;
	for( size_t ii=0; ii<(LENGTH( RallPars)*LENGTH(RallPars)); ii++){
		R_pointer[kount] = Hess.at( ii);
		kount++;
	}

	return( Rres);	//returning the logl



}}
