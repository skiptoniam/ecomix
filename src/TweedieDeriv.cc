#include"Tweedie.h"

//###################################
//####	Function Definitions		####
//###################################

//This function is for RCP code (probably very useful elsewhere too).
double dTweedieMu( const double &y, const double &mu, const double &phi, const double &p)
{
	double lambda, alpha, tau, muZ;
	
	lambda = ( pow( mu, (2-p))) / ( phi*(2-p));
	alpha = ( 2-p) / ( p-1);
	tau = phi*(p-1)*pow( mu, (p-1));
	muZ = alpha * tau;

	//from dTweedieDeriv (above)
	double Beta0, z1, z2;	//names are pretty obvious from supporting documentation
	vector<double> deri( 4, -9);
	vector<double> jmax( 4, -9);		//storing the js that maximise the log-desnity and its derivatives
	vector<double> jlims(8, -9);		//storing the lower and upper js for summation
	vector<double> logWjs;
	vector<double> tmpPt( 3, -9);

	if( y == 0){
		tmpPt.at(0) = -1;	//The derivatives for y=0
		tmpPt.at(1) = 0;
		tmpPt.at(2) = 0;
	}
	else{
		logWjs.clear();
		Beta0 = muZ / alpha;
		z1 = log( lambda) + alpha*log( y / muZ) + 1;
		z2 = 0.5*log( alpha) - log( 2*PI) + 1;
		findWDeriv( y, lambda, muZ, alpha, Beta0, z1, z2, jmax, jlims, deri);	//Calculate Bessel Function and its derivatives w.r.t. alpha, mu.Z and lambda. Note that logWjs and jmax, jlims and deri are all altered in this call
		tmpPt.at(0) = -1 + deri.at(1);			//updating the return object
		tmpPt.at(1) = y*alpha / pow( muZ, 2) + deri.at(2);
		tmpPt.at(2) = -y/muZ + deri.at(3);
	}
/*	for( size_t j=0; j<3; j++)
		tmpPt.at(j) *= -1;*/

	//putting it back on beta, phi, p scale
	vector<double> dBeta(3,-9);	//not worrying about p for now...
	double tmp = 0.0;

	//beta
	dBeta.at(0) = pow( mu, (1-p)) / phi;
	dBeta.at(1) = alpha*phi*( pow( ( p-1), 2))*( pow(mu,(p-2)));
	dBeta.at(2) = 0.0;
	tmp = 0.0;
	for( size_t j=0; j<3; j++)
		tmp += tmpPt.at(j) * dBeta.at(j);
//	tmp *= mu;
//	size_t kount=0;
//	for( size_t j=0; j<ncov; j++){
//		outDerivs.at(j) = tmp * X[MATREF(i,j,n)];
//		kount++;
//	}
	//phi
/*	dPhi.at(0) = -( pow( mu, (2-p)) / ( phi*phi*(2-p)));
	dPhi.at(1) = alpha*( p-1)*pow( mu, (p-1));
	dPhi.at(2) = 0.0;
	tmp = 0.0;
	for( size_t j=0; j<3; j++)
		tmp += tmpPt.at(j) * dPhi.at(j);
	outDerivs.at(kount) = tmp;
	kount++;
	outDerivs.at(kount) = -99999.99999; //a (hopefully) obviously stupid number.*/
	
	return( tmp);
}

//This function is for RCP code (probably very useful elsewhere too).
double dTweediePhi( const double &y, const double &mu, const double &phi, const double &p)
{
	double lambda, alpha, tau, muZ;
	
	lambda = ( pow( mu, (2-p))) / ( phi*(2-p));
	alpha = ( 2-p) / ( p-1);
	tau = phi*(p-1)*pow( mu, (p-1));
	muZ = alpha * tau;

	//from dTweedieDeriv (above)
	double Beta0, z1, z2;	//names are pretty obvious from supporting documentation
	vector<double> deri( 4, -9);
	vector<double> jmax( 4, -9);		//storing the js that maximise the log-desnity and its derivatives
	vector<double> jlims(8, -9);		//storing the lower and upper js for summation
	vector<double> logWjs;
	vector<double> tmpPt( 3, -9);

	if( y == 0){
		tmpPt.at(0) = -1;	//The derivatives for y=0
		tmpPt.at(1) = 0;
		tmpPt.at(2) = 0;
	}
	else{
		logWjs.clear();
		Beta0 = muZ / alpha;
		z1 = log( lambda) + alpha*log( y / muZ) + 1;
		z2 = 0.5*log( alpha) - log( 2*PI) + 1;
		findWDeriv( y, lambda, muZ, alpha, Beta0, z1, z2, jmax, jlims, deri);	//Calculate Bessel Function and its derivatives w.r.t. alpha, mu.Z and lambda. Note that logWjs and jmax, jlims and deri are all altered in this call
		tmpPt.at(0) = -1 + deri.at(1);			//updating the return object
		tmpPt.at(1) = y*alpha / pow( muZ, 2) + deri.at(2);
		tmpPt.at(2) = -y/muZ + deri.at(3);
	}
/*	for( size_t j=0; j<3; j++)
		tmpPt.at(j) *= -1;*/

	//putting it back on beta, phi, p scale
	vector<double> dPhi(3,-9);	//not worrying about p for now...
	double tmp = 0.0;

	//phi
	dPhi.at(0) = -( pow( mu, (2-p)) / ( phi*phi*(2-p)));
	dPhi.at(1) = alpha*( p-1)*pow( mu, (p-1));
	dPhi.at(2) = 0.0;
	tmp = 0.0;
	for( size_t j=0; j<3; j++)
		tmp += tmpPt.at(j) * dPhi.at(j);
	
	return( tmp);
}

void findWDeriv( double y2, double muN2, double muZ2, double alpha2, double beta2, double z12, double z22, vector<double> &jmax2, vector<double> &jlims2, vector<double> &derivsW)
//calculates the generalised Bessels function (normalising constant) via summing a finite series
{
  vector<double> logMaxs(4,-9);
  vector<double> logWjs, logdlambda, logdmuZ, logdalpha, signalpha;		//Storing the series terms for a particular y

  jmax2.at(0) = findjMax( y2, muN2, muZ2, alpha2, beta2, z12, z22, logMaxs.at(0));	//finding jmax for the bessel function
  findjMaxDerivs( y2, muN2, muZ2, alpha2, beta2, z12, z22, logMaxs, jmax2);		//finding jmax for derivatives
  findLogWjsForDeriv( y2, muN2, muZ2, alpha2, beta2, z12, z22, jmax2, jlims2, logMaxs, logWjs, logdlambda, logdmuZ, logdalpha, signalpha);	//find the summation series for derivatives and bessel - on log scale
  findEachDeriv( y2, muN2, muZ2, alpha2, beta2, z12, z22, logWjs, logdlambda, logdmuZ, logdalpha, signalpha, logMaxs, derivsW);		//find the derivatives of the bessel function w.r.t. parms
}

void findjMaxDerivs( double y3, double muN3, double muZ3, double alpha3, double beta3, double z13, double z23, vector<double> &logMaxs3, vector<double> &jmax3)
//Identify jmax using a NR algorithm for each derivative
{
  vector<double> jnr(2,1);	//putative js
  vector<double> deriv(2,0);	//derivatives at js
  vector<double> tmpMax(2,0);
  double dderiv, prev=-9, tmp;

  jnr.at(0) = jmax3.at(0);	//initialise j to be that which maximises Wj
  jnr.at(1) = jnr.at(0)+1;

  //maximum of derivative with lambda, same as derivative with muZ
  ddjOFlogdWjdLambda( jnr, z13, alpha3, deriv);		//calculate initial derivatives and stored in deriv
  while( sign( deriv.at(0))==sign( deriv.at( 1)) && jnr.at(0)!=prev) {	//NR loop
//    R_CheckUserInterrupt();
    prev = jnr.at(0);
    dderiv = deriv.at(1) - deriv.at( 0);		//approximate second derivs
    jnr.at(0) = fmax( floor( jnr.at(0) - deriv.at(0) / dderiv), 1.0);	//updating j to nearest lowest integer
    jnr.at(1) = jnr.at(0) + 1;					//updating j+1
    ddjOFlogdWjdLambda( jnr, z13, alpha3, deriv);		//update derivatives with new j
  }
  for( size_t i=0; i<2; i++)
    tmpMax.at(i) = log( jnr.at(i)) - log( muN3) + logWfun( jnr.at(i), y3, muN3, alpha3, beta3);		//derivative at j and j+1
  if( tmpMax.at(0) > tmpMax.at(1)){		//find which of j or j+1 is bigger and return it and update logWmax
    logMaxs3.at(1) = tmpMax.at(0);		//update lmabda derivative
    logMaxs3.at(2) = tmpMax.at(0) + log( muN3) - log( beta3);	//update mu.Z derivative
    jmax3.at(1) = jnr.at(0);			//update limits for lambda
    jmax3.at(2) = jnr.at(0);			//and for muZ
  }
  else{
    logMaxs3.at(1) = tmpMax.at(1);		//see `if' part of statement
    logMaxs3.at(2) = tmpMax.at(1) + log( muN3) - log( beta3);
    jmax3.at(1) = jnr.at(1);
    jmax3.at(2) = jnr.at(1);
  }

  jnr.at(0) = jmax3.at(0);	//initialise j to be that which maximises Wj
  jnr.at(1) = jnr.at(0) + 1;
  prev = -9;

  jmax3.at(3) = jmax3.at(0);	//assumes that deriv w.r.t. alpha is maximised at jmax
  tmp = 1 + log( y3 / beta3) - digamma( alpha3*jmax3.at(0));	//this could be either sign
  logMaxs3.at(3) = log( jmax3.at(0)) + logWfun( jmax3.at(0), y3, muN3, alpha3, beta3) + log( fabs( tmp));	//max abs derivative
}

void ddjOFlogdWjdLambda( const vector<double> &jnr4, double z14, double alpha4, vector<double> &deriv4)
//  d/dj (dlog( Wj / dlambda))  Function will alter deriv4
{
  for( size_t i=0; i<2; i++)
    deriv4.at(i) = 1/jnr4.at(i) + z14 - log( jnr4.at(i)+1) - alpha4*log( jnr4.at(i)) - ( 2*jnr4.at(i)+1)/( 2*( jnr4.at(i)+1)) + 1/( 2*jnr4.at(i));
}

void findLogWjsForDeriv( double y4, double muN4, double muZ4, double alpha4, double beta4, double z14, double z24, const vector<double> &jmax4, vector<double> &jlims4, const vector<double> &logMaxs4, vector<double> &logWjs4, vector<double> &logdlambda4, vector<double> &logdmuZ4, vector<double> &logdalpha4, vector<double> &signalpha4)
//finding the limits needed for the series evalulation and evaluating the series
{
  double eps = -37;	//Precision required to summing terms: corresponds to approx 10e-17
  double expEps = exp( eps);
  double jlow4, jupp4;
  double tmp;

  jlow4 = jmax4.at(0);		//starting at top of density's W and working down
  jupp4 = jmax4.at(0) + 1;	//starting at top and working out

  logWjs4.clear();	//Just in case there is residual

  logWjs4.push_back( logWfun( jlow4, y4, muN4, alpha4, beta4));			//initialising the first Wj
  logdlambda4.push_back( logWjs4.back() + log( jlow4) - log( muN4));		//first d / dlambda
  logdmuZ4.push_back( logWjs4.back() + log( jlow4) - log( beta4));		//this is actually the log negative deriv
  tmp = ( 1+log( y4/beta4)-digamma( jlow4*alpha4));
  logdalpha4.push_back( logWjs4.back() + log( jlow4) + log( fabs( tmp)));	//absolute derivative!
  signalpha4.push_back( sign( tmp));						//storing the sign
  while( checkTol( logMaxs4, logWjs4.back(), logdlambda4.back(), logdmuZ4.back(), logdalpha4.back(), eps, expEps) &&  jlow4 > 1) {		//not yet within tolerance and not at lower limit
//    R_CheckUserInterrupt();						//User interrupting from R?
    jlow4--;
    logWjs4.push_back( logWfun( jlow4, y4, muN4, alpha4, beta4));	//calculate next Wj
    logdlambda4.push_back( logWjs4.back() + log( jlow4) - log( muN4));	//calc next dlambda
    logdmuZ4.push_back( logWjs4.back() +log( jlow4) - log( beta4));		//this is actually the log negative deriv
    tmp = ( 1+log( y4/beta4)-digamma( jlow4*alpha4));
    logdalpha4.push_back( logWjs4.back() + log( jlow4) + log( fabs( tmp)));	//absolute derivative
    signalpha4.push_back( sign( tmp));						//store its sign
  }

  logWjs4.push_back( logWfun( jupp4, y4, muN4, alpha4, beta4));			//Same as previous loop but increasing j
  logdlambda4.push_back( logWjs4.back() + log( jupp4) - log( muN4));
  logdmuZ4.push_back( logWjs4.back() +log( jupp4) - log( beta4));		//this is actually the log negative deriv
  tmp = ( 1+log( y4/beta4)-digamma( jupp4*alpha4));
  logdalpha4.push_back( logWjs4.back() + log( jupp4) + log( fabs( tmp)));
  signalpha4.push_back( sign( tmp));
  while( checkTol( logMaxs4, logWjs4.back(), logdlambda4.back(), logdmuZ4.back(), logdalpha4.back(), eps, expEps)) {	//not yet within tolerance on upper side
//    R_CheckUserInterrupt();						//User interrrupt?
    jupp4++;
    logWjs4.push_back( logWfun( jupp4, y4, muN4, alpha4, beta4));	//Calculate next Wj
    logdlambda4.push_back( logWjs4.back() + log( jupp4) - log( muN4));
    logdmuZ4.push_back( logWjs4.back() + log( jupp4) - log( beta4));		//this is actually the log negative deriv
    tmp = ( 1+log( y4/beta4)-digamma( jupp4*alpha4));
    logdalpha4.push_back( logWjs4.back() + log( jupp4) + log( fabs( tmp)));
    signalpha4.push_back( sign( tmp));
  }
  jlims4.at(MATREF(0,0,4)) = jlims4.at(MATREF(1,0,4)) = jlims4.at(MATREF(2,0,4)) = jlims4.at(MATREF(3,0,4)) = jlow4;	//updating variables
  jlims4.at(MATREF(0,1,4)) = jlims4.at(MATREF(1,1,4)) = jlims4.at(MATREF(2,1,4)) = jlims4.at(MATREF(3,1,4)) = jupp4;
}

inline bool checkTol( const vector<double> &maxes, const double &currLogW, const double &currdlambda, const double &currdmuZ, const double &currdalpha, const double &eps1, const double &expeps1)
//check to see if all the conditions are met for stopping series expansion
//return TRUE if all not reached tolerance, FALSE if all have
{
  return( ( currLogW-maxes.at(0) > eps1) && ( currdlambda - maxes.at(1) > eps1) && ( currdmuZ - maxes.at(2) > eps1));
}

void findEachDeriv( const double y8, const double muN8, const double muZ8, const double alpha8, const double beta8, const double z18, const double z228, const vector<double> &logWjs8, const vector<double> &logdlambda8, const vector<double> &logdmuZ8, const vector<double> &logdalpha8, const vector<double> &signalpha8, const vector<double> &logMaxs8, vector<double> &derivsW8)
//adding up all the terms for the total term
{
  double tmpW=0, tmplambda=0, tmpmu=0, tmpalpha=0, signtmpalpha;

  for( size_t i=0; i<logWjs8.size(); i++){
    tmpW += exp( logWjs8.at(i) - logMaxs8.at(0));		//calc summands for W
    tmplambda += exp( logdlambda8.at(i) - logMaxs8.at(1));	//calc summands for d lambda
    tmpmu += exp( logdmuZ8.at(i) - logMaxs8.at(2));		//calc summands for d mu
    tmpalpha += signalpha8.at(i)*signalpha8.at(0) * exp( logdalpha8.at(i) - logMaxs8.at(3));	//calc summands for d alpha
  }
  derivsW8.at(0) = logMaxs8.at(0) + log( tmpW);		//calc log W
  derivsW8.at(1) = exp( logMaxs8.at(1) + log( tmplambda) - derivsW8.at(0));	//calc log dlambda
  derivsW8.at(2) = -exp( logMaxs8.at(2) + log( tmpmu) - derivsW8.at(0));	//calc log dmu
  signtmpalpha = sign( tmpalpha);
  derivsW8.at(3) = signtmpalpha*signalpha8.at(0) * exp( logMaxs8.at(3) + log( fabs( tmpalpha)) - derivsW8.at(0));	//calc d alpha
}
