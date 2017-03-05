#include"Tweedie.h"

//###################################
//####	Function Definitions		####
//###################################

double dTweedie( double y, double muN, double muZ, double alpha, int myLog)
//calculates the (log-)density for the parameters given
{
  double beta0, z1, z2, logW, logf;	//names are pretty obvious from supporting documentation

  if( y == 0)
    logf = -muN;
  else{
    beta0 = muZ / alpha;
    z1 = log( muN) + alpha*log( y / muZ) + 1;
    z2 = 0.5*log( alpha) - log( 2*PI) + 1;
    logW = findW( y, muN, muZ, alpha, beta0, z1, z2);	//Calculate Bessel Function
    logf = -y / beta0 - muN - log( y) + logW;
  }
  if( myLog !=1)
	 logf = exp( logf);
  return logf;
}


double findW( double y2, double muN2, double muZ2, double alpha2, double beta2, double z12, double z22)
//calculates the generalised Bessels function (normalising constant) via summing a finite series
{

  double jmax, jlow, jupp, logWmax, res;
  vector<double> logWjs2;

  jmax = findjMax( y2, muN2, muZ2, alpha2, beta2, z12, z22, logWmax);				//find the maximum j via NR.  Note that logWmax is updated during this procedure.
  findlogWjs( y2, muN2, muZ2, alpha2, beta2, z12, z22, jmax, jlow, jupp, logWmax, logWjs2);	//find the series of Wj elements. Note that logWjs2, jlow, and jupp are all altered by this function

  res = 0;
  for( size_t i=0; i<logWjs2.size(); i++)	//sum the series with alleviation for over-flow
    res += exp( logWjs2.at( i) - logWmax);
  res = logWmax + log( res);

  return( res);
}

double logWfun( double j1, double y1, double muN1, double alpha1, double beta1)
//Calculate individual Wj
{
  double tmp;
  tmp = j1*( log( muN1) + alpha1*log( y1/beta1)) - lgamma( j1+1) - lgamma( j1*alpha1);
  return tmp;
}

double logWderivApprox( double j6, double z16, double alpha6)
//Calculate the derivative of the approximate log Wj
{
  double tmp;
  tmp = z16 - log( j6+1) - alpha6*log( j6) - ( 2*j6+1) / ( 2*( j6+1)) + 1 / ( 2*j6);
  return( tmp);
}

double findjMax( double y3, double muN3, double muZ3, double alpha3, double beta3, double z13, double z23, double &logWmax3)
//Idenity jmax using a NR algorithm
{
  vector<double> jnr(2,1);	//putative js
  vector<double> deriv(2,0);	//derivatives at js
  vector<double> logW(2,0);
  double dderiv, prev=-9;

//not needed?  theta = ( 2+alpha3) / ( 1+alpha3);					//needed for Dunn and Smyth jmax
//Not needed?  phi = pow( muN3*muZ3, 2-theta) / ( muN3*( 2-theta));
//  jnr.at(0) = 1;
  jnr.at(0) = fmax( floor( exp( ( z13 - 1) / ( 1+ alpha3))), 1);
  jnr.at(1) = jnr.at(0)+1;

  deriv.at(0) = logWderivApprox( jnr.at(0), z13, alpha3);	//calculate initial derivatives
  deriv.at(1) = logWderivApprox( jnr.at(1), z13, alpha3);
//  while( sign( deriv.at(0))==sign( deriv.at( 1)) && jnr.at(0)!=prev) {	//NR loop
  while( deriv.at(0)*deriv.at(1)>0 && jnr.at(0)!=prev) {	//NR loop: while the sign of the derivatives are the same then max is not straddled.  Second condition for j=1 repeatedly.
//    R_CheckUserInterrupt();
    prev = jnr.at(0);
    dderiv = deriv.at(1) - deriv.at( 0);				//approximate second derivs
    jnr.at(0) = fmax( floor( jnr.at(0) - deriv.at(0) / dderiv), 1.0);	//updating j to nearest lowest integer
    jnr.at(1) = jnr.at(0) + 1;						//updating j+1
    deriv.at(0) = logWderivApprox( jnr.at(0), z13, alpha3);		//update derivatives with new j
    deriv.at(1) = logWderivApprox( jnr.at(1), z13, alpha3);
  }
  logW.at(0) = logWfun( jnr.at(0), y3, muN3, alpha3, beta3);		//logW at j and j+1
  logW.at(1) = logWfun( jnr.at(1), y3, muN3, alpha3, beta3);

  if( logW.at(0) > logW.at(1)){		//find which of j or j+1 is bigger and return it and update logWmax
    logWmax3 = logW.at(0);
    return jnr.at(0);
  }
  else{
    logWmax3 = logW.at(1);
    return jnr.at(1);
  }
}

void findlogWjs( double y4, double muN4, double muZ4, double alpha4, double beta4, double z14, double z24, double jmax4, double &jlow4, double &jupp4, double logWmax4, vector<double> &logWjs4)
//finding the limits needed for the series evalulation and evaluating the series
{
  double eps = -37;	//Precision required to summing terms: corresponds to 10e-17

  jlow4 = jmax4;	//starting at top and working down
  jupp4 = jmax4 + 1;	//starting at top and working out

  logWjs4.clear();	//Just in case there is residual

  logWjs4.push_back( logWfun( jlow4, y4, muN4, alpha4, beta4));
  while( logWjs4.back() - logWmax4 > eps &&  jlow4 > 1){		//not yet within tolerance and not at lower limit
//    R_CheckUserInterrupt();						//User interrupting from R?
    jlow4--;
    logWjs4.push_back( logWfun( jlow4, y4, muN4, alpha4, beta4));	//calculate next Wj
  }

  logWjs4.push_back( logWfun( jupp4, y4, muN4, alpha4, beta4));
  while( logWjs4.back() - logWmax4 > eps){				//not yet within tolerance on upper side
//    R_CheckUserInterrupt();						//User interrrupt?
    jupp4++;
    logWjs4.push_back( logWfun( jupp4, y4, muN4, alpha4, beta4));	//Calculate next Wj
  }
}
