#include <TMB.hpp>   //Links in the TMB libraries
template<class Type>
Type objective_function<Type>::operator() (){
  using namespace density;
  using namespace Eigen;

  //Read data from R
  DATA_MATRIX(Y);       //Responses
  DATA_MATRIX(X);       //X is the archetype design matrix
  DATA_MATRIC(W);       //W is the species design matrix
  DATA_VECTOR(offy);    //offy is the offset indexed by sites (i)
  DATA_MATRIX(wts);     //wts is a matrix indexed by sites, species (i,j).
  DATA_INTEGER(nObs);    //n sites.
  DATA_INTEGER(nG);     //n groups
  DATA_INTEGER(nS);     //n species
  DATA_INTEGER(npx);    // n coefs X
  DATA_INTEGER(npw);    // n coefs W

  // for doing the GAMy bits once glm version is working.
  // DATA_SPARSE_MATRIX(S);//Penalization matrix diag(S1,S2,S3,S4,S5) without storing off-diagonal zeros.
  // DATA_IVECTOR(Sdims);   //Dimensions of S1,S2,S3,S4 and S5
  // DATA_SPARSE_MATRIX(designMatrixForReport);//Design matrix for report of splines


  //Parameters
  // PARAMETER_VECTOR(alpha); //intercepts. // Could potentiak merge this into gamma for ease.
  PARAMETER_MATRIX(beta);  //archetype coefs.
  PARAMETER_MATRIX(gamma); //species specific coefs //species intercepts are in here.
  PARAMETER_VECTOR(pi);    //mixing coefs.
  PARAMETER_VECTOR(theta); //dispersion coefs.

  // intialise the negative loglike.
  Type logl=0;

  matrix<Type> sppEta(nObs,nS);
  matrix<Type> grpEta(nObs,nG);
  vector<Type> sppLogl(nG);
  vector<Type> summand(nG);

  sppEta = X*beta.transpose();//might need to transpose beta.
  grpEta = W*gamma.transpose();//might need to transpose gamma.

  for(int ss=0; ss<nS; ss++){
  // for( size_t ss=0; ss<nS1; ss++){
    // avSummand = 0.0;
    for(int gg=0; gg<nG; gg++){
      sppLogl.at(gg) = 0.0;
      for( size_t ii=0; ii<nn1; ii++){
        mu = exp(sppEta(ii,ss) + grpEta(ii,gg) + offy(ii));
        // sppLogl.at(gg) += ldnegbin( y1[MATREF(ii,ss,nn1)], mu, theta(ss));
        // sppLogl(gg) += dnbinom(Y(ii,ss), prob=mu, size=theta(ss),give_log=1);
        sppLogl(gg) += dnbinom(Y(ii,ss), prob=mu, size=theta(ss),give_log=1);
      }
      summand(gg) = log(pi(gg)) + sppLogl(gg);
    }
    tmp = 0.0;
    tau = 100;	//magic number
    for(int gg=0; gg<nG1; gg++){
      summandAlt(gg) = exp(summand(gg) / tau);
      tmp += summandAlt(gg);
    }
    for(int gg=0; gg<nG1; gg++)
      wt(gg) = summandAlt(gg) / tmp;
    avSummand = 0.0;
    for( int gg=0; gg<nG1; gg++)
      avSummand += wt.at(gg) * summand.at(gg);
    tmp = 0.0;
    for( int gg=0; gg<nG1; gg++)
      tmp += exp( summand(gg) - avSummand);
    sppContr = avSummand + log( tmp);
    logl += sppContr;
    penalty += dbeta( phi(ss), thetaRange1, penParm1, true);
  }
  return (logl + penalty);
}
