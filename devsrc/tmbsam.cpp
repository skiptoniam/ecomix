#include <TMB.hpp>   //Links in the TMB libraries

template <class Type>
Type invMultLogit( Type pi2, Type alpha2, int nG2){
  Type tmp(nG2, 0.0);
  Type sum=0;

  for( int gg=0; gg<(nG2-1); gg++){
    tmp( gg) = exp( alpha2(gg));
    sum += tmp(gg);
  }
  tmp( nG2-1) = 1.0;
  sum += tmp( nG2-1);

  for( int gg=0; gg<nG2; gg++)
    pi2.push_back(tmp(gg) / sum);
}


enum valid_family {
  bernoulli = 1,
  poisson  = 2,
  ippm  = 3,
  negative_binomial=4,
  tweedie = 5,
  normal  = 6
};

enum valid_link {
  identity_link = 0,
  log_link      = 1,
  logit_link    = 2
};

template <class Type>
Type InverseLink(Type eta, int link)
{
  Type out;
  switch (link) {
  case identity_link:
    out = eta;
    break;
  case log_link:
    out = exp(eta);
    break;
  case logit_link:
    out = invlogit(eta);
    break;
  default:
    error("Link not implemented.");
  }
  return out;
}

template<class Type>
Type objective_function<Type>::operator() (){
  using namespace density;
  using namespace Eigen;

  //Read data from R
  DATA_MATRIX(Y);       //Responses
  DATA_MATRIX(X);       //X is the archetype design matrix
  DATA_MATRIX(W);       //W is the species design matrix
  DATA_VECTOR(offy);    //offy is the offset indexed by sites (i)
  DATA_MATRIX(wts);     //wts is a matrix indexed by sites, species (i,j).
  DATA_INTEGER(nObs);    //n sites.
  DATA_INTEGER(nG);     //n groups
  DATA_INTEGER(nS);     //n species
  // DATA_VECTOR(thetaRange);

  // what distributiona and link function to use.
  DATA_INTEGER(family);
  DATA_INTEGER(link);
  // DATA_SCALAR(penParm1);

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
  Type mu, eta_i, logl= 0.0, tmp=0.0, avSummand = 0.0, sppContr= 0.0,tau = 100;

  matrix<Type> sppEta(nObs,nS);
  matrix<Type> grpEta(nObs,nG);
  vector<Type> sppLogl(nG);
  vector<Type> summand(nG);
  vector<Type> summandAlt(nG);
  vector<Type> wt(nG);
  // double tmp, tau, avSummand, sppContr;

  sppEta = X*beta;//might need to transpose beta.
  grpEta = W*gamma;//might need to transpose gamma.

  for(int ss=0; ss<nS; ss++){
    for(int gg=0; gg<nG; gg++){
      sppLogl(gg) = 0.0;
      for( size_t ii=0; ii<nObs; ii++){
        eta_i = sppEta(ii,ss) + grpEta(ii,gg) + offy(ii);
        mu = InverseLink(eta_i, link);
        sppLogl(gg) -= dnbinom(Y(ii,ss), theta(ss), mu, 1);
      }
      summand(gg) = log(pi(gg)) + sppLogl(gg);
    }
    // double tmp = 0.0;
    // double tau = 100;	//magic number
    for(int gg=0; gg<nG; gg++){
      summandAlt(gg) = exp(summand(gg) / tau);
      tmp += summandAlt(gg);
    }
    for(int gg=0; gg<nG; gg++)
      wt(gg) = summandAlt(gg) / tmp;

    for( int gg=0; gg<nG; gg++)
      avSummand += wt(gg) * summand(gg);
    Type tmp = 0.0;
    for( int gg=0; gg<nG; gg++)
      tmp += exp( summand(gg) - avSummand);
    sppContr = avSummand + log( tmp);
    logl -= sppContr;

// # pen.max <- theta.range[2]
// # pen.min <- theta.range[1]
// # shape1 <- shape2 <- 1.25
// # if(disty==4)  sppLogls <- sppLogls + dbeta( (theta-pen.min) / (pen.max-pen.min), shape1, shape2, log=TRUE)
    // theta_tmp = 0.0;
    // Type theta_tmp = (theta(ss)-thetaRange(0)) / (thetaRange(1)-thetaRange(0));
    // penalty += dbeta(theta_tmp, 1.25, 1.25, true);
  }
  return (logl);//# + penalty);
}
