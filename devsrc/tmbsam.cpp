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
  DATA_MATRIX(nObs);    //n sites.
  DATA_INTEGER(nG);     //n groups
  DATA_INTEGER(nS);     //n species
  DATA_INTEGER(npx);    // n coefs X
  DATA_INTEGER(npw);    // n coefs W

  //Parameters
  PARAMETER_VECTOR(alpha); //intercepts. // Could potentiak merge this into gamma for ease.
  PARAMETER_MATRIX(beta);  //archetype coefs.
  PARAMETER_MATRIX(gamma); //species specific coefs.
  PARAMETER_VECTOR(pi);    //mixing coefs.
  PARAMETER_VECTOR(theta); //dispersion coefs.

  // intialise the negative loglike.
  Type nll=0;

  matrix<Type> sppEta(nObs,nS);
  matrix<Type> grpEta(nObs,nG);
  vector<Type> sppLogl(nG);
  vector<Type> summand(nG);

  sppEta = X*beta;//might need to transpose beta.
  grpEta = W*gamma;//might need to transpose gamma.

  for(int ss=0; ss<nS; ss++){
  // for( size_t ss=0; ss<nS1; ss++){
    avSummand = 0.0;
    for(int gg=0; gg<nG; gg++){
      sppLogl.at(gg) = 0.0;
      for( size_t ii=0; ii<nn1; ii++){
        mu = exp(alpha(ss) + sppEta(ii,ss) + grpEta(ii,gg) + offy(ii));
        // sppLogl.at(gg) += ldnegbin( y1[MATREF(ii,ss,nn1)], mu, theta(ss));
        sppLogl(gg) += dnbinom(Y(ii,ss), prob=mu, size=theta(ss),give_log=1);
      }
      summand(gg) = log(pi(gg)) + sppLogl(gg);
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
