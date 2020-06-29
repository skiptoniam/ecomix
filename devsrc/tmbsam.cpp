#include <TMB.hpp>   //Links in the TMB libraries
#include "debug_print.hpp"

template<class Type>
vector<Type> invMultLogit( vector<Type> pi2, vector<Type> eta2, int nG2){
 
  Type sumTmp = 0.;

  vector<Type> expEta = exp(eta2);
  Type sumexpEta = sum(expEta);  
   
  for( int gg=0; gg<(nG2-1); gg++){
	  pi2(gg) = expEta(gg)/(1.0+sumexpEta);
	  sumTmp += pi2(gg);	  
  }
    
  pi2(nG2-1) = 1.0 - sumTmp;  
  return pi2;	
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
Type logit_inverse_linkfun(Type eta, int link) {
  Type ans;
  switch (link) {
  case logit_link:
    ans = eta;
    break;
  default:
    ans = logit( inverse_linkfun(eta, link) );
  } // End switch
  return ans;
}


//double log_ippm_sam(const double &y, const double &mu, const double &st_sp_wts){
	//double tmp, z;
	//z = y/st_sp_wts;
	//tmp = z * log(mu);
	//tmp -= mu;
	//tmp *= st_sp_wts;
	//return( tmp);
//}


template<class Type>
Type objective_function<Type>::operator() (){
  using namespace density;
  using namespace Eigen;

  //Read data from R
  DATA_MATRIX(Y);       //Responses
  DATA_MATRIX(y_is_na); //NA data in response
  DATA_MATRIX(X);       //X is the archetype design matrix
  DATA_MATRIX(W);       //W is the species design matrix
  DATA_VECTOR(size);    //Include a size in the tmb binomial model for nicole :)
  DATA_VECTOR(offy);    //offy is the offset indexed by sites (i)
  DATA_MATRIX(wts);     //wts is a matrix indexed by sites, species (i,j) this is for ippm.
  DATA_VECTOR(bb_wts);  //bb_wts is for doing bayesian bootstrap by species (j)
  DATA_INTEGER(nObs);   //n sites.
  DATA_INTEGER(nG);     //n groups
  DATA_INTEGER(nS);     //n species
  DATA_INTEGER(family); //What error distribution to fit.
  DATA_INTEGER(link);   //What link function to use.
  DATA_INTEGER(keep_mu);//logical 1 = return mus. 
  // DATA_VECTOR(thetaRange); penalties for overdispersion if needed.
  // DATA_SCALAR(penParm1);

  // for doing the GAMy bits once glm version is working.
  // DATA_SPARSE_MATRIX(S);//Penalization matrix diag(S1,S2,S3,S4,S5) without storing off-diagonal zeros.
  // DATA_IVECTOR(Sdims);   //Dimensions of S1,S2,S3,S4 and S5
  // DATA_SPARSE_MATRIX(designMatrixForReport);//Design matrix for report of splines

  //Parameters
  // PARAMETER_VECTOR(alpha); //intercepts. // Could potentiak merge this into gamma for ease.
  PARAMETER_MATRIX(beta);  //archetype coefs.
  PARAMETER_MATRIX(gamma); //species specific coefs //species intercepts are in here.
  PARAMETER_VECTOR(eta);    //mixing coefs.
  PARAMETER_VECTOR(theta); //dispersion coefs.

  // intialise the negative loglike.
  Type mu_i = 0., eta_i = 0., nll= 0.;

  array<Type> mus(nObs,nS,nG); //Array for storing mus
  matrix<Type> sppEta(nObs,nS); //Matrix of spp linear predictors
  matrix<Type> grpEta(nObs,nG); //Matrix of group linear predictors
  matrix<Type> loglGS(nG,nS); //loglike speceis grousps.
  vector<Type> pi2(nG);
  vector<Type> loglS(nS); 

  // additive transfrom.  
  vector<Type> pi = invMultLogit(pi2, eta, nG);

  
  //std::cout<<" exp(eta) "<< expEta <<"\n";//returns only one number
  //std::cout<<" sumEta "<< sumexpEta <<"\n";//returns only one number
  std::cout<<" pi "<< pi <<"\n";//returns only one number
  
  grpEta = X*beta; // Mixing coefs
  sppEta = W*gamma; // Species coefs
  std::cout<<"print grpEta"<< grpEta.head() <<"\n";
  std::cout<<"print sppEta"<< sppEta.head() <<"\n";
  // get the mus
  
  Type s1, s2 tmp_loglik=0.;
  
  for(int ss=0; ss<nS; ss++){
	     for(int gg=0; gg<nG; gg++){
			       for( int ii=0; ii<nObs; ii++){
					  if(y_is_na(ii,ss)>0){
					  eta_i = sppEta(ii,ss) + grpEta(ii,gg) + offy(ii);
                      mu_i = InverseLink(eta_i, link);
                      //std::cout<<"print mu"<< mu_i<<"\n";
                      if(keep_mu) mus(ii,ss,gg) = mu_i;
                      switch (family) {
					  case normal:
						loglGS(gg,ss) += dnorm(Y(ii,ss), mu_i, sqrt(theta(ss)), true);
						break;
					  case poisson:
						loglGS(gg,ss) += dpois(Y(ii,ss), mu_i, true);
						break;	
					  case ippm:
						loglGS(gg,ss) += dpois(Y(ii,ss), mu_i, true);
						break;
					  case bernoulli:
						s1 = logit_inverse_linkfun(eta_i, link); // logit(p)
						loglGS(gg,ss) += dbinom_robust(Y(ii,ss), size(ii), s1, true); //size(ii) if you want binomial with size.
						break;
					  case negative_binomial:
					    s1 = mu_i;
                        s2 = mu_i * (Type(1) + mu_i / theta(ss));
 					    loglGS(gg,ss) += dnbinom_robust(Y(ii,ss), s1, s2, true);
					    break;
					}
					////loglGS(gg,ss) *= wts(ii,ss);
				  }
			  }
            //loglGS(gg,ss) *= bb_wts(ss);// bayesian boostrap weights
	   }
    }

    // Sum up the species logls 
    for( int ss=0; ss<nS; ss++){
		
	Type eps=0.0, glogl=0.0, tloglike=0.0;

	////calcualte the G*S log conditional densities
    for( int gg=0; gg<nG; gg++){
		if(gg==0) eps = loglGS(gg,ss);
		if(loglGS(gg,ss) > eps) eps = loglGS(gg,ss);
	}

  // this will calculate the species specific loglikelihoods based on sum across Gs.
	for(int gg=0; gg<nG; gg++){
		glogl += pi(gg)*exp(loglGS(gg,ss) - eps);
	  }
	    
	 loglS(ss) = log(glogl) + eps;
	  nll -= tloglike;
	}
	
	

  return (nll);//# + penalty);
}
