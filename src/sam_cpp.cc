#include"sam_cpp.h"

/* Code for all SAMs distributions (except Tweedie).
 * I have tried to set this up like RCP, which makes more sense to me, in the future I can adapt Piers code into a single species_mix_cpp function.
 */

// this is the external C call which will be called by R using .Call.

extern "C" {
	SEXP species_mix_cpp(SEXP Ry, SEXP RX, SEXP Roffset, SEXP Rspp_wts, SEXP Rsite_spp_wts, SEXP Ry_not_na,
							  SEXP RnS, SEXP RnG, SEXP Rp, SEXP RnObs, SEXP Rdisty,
							  SEXP Ralpha, SEXP Rbeta, SEXP Reta, SEXP Rdisp, 
							  SEXP RderivsAlpha, SEXP RderivsBeta, SEXP RderivsEta, SEXP RderivsDisp, SEXP RgetScores, SEXP Rscores,
							  SEXP Rpis, SEXP Rmus, SEXP RlogliS, SEXP RlogliSG,
							  SEXP Rmaxit, SEXP Rtrace, SEXP RnReport, SEXP Rabstol, SEXP Rreltol, SEXP Rconv, SEXP Rprintparams,
							  SEXP Roptimise, SEXP RloglOnly, SEXP RderivsOnly, SEXP RoptiDisp){

	sam_cpp_all_classes all;

	//initialise the data structures -- they are mostly just pointers to REAL()s...
	all.data.setVals(Ry, RX, Roffset, Rspp_wts, Rsite_spp_wts, Ry_not_na, RnS, RnG, Rp, RnObs, Rdisty, RoptiDisp);	//read in the data
	all.params.setVals(all.data, Ralpha, Rbeta, Reta, Rdisp);	//read in the parameters
	all.derivs.setVals(all.data, RderivsAlpha, RderivsBeta, RderivsEta, RderivsDisp, RgetScores, Rscores);
	all.contr.setVals( Rmaxit, Rtrace, RnReport, Rabstol, Rreltol, Rconv, Rprintparams);
	all.fits.initialise(all.data.nObs, all.data.nG, all.data.nS, all.data.nP, 0);

	double logl = -999999;

	//doing the optimising
	if( *INTEGER(Roptimise) == 1)
		logl = SAM_optimise(all);
	//re-running to get pis and mus
	if( *INTEGER(RloglOnly) == 1)
	   	logl = sam_cpp_mix_loglike( all.data, all.params, all.fits);
	//and derivatives (inlcuding scores, for empirical info, if requested)
	if( *INTEGER(RderivsOnly) == 1)	    
	    sam_cpp_mix_gradient( all.data, all.params, all.derivs, all.fits);
	
	//bundling up things to return - will need to change these...
	//first the fitted pis
	//double *tmpPi = REAL( Reta);
	//additive_logistic_sam(all.fits.par_pis,1,all.data.nG);
	//for( int g=0; g<all.data.nG; g++)
			//tmpPi[g] = all.fits.par_pis.at(g);
	//the fitted expectations
	double *tmpMus = REAL( Rmus);
	for( size_t i=0; i<all.fits.allMus.size(); i++)
		tmpMus[i] = all.fits.allMus.at(i);
	//the log conditional densities
	double *tmplogliSG = REAL( RlogliSG);
	for( int g=0; g<all.data.nG;g++)
		for( int s=0; s<all.data.nS; s++)
				tmplogliSG[MATREF2D(s,g,all.data.nS)] = all.fits.log_like_species_group_contrib[MATREF2D(s,g,all.data.nS)];
	//the logl contributions
	double *tmplogliS = REAL( RlogliS);
	for( int s=0; s<all.data.nS; s++)
		tmplogliS[s] = all.fits.log_like_species_contrib.at(s);
	//Convergence code
	int *tmpconv = INTEGER( Rconv);
	tmpconv[0] = all.contr.ifail;
	//the logl
	//SEXP Rres;	//R object to return -- it is the logl!
	//Rres = PROTECT( allocVector(REALSXP,1));
    //REAL( Rres)[0] = logl;
	//UNPROTECT(1);
	//return( Rres);
	
	SEXP Rlogl = PROTECT(allocVector(REALSXP, 1));
    REAL(Rlogl)[0] = logl;
   	UNPROTECT(1);
    SEXP Ralpha_est = PROTECT(allocVector(REALSXP, all.data.nS));
    for( int s=0; s<all.data.nS; s++) REAL(Ralpha_est)[s] = all.params.Alpha[s];
	UNPROTECT(1);
	SEXP Rbeta_est = PROTECT(allocVector(REALSXP, all.data.nG*all.data.nP));
	for( int i=0; i<((all.data.nG*all.data.nP)); i++) REAL(Rbeta_est)[i] = all.params.Beta[i];
	UNPROTECT(1);
	SEXP Reta_est =PROTECT(allocVector(REALSXP, all.data.nG-1));
	for( int g=0; g<(all.data.nG-1);g++) REAL(Reta_est)[g] = all.params.Eta[g];
	UNPROTECT(1);
	SEXP Rdisp_est =PROTECT(allocVector(REALSXP, all.data.nS));
	for( int s=0; s<(all.data.nS);s++) REAL(Rdisp_est)[s] = all.params.Disp[s];
	UNPROTECT(1);
	
	
	const char *names[] = {"logl", "alpha", "beta", "eta", "disp",""};                   /* note the null string */
	SEXP Rres = PROTECT(mkNamed(VECSXP, names));  /* list of length 3 */
	SET_VECTOR_ELT(Rres, 0, Rlogl);        // loglike 
	SET_VECTOR_ELT(Rres, 1, Ralpha_est);   // species intercepts 
	SET_VECTOR_ELT(Rres, 2, Rbeta_est);    // mixing coefs
	SET_VECTOR_ELT(Rres, 3, Reta_est);     // etas - transformed pis.
	SET_VECTOR_ELT(Rres, 4, Rdisp_est);    // dispersion parameters
	UNPROTECT(1);
	return (Rres);

  }
}


double SAM_optimise( sam_cpp_all_classes &all){
	double *vmminGrad, *vmminParams, *oldParms;	//arrays to pass to vmmin
	double *vmminParamsIn;	//del
	vmminParams = (double *) R_alloc(all.params.nTot,sizeof(double));
	vmminParamsIn = (double *) R_alloc(all.params.nTot,sizeof(double));	//Del
	oldParms = (double *) R_alloc(all.params.nTot,sizeof(double));
	vmminGrad = (double *) R_alloc(all.params.nTot,sizeof(double));
	int *myMask;
	vector<int> vecMask(all.params.nTot, 1);
	double vmminLogl[1];

//	Rprintf( "Quasi-Newton iterations\n");
	all.params.getArray( oldParms, all.data);
	myMask = &vecMask[0];
	//optimise
	all.params.getArray( vmminParams, all.data);
	all.params.getArray( vmminParamsIn, all.data);	//del
	vmmin( all.params.nTot, vmminParams, vmminLogl, optimise_function_sam, gradient_function_sam, all.contr.maxitQN, 
	 all.contr.traceQN, myMask,  all.contr.abstol, all.contr.reltol,  all.contr.nReport, &all, &all.contr.fnKount,
	 &all.contr.grKount, &all.contr.ifail);

//	nmmin( all.params.nTot, vmminParamsIn, vmminParams, vmminLogl, optimise_function_rcp, &all.contr.ifail, all.contr.abstol, all.contr.reltol, &all, 1.0, 0.5, 2.0, all.contr.traceQN, &all.contr.fnKount, all.contr.maxitQN);

	//update parameters
	all.params.update( vmminParams, all.data);
	gradient_function_sam(all.params.nTot, vmminParams, vmminGrad, &all);
	all.derivs.update( vmminGrad, all.data);
    if(all.contr.printparams==1)all.params.printParms(all.data);

	return(vmminLogl[0]);
}

double optimise_function_sam(int n, double *par, void *ex){
	
	sam_cpp_all_classes *all = (sam_cpp_all_classes *) ex;
	double logl;

	all->params.update( par, all->data);
	logl = sam_cpp_mix_loglike( all->data, all->params, all->fits);
     
    //Rprintf("%f\n",logl);
    return((0.0-logl));
}

double sam_cpp_mix_loglike(const sam_data &dat, const sam_params &params, sam_fits &fits){

	double tloglike = 0.0, loglike = 0.0;
	vector<double> par_pi(dat.nG-1,0);
	
	fits.zero(0);
	for(int g=0; g<(dat.nG-1); g++) par_pi.at(g) = params.Eta[g];

	//transform additative pis to natural scale - need this to calc loglikes.	
	additive_logistic_sam(par_pi,1,dat.nG); // additive logistic transformation of pis.

	//calculate fitted values (constant over i)
	calc_mu_fits(fits.allMus, params, dat);

	//calculate the species/groups loglikes
	calc_sam_loglike_SG(fits.log_like_species_group_contrib, fits.allMus, dat, params);
	
	// calc loglike per species
	for( int s=0; s<dat.nS; s++){
		tloglike = calc_sam_loglike_S(fits.log_like_species_group_contrib, par_pi, dat, s); // this should give the mix loglike. ippm weights are calculated in this bit.
		fits.log_like_species_contrib.at(s) = tloglike;
		loglike += tloglike;
	}
	return(loglike);
}

// calculate mu fits for ippm this should calculate all the etas and mus for species and archetypes.
void calc_mu_fits(vector<double> &fits, const sam_params &params, const sam_data &dat){

	vector<double> lps(dat.nG*dat.nS, 0);	//the nG x nS intercepts
	double lp=0.0;	//the lin pred for the gth group, sth species and ith site

	//calcualte the G*S*n fits
	for( int g=0; g<dat.nG; g++){
		for( int s=0; s<dat.nS; s++){
			lps.at(MATREF2D(g,s,dat.nG)) = params.Alpha[s]; 
			for( int i=0; i<dat.nObs; i++){
				//std::cout << g << ' '<< s << ' '<<  i <<'\n'<<'\n';
				//Rprintf( " %i\n", dat.y_not_na[MATREF2D(i,s,dat.nObs)]); 
				// need logical flag which deals with NA data.
				if(dat.y_not_na[MATREF2D(i,s,dat.nObs)]>0){
				lp = lps.at(MATREF2D(g,s,dat.nG)) + dat.offset[i];
					for( int j=0;j<dat.nP; j++){
							lp += params.Beta[MATREF2D(g,j,(dat.nG))] * dat.X[MATREF2D(i,j,dat.nObs)];
						}
					if(dat.disty==1){//bernoulli
							fits.at( MATREF3D(i,s,g,dat.nObs,dat.nS)) = inverse_logit(lp);
						}							
					if(dat.disty==2 | dat.disty==3 | dat.disty==4 | dat.disty==5){//poisson, ippm, negative binomial, tweedie
							fits.at( MATREF3D(i,s,g,dat.nObs,dat.nS)) = exp(lp);
						}
					if(dat.disty==6){//normal
							fits.at( MATREF3D(i,s,g,dat.nObs,dat.nS)) = lp;
					}		
				}
			}
		}
	}
}



void calc_sam_loglike_SG(vector<double> &loglSG, vector<double> &fits, const sam_data &dat, const sam_params &params){
	
    //calcualte the G*S log conditional densities
	for(int s=0; s<dat.nS; s++){
		for( int g=0; g<dat.nG; g++){
			for(int i=0; i<dat.nObs; i++){
				if(dat.y_not_na[MATREF2D(i,s,dat.nObs)]>0){
					if(dat.disty==1){ // bernoulli
						loglSG.at(MATREF2D(g,s,dat.nG)) += log_bernoulli_sam(dat.y[MATREF2D(i,s,dat.nObs)], fits.at(MATREF3D(i,s,g,dat.nObs,dat.nS)));
						}
					if(dat.disty==2){ // poisson
						loglSG.at(MATREF2D(g,s,dat.nG)) += log_poisson_sam(dat.y[MATREF2D(i,s,dat.nObs)], fits.at(MATREF3D(i,s,g,dat.nObs,dat.nS)));
						}
					if(dat.disty==3){ // ippm
						loglSG.at(MATREF2D(g,s,dat.nG)) += log_ippm_sam(dat.y[MATREF2D(i,s,dat.nObs)], fits.at(MATREF3D(i,s,g,dat.nObs,dat.nS)), dat.site_spp_wts[MATREF2D(i,s,dat.nObs)]);	
						}
					if(dat.disty==4){ // negative binomial
						loglSG.at(MATREF2D(g,s,dat.nG)) += log_negative_binomial_sam(dat.y[MATREF2D(i,s,dat.nObs)], fits.at(MATREF3D(i,s,g,dat.nObs,dat.nS)), params.Disp[s]);	
						}
					//if(dat.disty==5){ // tweedie
						
					//}
					if(dat.disty==6){ // normal
						loglSG.at(MATREF2D(g,s,dat.nG)) += log_normal_sam(dat.y[MATREF2D(i,s,dat.nObs)], fits.at(MATREF3D(i,s,g,dat.nObs,dat.nS)), params.Disp[s]);	
						}
						
					}	
				}	
			// This is for the bayesian bootstrap. For ippm this will be set to default of 1 in R as BB is to hard to do with ippms.
			loglSG.at(MATREF2D(g,s,dat.nG)) = loglSG.at(MATREF2D(g,s,dat.nG))*dat.spp_wts[s]; //fix this up.
			}
		}
	}


// now we want to calculate the species specific likelihoods.
double calc_sam_loglike_S(vector<double> &fits, vector<double> const &pis, const sam_data &dat, int s){

	double eps=0.0, glogl=0.0;

	////calcualte the G*S log conditional densities
    for( int g=0; g<dat.nG; g++){
		if(g==0) eps = fits.at(MATREF2D(g,s,dat.nG));
		if(fits.at(MATREF2D(g,s,dat.nG)) > eps) eps = fits.at(MATREF2D(g,s,dat.nG));
  }

  // this will calculate the species specific loglikelihoods based on sum across Gs.
  for(int g=0; g<dat.nG; g++){
    glogl += pis.at(g)*exp(fits.at(MATREF2D(g,s,dat.nG)) - eps);
  }
  glogl = log(glogl) + eps;

  return(glogl);

}


double log_bernoulli_sam( const double &y, const double &mu){
	double tmp;
	if( y==1){
		tmp = log( mu);
		return( tmp);
	}
	tmp = log( 1-mu);
	return( tmp);
}

double log_bernoulli_deriv_sam(const double &y, const double &mu){
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

double log_poisson_sam( const double &y, const double &mu){
	double tmp;
	tmp = y * log( mu);
	tmp -= lgammafn(y+1);
	tmp -= mu;
	return( tmp);
}

double log_poisson_deriv_sam( const double &y, const double &mu){
	double tmp;
	tmp = y/mu;
	tmp -= 1;
	return( tmp);
}

// this should give the log-density for the ippm. fingers crossed...
double log_ippm_sam(const double &y, const double &mu, const double &st_sp_wts){
	double tmp, z;
	z = y/st_sp_wts;
	tmp = z * log(mu);
	tmp -= mu;
	tmp *= st_sp_wts;
	return( tmp);
}

double log_ippm_deriv_sam( const double &y, const double &mu, const double &st_sp_wts){
	double tmp = 0.0, z;
	z = y/st_sp_wts;
	tmp = z/mu;
	tmp -= 1;
	return( tmp);
}

double log_negative_binomial_sam( const double &y, const double &mu, const double &od){
	double tmp = 0.0, theta;
	theta = 1/exp( od);
	tmp = dnbinom_mu(y, theta, mu, 1);
	return( tmp);
}

double log_negative_binomial_deriv_disp_sam(const double &y, const double &mu, const double &od){
	
	double theta, tmp = 0.0;
	
    //sig = exp(od);
    theta = 1 / exp(od); 
    
    tmp = digamma( theta+y);
    tmp -= digamma( theta);
    tmp += 1 + log(theta) - log(theta+mu);
    tmp -= (theta+y)/(theta+mu);
    
    return(tmp);
    }
    
double log_negative_binomial_deriv_mu_sam( const double &y, const double &mu, const double &od){
	double tmp = 0.0, theta;
	theta = 1/exp( od);
	tmp = -(theta+y)/(theta+mu);
	tmp += y/mu;
	return( tmp);
}    

//double log_tweedie( const double &y, const double &mu, const double &phi, const double &p) {
 	//double lambda, alpha, tau, muZ, tmp, phi1;
 	//phi1 = exp( phi);
 	//lambda = R_pow( mu, (2-p)) / ( phi1*(2-p));
 	//alpha = ( 2-p) / ( p-1);
 	//tau = phi1*(p-1)*R_pow(mu,(p-1));
 	//muZ = alpha * tau;
 	//tmp = dTweedie( y, lambda, muZ, alpha, 1);
 	//return( tmp);
//}

//double log_tweedie_deriv_disp_sam( double y, double fit, double dispParm , double p){
 	//double phi, tmp;
 	//phi = exp( dispParm);
	//tmp = dTweediePhi( y, fit, phi, p);
	//tmp *= phi;
 	//return( tmp);
//}
    

double log_normal_sam( const double &y, const double &mu, const double &sig){
	double tmp = 0.0, sig1;
	sig1 = exp( sig);
	tmp = -log( sig1);
	tmp -= (y-mu)*(y-mu)/(2*sig1*sig1);
	return( tmp);
}    

double log_normal_deriv_mu_sam( const double &y, const double &mu, const double &sig){
	double tmp, sig1;

	sig1 = exp( sig);
	tmp = (y-mu) / (sig1*sig1);
	return( tmp);
}

double log_normal_deriv_disp_sam(const double &y, const double &mu, const double &sig){
	double sig1, res=0.0;

	sig1 = exp(sig);
	res = (y-mu)*(y-mu) / (sig1*sig1*sig1);
	res -= 1/sig1;
	res *= sig1; //for the change of variable dispParm --> sig

	return( res);
}

// this should do the additive transformation of pis (pi1/piN,pi2/piN,pi(N-1)/piN)
void additive_logistic_sam(vector< double > &x, int inv, int G){
  int i;
  // inv == 1 gives transfornmation
  // inv == 0 gives inverse transformation
  if(inv==0){
    for(i=0;i<G;i++) x.at(i) = log(x.at(i)/x.at(G-1));
    return;
  }

  vector< double > xt (x.size(),0);
  double sumx=0, sumxt=0;

  for(i=0;i<x.size();i++){
    xt.at(i) = exp(x.at(i));
    sumx+=xt.at(i);
  }
  for(i=0;i<x.size();i++){
    x.at(i) = xt.at(i)/(1+sumx);
    sumxt+=x.at(i);
  }
  x.push_back(1-sumxt);

}

// Gradient functions for IPPM.
// These are all the functions for the gradient function and hopefully they should run and help estimate the derivates.
// This is the wrapper for vmmin and it takes the sam_cpp_all_classes to
void gradient_function_sam(int n, double *par, double *gr, void *ex){
	sam_cpp_all_classes *all = (sam_cpp_all_classes *) ex;

    //Rprintf(all.derives) check rcp mod and ex.
	sam_cpp_mix_gradient( all->data, all->params, all->derivs, all->fits);

	all->derivs.getArray(gr, all->data);

	for( int i=0; i<n; i++){
	    gr[i] = 0.0 - gr[i];
	}
}


//// this function should work out the derivatives.
void sam_cpp_mix_gradient(const sam_data &dat, const sam_params &params, sam_derivs &derivs, sam_fits &fits){

	vector<double> parpi((dat.nG-1), 0);
	vector<double> eta_mu_derivs((dat.nS*dat.nG*dat.nObs), 0);
	vector<double> alphaDerivs(dat.nS, 0);//change to dat.NAN
	vector<double> betaDerivs((dat.nG*dat.nP), 0);
	vector<double> etaDerivs((dat.nG-1), 0); // check there should only be g pis
	vector<double> dispDerivs(dat.nS, 0);
	double logl;

    //calc loglike
    fits.zero(0);
    derivs.zeroDerivs(dat);
    
	logl = sam_cpp_mix_loglike(dat, params, fits);
	
	for(int g=0; g<(dat.nG-1); g++) parpi.at(g) = params.Eta[g];
	additive_logistic_sam(parpi,1,dat.nG);
	
	// derivate w.r.t the mean and link
	calc_mu_deriv(fits.all_derivs_mu, fits.allMus, dat, params);
	calc_eta_mu_deriv(eta_mu_derivs, dat, fits.all_derivs_mu, fits.allMus);
		
	//derivate w.r.t alpha
   	calc_dlog_dalpha(fits.dlogdalpha, eta_mu_derivs, dat);
	calc_alpha_deriv(alphaDerivs, fits.dlogdalpha, fits.log_like_species_group_contrib, fits.log_like_species_contrib, parpi, dat);

	//derivate w.r.t beta
	calc_dlog_dbeta(fits.dlogdbeta, eta_mu_derivs, dat);
	calc_beta_deriv(betaDerivs, fits.dlogdbeta, fits.log_like_species_group_contrib, fits.log_like_species_contrib, parpi, dat);
	
	//derivate w.r.t dispersions
	if( dat.isDispersion()){ // if no disperion?  move along please 
		calc_dlog_ddispersionS(fits.dlogddispersion, fits.allMus, dat, params);
		calc_dispersion_deriv(dispDerivs, fits.dlogddispersion, fits.log_like_species_group_contrib, fits.log_like_species_contrib, parpi, dat);	
    }
        
	//transform pis back to additative logistic scale to keep pi_dervis happy.
	additive_logistic_sam(parpi,0,dat.nG);
	
	//derivate w.r.t pi/eta
	calc_dlog_dpi(fits.dlogdpi, fits.log_like_species_group_contrib, fits.log_like_species_contrib, dat);
	calc_eta_deriv(etaDerivs, fits.dlogdpi, parpi, dat);

	//update the derivates.
	derivs.updateDerivs( dat, alphaDerivs, betaDerivs, etaDerivs, dispDerivs);	
	}

/* Ok I'm going to try and generalise the derivate function across all distributions */
/* firstly we are going to estimate DerivMu which will be used across all the derivates. 
 * This will replace the tmp_lpd or whatever I've called it */
void calc_mu_deriv( vector<double> &mu_derivs, const vector<double> &fits, const sam_data &dat, const sam_params &params){
	//derivatives of conditional density w.r.t. its mean
	//muDerivs is a GxS matrix of first derivatives

	for( int g=0; g<dat.nG; g++){
		for( int s=0; s<dat.nS; s++){
			for(int i=0; i<dat.nObs; i++){
			//switch( dat.disty){
				if(dat.disty==1){
					mu_derivs.at(MATREF3D(i,s,g,dat.nObs,dat.nS)) = log_bernoulli_deriv_sam( dat.y[MATREF2D(i,s,dat.nObs)], fits.at(MATREF3D(i,s,g,dat.nObs,dat.nS)));
				}
				if(dat.disty==2){
					mu_derivs.at(MATREF3D(i,s,g,dat.nObs,dat.nS)) = log_poisson_deriv_sam( dat.y[MATREF2D(i,s,dat.nObs)], fits.at(MATREF3D(i,s,g,dat.nObs,dat.nS)));
				}
				if(dat.disty==3){
					mu_derivs.at(MATREF3D(i,s,g,dat.nObs,dat.nS)) = log_ippm_deriv_sam(dat.y[MATREF2D(i,s,dat.nObs)], fits.at(MATREF3D(i,s,g,dat.nObs,dat.nS)), dat.site_spp_wts[MATREF2D(i,s,dat.nObs)]);
				}
				if(dat.disty==4){
					mu_derivs.at(MATREF3D(i,s,g,dat.nObs,dat.nS)) = log_negative_binomial_deriv_mu_sam( dat.y[MATREF2D(i,s,dat.nObs)], fits.at(MATREF3D(i,s,g,dat.nObs,dat.nS)), params.Disp[s]);
				}
				//if(dat.disty==5){
				 	//mu_derivs.at(MATREF3D(i,s,g,dat.nObs,dat.nS)) = log_tweedie_deriv_sam( dat.y[MATREF2D(i,s,dat.nObs)], fits.at(MATREF3D(i,s,g,dat.nObs,dat.nS)), exp( params.Disp[s]), params.Power[s]);
				//}
				if(dat.disty==6){
					mu_derivs.at(MATREF3D(i,s,g,dat.nObs,dat.nS)) = log_normal_deriv_mu_sam(dat.y[MATREF2D(i,s,dat.nObs)], fits.at(MATREF3D(i,s,g,dat.nObs,dat.nS)), params.Disp[s]);
				}
				//}	
			}
		}
	}
}

void calc_eta_mu_deriv( vector<double> &etaDerivs, const sam_data &dat, const vector<double> &muDerivs, const vector<double> &fits){
	//derivatives of conditional densities w.r.t. its linear predictor(s)
	for( int g=0; g<dat.nG; g++){
		for( int s=0; s<dat.nS; s++){
			for(int i=0; i<dat.nObs; i++){
				if(dat.y_not_na[MATREF2D(i,s,dat.nObs)]>0){
						if(dat.disty==1){ //bernoulli
							etaDerivs.at(MATREF3D(i,s,g,dat.nObs,dat.nS)) = fits.at(MATREF3D(i,s,g,dat.nObs,dat.nS)) * (1-fits.at(MATREF3D(i,s,g,dat.nObs,dat.nS))) * muDerivs.at(MATREF3D(i,s,g,dat.nObs,dat.nS));	//logit link
						}
						if(dat.disty==3){ // ippm
							etaDerivs.at(MATREF3D(i,s,g,dat.nObs,dat.nS)) = dat.site_spp_wts[MATREF2D(i,s,dat.nObs)] * fits.at(MATREF3D(i,s,g,dat.nObs,dat.nS)) * muDerivs.at(MATREF3D(i,s,g,dat.nObs,dat.nS)); // loglink + weights	
						}
						if(dat.disty==2 | dat.disty==4 | dat.disty==4){ // poisson, negative binomial, tweedie
							etaDerivs.at(MATREF3D(i,s,g,dat.nObs,dat.nS)) = fits.at(MATREF3D(i,s,g,dat.nObs,dat.nS)) * muDerivs.at(MATREF3D(i,s,g,dat.nObs,dat.nS));	//log link
						}
						if(dat.disty==6){ // normal 
							etaDerivs.at(MATREF3D(i,s,g,dat.nObs,dat.nS)) = muDerivs.at(MATREF3D(i,s,g,dat.nObs,dat.nS));	//identity link
						}	
				}
			}
		}
	}
}



void calc_dlog_dalpha(vector<double> &dlda, vector<double> const &mu_eta_derivs, const sam_data &dat){

	// dlda = dlogalpha passed as fits.dflogdalpha(dat.nG*dat.nS, dat.NAnum) from function call
	// mus = all the fitted values.
	//double tmp_lpd;

	for(int g=0; g<dat.nG; g++){
		for(int s=0;s<dat.nS; s++){
			for(int i=0; i<dat.nObs; i++){
			  if(dat.y_not_na[MATREF2D(i,s,dat.nObs)]>0){
					// this is the tmp log poisson derivative (lpd)
					//tmp_lpd = log_ippm_deriv_sam(dat.y[MATREF2D(i,s,dat.nObs)], mus.at(MATREF3D(i,s,g,dat.nObs,dat.nS)), dat.st_sp_wts[MATREF2D(i,s,dat.nObs)]);
					dlda.at(MATREF2D(g,s,dat.nG)) += mu_eta_derivs.at(MATREF3D(i,s,g,dat.nObs,dat.nS));//this is for the BB
				}
			}
			dlda.at(MATREF2D(g,s,dat.nG)) = dlda.at(MATREF2D(g,s,dat.nG))*dat.spp_wts[s]; //this is for the BB weights
		}
	}
}

void calc_dlog_dbeta(vector<double> &dldb, vector<double> const &mu_eta_derivs, const sam_data &dat){

	// dlda = dlogbeta passed as fits.dlogdbeta(dat.nG*dat.nS*dat.nP, dat.NAnum) from function call
	// mus = all the fitted values.
	
	//double tmp_lpd;

	for(int g=0; g<dat.nG; g++){
		for(int s=0;s<dat.nS; s++){
			for(int i=0; i<dat.nObs; i++){
				if(dat.y_not_na[MATREF2D(i,s,dat.nObs)]>0){
					//tmp_lpd = log_ippm_deriv_sam(dat.y[MATREF2D(i,s,dat.nObs)], mus.at(MATREF3D(i,s,g,dat.nObs,dat.nS)), dat.st_sp_wts[MATREF2D(i,s,dat.nObs)]);
						for(int j=0; j<dat.nP; j++){
							dldb.at(MATREF3D(g,j,s,dat.nG,dat.nP)) += mu_eta_derivs.at(MATREF3D(i,s,g,dat.nObs,dat.nS)) * dat.X[MATREF2D(i,j,dat.nObs)];
							//std::cout << dldb.at(MATREF3D(g,j,s,dat.nG,dat.nP)) << '\n';
							}
				}
			}
		for(int j=0; j<dat.nP; j++){
				dldb.at(MATREF3D(g,j,s,dat.nG,dat.nP)) = dldb.at(MATREF3D(g,j,s,dat.nG,dat.nP))*dat.spp_wts[s];
				//std::cout << dldb.at(MATREF3D(g,j,s,dat.nG,dat.nP)) << '\n';	
			}
		}
	}

}

void calc_dlog_ddispersionS(vector<double> &dldd, vector<double> const &mus, const sam_data &dat, const sam_params &params){

	// dlda = dlogalpha passed as fits.dflogdalpha(dat.nG*dat.nS, dat.NAnum) from function call
	// mus = all the fitted values.

	for(int g=0; g<dat.nG; g++){
		for(int s=0;s<dat.nS; s++){
			for(int i=0; i<dat.nObs; i++){
				//switch(dat.disty){
					if(dat.disty==4){ // negative binomial
						dldd.at(MATREF2D(g,s,dat.nG)) += log_negative_binomial_deriv_disp_sam(dat.y[MATREF2D(i,s,dat.nObs)], mus.at( MATREF3D(i,s,g,dat.nObs, dat.nS)), params.Disp[s]);
					}
					//case 5:	// tweedie
						//dldd.at(MATREF2D(g,s,dat.nG)) += log_tweedie_deriv_disp_sam(dat.y[MATREF2D(i,s,dat.nObs)], mus.at( MATREF3D(i,s,g,dat.nObs, dat.nS)), params.Disp[s]);
						//break;	
					if(dat.disty==6){ // normal
						dldd.at(MATREF2D(g,s,dat.nG)) += log_normal_deriv_disp_sam(dat.y[MATREF2D(i,s,dat.nObs)], mus.at( MATREF3D(i,s,g,dat.nObs, dat.nS)), params.Disp[s]);
					}
						//}
			}
		dldd.at(MATREF2D(g,s,dat.nG)) = dldd.at(MATREF2D(g,s,dat.nG))*dat.spp_wts[s];
		}
	}
}


void calc_dlog_dpi(vector<double> &dldpi, vector<double> const &llSG, vector<double> const &llS, const sam_data &dat){

	for(int g=0; g<(dat.nG); g++){
			//Rprintf( " %f\n", fits.dlogdpi.at(g)); 
			for(int s=0; s<(dat.nS); s++){
							//Rprintf( "logl_sg eta %f\n", fits.log_like_species_group_contrib.at(MATREF2D(g,s,dat.nG))); 
							//Rprintf( "logl_s eta %f\n", fits.log_like_species_contrib.at(s));
							dldpi.at(g) += exp(llSG.at(MATREF2D(g,s,dat.nG)) - llS.at(s));
			}
			//Rprintf( " %f\n", fits.dlogdpi.at(g)); 
	}
}

//// this should calculate the derivate w.r.t alpha.
void calc_alpha_deriv( vector<double> &alphaDerivs, vector<double> const &dlogdalpha, vector<double> const &llSG, vector<double> const &llS, vector<double> const &pis, const sam_data &dat){

	for(int g=0; g<(dat.nG); g++){
		for(int s=0;s<(dat.nS);s++){
			//calculate for alphas
    		alphaDerivs.at(s) +=  exp(llSG.at(MATREF2D(g,s,dat.nG)) - llS.at(s) + log(pis.at(g))) * dlogdalpha.at(MATREF2D(g,s,dat.nG));
			}
	}
 
}

// this should calculate the derivate w.r.t beta.
void calc_beta_deriv( vector<double> &betaDerivs, vector<double> const &dlogdbeta, vector<double> const &llSG, vector<double> const &llS, vector<double> const &pis, const sam_data &dat){
	
	for(int g=0; g<(dat.nG); g++){
	    for(int j=0; j<(dat.nP); j++){
			for(int s=0; s<(dat.nS); s++){
			// calculate for betas.
			betaDerivs.at(MATREF2D(g,j,dat.nG)) +=  exp(llSG.at(MATREF2D(g,s,dat.nG)) - llS.at(s) + log(pis.at(g))) * dlogdbeta.at(MATREF3D(g,j,s,dat.nG,dat.nP));
			}
		}
	}

}

//// this should calculate the derivate w.r.t dispersion parameter.
void calc_dispersion_deriv( vector<double> &dispDerivs, vector<double> const &dlogddispersionS, vector<double> const &llSG, vector<double> const &llS, vector<double> const &pis, const sam_data &dat){

	for(int g=0; g<(dat.nG); g++){
		for(int s=0;s<(dat.nS);s++){
			//calculate for dispersion (thetas)
    		dispDerivs.at(s) +=  exp(llSG.at(MATREF2D(g,s,dat.nG)) - llS.at(s) + log(pis.at(g))) * dlogddispersionS.at(MATREF2D(g,s,dat.nG));
			}
	}
 
}


// this should calculate the derivate w.r.t eta (transformed pi).
void calc_eta_deriv( vector<double> &etaDerivs, vector<double> const &dlogdpi, vector<double> const eta, const sam_data &dat){
	
  double add_log_trans=0;
  vector<double> pi_mat_deriv(dat.nG*(dat.nG-1),0);

  for(int g=0;g<(dat.nG-1);g++){
      add_log_trans+=exp(eta.at(g)); //add up transformed pi's
      //Rprintf( " %f\n", eta.at(g));
  }
  add_log_trans+=1;

	for(int g=0; g<(dat.nG); g++){
		for(int i=0;i<(dat.nG-1);i++){ // go through eta's
		   if(g<(dat.nG-1)){
			   	if(i==g){
					pi_mat_deriv.at(MATREF2D(i,g,(dat.nG-1))) = exp(eta.at(i))/add_log_trans - exp(2*eta.at(i))/(add_log_trans*add_log_trans);// diag
					pi_mat_deriv.at(MATREF2D(i,(dat.nG-1),(dat.nG-1))) += pi_mat_deriv.at(MATREF2D(i,g,(dat.nG-1)));
				}else{
					pi_mat_deriv.at(MATREF2D(i,g,(dat.nG-1))) = -exp(eta.at(i))*exp(eta.at(g)) / (add_log_trans*add_log_trans); //off-diag
					pi_mat_deriv.at(MATREF2D(i,(dat.nG-1),(dat.nG-1))) += pi_mat_deriv.at(MATREF2D(i,g,(dat.nG-1)));
				}
				//Rprintf( " %f\n",pi_mat_deriv.at(MATREF2D(i,(dat.nG-1),(dat.nG-1)))); 
			}
			//Rprintf( " %f\n",pi_mat_deriv.at(MATREF2D(i,g,(dat.nG-1))));
		}
	}
	for(int i=0;i<(dat.nG-1);i++) pi_mat_deriv.at(MATREF2D(i,(dat.nG-1),(dat.nG-1))) *= -1;


	for(int i=0; i<(dat.nG-1); i++){
			//Rprintf( "before dfdeta %f\n", etaDerivs.at(i));
		for(int g=0; g<(dat.nG); g++){
 
			etaDerivs.at(i) += dlogdpi.at(g) * pi_mat_deriv.at(MATREF2D(i,g,(dat.nG-1)));
			
		}
		//Rprintf( "after dfdeta %f\n", etaDerivs.at(i));
    }
}


bool converged_sam( double *oldP, double *newP, const sam_opt_contr &contr, int nTot){
	double tmp, eps=1e-8;

	for( int i=0; i<nTot; i++){
		tmp = fabs(newP[i]-oldP[i]);
		tmp /= fabs(oldP[i])+eps;
		if( tmp > contr.reltol)
			return( false);
	}
	return( true);

}

//inserve logistic link function
double inverse_logit(const double eta){
	double tmp;
	tmp = exp(eta);
	tmp = tmp / (1+tmp);
	return( tmp);
}
