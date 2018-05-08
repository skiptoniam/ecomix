#include"sam_bernoulli_sp_ints.h"

/* Code for inhomogenous poisson point process model.
 * I have tried to set this up like RCP, which makes more sense to me, in the future I can adapt Piers code into a single species_mix_cpp function.
 */

// this is the external C call which will be called by R using .Call.

extern "C" {
	SEXP species_mix_bernoulli_sp_ints(SEXP Ry, SEXP RX, SEXP Roffset, SEXP Rwts,
									   SEXP RnS, SEXP RnG, SEXP Rp, SEXP RnObs,
									   SEXP Ralpha, SEXP Rbeta, SEXP Reta, 
									   SEXP RderivsAlpha, SEXP RderivsBeta, SEXP RderivsEta, SEXP RgetScores, SEXP Rscores,
									   SEXP Rpis, SEXP Rmus, SEXP RlogliS, SEXP RlogliSG,
									   SEXP Rmaxit, SEXP Rtrace, SEXP RnReport, SEXP Rabstol, SEXP Rreltol, SEXP Rconv, SEXP Rprintparams,
									   SEXP Roptimise, SEXP RloglOnly, SEXP RderivsOnly){

	sam_bernoulli_sp_ints_all_classes all;

	//initialise the data structures -- they are mostly just pointers to REAL()s...
	all.data.setVals(Ry, RX, Roffset, Rwts, RnS, RnG, Rp, RnObs);	//read in the data
	all.params.setVals(all.data, Ralpha, Rbeta, Reta);	//read in the parameters
	all.derivs.setVals(all.data, RderivsAlpha, RderivsBeta, RderivsEta, RgetScores, Rscores);
	all.contr.setVals( Rmaxit, Rtrace, RnReport, Rabstol, Rreltol, Rconv, Rprintparams);
	all.fits.initialise(all.data.nObs, all.data.nG, all.data.nS, all.data.nP, 0);

	double logl = -999999;

	//doing the optimising
	if( *INTEGER(Roptimise) == 1)
		logl = sam_bernoulli_sp_ints_optimise(all);
	//re-running to get pis and mus
	if( *INTEGER(RloglOnly) == 1)
	   		logl = sam_bernoulli_sp_ints_mix_loglike( all.data, all.params, all.fits);
	//and derivatives (inlcuding scores, for empirical info, if requested)
	if( *INTEGER(RderivsOnly) == 1)	    
	    sam_bernoulli_sp_ints_mix_gradient_function( all.data, all.params, all.derivs, all.fits);
	
	//bundling up things to return - will need to change these...
	
	//first the fitted pis
	double *tmpPi = REAL( Reta);
	
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
    /* Construct named result list from variables containing the results */
    SEXP Rlogl = PROTECT(allocVector(REALSXP, 1));
    REAL(Rlogl)[0] = logl;
    SEXP Ralpha_est = PROTECT(allocVector(REALSXP, all.data.nS));
    for( int s=0; s<all.data.nS; s++) REAL(Ralpha_est)[s] = all.params.Alpha[s];
	SEXP Rbeta_est = PROTECT(allocVector(REALSXP, all.data.nG*all.data.nP));
	for( int i=0; i<((all.data.nG*all.data.nP)); i++) REAL(Rbeta_est)[i] = all.params.Beta[i];
	SEXP Reta_est =PROTECT(allocVector(REALSXP, all.data.nG-1));
	for( int g=0; g<(all.data.nG-1);g++) REAL(Reta_est)[g] = all.params.Eta[g];
	
	const char *names[] = {"logl", "alpha", "beta", "eta", ""};                   /* note the null string */
	SEXP Rres = PROTECT(mkNamed(VECSXP, names));  /* list of length 3 */
	SET_VECTOR_ELT(Rres, 0, Rlogl);       /* numeric(1) */ 
	SET_VECTOR_ELT(Rres, 1, Ralpha_est);   /* numeric(<some length>) */ 
	SET_VECTOR_ELT(Rres, 2, Rbeta_est);    /* integer(1) */
	SET_VECTOR_ELT(Rres, 3, Reta_est);
	UNPROTECT(1);
	return (Rres);
    


  }
}


double sam_bernoulli_sp_ints_optimise( sam_bernoulli_sp_ints_all_classes &all){
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
	vmmin( all.params.nTot, vmminParams, vmminLogl, optimise_function_bernoulli_sp_ints, gradient_function_bernoulli_sp_ints, all.contr.maxitQN, 
	 all.contr.traceQN, myMask,  all.contr.abstol, all.contr.reltol,  all.contr.nReport, &all, &all.contr.fnKount,
	 &all.contr.grKount, &all.contr.ifail);

//	nmmin( all.params.nTot, vmminParamsIn, vmminParams, vmminLogl, optimise_function_rcp, &all.contr.ifail, all.contr.abstol, all.contr.reltol, &all, 1.0, 0.5, 2.0, all.contr.traceQN, &all.contr.fnKount, all.contr.maxitQN);

	//update parameters
	all.params.update( vmminParams, all.data);
	gradient_function_bernoulli_sp_ints(all.params.nTot, vmminParams, vmminGrad, &all);
	all.derivs.update( vmminGrad, all.data);
    if(all.contr.printparams==1)all.params.printParms(all.data);

	return(vmminLogl[0]);
}

double optimise_function_bernoulli_sp_ints(int n, double *par, void *ex){
	
	sam_bernoulli_sp_ints_all_classes *all = (sam_bernoulli_sp_ints_all_classes *) ex;
	double logl;

	all->params.update( par, all->data);
	logl = sam_bernoulli_sp_ints_mix_loglike( all->data, all->params, all->fits);
     
    return((0.0-logl));
}

double sam_bernoulli_sp_ints_mix_loglike(const sam_bernoulli_sp_ints_data &dat, const sam_bernoulli_sp_ints_params &params, sam_bernoulli_sp_ints_fits &fits){

	double tloglike = 0.0, loglike = 0.0;
	vector<double> par_pi(dat.nG-1,0);
	
	fits.zero(0);
	for(int g=0; g<(dat.nG-1); g++) par_pi.at(g) = params.Eta[g];

	//transform additative pis to natural scale - need this to calc loglikes.	
	additive_logistic_bernoulli_sp_ints(par_pi,1,dat.nG); // additive logistic transformation of pis.

	//calculate fitted values (constant over i)
	calc_mu_fits(fits.allMus, params, dat);

	//calculate the species/groups loglikes
	calc_bernoulli_loglike_SG(fits.log_like_species_group_contrib, fits.allMus, dat);
	
	// calc loglike per species
	for( int s=0; s<dat.nS; s++){
		tloglike = calc_bernoulli_loglike_S(fits.log_like_species_group_contrib, par_pi, dat, s); // this should give the mix loglike. bernoulli_sp_ints weights are calculated in this bit.
		fits.log_like_species_contrib.at(s) = tloglike;
		loglike += tloglike;
	}
	return(loglike);
}

// calculate mu fits for bernoulli_sp_ints this should calculate all the etas and mus for species and archetypes.
void calc_mu_fits(vector<double> &fits, const sam_bernoulli_sp_ints_params &params, const sam_bernoulli_sp_ints_data &dat){

	vector<double> lps(dat.nG*dat.nS, 0);	//the nG x nS intercepts
	double lp=0.0;	//the lin pred for the gth group, sth species and ith site

	//calcualte the G*S*n fits
	for( int g=0; g<dat.nG; g++){
		for( int s=0; s<dat.nS; s++){
			lps.at(MATREF2D(g,s,dat.nG)) = params.Alpha[s]; 
			for( int i=0; i<dat.nObs; i++){
				lp = lps.at(MATREF2D(g,s,dat.nG)) + dat.offset[i];
					for( int j=0;j<dat.nP; j++){
						lp += params.Beta[MATREF2D(g,j,(dat.nG))] * dat.X[MATREF2D(i,j,dat.nObs)];
				   	}
				fits.at( MATREF3D(i,s,g,dat.nObs,dat.nS)) = invLogit_bern(lp);
			}
		}
	}

}


void calc_bernoulli_loglike_SG(vector<double> &loglSG, vector<double> &fits, const sam_bernoulli_sp_ints_data &dat){
	
    //calcualte the G*S log conditional densities
	for(int s=0; s<dat.nS; s++){
		for( int g=0; g<dat.nG; g++){
			for(int i=0; i<dat.nObs; i++){
				loglSG.at(MATREF2D(g,s,dat.nG)) += log_bernoulli(dat.y[MATREF2D(i,s,dat.nObs)], fits.at(MATREF3D(i,s,g,dat.nObs,dat.nS)));//*dat.wts[i];
				}
				loglSG.at(MATREF2D(g,s,dat.nG)) = loglSG.at(MATREF2D(g,s,dat.nG))*dat.wts[s]; //fix this up.
			}
		}
}

// now we want to calculate the species specific likelihoods.
double calc_bernoulli_loglike_S(vector<double> &fits, vector<double> const &pis, const sam_bernoulli_sp_ints_data &dat, int s){

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

// this should give the log-density for the bernoulli.
double log_bernoulli( const double &y, const double &mu){
	double tmp;
	if( y==1){
		tmp = log( mu);
		return( tmp);
	} 
	tmp = log( 1-mu);
	return( tmp);
	
}

//inserve logistic link function
double invLogit_bern(const double eta){
	double tmp;
	tmp = exp(eta);
	tmp = tmp / (1+tmp);
	return( tmp);
}

// this should do the additive transformation of pis (pi1/piN,pi2/piN,pi(N-1)/piN)
void additive_logistic_bernoulli_sp_ints(vector< double > &x, int inv, int G){
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

// Gradient functions for bernoulli.
// These are all the functions for the gradient function and hopefully they should run and help estimate the derivates.
// This is the wrapper for vmmin and it takes the sam_bernoulli_sp_ints_all_classes to
void gradient_function_bernoulli_sp_ints(int n, double *par, double *gr, void *ex){
	sam_bernoulli_sp_ints_all_classes *all = (sam_bernoulli_sp_ints_all_classes *) ex;

    //Rprintf(all.derives) check rcp mod and ex.
	sam_bernoulli_sp_ints_mix_gradient_function( all->data, all->params, all->derivs, all->fits);

	all->derivs.getArray(gr, all->data);

	for( int i=0; i<n; i++){
	    gr[i] = 0.0 - gr[i];
	}
}


//// this function should work out the derivatives.
void sam_bernoulli_sp_ints_mix_gradient_function(const sam_bernoulli_sp_ints_data &dat, const sam_bernoulli_sp_ints_params &params, sam_bernoulli_sp_ints_derivs &derivs, sam_bernoulli_sp_ints_fits &fits){

	vector<double> parpi((dat.nG-1), 0);
	vector<double> alphaDerivs(dat.nS, 0);//change to dat.NAN
	vector<double> betaDerivs((dat.nG*dat.nP), 0);
	vector<double> etaDerivs((dat.nG-1), 0); // check there should only be g pis
	double logl;

    //calc loglike
    fits.zero(0);
    derivs.zeroDerivs(dat);
    
	logl = sam_bernoulli_sp_ints_mix_loglike(dat, params, fits);
	
	for(int g=0; g<(dat.nG-1); g++) parpi.at(g) = params.Eta[g];
	additive_logistic_bernoulli_sp_ints(parpi,1,dat.nG);
	
	//derivate w.r.t alpha
   	calc_dlog_dalpha(fits.dlogdalpha, fits.allMus, dat);
	calc_alpha_deriv(alphaDerivs, fits.dlogdalpha, fits.log_like_species_group_contrib, fits.log_like_species_contrib, parpi, dat);

	//derivate w.r.t beta
	calc_dlog_dbeta(fits.dlogdbeta, fits.allMus, dat);
	calc_beta_deriv(betaDerivs, fits.dlogdbeta, fits.log_like_species_group_contrib, fits.log_like_species_contrib, parpi, dat);
	
	//transform pis back to additative logistic scale to keep pi_dervis happy.
	additive_logistic_bernoulli_sp_ints(parpi,0,dat.nG);
	
	//derivate w.r.t pi/eta
	calc_dlog_dpi(fits.dlogdpi, fits.log_like_species_group_contrib, fits.log_like_species_contrib, dat);
	calc_eta_deriv(etaDerivs, fits.dlogdpi, parpi, dat);

	//update the derivates.
	derivs.updateDerivs( dat, alphaDerivs, betaDerivs, etaDerivs);	
	}

double log_bernoulli_deriv(const double &y, const double &mu){
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

double dmu_deta_bernoulli(const double &mu){
	
	double tmp;
	tmp = mu;
	tmp *= (1-mu);
	return(tmp);
	
}


void calc_dlog_dalpha(vector<double> &dlda, vector<double> const &mus, const sam_bernoulli_sp_ints_data &dat){

	// dlda = dlogalpha passed as fits.dflogdalpha(dat.nG*dat.nS, dat.NAnum) from function call
	// mus = all the fitted values.
	double tmp_lbd, tmp_dmde;

	for(int g=0; g<dat.nG; g++){
		for(int s=0;s<dat.nS; s++){
			for(int i=0; i<dat.nObs; i++){
				// this is the tmp log bernoulli derivative (lbd)
				tmp_lbd = log_bernoulli_deriv(dat.y[MATREF2D(i,s,dat.nObs)], mus.at(MATREF3D(i,s,g,dat.nObs,dat.nS)));//*dat.wts[i];
				tmp_dmde = dmu_deta_bernoulli(mus.at(MATREF3D(i,s,g,dat.nObs,dat.nS)));
				dlda.at(MATREF2D(g,s,dat.nG)) += (tmp_lbd * tmp_dmde * 1);
			}
			dlda.at(MATREF2D(g,s,dat.nG)) = dlda.at(MATREF2D(g,s,dat.nG))*dat.wts[s];
		}
	}
}

void calc_dlog_dbeta(vector<double> &dldb, vector<double> const &mus, const sam_bernoulli_sp_ints_data &dat){

	double tmp_lbd, tmp_dmde;

	for(int g=0; g<dat.nG; g++){
		for(int s=0;s<dat.nS; s++){
			for(int i=0; i<dat.nObs; i++){
					// calc the log bernoulli deriv
					tmp_lbd = log_bernoulli_deriv(dat.y[MATREF2D(i,s,dat.nObs)], mus.at(MATREF3D(i,s,g,dat.nObs,dat.nS)));//*dat.wts[i];
					tmp_dmde = dmu_deta_bernoulli(mus.at(MATREF3D(i,s,g,dat.nObs,dat.nS)));
						for(int j=0; j<dat.nP; j++){
							dldb.at(MATREF3D(g,j,s,dat.nG,dat.nP)) += (tmp_lbd * tmp_dmde * dat.X[MATREF2D(i,j,dat.nObs)]);
				}
			}
		for(int j=0; j<dat.nP; j++){
				dldb.at(MATREF3D(g,j,s,dat.nG,dat.nP)) = dldb.at(MATREF3D(g,j,s,dat.nG,dat.nP))*dat.wts[s];	
			}		
		}
	}

}

void calc_dlog_dpi(vector<double> &dldpi, vector<double> const &llSG, vector<double> const &llS, const sam_bernoulli_sp_ints_data &dat){

	for(int g=0; g<(dat.nG); g++){
			for(int s=0; s<(dat.nS); s++){
							dldpi.at(g) += exp(llSG.at(MATREF2D(g,s,dat.nG)) - llS.at(s));
			}
	}
}

//// this should calculate the derivate w.r.t alpha.
void calc_alpha_deriv( vector<double> &alphaDerivs, vector<double> const &dlogdalpha, vector<double> const &llSG, vector<double> const &llS, vector<double> const &pis, const sam_bernoulli_sp_ints_data &dat){

	for(int g=0; g<(dat.nG); g++){
		for(int s=0;s<(dat.nS);s++){
			//calculate for alphas
    		alphaDerivs.at(s) +=  exp(llSG.at(MATREF2D(g,s,dat.nG)) - llS.at(s) + log(pis.at(g))) * dlogdalpha.at(MATREF2D(g,s,dat.nG));
			}
	}
 
}

// this should calculate the derivate w.r.t beta.
void calc_beta_deriv( vector<double> &betaDerivs, vector<double> const &dlogdbeta, vector<double> const &llSG, vector<double> const &llS, vector<double> const &pis, const sam_bernoulli_sp_ints_data &dat){
	
	for(int g=0; g<(dat.nG); g++){
	    for(int j=0; j<(dat.nP); j++){
			for(int s=0; s<(dat.nS); s++){
			// calculate for betas.
			betaDerivs.at(MATREF2D(g,j,dat.nG)) +=  exp(llSG.at(MATREF2D(g,s,dat.nG)) - llS.at(s) + log(pis.at(g))) * dlogdbeta.at(MATREF3D(g,j,s,dat.nG,dat.nP));
			}
		}
	}

}

// this should calculate the derivate w.r.t eta (transformed pi).
void calc_eta_deriv( vector<double> &etaDerivs, vector<double> const &dlogdpi, vector<double> const eta, const sam_bernoulli_sp_ints_data &dat){
	
  double add_log_trans=0;
  vector<double> pi_mat_deriv(dat.nG*(dat.nG-1),0);

  for(int g=0;g<(dat.nG-1);g++){
      add_log_trans+=exp(eta.at(g)); //add up transformed pi's

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
			}
		}
	}
	for(int i=0;i<(dat.nG-1);i++) pi_mat_deriv.at(MATREF2D(i,(dat.nG-1),(dat.nG-1))) *= -1;


	for(int i=0; i<(dat.nG-1); i++){
		for(int g=0; g<(dat.nG); g++){
 			etaDerivs.at(i) += dlogdpi.at(g) * pi_mat_deriv.at(MATREF2D(i,g,(dat.nG-1)));
		}
    }
}


bool converged_bernoulli( double *oldP, double *newP, const sam_bernoulli_sp_ints_opt_contr &contr, int nTot){
	double tmp, eps=1e-8;

	for( int i=0; i<nTot; i++){
		tmp = fabs(newP[i]-oldP[i]);
		tmp /= fabs(oldP[i])+eps;
		if( tmp > contr.reltol)
			return( false);
	}
	return( true);

}
