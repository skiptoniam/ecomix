#include"sam_ippm.h"

/* Code for inhomogenous poisson point process model.
 * This should be very similar to the negative binomial with y_is_na added in and z and setups data classes for negative binomial.
 *
 * I have tried to set this up like RCP, which makes more sense to me, in the future I can adapt Piers code into a single species_mix_cpp function.
 */

// this is the external C call which will be called by R using .Call.

extern "C" {
	SEXP species_mix_ippm_cpp(SEXP Ry, SEXP RX, SEXP Roffset, SEXP Rwts, SEXP Ry_not_na,
							  SEXP RnS, SEXP RnG, SEXP Rp, SEXP RnObs,
							  SEXP Ralpha, SEXP Rbeta, SEXP Rpi,
							  SEXP RderivsAlpha, SEXP RderivsBeta, SEXP RderivsPi, SEXP RgetScores, SEXP Rscores,
							  SEXP Rpis, SEXP Rmus, SEXP RlogliS, SEXP RlogliSG,
							  SEXP Rmaxit, SEXP Rtrace, SEXP RnReport, SEXP Rabstol, SEXP Rreltol, SEXP Rconv,
							  SEXP Roptimise, SEXP RloglOnly, SEXP RderivsOnly){

	sam_ippm_all_classes all;

	//initialise the data structures -- they are mostly just pointers to REAL()s...
	all.data.setVals(Ry, RX, Roffset, Rwts, Ry_not_na, RnS, RnG, Rp, RnObs);	//read in the data
	all.params.setVals(all.data, Ralpha, Rbeta, Rpi);	//read in the parameters
	all.derivs.setVals(all.data, RderivsAlpha, RderivsBeta, RderivsPi, RgetScores, Rscores);
	all.contr.setVals( Rmaxit, Rtrace, RnReport, Rabstol, Rreltol, Rconv);
	all.fits.initialise(all.data.nObs, all.data.nG, all.data.nS, all.data.nP, all.data.NAnum);

	double logl = -999999;

	//doing the optimising
	if( *INTEGER(Roptimise) == 1)
		logl = sam_ippm_optimise(all);
	//re-running to get pis and mus
	if( *INTEGER(RloglOnly) == 1)
		logl = sam_ippm_mix_loglike( all.data, all.params, all.fits);
	//and derivatives (inlcuding scores, for empirical info, if requested)
	if( *INTEGER(RderivsOnly) == 1)
		sam_ippm_mix_gradient_function( all.data, all.params, all.derivs, all.fits);
		
	//for( int i=0; i<(all.params.nTot); i++)
		//Rprintf( " %f", (REAL(Rscores))[i]);	

	//bundling up things to return - will need to change these...
	//first the fitted pis
	double *tmpPi = REAL( Rpis);
	for( int g=0; g<all.data.nG; g++)
			tmpPi[g] = all.fits.estpi.at(g);
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
	SEXP Rres;	//R object to return -- it is the logl!
	Rres = PROTECT( allocVector(REALSXP,1));
    REAL( Rres)[0] = logl;
	UNPROTECT(1);
	return( Rres);

  }
}


double sam_ippm_optimise( sam_ippm_all_classes &all){
	double *vmminGrad, *vmminParms, *oldParms;	//arrays to pass to vmmin
	double *vmminParmsIn;	//del
	vmminParms = (double *) R_alloc(all.params.nTot,sizeof(double));
	vmminParmsIn = (double *) R_alloc(all.params.nTot,sizeof(double));	//Del
	oldParms = (double *) R_alloc(all.params.nTot,sizeof(double));
	vmminGrad = (double *) R_alloc(all.params.nTot,sizeof(double));
	int *myMask;
	vector<int> vecMask(all.params.nTot, 1);
	double vmminLogl[1];

//	Rprintf( "Quasi-Newton iterations\n");
	all.params.getArray( oldParms, all.data);
	myMask = &vecMask[0];
	//optimise
	all.params.getArray( vmminParms, all.data);
	all.params.getArray( vmminParmsIn, all.data);	//del
	vmmin( all.params.nTot, vmminParms, vmminLogl, optimise_function_ippm, gradient_function_ippm, all.contr.maxitQN, 
	all.contr.traceQN, myMask,  all.contr.abstol, all.contr.reltol,  all.contr.nReport, &all, &all.contr.fnKount, &all.contr.grKount, &all.contr.ifail);

//	nmmin( all.params.nTot, vmminParmsIn, vmminParms, vmminLogl, optimise_function_rcp, &all.contr.ifail, all.contr.abstol, all.contr.reltol, &all, 1.0, 0.5, 2.0, all.contr.traceQN, &all.contr.fnKount, all.contr.maxitQN);

	//update parameters
	all.params.update( vmminParms, all.data);
	gradient_function_ippm(all.params.nTot, vmminParms, vmminGrad, &all);
	all.derivs.update( vmminGrad, all.data);

	return( vmminLogl[0]);
}

double optimise_function_ippm(int n, double *par, void *ex){
	sam_ippm_all_classes *all = (sam_ippm_all_classes *) ex;
	double logl;

	all->params.update( par, all->data);
	logl = sam_ippm_mix_loglike( all->data, all->params, all->fits);

	return( (0.0-logl));
}

double sam_ippm_mix_loglike(sam_ippm_data &dat, sam_ippm_params &params, sam_ippm_fits &fits){

	vector<double> log_density(dat.nS*dat.nG, dat.NAnum);//, pis( dat.nG, dat.NAnum);
	double tloglike, loglike;	//for coding purposes really -- note that this is NOT related to the spp design matrix W
	fits.zero( 0);

	//calculate fitted values (constant over i)
	calc_mu_fits(fits.allMus, dat, params);// this should return the fitted valyes for all species, sites and groups.]
	for(int g=0; g<(dat.nG); g++) fits.estpi.at(g) = params.Pi[g];
	
	for(int g=0; g<(dat.nG); g++) //Rprintf( " %f", fits.estpi.at(g),"\n");
	//Rprintf("\n");
	additive_logistic_ippm(fits.estpi,1); // additive logistic transformation of pis.
	
	for(int g=0; g<(dat.nG); g++) //Rprintf( " %f", fits.estpi.at(g));

	for( int s=0; s<dat.nS; s++){
		tloglike = calc_ippm_loglike_per_species(dat, params, fits, s); // this should give the mix loglike. ippm weights are calculated in this bit.
		fits.log_like_species_contrib.at(s) = tloglike;
		loglike += tloglike;
	}
	return(loglike);
}

// calculate mu fits for ippm this should calculate all the etas and mus for species and archetypes.
void calc_mu_fits( vector<double> &fits, const sam_ippm_data &dat, const sam_ippm_params &params){
	//fits is a G*S matrix of the fitted values if dat.npw==0 and a G*S*nObs array if npw>0
	//vector<double> newTau( dat.nG*dat.nS, dat.NAnum);
	vector<double> lps(dat.nG*dat.nS, 0);	//the nG x nS intercepts
	double lp=0.0;	//the lin pred for the gth group, sth species and ith site

	//calcualte the G*S*n fits
	for( int g=0; g<dat.nG; g++){
		for( int s=0; s<dat.nS; s++){
			lps.at(MATREF2D(g,s,dat.nG)) = params.Alpha[s];
			for( int i=0; i<dat.nObs; i++){
				// need logical flag which deals with NA data.
				if(dat.y_not_na[MATREF2D(i,s,dat.nObs)]>0){
				lp = lps.at(MATREF2D(g,s,dat.nG)) + dat.offset[i];
					for( int j=0;j<dat.nP; j++){
						lp += params.Beta[MATREF2D(g,j,dat.nG)] * dat.X[MATREF2D(i,j,dat.nObs)];
				   	}
				fits.at( MATREF3D(i,s,g,dat.nObs,dat.nS)) = exp(lp);
				}
			}
		}
	}

}

// now we want to calculate the species specific likelihoods.
double calc_ippm_loglike_per_species(sam_ippm_data &dat, sam_ippm_params &params, sam_ippm_fits &fits, int s){

	double eps=0, glogl=0;

	//calcualte the G*S log conditional densities
    for( int g=0; g<dat.nG; g++){
		for(int i=0; i<dat.nObs; i++){
			if(dat.y_not_na[MATREF2D(i,s,dat.nObs)]>0){
		    fits.log_like_species_group_contrib.at(MATREF2D(g,s,dat.nG)) += log_ippm(dat.y[MATREF2D(i,s,dat.nObs)], fits.allMus.at(MATREF3D(i,s,g,dat.nObs,dat.nS)), dat.wts[MATREF2D(i,s,dat.nObs)]);
		    }
		}

	if(g==0) eps = fits.log_like_species_group_contrib.at(MATREF2D(g,s,dat.nG));
    if(fits.log_like_species_group_contrib.at(MATREF2D(g,s,dat.nG)) > eps) eps = fits.log_like_species_group_contrib.at(MATREF2D(g,s,dat.nG));

  }

  // this will calculate the species specific loglikelihoods based on sum across Gs.
  for(int g=0; g<dat.nG; g++){
    fits.log_like_species_group_contrib.at(MATREF2D(g,s,dat.nG)) = fits.log_like_species_group_contrib.at(MATREF2D(g,s,dat.nG));
    glogl += fits.estpi.at(g)*exp(fits.log_like_species_group_contrib.at(MATREF2D(g,s,dat.nG)) - eps);
  }
  glogl = log(glogl) + eps;

  return(glogl);

}

// this should give the log-density for the ippm. fingers crossed...
double log_ippm(const double &y, const double &mu, const double &wts){
	double tmp, z;
	z = y/wts;
	tmp = z * log(mu);
	tmp -= mu;
	tmp *= wts;
	return( tmp);
}

// this should do the additive transformation of pis (pi1/piN,pi2/piN,pi(N-1)/piN)
void additive_logistic_ippm(vector< double > &x,int inv){
  int i;
  // inv == 1 gives transfornmation
  // inv == 0 gives inverse transformation
  if(inv==0){
    for(i=0;i<x.size();i++) x.at(i) = log(x.at(i)/x.back());
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
// This is the wrapper for vmmin and it takes the sam_ippm_all_classes to
void gradient_function_ippm(int n, double *par, double *gr, void *ex){
	sam_ippm_all_classes *all = (sam_ippm_all_classes *) ex;

    //Rprintf(all.derives) check rcp mod and ex.
	sam_ippm_mix_gradient_function( all->data, all->params, all->derivs, all->fits);

	all->derivs.getArray(gr, all->data);

	for( int i=0; i<n; i++)
		gr[i] = 0-gr[i];
		//Rprintf(gr[i])
}


//// this function should work out the derivatives.
void sam_ippm_mix_gradient_function( const sam_ippm_data &dat, const sam_ippm_params &params, sam_ippm_derivs &derivs, sam_ippm_fits &fits){

	double tmp1;
	vector<double> alphaDerivs(dat.nS, 0);
	vector<double> betaDerivs((dat.nG*dat.nP), 0);
	vector<double> piDerivs(dat.nG*dat.nG, 0); // check there should only be g pis


	//calculate fitted values
	calc_mu_fits(fits.allMus, dat, params);// this should return the fitted valyes for all species, sites and groups.
	//additive_logistic_ippm(fits.estpi,1); // additive logistic transformation of pis.
	
	//derivate w.r.t alpha
   	calc_dlog_dalpha(dat, fits);
	calc_alpha_deriv(alphaDerivs, dat, fits);
	for (int s=0; s<dat.nS;s++)
	Rprintf( " %f\n", alphaDerivs[s]);

	//derivate w.r.t beta
	calc_dlog_dbeta(dat, fits);
	calc_beta_deriv(betaDerivs, dat, fits);
	
	for(int g=0; g<(dat.nG); g++){
	    for(int j=0; j<(dat.nP); j++){
			Rprintf( "%f\n", betaDerivs.at(MATREF2D(g,j,dat.nG)));
		}
	}
	//derivate w.r.t pi
	calc_pi_deriv(piDerivs, dat, derivs, fits);
	
	for(int g=0; g<(dat.nG); g++)
		Rprintf( " %f\n", piDerivs[g]);

	//update the derivates.
	derivs.updateDerivs( dat, alphaDerivs, betaDerivs, piDerivs);	//if dat.disty specifies no dispersion then no place for disp derivs (and are not updated, of course)

	//(void)tmp1;	//tmp1 is not used again, this is just a little trick to avoid a compile warning.
}

double log_poisson_deriv( const double &y, const double &mu, const double &wts){
	double tmp, z;
	z = y/wts;
	tmp = z/mu;
	tmp -= 1;
	return( tmp);
}


void calc_dlog_dalpha(const sam_ippm_data &dat, sam_ippm_fits &fits){

	//vector<double> dlogdalpha(dat.nG*dat.nS, dat.NAnum);
	//vector<double> muGS(dat.nG*dat.nS, dat.NAnum);
	double tmp_lpd;

	for(int g=0; g<dat.nG; g++){
		for(int s=0;s<dat.nS; s++){
			for(int i=0; i<dat.nObs; i++){
				if(dat.y_not_na[MATREF2D(i,s,dat.nObs)]>0){
					tmp_lpd = log_poisson_deriv(dat.y[MATREF2D(i,s,dat.nObs)], fits.allMus.at(MATREF3D(i,s,g,dat.nObs,dat.nS)), dat.wts[MATREF2D(i,s,dat.nObs)]);
					//muGS.at(MATREF2D(s,g,dat.nG)) =+ fits.at(MATREF3D(i,s,g,dat.nObs,dat.nS));	// why did i have this? maybe for pis?
					fits.dlogdalpha.at(MATREF2D(g,s,dat.nG)) +=	(dat.wts[MATREF2D(i,s,dat.nObs)] * tmp_lpd * fits.allMus.at(MATREF3D(i,s,g,dat.nObs,dat.nS)) * 1);
				}
			}
		}
	}
}

void calc_dlog_dbeta(const sam_ippm_data &dat, sam_ippm_fits &fits){

	//vector<double> dlogdbeta(dat.nG*dat.nS*dat.nP, dat.NAnum);
	//vector<double> muGS(dat.nG*dat.nS, dat.NAnum);
	double tmp_lpd;

	for(int g=0; g<dat.nG; g++){
		for(int s=0;s<dat.nS; s++){
			for(int i=0; i<dat.nObs; i++){
				if(dat.y_not_na[MATREF2D(i,s,dat.nObs)]>0){
					tmp_lpd = log_poisson_deriv(dat.y[MATREF2D(i,s,dat.nObs)], fits.allMus.at(MATREF3D(i,s,g,dat.nObs,dat.nS)), dat.wts[MATREF2D(i,s,dat.nObs)]);
						for(int j=0; j<dat.nP; j++){
							fits.dlogdbeta.at(MATREF2D(g,j,dat.nG)) +=	(dat.wts[MATREF2D(i,s,dat.nObs)] * tmp_lpd * fits.allMus.at(MATREF3D(i,s,g,dat.nObs,dat.nS)) * dat.X[MATREF2D(i,j,dat.nObs)]);
					}
				}
			}
		}
	}

}



//// this should calculate the derivate w.r.t alpha.
void calc_alpha_deriv( vector<double> &alphaDerivs, const sam_ippm_data &dat, sam_ippm_fits &fits){

	for(int g=0; g<(dat.nG); g++){
		for(int s=0;s<(dat.nS);s++){
			//calculate for alphas
			////	for( int i=0; i<(all.data.nObs*all.parms.nTot); i++)
    		alphaDerivs.at(s) +=  exp(fits.log_like_species_group_contrib.at(MATREF2D(g,s,dat.nG)) - fits.log_like_species_contrib.at(s) + log(fits.estpi.at(g))) * fits.dlogdalpha.at(MATREF2D(g,s,dat.nG));
    		
			}
	}
	//Rprintf( "\n");
 
}

// this should calculate the derivate w.r.t beta.
void calc_beta_deriv( vector<double> &betaDerivs, const sam_ippm_data &dat, sam_ippm_fits &fits){

    //vector<double> dlogdpi(dat.nG,0);

	for(int g=0; g<(dat.nG); g++){
	    for(int j=0; j<(dat.nP); j++){
			for(int s=0; s<(dat.nS); s++){
			// calculate for betas.
			betaDerivs.at(MATREF2D(g,j,dat.nG)) +=  exp(fits.log_like_species_group_contrib.at(MATREF2D(g,s,dat.nG)) - fits.log_like_species_contrib.at(s) + log(fits.estpi.at(g))) * fits.dlogdbeta.at(MATREF2D(g,j,dat.nG));
			//Rprintf( " %f\n", betaDerivs.at(MATREF2D(g,j,dat.nG)));
			if(j==0) fits.dlogdpi.at(g) += exp(fits.log_like_species_group_contrib.at(MATREF2D(g,s,dat.nG)) - fits.log_like_species_contrib.at(s) + log(fits.estpi.at(g)));
			}
		}
	}
	
	//Rprintf( "\n");
}

// this should calculate the derivate w.r.t tau.
void calc_pi_deriv( vector<double> &piDerivs, const sam_ippm_data &dat, sam_ippm_derivs &derivs, sam_ippm_fits &fits){

  double add_log_trans=0;
  vector<double> grad_pis(dat.nG,0);
  vector<double> pi_mat_deriv(dat.nG*(dat.nG-1),0);

  for(int i=0;i<(dat.nG-1);i++){
	grad_pis.at(i) = fits.estpi.at(i);
	//Rprintf( " %f", grad_pis.at(i));
	}

  for(int g=0;g<(dat.nG-1);g++){
      add_log_trans+=exp(grad_pis.at(g)); //add up transformed pi's
  }
  add_log_trans+=1;

	for(int g=0; g<(dat.nG); g++){
		for(int i=0;i<(dat.nG-1);i++){ // go through eta's
		   if(g<(dat.nG-1)){
				if(i==g){
					pi_mat_deriv.at(MATREF2D(i,g,(dat.nG-1))) = exp(grad_pis.at(i))/add_log_trans - exp(2*grad_pis.at(i))/(add_log_trans*add_log_trans);// diag
					pi_mat_deriv.at(MATREF2D(i,(dat.nG-1),(dat.nG-1))) += pi_mat_deriv.at(MATREF2D(i,g,(dat.nG-1)));
				}else{
					pi_mat_deriv.at(MATREF2D(i,g,(dat.nG-1))) = -exp(grad_pis.at(i))*exp(grad_pis.at(g)) / (add_log_trans*add_log_trans); //off-diag
					pi_mat_deriv.at(MATREF2D(i,(dat.nG-1),(dat.nG-1))) += pi_mat_deriv.at(MATREF2D(i,g,(dat.nG-1)));
				}
			}
		}
	}
	for(int i=0;i<(dat.nG-1);i++) pi_mat_deriv.at(MATREF2D(i,(dat.nG-1),(dat.nG-1))) *= -1;

	for(int i=0; i<(dat.nG-1); i++){
		for(int g=0; g<dat.nG; g++){
			piDerivs.at(i) += fits.dlogdpi.at(g)* pi_mat_deriv.at(MATREF2D(i,g,(dat.nG-1)));
			//Rprintf( " %f", piDerivs.at(i));
		}
    }
}


bool converged_ippm( double *oldP, double *newP, const sam_ippm_opt_contr &contr, int nTot)
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
