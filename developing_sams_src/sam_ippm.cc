#include"sam_ippm.h"

/* Code for inhomogenous poisson point process model.
 * This should be very similar to the negative binomial with y_is_na added in and z and setups data classes for negative binomial.   
 *
 * I have tried to set this up like RCP, which makes more sense to me, in the future I can adapt Piers code into a single species_mix_cpp function.
 */

// this is the external C call which will be called by R using .Call.

extern "C" { 
		SEXP species_mix_ippm_cpp(SEXP Ry, SEXP RX, SEXP Roffset, SEXP Rwts, SEXP Ry_not_na,
								  SEXP RnS, SEXP RnG, SEXP RnP, SEXP RnObs,
								  SEXP Ralpha, SEXP Rbeta, SEXP Rtau, 
								  SEXP RderivsAlpha, SEXP RderivsBeta, SEXP RderivsTau, SEXP Rscores,
								  SEXP Rpis, SEXP Rmus, SEXP RlogDens, SEXP Rloglike,
								  SEXP Rmaxit, SEXP Rtrace, SEXP RnReport, SEXP Rabstol, SEXP Rreltol, SEXP Rconv,
								  SEXP Roptimise, SEXP RloglOnly, SEXP RderivsOnly, SEXP RoptiDisp, SEXP RgetScores){
	
	sam_ippm_all_classes all;

	//initialise the data structures -- they are mostly just pointers to REAL()s...
	all.data.setVals(Ry, RX, Roffset, Rwts, Ry_not_na, RnS, RnG, RnP, RnObs);	//read in the data
	all.parms.setVals(all.data, Ralpha, Rbeta, Rtau, Rgamma, Rdisps, Rpowers, Rconc, Rsd, RsdGamma, RdispLocat, RdispScale);	//read in the parameters
	all.derivs.setVals(all.data, RderivsAlpha, RderivsTau, RderivsBeta, RderivsGamma, RderivsDisps, RgetScores, Rscores);
	all.contr.setVals( Rmaxit, Rtrace, RnReport, Rabstol, Rreltol, Rconv);
	all.fits.initialise(all.data.nObs, all.data.nG, all.data.nS, all.data.NAnum);

	double logl = -999999;
	
	//doing the optimising
	if( *INTEGER(Roptimise) == 1)
		logl = sam_optimise(all);
	//re-running to get pis and mus
	if( *INTEGER(RloglOnly) == 1)
		logl = sam_ippm_mix_loglike( all.data, all.parms, all.fits);
	//and derivatives (inlcuding scores, for empirical info, if requested)
	if( *INTEGER(RderivsOnly) == 1)
		sam_ippm_mix_logl_derivs( all.data, all.parms, all.derivs, all.fits);	
		
		//bundling up things to return
	//first the fitted pis
	double *tmpPi = REAL( Rpis);
	for( int i=0; i<all.data.nObs; i++)
		for( int g=0; g<all.data.nG; g++)
			tmpPi[MATREF2D(i,g,all.data.nObs)] = all.fits.allPis.at(i).at(g);
	//the fitted expectations
	double *tmpMus = REAL( Rmus);
	for( size_t i=0; i<all.fits.allMus.size(); i++)
		tmpMus[i] = all.fits.allMus.at(i);
	//the log conditional densities
	double *tmpDens = REAL( RlogDens);
	for( int g=0; g<all.data.nG;g++)
		for( int i=0; i<all.data.nObs; i++)
				tmpDens[MATREF2D(i,g,all.data.nObs)] = all.fits.allLogDens.at(i).at(g);
	//the logl contributions
	double *tmplogls = REAL( Rlogli);
	for( int i=0; i<all.data.nObs; i++)
		tmplogls[i] = all.fits.allLogls.at(i);
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


double sam_optimise( allClasses &all){
	double *vmminGrad, *vmminParms, *oldParms;	//arrays to pass to vmmin
	double *vmminParmsIn;	//del
	vmminParms = (double *) R_alloc(all.parms.nTot,sizeof(double));
	vmminParmsIn = (double *) R_alloc(all.parms.nTot,sizeof(double));	//Del
	oldParms = (double *) R_alloc(all.parms.nTot,sizeof(double));
	vmminGrad = (double *) R_alloc(all.parms.nTot,sizeof(double));
	int *myMask;
	vector<int> vecMask(all.parms.nTot, 1);
	double vmminLogl[1];

//	Rprintf( "Quasi-Newton iterations\n");
	all.parms.getArray( oldParms, all.data);
	myMask = &vecMask[0];
	//optimise
	all.parms.getArray( vmminParms, all.data);
	all.parms.getArray( vmminParmsIn, all.data);	//del
	vmmin( all.parms.nTot, vmminParms, vmminLogl, optimise_function_ippm_sam, gradient_function_ippm_sam, all.contr.maxitQN, all.contr.traceQN, myMask, all.contr.abstol, all.contr.reltol, all.contr.nReport, &all, &all.contr.fnKount, &all.contr.grKount, &all.contr.ifail);
//	nmmin( all.parms.nTot, vmminParmsIn, vmminParms, vmminLogl, optimise_function_rcp, &all.contr.ifail, all.contr.abstol, all.contr.reltol, &all, 1.0, 0.5, 2.0, all.contr.traceQN, &all.contr.fnKount, all.contr.maxitQN);
	//update parameters
	all.parms.update( vmminParms, all.data);
	gradient_function_ippm_sam(all.parms.nTot, vmminParms, vmminGrad, &all);
	all.derivs.update( vmminGrad, all.data);

	return( vmminLogl[0]);
}

double optimise_function_sam_ippm(int n, double *par, void *ex)
{
	allClasses *all = (allClasses *) ex;
	double logl;

	all->parms.update( par, all->data);
	logl = sam_ippm_mix_loglike( all->data, all->parms, all->fits);

	return( (0.0-logl));
}

double sam_ippm_mix_loglike( const sam_ippm_data &dat, const sam_ippm_params &params, sam_ippm_fits &fits){
	
	vector<double> log_density(dat.nS*dat.nG, dat.NAnum);//, pis( dat.nG, dat.NAnum);
	double tloglike, loglike;	//for coding purposes really -- note that this is NOT related to the spp design matrix W
	fits.zero( dat.NAnum);

	//calculate fitted values (constant over i)
	calc_mu_fits(fits.allMus, dat, parms);// this should return the fitted valyes for all species, sites and groups.
	additive_logistic(fits.estpi,1) // additive logistic transformation of pis.

	for( int s=0; s<dat.nS; s++){
		tloglike.at(0) = calc_ippmsam_loglike_per_species(estpi, log_densitySG, fits.allMus, dat, parms, s); // this should give the mix loglike. ippm weights are calculated in this bit.
		fits.allSpecies_loglike_contribution.at(s) = tloglike.at(0);
		loglike.at(0) += tloglike.at(0);
	}
	return(loglike);
}

// calculate mu fits for ippm this should calculate all the etas and mus for species and archetypes. 
void calc_mu_fits( vector<double> &fits, const sam_ipp_data &dat, const sam_ippm_params &parms){
	//fits is a G*S matrix of the fitted values if dat.npw==0 and a G*S*nObs array if npw>0
	//vector<double> newTau( dat.nG*dat.nS, dat.NAnum);
	vector<double> lps(dat.nG*dat.nS, dat.NAnum);	//the nG x nS intercepts
	double lp=0.0;	//the lin pred for the gth group, sth species and ith site

	//calcualte the G*S*n fits
	for( int g=0; g<dat.nG; g++){
		for( int s=0; s<dat.nS; s++){
			lps.at(MATREF2D(g,s,dat.nG)) = parms.Alpha[s];
			for( int i=0; i<dat.nObs; i++){
				// need logical flag which deals with NA data.
				if(dat.y_not_na[MATREF2D(i,s,dat.nObs)]>0){
				lp = lps.at(MATREF2D(g,s,dat.nG)) + dat.offset[i];
					for( int j=0;j<dat.np; j++){
						lp += parms.Beta[MATREF2D(g,j,(dat.nG-1))] * dat.X[MATREF2D(i,j,dat.nObs)];
				   	}
				fits.at( MATREF3D(i,s,g,dat.nObs,dat.nS)) = exp(lp);
				}
			}
		}
	}	

}

// now we want to calculate the species specific likelihoods.
double calc_ippmsam_loglike_per_species(vector<double> &estpi, vector<double> &log_densitySG, const vector<double> &fits, const sam_ippm_data &dat, const sam_ippm_params &parms, int s){
	
	double eps=0, glogl=0;
    
	//calcualte the G*S log conditional densities
    for( int g=0; g<dat.nG; g++){
		for(int i=0; i<dat.nObs; i++){
			if(dat.y_not_na[MATREF2D(i,s,dat.nObs)]>0){
		    log_densitySG.at(MATREF2D(g,s,dat.nG)) += log_ippm(dat.y[MATREF2D(i,s,dat.nObs)], fits.at(MATREF3D(i,s,g,dat.nObs,dat.nS)), dat.wts[MATREF2D(i,s,dat.nObs)]);
		    }
		}
	
	if(g==0) eps = log_densitySG.at(MATREF2D(g,s,dat.nG));
    if(log_densitySG.at(MATREF2D(g,s,dat.nG)) > eps) eps = log_densitySG.at(MATREF2D(g,s,dat.nG));
    
  }
  
  // this will calculate the species specific loglikelihoods based on sum across Gs.
  for(g=0;g<G;g++){
    fits.allSum_f_species.at(MAT_RF(g,s,G)) = log_densitySG.at(MATREF2D(g,s,dat.nG);
    glogl+= estpi.at(g)*exp(log_densitySG.at(MATREF2D(g,s,dat.nG) - eps);
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
void additive_logistic(vector< double > &x,int inv){
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
// This is the wrapper for vmmin and it takes the sam_ippm_all_classes to 
void gradient_function_sam_ippm(int n, double *par, double *gr, void *ex){
	sam_ippm_all_classes *all = (sam_ippm_all_classes *) ex;

	sam_ippm_mix_gradient_function( all->data, all->parms, all->derivs, all->fits);

	all->derivs.getArray(gr, all->data);

	for( int i=0; i<n; i++)
		gr[i] = 0-gr[i];
}


//// this function should work out the derivatives.
void sam_ippm_mix_gradient_function( const sam_ippm_data &dat, const sam_ippm_params &parms, sam_ippm_derivs &derivs, sam_ippm_fits &fits){
	
	vector<double> alphaDerivsI( dat.nS, dat.NAnum);
	vector<double> betaDerivsI( (dat.nG-1)*dat.np, dat.NAnum);
	vector<double> piDerivsI( dat.nG*dat.nG, dat.NAnum);


	//calculate fitted values (constant over i)
	calc_mu_fits(fits.allMus, dat, parms);// this should return the fitted valyes for all species, sites and groups.
	additive_logistic(fits.estpi,1) // additive logistic transformation of pis.
    
   	//derivate w.r.t alpha 
	calc_alpha_deriv(alphaDerivsI, dat, fits.allMus)
		
	//derivate w.r.t beta
	calc_beta_deriv(betaDerivI, dat, fits.allMus)		
				
	//derivate w.r.t pi 
	calc_pi_deriv(piDerivsI, dat, fits.allMus)

	//update the derivates for each g.
	derivs.updateDerivs( dat, alphaDerivsI, betaDerivsI, piDerivsI);	//if dat.disty specifies no dispersion then no place for disp derivs (and are not updated, of course)
	     	  
	(void)tmp1;	//tmp1 is not used again, this is just a little trick to avoid a compile warning.
}

double log_poisson_deriv( const double &y, const double &mu){
	double tmp, z;
	z = y/wts;
	tmp = z/mu;
	tmp -= 1;
	return( tmp);
}

// this should calculate the derivate w.r.t alpha.
void calc_dlog_dalpha(){
	
	vector<double> dlogdalpha(dat.nG*dat.nS, dat.NAnum);
	double tmp_lpd;
	
	for(int g=0; g<dat.nG; g++){	
		for(int s=0;s<dat.nS; s++){
			for(int i=0; i<dat.nObs; i++){
				if(dat.y_not_na[MATREF2D()]>0){
					tmp_lpd = log_poisson_deriv(dat.y[MATREF2D(i,s,dat.nObs)], fits.at(MATREF3D(i,s,g,dat.nObs,dat.nS));	
					dlogdalpha.at(MATREF2D(s,g,(dat.nG)) +=	(dat.wts[MATREF2D(i,s,dat.nObs)] * tmp_ldp * fits.at(MATREF3D(i,s,g,dat.nObs,dat.nS) * 1);
				}
			}
		}
	}				
}


//

void calc_alpha_deriv( vector<double> &alphaDerivsI, const sam_ippm_data &dat, const double &mu, int g){



	
	for(g=0; g<dat.nG; g++){	
		for(s=0;s<S;s++){
			//calculate for alphas
			alphaDerivsI.at(s) += exp(data->sum_f_species.at(MAT_RF(g,s,G)) - data->species_l_contrib.at(s)+ log(data->parpi.at(g)))* data->deriv_f_alphaS.at(MAT_RF(g,s,G));
			}
	}			
				
}

// this should calculate the derivate w.r.t beta.

void calc_beta_deriv( vector<double> &betaDerivsI, const sam_ippm_data &dat, const double &mu, int g){

	for(g=0; g<dat.nG; g++){	
	    for(j=0; j<dat.np; j++){
			for(s=0; s<dat.nS; s++){
			// calculate for betas. 
			ad_g.at((G-1)+ MAT_RF(g,j,G)) += exp( -1*data->species_l_contrib.at(s) + log(data->parpi.at(g)) + data->sum_f_species.at(MAT_RF(g,s,G))) *data->deriv_f_B.at(MAT_3D(g,j,s,G,Xc));
			if(j==0) dl_dpi.at(g) += exp(  -1*data->species_l_contrib.at(s)+ data->sum_f_species.at(MAT_RF(g,s,G)));
			}
		}
	}
}

// this should calculate the derivate w.r.t tau.
void calc_pi_deriv( vector<double> &piDerivsI, const sam_ippm_data &dat, const double &mu, int g){
	
  for(g=0;g<(G-1);g++){ 
      add_log_trans+=exp(x.at(g)); //add up transformed pi's
  }
  add_log_trans+=1;
	
	for(g=0; g<dat.nG; g++){
		for(i=0;i<(G-1);i++){ // go through eta's
		   if(g<(G-1)){
				if(i==g){
					pi_mat_deriv.at(MAT_RF(i,g,(G-1))) = exp(x.at(i))/add_log_trans - exp(2*x.at(i))/(add_log_trans*add_log_trans);// diag
					pi_mat_deriv.at(MAT_RF(i,(G-1),(G-1))) += pi_mat_deriv.at(MAT_RF(i,g,(G-1)));
				}else{
					pi_mat_deriv.at(MAT_RF(i,g,(G-1))) = -exp(x.at(i))*exp(x.at(g)) / (add_log_trans*add_log_trans); //off-diag
					pi_mat_deriv.at(MAT_RF(i,(G-1),(G-1))) += pi_mat_deriv.at(MAT_RF(i,g,(G-1)));
				}
			}
		}
   }
}	
