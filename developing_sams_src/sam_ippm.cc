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

double sam_ippm_mix_loglike( const sam_ippm_data &dat, const sam_ippm_params &parms, sam_ippm_fits &fits){
	
	vector<double> logPis(dat.nG, dat.NAnum);//, pis( dat.nG, dat.NAnum);
	double loglike;	//for coding purposes really -- note that this is NOT related to the spp design matrix W
	vector<double> wij( dat.nG, dat.NAnum);	//for coding purposes -- note that this is NOT related to the spp design matrix W
	int m;	//for coding purposes

	fits.zero( dat.NAnum);

	//calculate fitted values (constant over i)
	calc_mu_fits(fits.allMus, dat, parms);// this should return the fitted valyes for all species, sites and groups.
	
	for( int i=0; i<dat.nObs; i++){
		calc_log_pis( logPis, fits.allPis.at(i), dat, parms, i); // this should estimate the pis per group
		calc_log_cond_dens( fits.allLogDens.at(i), fits.allMus, dat, parms, i); // this should give the log-densities for ippm. 
		loglike += calcMixSum( logPis, fits.allLogDens.at(i), wi, wij, m); // this should give the mix loglike. ippm weights are calculated in this bit.
	}
	return(loglike);
}

// calculate mu fits for ippm
// this should calculate all the etas and mus for species and archetypes. 
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

// this should calculate the log_pis for
calc_log_pis( logPis, fits.allPis.at(i), dat, parms, i); // this should estimate the pis per group

void calcLogPis( vector<double> &logPis, vector<double> &pis, const myData &dat, const myParms &parms, int i)
{
	vector<double> lp((dat.nG-1),0.0);
	double sumlp=0.0, sumpi=0.0;

	lp.assign( (dat.nG-1), 0.0);
	for( int k=0; k<(dat.nG-1); k++){
		for( int p=0; p<dat.np; p++)
			lp.at(k) += parms.Beta[MATREF2D(k,p,(dat.nG-1))] * dat.X[MATREF2D(i,p,dat.nObs)];
		lp.at(k) = exp( lp.at(k));
		sumlp += lp.at(k);
	}
	for( int k=0; k<(dat.nG-1); k++){
		pis.at(k) = lp.at(k) / ( 1+sumlp);
		sumpi += pis.at(k);
	}
	pis.at(dat.nG-1) = 1-sumpi;
	for( int k=0; k<dat.nG; k++)
		logPis.at(k) = log( pis.at(k));

	for( int k=0; k<dat.nG; k++){
		if( logPis.at(k)>=0)
			logPis.at(k) = -DBL_MIN;	//Smallest (absolute) non-zero number on your machine
		if( !R_FINITE(logPis.at(k)))
			logPis.at(k) = -DBL_MAX;	//Smallest (most negative) number on your machine
	}

}

get_taus_ippm <- function(pi, logls, G, S){
  fullLogPis <- matrix( rep( log(pi), each=S), nrow=S, ncol=G)
  a_k <- fullLogPis + logls
  a_m <- apply( a_k, 1, max)
  tmp <- exp( a_k - rep( a_m, times=G))
  log_denom <- a_m + log( rowSums( tmp))
  return( exp( a_k - log_denom))
}

colSums(taus)/S



// now we want to calculate the log densities. ippm has a slightly different derivation to logPoisson.

void calc_log_cond_dens( vector<double> &condDens, const vector<double> &fits, const sam_ippm_data &dat, const sam_ippm_params &parms, int i){
	vector<double> condDensSG( dat.nG*dat.nS, dat.NAnum);

	//calcualte the G*S log conditional densities
	for( int g=0; g<dat.nG; g++){
		for( int s=0; s<dat.nS; s++){
			if(dat.y_not_na[MATREF2D(i,s,dat.nObs)]>0){
		    condDensSG.at(MATREF2D(g,s,dat.nG)) = log_ippm(dat.y[MATREF2D(i,s,dat.nObs)], fits.at(MATREF3D(i,s,g,dat.nObs,dat.nS)), dat.wts[MATREF2D(i,s,dat.nObs)]);
		    }
		    // here are all the other options for other distributions. 
		}
//		Rprintf("\n");
	}
	//calculate the G log conditional densities (under independence)
	for( int g=0; g<dat.nG; g++){
		condDens.at(g) = 0.0;
		for( int s=0; s<dat.nS; s++)
			condDens.at(g) += condDensSG.at(MATREF2D(g,s,dat.nG));
	}
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

// Gradient functions for IPPM. 
// This is the wrapper for vmmin and it takes the sam_ippm_all_classes to 
void gradient_function_sam_ippm(int n, double *par, double *gr, void *ex){
	sam_ippm_all_classes *all = (sam_ippm_all_classes *) ex;

	logl_derivs_ippm( all->data, all->parms, all->derivs, all->fits);

	all->derivs.getArray(gr, all->data);

	for( int i=0; i<n; i++)
		gr[i] = 0-gr[i];
}


// this function should work out the derivatives.
void sam_ippm_mix_logl_derivs( const sam_ippm_data &dat, const sam_ippm_params &parms, sam_ippm_derivs &derivs, sam_ippm_fits &fits)
{
	vector<double> logPis(dat.nG, dat.NAnum), pis( dat.nG, dat.NAnum);
	vector<double> logCondDens( dat.nG, dat.NAnum);
	double wi, tmp1;
	vector<double> wij( dat.nG, dat.NAnum);
	int m;	//location of the maximum group contribution

	vector<double> muDerivsI( dat.nG*dat.nS, dat.NAnum);
	vector<double> etaDerivsI( dat.nG*dat.nS, dat.NAnum);
	vector<double> alphaDerivsI( dat.nS, dat.NAnum);
	vector<double> tauDerivsI( (dat.nG-1)*dat.nS, dat.NAnum);
	vector<double> piDerivsI( dat.nG, dat.NAnum);
	vector<double> betaDerivsI( (dat.nG-1)*dat.np, dat.NAnum);
	//vector<double> gammaDerivsI( dat.nS*dat.npw, dat.NAnum);
	//vector<double> dispDerivsI( dat.nS, dat.NAnum);

	vector<double> tmpPiDerivs( dat.nG*dat.nS, 0.0);

	calc_mu_fits( fits.allMus, dat, parms);
	derivs.zeroDerivs( dat);
	for( int i=0; i<dat.nObs; i++){
		//calculating the w_i and the {w_ig}
		calcLogPis( logPis, pis, dat, parms, i);
		calcLogCondDens( logCondDens, fits.allMus, dat, parms, i);
		tmp1 = calcMixSum( logPis, logCondDens, wi, wij, m);
		
		//calc deriv w.r.t. mu and the eta (all of them)
		calcDerivMu( muDerivsI, fits.allMus, dat, parms, wi, wij, m, i);
		calcDerivEtaMu( etaDerivsI, dat, muDerivsI, fits.allMus, i);
		
		//calc deriv w.r.t. alpha and then tau and then gamma
		calcAlphaDeriv( alphaDerivsI, etaDerivsI, dat);
		calcTauDeriv( tauDerivsI, etaDerivsI, dat, parms);
		
		//calc deriv w.r.t. beta
		calcBetaDeriv( betaDerivsI, piDerivsI, pis, dat, i);
		
		//update derivatives
		derivs.updateDerivs( dat, alphaDerivsI, tauDerivsI, betaDerivsI, gammaDerivsI, dispDerivsI, i);	//if dat.disty specifies no dispersion then no place for disp derivs (and are not updated, of course)
	}
	calcTauPenDeriv( tauDerivsI, dat, parms);
	calcGammaPenDeriv( gammaDerivsI, dat, parms);
	calcDispPenDeriv( dispDerivsI, dat, parms);

	alphaDerivsI.assign(alphaDerivsI.size(), 0.0);
	betaDerivsI.assign(betaDerivsI.size(), 0.0);

	derivs.updateDerivs( dat, alphaDerivsI, tauDerivsI, betaDerivsI, gammaDerivsI, dispDerivsI, -1);
	(void)tmp1;	//tmp1 is not used again, this is just a little trick to avoid a compile warning.
}

void additive_logistic(vector< double > &x,int inv){
  int i; 
  // inv == 1 gives transfornmation
  //inv = 0 give inverse transformatio
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

// this should calculate the derivate w.r.t alpha.

void calcAlphaDeriv( vector<double> &alphaDerivsI, const sam_ippm_data &dat, const double &mu){
	
	
	
	
	
	
	
	}

// this should calculate the derivate w.r.t beta.

void calcBetaDeriv( vector<double> &betaDerivsI, const sam_ippm_data &dat, const double &mu){
	
	
	
	
	
	}

// this should calculate the derivate w.r.t tau.

void calcTauDeriv( vector<double> &tauDerivsI, const vector<double> &etaDerivs, const sam_ippm_data &dat, const sam_ippm_data &parms){
	
	
	
	
	
	}

//void gradient_ippm(int n, double *pars, double *gr, void *ex ){
  //Optimise_data_ippm *data = (Optimise_data_ippm *) ex;
  //int Xr,Xc,i,j;
  //double d1;
  //Xr=data->Xr;
  //Xc=data->Xc;
  ////pars[0]=1;
  //for(j=0;j<n;j++) gr[j]=0;

  
  //// need to setup wts based on MATREF2D so w[MATREF2D(i,(j-1),Xr] should do the trick...
  //// need to setup y_is_not_na based on MATREF2D so y_is_not_na[MATREF2D(i,(j-1),Xr] should also do the trick...
  
  ////This is estimates the gradient for betas.  
  //for(i=0;i<Xr;i++){
	  
    //d1=data->y[i]/data->lp.at(i)-(pars[0]+data->y[i])/(data->lp.at(i)+pars[0]); //But where is the alphas? maybe pars[0]
    //for(j=1;j<=Xc;j++){
      //gr[j]+=(d1*data->lp.at(i)*data->X[MATREF2D(i,(j-1),Xr)])*data->w[i];

    //}
    //// gradient for theta
    //gr[0] += (digamma(pars[0]+data->y[i]) - digamma(pars[0]) + log(pars[0]) + 1 - log( data->lp.at(i) + pars[0]) - (pars[0]+data->y[i])/(data->lp.at(i) + pars[0]))*data->w[i];
  //}
  ////gr[0]=0;
  //for(i=0;i<n;i++) gr[i]=0-gr[i];

//}


//void gradient_mix_ippm_function(int n, double *pars, double *gr, void *ex ){
  //Optimise_data_ippm *data = (Optimise_data_ippm *) ex;
  ////Optimise_data data = *(Optimise_data *) ex;
  //vector<double> ad_g(n,0);
  //vector<double> x (n,0);
  //int i,G,g,j,Xc,s,S;
  //double add_log_trans=0;
  //vector<double> logl(1,0);

  //G=data->G;
  //S=data->S;
  //vector<double> pi_mat_deriv(G*(G-1),0);
  //vector<double> dl_dpi(G,0);

  //Xc=data->Xc;
  //for(i=0;i<n;i++) x.at(i) = pars[i];
 
  //logl = calc_mix_ippm_logl(x,*data);

  //for(g=0;g<(G-1);g++){ 
      //add_log_trans+=exp(x.at(g)); //add up transformed pi's
  //}
  //add_log_trans+=1;
 
  //for(g=0;g<G;g++){ //GO through all pi's

   //for(j=0;j<Xc;j++){
	   
	//for(s=0;s<S;s++){
	  //ad_g.at((G-1)+ MATREF2D(g,j,G)) += exp( -1*data->species_l_contrib.at(s) + log(data->parpi.at(g)) + data->sum_f_species.at(MATREF2D(g,s,G))) *data->deriv_f_B.at(MATREF3D(g,j,s,G,Xc));
	  //if(j==0) dl_dpi.at(g) += exp(  -1*data->species_l_contrib.at(s)+ data->sum_f_species.at(MATREF2D(g,s,G)));
	//}
      //}

      //for(i=0;i<(G-1);i++){ // go through eta's
	//if(g<(G-1)){
	  //if(i==g){
	    //pi_mat_deriv.at(MATREF2D(i,g,(G-1))) = exp(x.at(i))/add_log_trans - exp(2*x.at(i))/(add_log_trans*add_log_trans);// diag
	    //pi_mat_deriv.at(MATREF2D(i,(G-1),(G-1))) += pi_mat_deriv.at(MATREF2D(i,g,(G-1)));
	  //}else{
	    //pi_mat_deriv.at(MATREF2D(i,g,(G-1))) = -exp(x.at(i))*exp(x.at(g)) / (add_log_trans*add_log_trans); //off-diag
	    //pi_mat_deriv.at(MATREF2D(i,(G-1),(G-1))) += pi_mat_deriv.at(MATREF2D(i,g,(G-1)));
	  //}
	//}
      //}

      //for(s=0;s<S;s++){
      ////calculate for alphas
	//ad_g.at((G*Xc + G-1)+s)+= exp(data->sum_f_species.at(MATREF2D(g,s,G)) - data->species_l_contrib.at(s)+ log(data->parpi.at(g)))* data->deriv_f_alphaS.at(MATREF2D(g,s,G));
	
	  //// remove the dispersion parameters. We don't need this for ippm
      //////calculate for dispersion
	////ad_g.at((G*Xc + G-1 + S)+s)+= exp(data->sum_f_species.at(MATREF2D(g,s,G)) - data->species_l_contrib.at(s)+ log(data->parpi.at(g)))* data->deriv_f_dispersionS.at(MATREF2D(g,s,G));
      //}
  //}
  //for(i=0;i<(G-1);i++) pi_mat_deriv.at(MATREF2D(i,(G-1),(G-1))) *= -1;



    //for(i=0;i<(G-1);i++){
      //for(g=0;g<G;g++){
	//ad_g.at(i) += dl_dpi.at(g)* pi_mat_deriv.at(MATREF2D(i,g,(G-1)));
      //}
    //}

  //for(i=0;i<n;i++) gr[i] = 0-ad_g.at(i); 

//}

// going to have to fix this up to exclude NA data...

 //vector <double> calc_mix_ippm_logl(const vector<double> &pars, Optimise_data &data){
  //int G,S,Xr,Xc;
  //G=data.G;
  //S=data.S;
  //Xr=data.Xr;
  //Xc=data.Xc;
  //vector< double > estpi(G-1,0); //vector to hold pi's
  //vector< double > coef(Xc*G,0); //vector to hold all coefficents *DOES NOT INCLUDE INTERCEPT*
  //vector< double > sp_int(S,0); // vector for species intercepts
  ////vector< double > sp_dispersion(S,0); //vector for species specific dispersion parameter
  //vector< double > logl(1 ,0), tlog(1 ,0); //output log likelihood
  ////  double sumlogl=0;
  //int i,s,g,j;
  //int start, end;
  //double tmp;
  
  
  //for(i=0;i<(G-1);i++){
    //estpi.at(i) = pars.at(i); //G-1 values for pi
  //}

  //additive_logistic(estpi,1); // additive logistic transfor on pi;
  ////for(i=0;i<G;i++) Rprintf("%f,",estpi.at(i));
  ////Rprintf("\n");
  //for(i=0;i<G;i++){
    //data.parpi.at(i) = estpi.at(i);
    //for(s=0;s<S;s++){ 
      //data.sum_f_species.at(MATREF2D(i,s,G))=0;
      ////  data.deriv_f_dispersionS.at(MATREF2D(i,s,G))=0;
      //data.deriv_f_alphaS.at(MATREF2D(i,s,G))=0;
      //for(j=0;j<Xc;j++) data.deriv_f_B.at(MATREF3D(i,j,s,G,Xc)) =0;
    //}
  //}

 

  //for(i=(G-1);i<(G*Xc + G-1);i++){coef.at(i-(G-1))= pars.at(i);
    ////Rprintf("%f,",coef.at(i-(G-1)));
  //}
  //for(i= (G*Xc + G-1); i<(G*Xc + G-1 + S) ; i++) {sp_int.at(i-(G*Xc + G-1)) = pars.at(i);
    //// Rprintf("%f,",sp_int.at(i-(G*Xc + G-1)) );
//}

  ////Rprintf("\n");
  //// for(g=0;g<G;g++){
    //for(j=0;j<Xc;j++)
      //Rprintf("%f, ",coef[MATREF2D(g,j,G)]);
    //Rprintf("| %d\n",g);
    //}*/
  //logl.at(0) = 0;

  //for(s=0;s<S;s++){
      ////start = data.StartEndPos.at(s*2); // y is now a matrix so this is obsolte.
      ////end = data.StartEndPos.at(s*2+1);

      //tlog.at(0) = like_mix_ippm_function(estpi, coef, sp_int, data.y, data.X, Xr, Xc, data.tau, s, data.sum_f_species, data.deriv_f_B, data.deriv_f_alphaS, data.log_y_factorial, data.offset, data.y_is_not_na);
      //logl.at(0)+= tlog.at(0);
  
      //// this is the species loglike contribution. 
      //data.species_l_contrib.at(s) = tlog.at(0);
  //}
  //return(logl);

//}

// log-ippm-derivative - this will be used for dfdlogalpha and dfdlogbeta.
double log_df_dalpha_ippm( const double &y, const double &mu, const double &wts){
	double tmp, z;
	z = y/wts;
	tmp = z/mu;
	tmp -= 1;
	tmp *= wts;
	return( tmp);
}

//double like_mix_ippm_function(vector< double > &estpi, vector < double > &coef, vector < double > &sp_int, const double *y, const double *X, const double *wts, int Xr, int Xc, int start, int end, double *tau, int s, vector<double> &sum_f_species, vector<double> &deriv_f_B, vector<double> &deriv_f_alphaS, vector<double> &deriv_f_dispersionS, vector<double> &log_y_factorial, const double *offset, const int *y_is_not_na){
  //int len,i,j,G,g;
  //len=end-start+1;
  ////vector< AD<double > > p(len,0);
  //vector< double  > p(1,0);
  //double lpre=0, eps=0,glogl=0, d1=0;
  //G = estpi.size();
  //vector< double  > sump(G,0);

 
  //for(g=0;g<G;g++){
    //for(i=0;i<len;i++){
	  //if(y_is_not_na[MATREF2D(i,s,Xr)]>0){ //hopefully this will result in only selecting sites which have data. 
      //lpre=offset[i]+sp_int.at(s);
      //for(j=0;j<Xc;j++){ 
		//lpre +=  X[MATREF2D(i,j,Xr)] * coef[MATREF2D(g,j,G)];
      //}
    
      //p.at(0)=exp(lpre);  // mu for each archetype at site k
      //// datalp.at(i) = p.at(0);
      ////Rprintf("%f,",lpre);
      //dl = log_ippm_derivative(y[MATREF2D(i,s,Xr)], p.at(0), wts[MATREF2D(i,s,Xr)]);
      ////d1=y[start+i]/p.at(0)-(y[start+i])/(p.at(0));
      //for(j=0;j<Xc;j++) {deriv_f_B.at(MATREF3D(g,j,s,G,Xc)) +=  (d1*p.at(0)*X[MATREF2D(i,j,Xr)]);} 
      //deriv_f_alphaS.at(MATREF2D(g,s,G))+= (d1*p.at(0)*1);

      ////if(y[start+i]==0) p.at(0) = 1-p.at(0);

      //sump.at(g) += (lgammafn(y[MATREF2D(i,s,Xr)]) + y[MATREF2D(i,s,Xr)]*log(p.at(0)) - (y[MATREF2D(i,s,Xr)])*log(p.at(0)));
	  //}
    //}
     
 
    //for(j=0;j<Xc;j++)
	//Rprintf("%f, ",coef[MATREF2D(g,j,G)]);
      //Rprintf("\n");
    //*/

    //if(g==0) eps = sump.at(g);
    //if(sump.at(g) > eps) eps = sump.at(g);

    
  //}

  //for(g=0;g<G;g++){
    //sum_f_species.at(MATREF2D(g,s,G)) = sump.at(g);
    //glogl+= estpi.at(g)*exp(sump.at(g) - eps);
  //}
  //// code for taus can go here
  //glogl = log(glogl) + eps;
 
  //return(glogl);

//}
