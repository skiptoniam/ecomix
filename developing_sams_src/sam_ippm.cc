#include"sam.h"

/* Code for inhomogenous poisson point process model.
 * This should be very similar to the negative binomial with y_is_na added in and z and setups data classes for negative binomial.   
 *
 * I have tried to set this up like RCP, which makes more sense to me, in the future I can adapt Piers code into a single species_mix_cpp function.
 */

extern "C" { SEXP species_mix_ippm_cpp(SEXP Ry, SEXP RX, SEXP Roffset, SEXP Rwts, SEXP Ry_not_na,
								  SEXP RS, SEXP RG, SEXP Rp, SEXP RnObs,
								  SEXP Ralpha, SEXP Rbeta, SEXP Rtau, 
								  SEXP RderivsAlpha, SEXP RderivsBeta, SEXP RderivsTau, SEXP Rscores,
								  SEXP Rpis, SEXP Rmus, SEXP RlogDens, SEXP Rloglike,
								  SEXP Rmaxit, SEXP Rtrace, SEXP RnReport, SEXP Rabstol, SEXP Rreltol, SEXP Rconv,
								  SEXP Roptimise, SEXP RloglOnly, SEXP RderivsOnly, SEXP RoptiDisp, SEXP RgetScores){
	
	sam_ippm_all_classes all;

	//initialise the data structures -- they are mostly just pointers to REAL()s...
	all.data.setVals(Ry, RX, Roffset, RS, RG, Rp, RnObs, Rwts);	//read in the data
	all.parms.setVals(all.data, Ralpha, Rbeta, Rtau, Rgamma, Rdisps, Rpowers, Rconc, Rsd, RsdGamma, RdispLocat, RdispScale);	//read in the parameters
	all.derivs.setVals(all.data, RderivsAlpha, RderivsTau, RderivsBeta, RderivsGamma, RderivsDisps, RgetScores, Rscores);
	all.contr.setVals( Rmaxit, Rtrace, RnReport, Rabstol, Rreltol, Rconv);
	all.fits.initialise(all.data.nObs, all.data.nG, all.data.nS, all.data.NAnum);

	double logl = -999999;
	
	//doing the optimising
	if( *INTEGER(Roptimise) == 1)
		logl = ALLoptimise( all);

	//re-running to get pis and mus
	if( *INTEGER(RloglOnly) == 1)
		logl = mixLogl( all.data, all.parms, all.fits);
	//and derivatives (inlcuding scores, for empirical info, if requested)
	if( *INTEGER(RderivsOnly) == 1)
		loglDerivs( all.data, all.parms, all.derivs, all.fits);								  
									  
									  
									  }	
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


extern "C" {  SEXP species_mix_ippm_cpp(SEXP R_pars, SEXP R_X, SEXP R_y, SEXP R_w, SEXP R_offset, SEXP R_y_is_not_na,
										SEXP R_gradient, SEXP R_fitted_values,
										SEXP Rmaxit, SEXP Rtrace, SEXP RnReport, SEXP Rabstol, SEXP Rreltol, SEXP Rconv){
											
	allClasses all;
	
	all.contr.setVals( Rmaxit, Rtrace, RnReport, Rabstol, Rreltol, Rconv);										
	
											
    // tau is the parameter vector, 
    // od is the overdispersion parameter, 
    // X is the design matrix, 
    // N is the NB distributed variable
    int *y_is_not_na=NULL, Xr, Xc;
    double *pars=NULL, *y=NULL, *X=NULL, *w=NULL,*offset=NULL,*r_pointer=NULL, *gradient=NULL, *fitted_values=NULL;
    double logl;
    int lpar,i;
    double abstol,reltol;
    int fncount, grcount, ifail,trace,nREPORT;
    SEXP dimensions, R_logl;
    lpar = LENGTH(R_pars);

    //    vector<double> prams( lpar );
    vector< double > params(lpar);
    pars=REAL(R_pars);
    y=REAL(R_y);
    X=REAL(R_X);
    w=REAL(R_w);
    //ID=INTEGER(R_ID);
    offset=REAL(R_offset);
    y_is_not_na=INTEGER(R_y_is_not_na);
    gradient=REAL(R_gradient);
    fitted_values=REAL(R_fitted_values);
  
    dimensions=getAttrib(R_X, R_DimSymbol);
    Xr=INTEGER(dimensions)[0];
    Xc=INTEGER(dimensions)[1];

    dimensions=getAttrib(R_tau, R_DimSymbol);
    S=INTEGER(dimensions)[0];
    G=INTEGER(dimensions)[1];

    Optimise_data_ippm data;
    
    data.set_vars_ippm(y, X, w, *offset, y_is_not_na, ID, S, G, Xr, Xc, lpar, ly, tau)
    //data.SetVars(y, X, w, *offset, y_is_not_na, Xr, Xc);


    //all.contr.setVals( Rmaxit, Rtrace, RnReport, Rabstol, Rreltol, Rconv);
    abstol = 1e-8;
    reltol = 1e-8;
    nREPORT = 50;//how often to report
    fncount=0;
    grcount=0;
    ifail=0;
    trace=0;
    vector <int> mask (lpar,1); 
    vector < double > logl_out(lpar,0);

    vmmin(lpar, pars, &logl_out.at(0), optimise_ippm, gradient_ippm, 1000, trace, &mask.at(0),  abstol,  reltol,  nREPORT, &data, &fncount, &grcount, &ifail);
    
    //std::cout << "ifail = " << ifail << "\n" ;
    logl=1;
    logl = optimise_ippm(lpar, pars, &data);
    gradient_ippm(lpar,pars,gradient,&data);
    for(i=0;i<Xr;i++) fitted_values[i] = data.lp.at(i);

    R_logl = allocVector(REALSXP,1);
    r_pointer = REAL(R_logl);
    *r_pointer = logl;
    return(R_logl);
  }
}
extern "C" { SEXP species_mix_ippm_gradient_cpp(SEXP R_pars, SEXP R_X, SEXP R_y, SEXP R_w, SEXP R_offset, SEXP R_y_is_not_na, SEXP R_gradient){
    // tau is the parameter vector, 
    // od is the overdispersion parameter, 
    // X is the design matrix, 
    // N is the NB distributed variable
    int *y_is_not_na=NULL, Xr, Xc;
    double *pars=NULL, *y=NULL, *X=NULL, *w=NULL ,*offset=NULL,*r_pointer=NULL, *gradient=NULL;
    double logl;
    int lpar;
    double abstol,reltol;
    int fncount, grcount, ifail,trace,nREPORT;
    SEXP dimensions, R_logl;
    lpar = LENGTH(R_pars);

    //    vector<double> prams( lpar );
    vector< double > params(lpar);
    pars=REAL(R_pars);
    y=REAL(R_y);
    X=REAL(R_X);
    w=REAL(R_w);
    offset=REAL(R_offset);
    y_is_not_na = INTEGER(R_y_is_not_na);
    gradient=REAL(R_gradient);
  
    dimensions=getAttrib(R_X, R_DimSymbol);
    Xr=INTEGER(dimensions)[0];
    Xc=INTEGER(dimensions)[1];

    Optimise_data_ippm data;
    data.SetVars(y, X, w, *offset, y_is_not_na, Xr, Xc);

    abstol = 1e-8;
    reltol = 1e-8;
    nREPORT = 5;//how often to report
    fncount=0;
    grcount=0;
    ifail=0;
    trace=1;
    vector <int> mask (lpar,1); 
    vector < double > logl_out(lpar,0);

      
    logl=1;
    logl = optimise_ippm(lpar, pars, &data);
    gradient_ippm(lpar,pars,gradient,&data);

    R_logl = allocVector(REALSXP,1);
    r_pointer = REAL(R_logl);
    *r_pointer = logl;
    return(R_logl);
  }
}


void gradient_ippm(int n, double *pars, double *gr, void *ex ){
  Optimise_data_ippm *data = (Optimise_data_ippm *) ex;
  int Xr,Xc,i,j;
  double d1;
  Xr=data->Xr;
  Xc=data->Xc;
  //pars[0]=1;
  for(j=0;j<n;j++) gr[j]=0;

  
  // need to setup wts based on MATREF2D so w[MATREF2D(i,(j-1),Xr] should do the trick...
  // need to setup y_is_not_na based on MATREF2D so y_is_not_na[MATREF2D(i,(j-1),Xr] should also do the trick...
  
  //This is estimates the gradient for betas.  
  for(i=0;i<Xr;i++){
	  
    d1=data->y[i]/data->lp.at(i)-(pars[0]+data->y[i])/(data->lp.at(i)+pars[0]); //But where is the alphas? maybe pars[0]
    for(j=1;j<=Xc;j++){
      gr[j]+=(d1*data->lp.at(i)*data->X[MATREF2D(i,(j-1),Xr)])*data->w[i];

    }
    // gradient for theta
    gr[0] += (digamma(pars[0]+data->y[i]) - digamma(pars[0]) + log(pars[0]) + 1 - log( data->lp.at(i) + pars[0]) - (pars[0]+data->y[i])/(data->lp.at(i) + pars[0]))*data->w[i];
  }
  //gr[0]=0;
  for(i=0;i<n;i++) gr[i]=0-gr[i];

}


double ippm_logl(vector<double> &pars, Optimise_data_ippm &data ){//vector<double> &tau, double theta, vector<double> &X, vector<double> &y, double offset){
  vector<double> lp(data.Xr,0);
  int Xr,Xc,i,j;
  double offset=0,logl=0;
  //int *y_is_not_na=NULL;
  Xr=data.Xr;
  Xc=data.Xc;
  offset=data.offset;
  //y_is_not_na = data.y_is_not_na;
  
  for(i=0;i<Xr;i++){
	 for(j=0;j<Xc;j++){
	   if(data.y_is_not_na[MATREF2D(i,j,Xr)]>0){
	   lp.at(i)+=data.X[MATREF2D(i,j,Xr)]*pars.at(j+1);
	   lp.at(i)+=offset;
	   }
    }
    lp.at(i)=exp(lp.at(i));
    data.lp.at(i)=lp.at(i); //fitted values
   
   //loglikespp sum(first_fit$weights[sp_idx,ss]*(((first_fit$y[sp_idx,ss]/first_fit$weights[sp_idx,ss])*eta) -exp(eta)))
   
    logl+=(lgammafn(pars.at(0)+data.y[i]) - lgammafn(pars.at(0)) + data.y[i]*log(lp.at(i)) + pars.at(0)*log(pars.at(0)) - (pars.at(0) + data.y[i])*log(lp.at(i)+pars.at(0)) - data.log_y_factorial.at(i))*(data.w[i]) ;
  }
  return( logl);

} 
 //Optimise_data_ippm::Optimise_data_ippm(){}
//Optimise_data_ippm::~Optimise_data_ippm(){}

//void Optimise_data_ippm::SetVars(double *ty, double *tX, double *tw, double toffset, int *ty_is_not_na, int tXr, int tXc){
  //int i,j,yt;
  //double tmp;
  //y=ty;
  //X=tX;
  //Xr=tXr;
  //Xc=tXc;
  //w=tw;
  //offset=toffset;
  //y_is_not_na=ty_is_not_na;
  //for(i=0;i<Xr;i++){
    //lp.push_back(0); //fitted values
    //tmp=0;
    //yt = (int)y[i];
    //for(j=1;j<=yt;j++){tmp+=log((double) j); }
    //log_y_factorial.push_back(tmp); 
  //}
//}


Optimise_data_ippm::Optimise_data_ippm(){}
Optimise_data_ippm::~Optimise_data_ippm(){}

void Optimise_data_ippm::set_vars_ippm(double *ty, double *tX, double *tw, double *toffset, int *ty_is_not_na, int *tID, int tS, int tG, int tXr, int tXc, int tlpar, int tly, double *ttau){
  
  int i, s, j;
  double tmp,yt;
  tau = ttau;
  y=ty;
  X=tX;
  ly = tly;
  w=tw;
  offset=toffset;
  y_is_not_na=ty_is_not_na;
  ID=tID;
  S=tS;
  G=tG;
  Xr=tXr;
  Xc=tXc;
  
  lpar=tlpar;
  ly = tly;
 
  //s=1;
  for(i=0;i<G;i++){ // vectors of length G
    parpi.push_back(0);
    for(s=0;s<S;s++){ 
      sum_f_species.push_back(0);
      species_group_l_contrib.push_back(0);
      deriv_f_alphaS.push_back(0);
      for(j=0;j<Xc;j++) deriv_f_B.push_back(0);
    }
  }

  for(s=0;s<S;s++) {
    species_l_contrib.push_back(0);
  }

  s=ID[0];
 
}
 

double optimise_ippm(int n, double *pars, void *ex){

  Optimise_data_ippm *data = (Optimise_data_ippm *) ex;
  //Optimise_data data =  * (Optimise_data *) ex;
  
  vector<double> logl(1,0);
  int i;
  vector<double> x (n,0);
  //pars[0]=1;
  
  for(i=0;i<n;i++) x.at(i) = pars[i];


  logl.at(0) = ippm_logl(x,*data);
 
  //logl = data->F.Forward(0,x);

  return(0-logl.at(0));

}

// Not sure what the difference between these two are? Anyway - I've followed Piers' wacky method.

double optimise_mix_ippm_function(int n, double *pars, void *ex){

  Optimise_data_ippm *data = (Optimise_data_ippm *) ex;
  //Optimise_data data =  * (Optimise_data *) ex;
  
  vector<double> logl(1,0);
  int i;
  vector<double> x (n,0);
  
  for(i=0;i<n;i++) x.at(i) = pars[i];

  logl = calc_mix_ippm_logl(x,*data);

  //logl = data->F.Forward(0,x);

  return(0-logl.at(0));

}

void gradient_mix_ippm_function(int n, double *pars, double *gr, void *ex ){
  Optimise_data_ippm *data = (Optimise_data_ippm *) ex;
  //Optimise_data data = *(Optimise_data *) ex;
  vector<double> ad_g(n,0);
  vector<double> x (n,0);
  int i,G,g,j,Xc,s,S;
  double add_log_trans=0;
  vector<double> logl(1,0);

  G=data->G;
  S=data->S;
  vector<double> pi_mat_deriv(G*(G-1),0);
  vector<double> dl_dpi(G,0);

  Xc=data->Xc;
  for(i=0;i<n;i++) x.at(i) = pars[i];
 
  logl = calc_mix_ippm_logl(x,*data);

  for(g=0;g<(G-1);g++){ 
      add_log_trans+=exp(x.at(g)); //add up transformed pi's
  }
  add_log_trans+=1;
 
  for(g=0;g<G;g++){ //GO through all pi's

   for(j=0;j<Xc;j++){
	   
	for(s=0;s<S;s++){
	  ad_g.at((G-1)+ MATREF2D(g,j,G)) += exp( -1*data->species_l_contrib.at(s) + log(data->parpi.at(g)) + data->sum_f_species.at(MATREF2D(g,s,G))) *data->deriv_f_B.at(MATREF3D(g,j,s,G,Xc));
	  if(j==0) dl_dpi.at(g) += exp(  -1*data->species_l_contrib.at(s)+ data->sum_f_species.at(MATREF2D(g,s,G)));
	}
      }

      for(i=0;i<(G-1);i++){ // go through eta's
	if(g<(G-1)){
	  if(i==g){
	    pi_mat_deriv.at(MATREF2D(i,g,(G-1))) = exp(x.at(i))/add_log_trans - exp(2*x.at(i))/(add_log_trans*add_log_trans);// diag
	    pi_mat_deriv.at(MATREF2D(i,(G-1),(G-1))) += pi_mat_deriv.at(MATREF2D(i,g,(G-1)));
	  }else{
	    pi_mat_deriv.at(MATREF2D(i,g,(G-1))) = -exp(x.at(i))*exp(x.at(g)) / (add_log_trans*add_log_trans); //off-diag
	    pi_mat_deriv.at(MATREF2D(i,(G-1),(G-1))) += pi_mat_deriv.at(MATREF2D(i,g,(G-1)));
	  }
	}
      }

      for(s=0;s<S;s++){
      //calculate for alphas
	ad_g.at((G*Xc + G-1)+s)+= exp(data->sum_f_species.at(MATREF2D(g,s,G)) - data->species_l_contrib.at(s)+ log(data->parpi.at(g)))* data->deriv_f_alphaS.at(MATREF2D(g,s,G));
	
	  // remove the dispersion parameters. We don't need this for ippm
      ////calculate for dispersion
	//ad_g.at((G*Xc + G-1 + S)+s)+= exp(data->sum_f_species.at(MATREF2D(g,s,G)) - data->species_l_contrib.at(s)+ log(data->parpi.at(g)))* data->deriv_f_dispersionS.at(MATREF2D(g,s,G));
      }
  }
  for(i=0;i<(G-1);i++) pi_mat_deriv.at(MATREF2D(i,(G-1),(G-1))) *= -1;



    for(i=0;i<(G-1);i++){
      for(g=0;g<G;g++){
	ad_g.at(i) += dl_dpi.at(g)* pi_mat_deriv.at(MATREF2D(i,g,(G-1)));
      }
    }

  for(i=0;i<n;i++) gr[i] = 0-ad_g.at(i); 

}

// going to have to fix this up to exclude NA data...

 vector <double> calc_mix_ippm_logl(const vector<double> &pars, Optimise_data &data){
  int G,S,Xr,Xc;
  G=data.G;
  S=data.S;
  Xr=data.Xr;
  Xc=data.Xc;
  vector< double > estpi(G-1,0); //vector to hold pi's
  vector< double > coef(Xc*G,0); //vector to hold all coefficents *DOES NOT INCLUDE INTERCEPT*
  vector< double > sp_int(S,0); // vector for species intercepts
  //vector< double > sp_dispersion(S,0); //vector for species specific dispersion parameter
  vector< double > logl(1 ,0), tlog(1 ,0); //output log likelihood
  //  double sumlogl=0;
  int i,s,g,j;
  int start, end;
  double tmp;
  
  
  for(i=0;i<(G-1);i++){
    estpi.at(i) = pars.at(i); //G-1 values for pi
  }

  additive_logistic(estpi,1); // additive logistic transfor on pi;
  //for(i=0;i<G;i++) Rprintf("%f,",estpi.at(i));
  //Rprintf("\n");
  for(i=0;i<G;i++){
    data.parpi.at(i) = estpi.at(i);
    for(s=0;s<S;s++){ 
      data.sum_f_species.at(MATREF2D(i,s,G))=0;
      //  data.deriv_f_dispersionS.at(MATREF2D(i,s,G))=0;
      data.deriv_f_alphaS.at(MATREF2D(i,s,G))=0;
      for(j=0;j<Xc;j++) data.deriv_f_B.at(MATREF3D(i,j,s,G,Xc)) =0;
    }
  }

 

  for(i=(G-1);i<(G*Xc + G-1);i++){coef.at(i-(G-1))= pars.at(i);
    //Rprintf("%f,",coef.at(i-(G-1)));
  }
  for(i= (G*Xc + G-1); i<(G*Xc + G-1 + S) ; i++) {sp_int.at(i-(G*Xc + G-1)) = pars.at(i);
    // Rprintf("%f,",sp_int.at(i-(G*Xc + G-1)) );
}

  //Rprintf("\n");
  /*  for(g=0;g<G;g++){
    for(j=0;j<Xc;j++)
      Rprintf("%f, ",coef[MATREF2D(g,j,G)]);
    Rprintf("| %d\n",g);
    }*/
  logl.at(0) = 0;

  for(s=0;s<S;s++){
      //start = data.StartEndPos.at(s*2); // y is now a matrix so this is obsolte.
      //end = data.StartEndPos.at(s*2+1);

      tlog.at(0) = like_mix_ippm_function(estpi, coef, sp_int, data.y, data.X, Xr, Xc, data.tau, s, data.sum_f_species, data.deriv_f_B, data.deriv_f_alphaS, data.log_y_factorial, data.offset, data.y_is_not_na);
      logl.at(0)+= tlog.at(0);
  
      // this is the species loglike contribution. 
      data.species_l_contrib.at(s) = tlog.at(0);
  }
  return(logl);

}

// log-ippm-derivative - this will be used for dfdlogalpha and dfdlogbeta.
double log_ippm_derivative( const double &y, const double &mu, const double &wts){
	double tmp, z;
	z = y/wts;
	tmp = z/mu;
	tmp -= 1;
	tmp *= wts;
	return( tmp);
}

double like_mix_ippm_function(vector< double > &estpi, vector < double > &coef, vector < double > &sp_int, const double *y, const double *X, const double *wts, int Xr, int Xc, int start, int end, double *tau, int s, vector<double> &sum_f_species, vector<double> &deriv_f_B, vector<double> &deriv_f_alphaS, vector<double> &deriv_f_dispersionS, vector<double> &log_y_factorial, const double *offset, const int *y_is_not_na){
  int len,i,j,G,g;
  len=end-start+1;
  //vector< AD<double > > p(len,0);
  vector< double  > p(1,0);
  double lpre=0, eps=0,glogl=0, d1=0;
  G = estpi.size();
  vector< double  > sump(G,0);

 
  for(g=0;g<G;g++){
    for(i=0;i<len;i++){
	  if(y_is_not_na[MATREF2D(i,s,Xr)]>0){ //hopefully this will result in only selecting sites which have data. 
      lpre=offset[i]+sp_int.at(s);
      for(j=0;j<Xc;j++){ 
		lpre +=  X[MATREF2D(i,j,Xr)] * coef[MATREF2D(g,j,G)];
      }
    
      p.at(0)=exp(lpre);  // mu for each archetype at site k
      // datalp.at(i) = p.at(0);
      //Rprintf("%f,",lpre);
      dl = log_ippm_derivative(y[MATREF2D(i,s,Xr)], p.at(0), wts[MATREF2D(i,s,Xr)]);
      //d1=y[start+i]/p.at(0)-(y[start+i])/(p.at(0));
      for(j=0;j<Xc;j++) {deriv_f_B.at(MATREF3D(g,j,s,G,Xc)) +=  (d1*p.at(0)*X[MATREF2D(i,j,Xr)]);} 
      deriv_f_alphaS.at(MATREF2D(g,s,G))+= (d1*p.at(0)*1);

      //if(y[start+i]==0) p.at(0) = 1-p.at(0);

      sump.at(g) += (lgammafn(y[MATREF2D(i,s,Xr)]) + y[MATREF2D(i,s,Xr)]*log(p.at(0)) - (y[MATREF2D(i,s,Xr)])*log(p.at(0)));
	  }
    }
     
 
    /*Rprintf("%d, %d, %f, %f, %f \n",s,g, sump.at(g), sp_int.at(s),sp_dispersion.at(s));
      for(j=0;j<Xc;j++)
	Rprintf("%f, ",coef[MATREF2D(g,j,G)]);
      Rprintf("\n");
    */

    if(g==0) eps = sump.at(g);
    if(sump.at(g) > eps) eps = sump.at(g);

    
  }

  for(g=0;g<G;g++){
    sum_f_species.at(MATREF2D(g,s,G)) = sump.at(g);
    glogl+= estpi.at(g)*exp(sump.at(g) - eps);
  }
  // code for taus can go here
  glogl = log(glogl) + eps;
 
  return(glogl);

}
