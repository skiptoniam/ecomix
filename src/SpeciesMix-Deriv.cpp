#include"SpeciesMix-Deriv.h"

extern "C" 
{ 
  SEXP SpeciesMix(SEXP R_pars, SEXP R_y, SEXP R_X, SEXP R_ID,SEXP R_tau, SEXP R_gradient, SEXP R_offset, SEXP R_model_type){

    // y is response
    // X is design matrix
    // ID is vector length(y) with species names for each observation
    // ID_names is vector length S of unique ID
    // tau is easy pass out of matrix of tau
    // estpi is easy pass out of pi
    // model_type selects model, 1 = binomial, 2 = negative binomial, 3 = tweedie, 4 = ipp
    // pars is pi, coef matrix (including intercept in the binomial case), species intercept (in case of neg bin, tweedie),  dispersion (in case of neg bin, tweedie)
    double *pars=NULL, *y=NULL, *X=NULL, *estpi=NULL, *tau=NULL, *pi=NULL, *r_pointer=NULL, *gradient=NULL, *offset=NULL;
    int *ID=NULL;
    double logl=0;
    int S,G, Xr, Xc, lpar,s,i ,ly;
    SEXP dimensions, R_logl=NULL;
    double abstol,reltol;
    int fncount, grcount, ifail,trace,nREPORT, model_type;
  
  
    Optimise_data data; // this is were everything is set up

    //R_pars has length = length(pi) + length( G* (number of covariates + intercept))
    // ordered so that first G elements are pi and remaining elements are coefficents. all coefficents from same group are together
    lpar = LENGTH(R_pars);
    ly = LENGTH(R_y);
    vector<double> prams( lpar );
    vector< double > params(lpar,0);
    pars=REAL(R_pars);
    y=REAL(R_y);
    X=REAL(R_X);
    tau=REAL(R_tau);
    ID=INTEGER(R_ID);
    gradient=REAL(R_gradient);
    model_type=*INTEGER(R_model_type);
    offset=REAL(R_offset);

 
    for(i=0;i<lpar;i++){ params.at(i) = pars[i]; }
    // Rprintf("\n");


    dimensions=getAttrib(R_X, R_DimSymbol);
    Xr=INTEGER(dimensions)[0];
    Xc=INTEGER(dimensions)[1];34

    dimensions=getAttrib(R_tau, R_DimSymbol);
    S=INTEGER(dimensions)[0];
    G=INTEGER(dimensions)[1];

    // this is were the data is added to vmin function. This is were I need to define the data classes - write a new class. 
    data.SetVars(y, X, ID, S, G,  Xr, Xc,  lpar, ly, tau, offset);
    //    Rprintf("%d,%d,%d,%d,%d\n",ly,Xr,Xc,S,G);
    vector< double > pll( 1, 0);
    vector< double > grd(lpar,0);
   
     // set up variables for optimisation
    abstol = 1e-8;
    reltol = 1e-8;
    nREPORT = 10;//how often to report
    fncount=0;
    grcount=0;
    ifail=0;
    trace=1;
    vector <int> mask (lpar,1); 
    vector < double > logl_out(lpar,0);
    
    
    if(model_type==1){
      vmmin(lpar, pars, &logl_out.at(0), optimise_function, gradient_function, 1000, trace, &mask.at(0),  abstol,  reltol,  nREPORT, &data, &fncount, &grcount, &ifail);
    }
    if(model_type==2){ // this is where the function data will be optimised and &data is where the all data comes from this will be from a new ppm class. 
      vmmin(lpar, pars, &logl_out.at(0), optimise_mixnbinom_function, gradient_mixnbinom_function, 1000, trace, &mask.at(0),  abstol,  reltol,  nREPORT, &data, &fncount, &grcount, &ifail);
    }
    if(model_type==3){
      //for(i=(G*Xc + G-1 + S) ; i< data.lpar; i++) mask.at(i)=0;
      //for(i=(G*Xc + G-1 + S+G) ; i< data.lpar; i++) mask.at(i)=0;
      vmmin(lpar, pars, &logl_out.at(0), optimise_tweedie_function, gradient_tweedie_function, 1000, trace, &mask.at(0),  abstol,  reltol,  nREPORT, &data, &fncount, &grcount, &ifail);
    }
     // will need to develop a optimise_ipp_function & gradient_ipp_function
     //if(model_type==4){
      //vmmin(lpar, pars, &logl_out.at(0), optimise_ipp_function, gradient_ipp_function, 1000, trace, &mask.at(0),  abstol,  reltol,  nREPORT, &data, &fncount, &grcount, &ifail);
    //}
    
     logl=1;

    if(model_type==1){
      logl = optimise_function(lpar, pars, &data);
      gradient_function(lpar, pars, &grd.at(0),&data);
    }


    if(model_type==2){
      logl = optimise_mixnbinom_function(lpar, pars, &data);
      gradient_mixnbinom_function(lpar, pars, &grd.at(0),&data);
    }
    if(model_type==3){
      logl = optimise_tweedie_function(lpar, pars, &data);
      gradient_tweedie_function(lpar, pars, &grd.at(0),&data);
    }
  
	// get the log-likelihood from this call.
    //if(model_type==4){
      //logl = optimise_ipp_function(lpar, pars, &data);
      //gradient_ipp_function(lpar, pars, &grd.at(0),&data);
    //}
  
    for(i=0;i<lpar;i++){
      gradient[i] = grd.at(i);
    }

    R_logl = allocVector(REALSXP,1);
    r_pointer = REAL(R_logl);
    *r_pointer = logl;
    return(R_logl);
   }
}

double optimise_function(int n, double *pars, void *ex){

  Optimise_data *data = (Optimise_data *) ex;
  //Optimise_data data =  * (Optimise_data *) ex;
  
  vector<double> logl(1,0);
  int i;
  vector<double> x (n,0);
  
  for(i=0;i<n;i++) x.at(i) = pars[i];

  logl = calc_logl(x,*data);

  //logl = data->F.Forward(0,x);

  return(0-logl.at(0));

}

double optimise_ipp_function(int n, double *pars, void *ex){

  Optimise_data *data = (Optimise_data *) ex;
  //Optimise_data data =  * (Optimise_data *) ex;
  
  vector<double> logl(1,0);
  int i;
  vector<double> x (n,0);
  
  for(i=0;i<n;i++) x.at(i) = pars[i];

  logl = calc_ipp_logl(x,*data);

  //logl = data->F.Forward(0,x);

  return(0-logl.at(0));

}

void gradient_function(int n, double *pars, double *gr, void *ex ){
  Optimise_data *data = (Optimise_data *) ex;
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
 
  logl = calc_logl(x,*data);

  for(g=0;g<(G-1);g++){ 
      add_log_trans+=exp(x.at(g)); //add up transformed pi's
  }
  add_log_trans+=1;
 
  for(g=0;g<G;g++){ //GO through all pi's

      for(j=0;j<Xc;j++){
	for(s=0;s<S;s++){
	  ad_g.at((G-1)+ MAT_RF(g,j,G)) += exp( -1*data->species_l_contrib.at(s) + log(data->parpi.at(g)) + data->sum_f_species.at(MAT_RF(g,s,G))) *data->deriv_f_B.at(MAT_3D(g,j,s,G,Xc));
	  if(j==0) dl_dpi.at(g) += exp(  -1*data->species_l_contrib.at(s)+ data->sum_f_species.at(MAT_RF(g,s,G)));
	}
      }

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
  for(i=0;i<(G-1);i++) pi_mat_deriv.at(MAT_RF(i,(G-1),(G-1))) *= -1;



    for(i=0;i<(G-1);i++){
      for(g=0;g<G;g++){
	ad_g.at(i) += dl_dpi.at(g)* pi_mat_deriv.at(MAT_RF(i,g,(G-1)));
      }
    }

  for(i=0;i<n;i++) gr[i] = 0-ad_g.at(i); 

}


 vector <double> calc_logl(const vector<double> &pars, Optimise_data &data){
  int G,S,Xr,Xc;
  G=data.G;
  S=data.S;
  Xr=data.Xr;
  Xc=data.Xc;
  vector< double > estpi(G-1,0); //vector to hold pi's
  vector< double > coef(Xc*G,0); //vector to hold all coefficents
  //vector< double > logl( S * G ,0); //output log likelihood
  vector< double > logl(1 ,0), tlog(1 ,0); //output log likelihood
  //  double sumlogl=0;
  int i,s,g,j;
  int start, end;
  double tmp;
  
  for(i=0;i<(G-1);i++){
    estpi.at(i) = pars.at(i); //G-1 values for pi
  }

  additive_logistic(estpi,1); // additive logistic transfor on pi;

  for(i=0;i<G;i++){
    data.parpi.at(i) = estpi.at(i);
    for(s=0;s<S;s++){ 
      data.sum_f_species.at(MAT_RF(i,s,G))=0;
      for(j=0;j<Xc;j++) data.deriv_f_B.at(MAT_3D(i,j,s,G,Xc)) =0;
    }
  }
 
   
  for(i=(G-1);i<data.lpar;i++){coef.at(i-(G-1))= pars.at(i);}

  logl.at(0) = 0;

  for(s=0;s<S;s++){
      start = data.StartEndPos.at(s*2);
      end = data.StartEndPos.at(s*2+1);

      tlog.at(0) = like_function(estpi, coef,data.y,data.X,Xr,Xc,start,end, data.tau, s ,data.sum_f_species,data.deriv_f_B);
      logl.at(0)+= tlog.at(0);
      data.species_l_contrib.at(s) = tlog.at(0);
  }
  return(logl);

}

 double like_function(vector< double > &estpi, vector < double > &coef, const double *y, const double *X, int Xr, int Xc, int start, int end, double *tau, int s, vector<double> &sum_f_species, vector<double> &deriv_f_B){
  int len,i,j,G,g;
  len=end-start+1;
  //vector< AD<double > > p(len,0);
  vector< double  > p(1,0);
  double lpre=0, eps=0,glogl=0;
  G = estpi.size();
  vector< double  > sump(G,0);
 
  for(g=0;g<G;g++){
    for(i=0;i<len;i++){
      lpre=0;
      for(j=0;j<Xc;j++){ 
	lpre +=  X[MAT_RF(i,j,Xr)] * coef[MAT_RF(g,j,G)];
      }
      p.at(0)=inv_link_function(lpre,0);  // mu for each archetype at site k
 
      for(j=0;j<Xc;j++) {deriv_f_B.at(MAT_3D(g,j,s,G,Xc)) +=   (y[start+i] - p.at(0)) * X[MAT_RF(i,j,Xr)] ;} 

      if(y[start+i]==0) p.at(0) = 1-p.at(0);

      sump.at(g) += log(p.at(0));
  
    }
    if(g==0) eps = sump.at(g);
    if(sump.at(g) > eps) eps = sump.at(g);
    
  }

  for(g=0;g<G;g++){
    sum_f_species.at(MAT_RF(g,s,G)) = sump.at(g);
    glogl+= estpi.at(g)*exp(sump.at(g) - eps);
  }
  // code for taus can go here
  glogl = log(glogl) + eps;
 
  return(glogl);

}




double link_function(double p, int link_type){
  // using logit link function
  if(link_type==0) return(log(p)-log(1-p));
  return(log(p));
}
double inv_link_function(double lpre, int link_type){
  // using logit link function
  if(link_type==0) return(exp(lpre)/(1+exp(lpre)));
  return(exp(lpre));
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

Optimise_data::Optimise_data(){}
Optimise_data::~Optimise_data(){}

// this is where the data is setup. 
void Optimise_data::SetVars(double *ty, double *tX, int *tID, int tS, int tG, int tXr, int tXc, int tlpar, int tly, double *ttau, double *toffset){
  int i, s, j;

  double tmp,yt;
  tau = ttau;
  y=ty;
  X=tX;
  ID=tID;
  S=tS;
  G=tG;
  Xr=tXr;
  Xc=tXc;
  lpar=tlpar;
  ly = tly;
  offset=toffset;
  //vector<int> StartEndPos(S*2);

  StartEndPos.push_back(0);
  //s=1;
  for(i=0;i<G;i++){ // vectors of length G
    parpi.push_back(0);
    for(s=0;s<S;s++){ 
      sum_f_species.push_back(0);
      species_group_l_contrib.push_back(0);
      deriv_f_alphaS.push_back(0);
      deriv_f_dispersionS.push_back(0);
      for(j=0;j<Xc;j++) deriv_f_B.push_back(0);
    }
  }

  for(s=0;s<S;s++) {
    species_l_contrib.push_back(0);
  
  }

  s=ID[0];

  for(i=0;i<ly;i++) // find start & end positions of each species data
    if(ID[i]!=s){
      StartEndPos.push_back(i-1);
      //s++; // index to next species
      //std::cout << s << "," << ID[i] << "," << i <<"\n" ; 
      s=ID[i];
      StartEndPos.push_back(i); // next species start pos
    }
  
  StartEndPos.push_back(ly-1);  //add final position


  for(i=0;i<ly;i++){
   
    tmp=0;
    yt = y[i];
    for(j=1;j<=yt;j++){tmp+=log((double) j); }
    log_y_factorial.push_back(tmp); 
  }
 
}


extern "C" 
{ 
  SEXP Calculate_Gradient(SEXP R_pars, SEXP R_y, SEXP R_X, SEXP R_ID,SEXP R_tau, SEXP R_gradient, SEXP R_offset, SEXP R_model_type){
    // y is response
    //X is design matrix
    // ID is vector length(y) with species names for each observation
    // ID_names is vector length S of unique ID
    // tau is easy pass out of matrix of tau
    // estpi is easy pass out of pi
    double *pars=NULL, *y=NULL, *X=NULL, *estpi=NULL, *tau=NULL, *pi=NULL, *r_pointer=NULL, *gradient=NULL,*offset=NULL;
    int *ID=NULL;
    double logl=0;
    int S,G, Xr, Xc, lpar,s,i ,ly;
    SEXP dimensions, R_logl;
    double abstol,reltol;
    int fncount, grcount, ifail,trace,nREPORT,model_type=0;
  
  
     Optimise_data data;

    //R_pars has length = length(pi) + length( G* (number of covariates + intercept))
    // ordered so that first G elements are pi and remaining elements are coefficents. all coefficents from same group are together
    lpar = LENGTH(R_pars);
    ly = LENGTH(R_y);
    //    vector<double> prams( lpar );
     vector< double > params(lpar);
    pars=REAL(R_pars);
    y=REAL(R_y);
    X=REAL(R_X);
    tau=REAL(R_tau);
    ID=INTEGER(R_ID);
    gradient=REAL(R_gradient);
     offset=REAL(R_offset);
    model_type=*INTEGER(R_model_type);

      for(i=0;i<lpar;i++){ params.at(i) = pars[i]; }



    dimensions=getAttrib(R_X, R_DimSymbol);
    Xr=INTEGER(dimensions)[0];
    Xc=INTEGER(dimensions)[1];

    dimensions=getAttrib(R_tau, R_DimSymbol);
    S=INTEGER(dimensions)[0];
    G=INTEGER(dimensions)[1];

    
     data.SetVars(y, X, ID, S, G,  Xr, Xc,  lpar, ly, tau, offset);
    //    Rprintf("%d,%d,%d,%d,%d\n",ly,Xr,Xc,S,G);
    //vector< double > pll( 1, 0);
    vector< double > grd(lpar,0);
   
     // set up variables for optimisation
    abstol = 1e-8;
    reltol = 1e-8;
    nREPORT = 5;//how often to report
    fncount=0;
    grcount=0;
    ifail=0;
    trace=1;
    // vector <int> mask (lpar,1); 
    // vector < double > logl_out(lpar,0);

    logl=1;
    if(model_type==1){
      logl = optimise_function(lpar, pars, &data);

    //std::cout << "Optimise = " << logl << "\n";
       gradient_function(lpar, pars, &grd.at(0),&data);

    }
    if(model_type==2){
      logl = optimise_mixnbinom_function(lpar, pars, &data);

    //std::cout << "Optimise = " << logl << "\n";
      gradient_mixnbinom_function(lpar, pars, &grd.at(0),&data);
    }
 if(model_type==3){
      logl = optimise_tweedie_function(lpar, pars, &data);

    //std::cout << "Optimise = " << logl << "\n";
      gradient_tweedie_function(lpar, pars, &grd.at(0),&data);
    }

 //if(model_type==4){
      //logl = optimise_ipp_function(lpar, pars, &data);

    ////std::cout << "Optimise = " << logl << "\n";
      //gradient_ipp_function(lpar, pars, &grd.at(0),&data);
    //}

    // std::cout << "Gradient  = ";
    for(i=0;i<lpar;i++){
      //  std::cout << i<< " : " << grd.at(i) << ", ";
       gradient[i] = grd.at(i);
    }

    // std::cout << "\n";


    R_logl = allocVector(REALSXP,1);
    r_pointer = REAL(R_logl);
    *r_pointer = logl;
    return(R_logl);
  }

}

extern "C" 
{
  SEXP Neg_Bin(SEXP R_pars, SEXP R_X, SEXP R_y, SEXP R_w, SEXP R_offset, SEXP R_gradient, SEXP R_fitted_values){
    // tau is the parameter vector, 
    // od is the overdispersion parameter, 
    // X is the design matrix, 
    // N is the NB distributed variable
    int Xr, Xc;
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
    offset=REAL(R_offset);
    gradient=REAL(R_gradient);
    fitted_values=REAL(R_fitted_values);
  


    dimensions=getAttrib(R_X, R_DimSymbol);
    Xr=INTEGER(dimensions)[0];
    Xc=INTEGER(dimensions)[1];

    Optimise_data_nbinom data;
    data.SetVars(y,X, w, *offset,Xr,Xc);

    abstol = 1e-8;
    reltol = 1e-8;
    nREPORT = 50;//how often to report
    fncount=0;
    grcount=0;
    ifail=0;
    trace=0;
    vector <int> mask (lpar,1); 
    vector < double > logl_out(lpar,0);

    vmmin(lpar, pars, &logl_out.at(0), optimise_nbinom, gradient_nbinom, 1000, trace, &mask.at(0),  abstol,  reltol,  nREPORT, &data, &fncount, &grcount, &ifail);
    
    //std::cout << "ifail = " << ifail << "\n" ;
    logl=1;
    logl = optimise_nbinom(lpar, pars, &data);
    gradient_nbinom(lpar,pars,gradient,&data);
    for(i=0;i<Xr;i++) fitted_values[i] = data.lp.at(i);

    R_logl = allocVector(REALSXP,1);
    r_pointer = REAL(R_logl);
    *r_pointer = logl;
    return(R_logl);
  }
}
extern "C" 
{
  SEXP Neg_Bin_Gradient(SEXP R_pars, SEXP R_X, SEXP R_y, SEXP R_w, SEXP R_offset, SEXP R_gradient){
    // tau is the parameter vector, 
    // od is the overdispersion parameter, 
    // X is the design matrix, 
    // N is the NB distributed variable
    int Xr, Xc;
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
    gradient=REAL(R_gradient);
  


    dimensions=getAttrib(R_X, R_DimSymbol);
    Xr=INTEGER(dimensions)[0];
    Xc=INTEGER(dimensions)[1];

    Optimise_data_nbinom data;
    data.SetVars(y,X,w,*offset,Xr,Xc);

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
    logl = optimise_nbinom(lpar, pars, &data);
    gradient_nbinom(lpar,pars,gradient,&data);

    R_logl = allocVector(REALSXP,1);
    r_pointer = REAL(R_logl);
    *r_pointer = logl;
    return(R_logl);
  }
}


double optimise_nbinom(int n, double *pars, void *ex){

  Optimise_data_nbinom *data = (Optimise_data_nbinom *) ex;
  //Optimise_data data =  * (Optimise_data *) ex;
  
  vector<double> logl(1,0);
  int i;
  vector<double> x (n,0);
  //pars[0]=1;
  
  for(i=0;i<n;i++) x.at(i) = pars[i];


  logl.at(0) = NBlogl(x,*data);
 
  //logl = data->F.Forward(0,x);

  return(0-logl.at(0));

}

void gradient_nbinom(int n, double *pars, double *gr, void *ex ){
  Optimise_data_nbinom *data = (Optimise_data_nbinom *) ex;
  int Xr,Xc,i,j;
  double d1;
  Xr=data->Xr;
  Xc=data->Xc;
  //pars[0]=1;
  for(j=0;j<n;j++) gr[j]=0;

  //gradient for betas

  for(i=0;i<Xr;i++){
    d1=data->y[i]/data->lp.at(i)-(pars[0]+data->y[i])/(data->lp.at(i)+pars[0]);
    for(j=1;j<=Xc;j++){
      gr[j]+=(d1*data->lp.at(i)*data->X[MAT_RF(i,(j-1),Xr)])*data->w[i];

    }
    // gradient for theta
    gr[0] += (digamma(pars[0]+data->y[i]) - digamma(pars[0]) + log(pars[0]) + 1 - log( data->lp.at(i) + pars[0]) - (pars[0]+data->y[i])/(data->lp.at(i) + pars[0]))*data->w[i];
  }
  //gr[0]=0;
  for(i=0;i<n;i++) gr[i]=0-gr[i];

}


double NBlogl(vector<double> &pars, Optimise_data_nbinom &data ){//vector<double> &tau, double theta, vector<double> &X, vector<double> &y, double offset){
  vector<double> lp(data.Xr,0);
  int Xr,Xc,i,j;
  double offset=0,logl=0;
  Xr=data.Xr;
  Xc=data.Xc;
  offset=data.offset;

  for(i=0;i<Xr;i++){
    lp.at(i)=offset;
    for(j=0;j<Xc;j++){
       lp.at(i)+=data.X[MAT_RF(i,j,Xr)]*pars.at(j+1);
    }
    lp.at(i)=exp(lp.at(i));
    data.lp.at(i)=lp.at(i); //fitted values

    logl+=(lgammafn(pars.at(0)+data.y[i]) - lgammafn(pars.at(0)) + data.y[i]*log(lp.at(i)) + pars.at(0)*log(pars.at(0)) - (pars.at(0) + data.y[i])*log(lp.at(i)+pars.at(0)) - data.log_y_factorial.at(i))*(data.w[i]) ;
  }
  return( logl);

}

double ipplogl(vector<double> &pars, Optimise_data_ipp &data ){//vector<double> &tau, double theta, vector<double> &X, vector<double> &y, double offset, vector<double> &weights, vector<int> &y_id){
  vector<double> lp(data.Xr,0);
  int Xr,Xc,i,j;
  double offset=0,logl=0;
  Xr=data.Xr;
  Xc=data.Xc;
  offset=data.offset;

  for(i=0;i<Xr;i++){
    lp.at(i)=offset;
    for(j=0;j<Xc;j++){
       lp.at(i)+=data.X[MAT_RF(i,j,Xr)]*pars.at(j+1);
    }
    lp.at(i)=exp(lp.at(i));
    data.lp.at(i)=lp.at(i); //fitted values

    logl+=(lgammafn(pars.at(0)+data.y[i]) - lgammafn(pars.at(0)) + data.y[i]*log(lp.at(i)) + pars.at(0)*log(pars.at(0)) - (pars.at(0) + data.y[i])*log(lp.at(i)+pars.at(0)) - data.log_y_factorial.at(i))*(data.w[i]) ;
  }
  return( logl);

}

Optimise_data_nbinom::Optimise_data_nbinom(){}
Optimise_data_nbinom::~Optimise_data_nbinom(){}

void Optimise_data_nbinom::SetVars(double *ty, double *tX, double *tw, double toffset, int tXr, int tXc){
  int i,j,yt;
  double tmp;
  y=ty;
  X=tX;
  Xr=tXr;
  Xc=tXc;
  w=tw;
  offset=toffset;
  for(i=0;i<Xr;i++){
    lp.push_back(0); //fitted values
    tmp=0;
    yt = (int)y[i];
    for(j=1;j<=yt;j++){tmp+=log((double) j); }
    log_y_factorial.push_back(tmp); 
  }
}


/*
 *
 * code for negative binomial mixture
 *
 *
 */

// this is wrapper for the all the other stuff that needs to be calculed. 

double optimise_mixnbinom_function(int n, double *pars, void *ex){

  Optimise_data *data = (Optimise_data *) ex;
  //Optimise_data data =  * (Optimise_data *) ex;
  
  vector<double> logl(1,0);
  int i;
  vector<double> x (n,0);
  
  for(i=0;i<n;i++) x.at(i) = pars[i];

  logl = calc_mixnbinom_logl(x,*data);

  //logl = data->F.Forward(0,x);

  return(0-logl.at(0));

}

void gradient_mix_ipp_function(int n, double *pars, double *gr, void *ex ){
  Optimise_data *data = (Optimise_data *) ex;
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
 
  logl = calc_mixnbinom_logl(x,*data);

  for(g=0;g<(G-1);g++){ 
      add_log_trans+=exp(x.at(g)); //add up transformed pi's
  }
  add_log_trans+=1;
 
  for(g=0;g<G;g++){ //GO through all pi's

      for(j=0;j<Xc;j++){
	for(s=0;s<S;s++){
	  ad_g.at((G-1)+ MAT_RF(g,j,G)) += exp( -1*data->species_l_contrib.at(s) + log(data->parpi.at(g)) + data->sum_f_species.at(MAT_RF(g,s,G))) *data->deriv_f_B.at(MAT_3D(g,j,s,G,Xc));
	  if(j==0) dl_dpi.at(g) += exp(  -1*data->species_l_contrib.at(s)+ data->sum_f_species.at(MAT_RF(g,s,G)));
	}
      }

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

      for(s=0;s<S;s++){
      //calculate for alphas
	ad_g.at((G*Xc + G-1)+s)+= exp(data->sum_f_species.at(MAT_RF(g,s,G)) - data->species_l_contrib.at(s)+ log(data->parpi.at(g)))* data->deriv_f_alphaS.at(MAT_RF(g,s,G));
	
      //calculate for dispersion
	ad_g.at((G*Xc + G-1 + S)+s)+= exp(data->sum_f_species.at(MAT_RF(g,s,G)) - data->species_l_contrib.at(s)+ log(data->parpi.at(g)))* data->deriv_f_dispersionS.at(MAT_RF(g,s,G));
      }
  }
  for(i=0;i<(G-1);i++) pi_mat_deriv.at(MAT_RF(i,(G-1),(G-1))) *= -1;



    for(i=0;i<(G-1);i++){
      for(g=0;g<G;g++){
	ad_g.at(i) += dl_dpi.at(g)* pi_mat_deriv.at(MAT_RF(i,g,(G-1)));
      }
    }

  for(i=0;i<n;i++) gr[i] = 0-ad_g.at(i); 

}

// this is main gradient function. This is the expression of the derivates what we have. 
void gradient_mixnbinom_function(int n, double *pars, double *gr, void *ex ){
  Optimise_data *data = (Optimise_data *) ex;
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
 
  logl = calc_mixnbinom_logl(x,*data);

  for(g=0;g<(G-1);g++){ 
      add_log_trans+=exp(x.at(g)); //add up transformed pi's
  }
  add_log_trans+=1;
 
  for(g=0;g<G;g++){ //GO through all pi's

    for(j=0;j<Xc;j++){
		for(s=0;s<S;s++){
			// calculate for betas. 
			ad_g.at((G-1)+ MAT_RF(g,j,G)) += exp( -1*data->species_l_contrib.at(s) + log(data->parpi.at(g)) + data->sum_f_species.at(MAT_RF(g,s,G))) *data->deriv_f_B.at(MAT_3D(g,j,s,G,Xc));
			if(j==0) dl_dpi.at(g) += exp(  -1*data->species_l_contrib.at(s)+ data->sum_f_species.at(MAT_RF(g,s,G)));
			}
    }

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

      for(s=0;s<S;s++){
      //calculate for alphas
	ad_g.at((G*Xc + G-1)+s)+= exp(data->sum_f_species.at(MAT_RF(g,s,G)) - data->species_l_contrib.at(s)+ log(data->parpi.at(g)))* data->deriv_f_alphaS.at(MAT_RF(g,s,G));
	
      //calculate for dispersion
	ad_g.at((G*Xc + G-1 + S)+s)+= exp(data->sum_f_species.at(MAT_RF(g,s,G)) - data->species_l_contrib.at(s)+ log(data->parpi.at(g)))* data->deriv_f_dispersionS.at(MAT_RF(g,s,G));
      }
  }
  for(i=0;i<(G-1);i++) pi_mat_deriv.at(MAT_RF(i,(G-1),(G-1))) *= -1;



    for(i=0;i<(G-1);i++){
      for(g=0;g<G;g++){
	ad_g.at(i) += dl_dpi.at(g)* pi_mat_deriv.at(MAT_RF(i,g,(G-1)));
      }
    }
  
  // last line is the negative sum.
  for(i=0;i<n;i++) gr[i] = 0-ad_g.at(i); 

}

//Calculate the ipp logl
 vector <double> calc_mix_ipp_logl(const vector<double> &pars, Optimise_data &data){
  int G,S,Xr,Xc;
  G=data.G; //groups
  S=data.S; //species
  Xr=data.Xr; // design matrix rows
  Xc=data.Xc; // design matrix cols
  vector< double > estpi(G-1,0); //vector to hold pi's
  vector< double > coef(Xc*G,0); //vector to hold all coefficents *DOES NOT INCLUDE INTERCEPT*
  vector< double > sp_int(S,0); // vector for species intercepts
  vector< double > sp_dispersion(S,0); //vector for species specific dispersion parameter
  //vector< double > logl( S * G ,0); //output log likelihood
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
      data.sum_f_species.at(MAT_RF(i,s,G))=0;
      data.deriv_f_dispersionS.at(MAT_RF(i,s,G))=0;
      data.deriv_f_alphaS.at(MAT_RF(i,s,G))=0;
      for(j=0;j<Xc;j++) data.deriv_f_B.at(MAT_3D(i,j,s,G,Xc)) =0;
    }
  }

 

  for(i=(G-1);i<(G*Xc + G-1);i++){coef.at(i-(G-1))= pars.at(i);
    //Rprintf("%f,",coef.at(i-(G-1)));
  }
  for(i= (G*Xc + G-1); i<(G*Xc + G-1 + S) ; i++) {sp_int.at(i-(G*Xc + G-1)) = pars.at(i);
    // Rprintf("%f,",sp_int.at(i-(G*Xc + G-1)) );
}
  for(i=(G*Xc + G-1 + S) ; i< data.lpar; i++) {sp_dispersion.at(i-(G*Xc + G-1 + S)) = pars.at(i);
    // Rprintf("%f,",sp_dispersion.at(i-(G*Xc + G-1 + S)) );
}

  //Rprintf("\n");
  /*  for(g=0;g<G;g++){
    for(j=0;j<Xc;j++)
      Rprintf("%f, ",coef[MAT_RF(g,j,G)]);
    Rprintf("| %d\n",g);
    }*/
  logl.at(0) = 0;

  for(s=0;s<S;s++){
      start = data.StartEndPos.at(s*2);
      end = data.StartEndPos.at(s*2+1);

      tlog.at(0) = like_mixnbinom_function(estpi, coef,sp_int,sp_dispersion,data.y,data.X,Xr,Xc,start,end, data.tau, s ,data.sum_f_species,data.deriv_f_B,data.deriv_f_alphaS,data.deriv_f_dispersionS,data.log_y_factorial,data.offset);
      logl.at(0)+= tlog.at(0);
  
      data.species_l_contrib.at(s) = tlog.at(0);
  }
  return(logl);

}


 vector <double> calc_mixnbinom_logl(const vector<double> &pars, Optimise_data &data){
  int G,S,Xr,Xc;
  G=data.G;
  S=data.S;
  Xr=data.Xr;
  Xc=data.Xc;
  vector< double > estpi(G-1,0); //vector to hold pi's
  vector< double > coef(Xc*G,0); //vector to hold all coefficents *DOES NOT INCLUDE INTERCEPT*
  vector< double > sp_int(S,0); // vector for species intercepts
  vector< double > sp_dispersion(S,0); //vector for species specific dispersion parameter
  //vector< double > logl( S * G ,0); //output log likelihood
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
      data.sum_f_species.at(MAT_RF(i,s,G))=0;
      data.deriv_f_dispersionS.at(MAT_RF(i,s,G))=0;
      data.deriv_f_alphaS.at(MAT_RF(i,s,G))=0;
      for(j=0;j<Xc;j++) data.deriv_f_B.at(MAT_3D(i,j,s,G,Xc)) =0;
    }
  }

 

  for(i=(G-1);i<(G*Xc + G-1);i++){coef.at(i-(G-1))= pars.at(i);
    //Rprintf("%f,",coef.at(i-(G-1)));
  }
  for(i= (G*Xc + G-1); i<(G*Xc + G-1 + S) ; i++) {sp_int.at(i-(G*Xc + G-1)) = pars.at(i);
    // Rprintf("%f,",sp_int.at(i-(G*Xc + G-1)) );
}
  for(i=(G*Xc + G-1 + S) ; i< data.lpar; i++) {sp_dispersion.at(i-(G*Xc + G-1 + S)) = pars.at(i);
    // Rprintf("%f,",sp_dispersion.at(i-(G*Xc + G-1 + S)) );
}

  //Rprintf("\n");
  /*  for(g=0;g<G;g++){
    for(j=0;j<Xc;j++)
      Rprintf("%f, ",coef[MAT_RF(g,j,G)]);
    Rprintf("| %d\n",g);
    }*/
  logl.at(0) = 0;

  for(s=0;s<S;s++){
      start = data.StartEndPos.at(s*2);
      end = data.StartEndPos.at(s*2+1);

      tlog.at(0) = like_mixnbinom_function(estpi, coef,sp_int,sp_dispersion,data.y,data.X,Xr,Xc,start,end, data.tau, s ,data.sum_f_species,data.deriv_f_B,data.deriv_f_alphaS,data.deriv_f_dispersionS,data.log_y_factorial,data.offset);
      logl.at(0)+= tlog.at(0);
  
      data.species_l_contrib.at(s) = tlog.at(0);
  }
  return(logl);

}

//likelihood function for ipp

double like_mix_ipp_function(vector< double > &estpi, vector < double > &coef, vector < double > &sp_int, vector < double > &sp_dispersion, const double *y, const double *X, int Xr, int Xc, int start, int end, double *tau, int s, vector<double> &sum_f_species, vector<double> &deriv_f_B, vector<double> &deriv_f_alphaS, vector<double> &deriv_f_dispersionS, vector<double> &log_y_factorial, const double *offset){
  int len,i,j,G,g;
  len=end-start+1;
  //vector< AD<double > > p(len,0);
  vector< double  > p(1,0);
  double lpre=0, eps=0,glogl=0, d1=0;
  G = estpi.size();
  vector< double  > sump(G,0);

 
  for(g=0;g<G;g++){
    for(i=0;i<len;i++){
      lpre=offset[i]+sp_int.at(s);
      for(j=0;j<Xc;j++){ 
	lpre +=  X[MAT_RF(i,j,Xr)] * coef[MAT_RF(g,j,G)];

      }
    
      p.at(0)=exp(lpre);  // mu for each archetype at site k
      // datalp.at(i) = p.at(0);
      //Rprintf("%f,",lpre);

      deriv_f_dispersionS.at(MAT_RF(g,s,G)) += (digamma(sp_dispersion.at(s)+y[start+i]) - digamma(sp_dispersion.at(s)) + log(sp_dispersion.at(s)) + 1 - log( p.at(0) + sp_dispersion.at(s)) - (sp_dispersion.at(s)+y[start+i])/(p.at(0) + sp_dispersion.at(s))); // df/dtheta for theta0
      
      //change dl to dfdm which will be (y[start]/p.at(0))-1
      d1=y[start+i]/p.at(0)-(sp_dispersion.at(s)+y[start+i])/(p.at(0)+sp_dispersion.at(s));
      for(j=0;j<Xc;j++) {deriv_f_B.at(MAT_3D(g,j,s,G,Xc)) +=  (d1*p.at(0)*X[MAT_RF(i,j,Xr)]);} 
      deriv_f_alphaS.at(MAT_RF(g,s,G))+= (d1*p.at(0)*1); //This is the chain rule for the df/dalpah

      //if(y[start+i]==0) p.at(0) = 1-p.at(0);

      sump.at(g) += (lgammafn(sp_dispersion.at(s)+y[start+i]) - lgammafn(sp_dispersion.at(s)) + y[start+i]*log(p.at(0)) + sp_dispersion.at(s)*log(sp_dispersion.at(s)) - (sp_dispersion.at(s) + y[start+i])*log(p.at(0)+sp_dispersion.at(s)) - log_y_factorial.at(start+i));

    }
     
 
    /*Rprintf("%d, %d, %f, %f, %f \n",s,g, sump.at(g), sp_int.at(s),sp_dispersion.at(s));
      for(j=0;j<Xc;j++)
	Rprintf("%f, ",coef[MAT_RF(g,j,G)]);
      Rprintf("\n");
    */

    if(g==0) eps = sump.at(g);
    if(sump.at(g) > eps) eps = sump.at(g);

    
  }

  for(g=0;g<G;g++){
    sum_f_species.at(MAT_RF(g,s,G)) = sump.at(g);
    glogl+= estpi.at(g)*exp(sump.at(g) - eps);
  }
  // code for taus can go here
  glogl = log(glogl) + eps;
 
  return(glogl);

}

double like_mixnbinom_function(vector< double > &estpi, vector < double > &coef, vector < double > &sp_int, vector < double > &sp_dispersion, const double *y, const double *X, int Xr, int Xc, int start, int end, double *tau, int s, vector<double> &sum_f_species, vector<double> &deriv_f_B, vector<double> &deriv_f_alphaS, vector<double> &deriv_f_dispersionS, vector<double> &log_y_factorial, const double *offset){
  int len,i,j,G,g;
  len=end-start+1;
  //vector< AD<double > > p(len,0);
  vector< double  > p(1,0);
  double lpre=0, eps=0,glogl=0, d1=0;
  G = estpi.size();
  vector< double  > sump(G,0);

 
  for(g=0;g<G;g++){
    for(i=0;i<len;i++){
      lpre=offset[i]+sp_int.at(s);
      for(j=0;j<Xc;j++){ 
	lpre +=  X[MAT_RF(i,j,Xr)] * coef[MAT_RF(g,j,G)];

      }
    
      p.at(0)=exp(lpre);  // mu for each archetype at site k
      // datalp.at(i) = p.at(0);
      //Rprintf("%f,",lpre);

      deriv_f_dispersionS.at(MAT_RF(g,s,G)) += (digamma(sp_dispersion.at(s)+y[start+i]) - digamma(sp_dispersion.at(s)) + log(sp_dispersion.at(s)) + 1 - log( p.at(0) + sp_dispersion.at(s)) - (sp_dispersion.at(s)+y[start+i])/(p.at(0) + sp_dispersion.at(s))); // df/dtheta for theta0
      
      d1=y[start+i]/p.at(0)-(sp_dispersion.at(s)+y[start+i])/(p.at(0)+sp_dispersion.at(s));
      for(j=0;j<Xc;j++) {deriv_f_B.at(MAT_3D(g,j,s,G,Xc)) +=  (d1*p.at(0)*X[MAT_RF(i,j,Xr)]);} 
      deriv_f_alphaS.at(MAT_RF(g,s,G))+= (d1*p.at(0)*1);

      //if(y[start+i]==0) p.at(0) = 1-p.at(0);

      sump.at(g) += (lgammafn(sp_dispersion.at(s)+y[start+i]) - lgammafn(sp_dispersion.at(s)) + y[start+i]*log(p.at(0)) + sp_dispersion.at(s)*log(sp_dispersion.at(s)) - (sp_dispersion.at(s) + y[start+i])*log(p.at(0)+sp_dispersion.at(s)) - log_y_factorial.at(start+i));

    }
     
 
    /*Rprintf("%d, %d, %f, %f, %f \n",s,g, sump.at(g), sp_int.at(s),sp_dispersion.at(s));
      for(j=0;j<Xc;j++)
	Rprintf("%f, ",coef[MAT_RF(g,j,G)]);
      Rprintf("\n");
    */

    if(g==0) eps = sump.at(g);
    if(sump.at(g) > eps) eps = sump.at(g);

    
  }

  for(g=0;g<G;g++){
    sum_f_species.at(MAT_RF(g,s,G)) = sump.at(g);
    glogl+= estpi.at(g)*exp(sump.at(g) - eps);
  }
  // code for taus can go here
  glogl = log(glogl) + eps;
 
  return(glogl);

}



/*
 *
 * code for tweedie mixture
 *
 *
 */

double optimise_tweedie_function(int n, double *pars, void *ex){

  Optimise_data *data = (Optimise_data *) ex;
  //Optimise_data data =  * (Optimise_data *) ex;
  
  vector<double> logl(1,0);
  int i;
  vector<double> x (n,0);
  
  for(i=0;i<n;i++) x.at(i) = pars[i];

  logl = calc_tweedie_logl(x,*data);

  //logl = data->F.Forward(0,x);

  return(0-logl.at(0));

}

void gradient_tweedie_function(int n, double *pars, double *gr, void *ex ){
  Optimise_data *data = (Optimise_data *) ex;
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
 
  logl = calc_tweedie_logl(x,*data);

  for(g=0;g<(G-1);g++){ 
      add_log_trans+=exp(x.at(g)); //add up transformed pi's
  }
  add_log_trans+=1;
 
  for(g=0;g<G;g++){ //GO through all pi's

      for(j=0;j<Xc;j++){
	for(s=0;s<S;s++){
	  ad_g.at((G-1)+ MAT_RF(g,j,G)) += exp( -1*data->species_l_contrib.at(s) + log(data->parpi.at(g)) + data->sum_f_species.at(MAT_RF(g,s,G))) *data->deriv_f_B.at(MAT_3D(g,j,s,G,Xc));
	  if(j==0) dl_dpi.at(g) += exp(  -1*data->species_l_contrib.at(s)+ data->sum_f_species.at(MAT_RF(g,s,G)));
	}
      }

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

      for(s=0;s<S;s++){
      //calculate for alphas
	ad_g.at((G*Xc + G-1)+s)+= exp(data->sum_f_species.at(MAT_RF(g,s,G)) - data->species_l_contrib.at(s)+ log(data->parpi.at(g)))* data->deriv_f_alphaS.at(MAT_RF(g,s,G));
	
      //calculate for dispersion
	ad_g.at((G*Xc + G-1 + S)+s)+= exp(data->sum_f_species.at(MAT_RF(g,s,G)) - data->species_l_contrib.at(s)+ log(data->parpi.at(g)))* data->deriv_f_dispersionS.at(MAT_RF(g,s,G));
      }
  }
  for(i=0;i<(G-1);i++) pi_mat_deriv.at(MAT_RF(i,(G-1),(G-1))) *= -1;



    for(i=0;i<(G-1);i++){
      for(g=0;g<G;g++){
	ad_g.at(i) += dl_dpi.at(g)* pi_mat_deriv.at(MAT_RF(i,g,(G-1)));
      }
    }

  for(i=0;i<n;i++) gr[i] = 0-ad_g.at(i); 

}


 vector <double> calc_tweedie_logl(const vector<double> &pars, Optimise_data &data){
  int G,S,Xr,Xc;
  G=data.G;
  S=data.S;
  Xr=data.Xr;
  Xc=data.Xc;
  vector< double > estpi(G-1,0); //vector to hold pi's
  vector< double > coef(Xc*G,0); //vector to hold all coefficents *DOES NOT INCLUDE INTERCEPT*
  vector< double > sp_int(S,0); // vector for species intercepts
  vector< double > sp_dispersion(S,0); //vector for species specific dispersion parameter
  //vector< double > logl( S * G ,0); //output log likelihood
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
      data.sum_f_species.at(MAT_RF(i,s,G))=0;
      data.deriv_f_dispersionS.at(MAT_RF(i,s,G))=0;
      data.deriv_f_alphaS.at(MAT_RF(i,s,G))=0;
      for(j=0;j<Xc;j++) data.deriv_f_B.at(MAT_3D(i,j,s,G,Xc)) =0;
    }
  }



    for(i=(G-1);i<(G*Xc + G-1);i++){coef.at(i-(G-1))= pars.at(i);
      //     Rprintf("%f,",coef.at(i-(G-1)));
  }

  for(i= (G*Xc + G-1); i<(G*Xc + G-1 + S) ; i++) {sp_int.at(i-(G*Xc + G-1)) = pars.at(i);
    // Rprintf("%f,",sp_int.at(i-(G*Xc + G-1)) );
}

  for(i=(G*Xc + G-1 + S) ; i< data.lpar; i++) {sp_dispersion.at(i-(G*Xc + G-1 + S)) = pars.at(i);
    //     Rprintf("%f,",sp_dispersion.at(i-(G*Xc + G-1 + S)) );
}

  //Rprintf("\n");
  /*  for(g=0;g<G;g++){
    for(j=0;j<Xc;j++)
      Rprintf("%f, ",coef[MAT_RF(g,j,G)]);
    Rprintf("| %d\n",g);
    }*/
  logl.at(0) = 0;

  for(s=0;s<S;s++){
      start = data.StartEndPos.at(s*2);
      end = data.StartEndPos.at(s*2+1);

      tlog.at(0) = like_tweedie_function(estpi, coef,sp_int,sp_dispersion,data.y,data.X,Xr,Xc,start,end, data.tau, s ,data.sum_f_species,data.deriv_f_B,data.deriv_f_alphaS,data.deriv_f_dispersionS,data.log_y_factorial,data.offset);
      logl.at(0)+= tlog.at(0);
  
      data.species_l_contrib.at(s) = tlog.at(0);
  }
  return(logl);

}


double like_tweedie_function(vector< double > &estpi, vector < double > &coef, vector < double > &sp_int, vector < double > &sp_dispersion, const double *y, const double *X, int Xr, int Xc, int start, int end, double *tau, int s, vector<double> &sum_f_species, vector<double> &deriv_f_B, vector<double> &deriv_f_alphaS, vector<double> &deriv_f_dispersionS, vector<double> &log_y_factorial, const double *offset){
  int len,i,j,G,g;
  len=end-start+1;
  //vector< AD<double > > p(len,0);
  vector< double  > p(1,0);
  double lpre=0, eps=0,glogl=0, d1=0,mu,alpha,muN,muZ,power=1.6;
  G = estpi.size();
  vector< double  > sump(G,0);
  vector<double> beta_for_deriv(Xc+1,0), outDerivs(Xc+2,0),X_for_deriv(Xc+1,0);


  //  lambda <- ( mu^( 2-p)) / ( phi*(2-p))
  //alpha <- ( 2-p) / ( p-1)
  //tau <- phi*(p-1)*mu^(p-1)
  //mu.Z <- alpha * tau

 
  for(g=0;g<G;g++){
    for(i=0;i<len;i++){
      mu=offset[i]+sp_int.at(s);
      beta_for_deriv.at(0) = sp_int.at(s);
      X_for_deriv.at(0) = 1;
      outDerivs.at(0)=0;
      for(j=0;j<Xc;j++){ 
	mu +=  X[MAT_RF(i,j,Xr)] * coef[MAT_RF(g,j,G)];
	beta_for_deriv.at(j+1) = coef[MAT_RF(g,j,G)];
	X_for_deriv.at(j+1) = X[MAT_RF(i,j,Xr)];
	outDerivs.at(j+1)=0;
      }
      outDerivs.at(Xc+1)=0;

      mu = exp(mu);
      muN=pow(mu,(2-power))/(sp_dispersion.at(s)*(2-power)); //changed
      alpha=((2-power)/(power-1));
      muZ=alpha*sp_dispersion.at(s)*(power-1)*pow(mu,(power-1));//changed
      p.at(0)=dTweedie( y[start+i], muN, muZ, alpha, 1);  // mu for each archetype at site k
      // datalp.at(i) = p.at(0);
      //Rprintf("%f,",lpre);
      dTGLM(outDerivs, X_for_deriv, y[start+i], offset[i], beta_for_deriv, sp_dispersion.at(s), power); //changed
      deriv_f_dispersionS.at(MAT_RF(g,s,G)) += 0-outDerivs.at(Xc+1);
            
      for(j=0;j<Xc;j++) {deriv_f_B.at(MAT_3D(g,j,s,G,Xc)) +=  0-outDerivs.at(j+1);}
      deriv_f_alphaS.at(MAT_RF(g,s,G))+= 0-outDerivs.at(0);
      //if(y[start+i]==0) p.at(0) = 1-p.at(0);

      sump.at(g) += p.at(0);

    }
     
 
    //Rprintf("%d, %d, %f, %f, %f \n",s,g, sump.at(g), sp_int.at(s),sp_dispersion.at(s));
      /*   for(j=0;j<Xc;j++)
	Rprintf("%f, ",coef[MAT_RF(g,j,G)]);
      Rprintf("\n");
    */

    if(g==0) eps = sump.at(g);
    if(sump.at(g) > eps) eps = sump.at(g);

    
  }

  for(g=0;g<G;g++){
    sum_f_species.at(MAT_RF(g,s,G)) = sump.at(g);
    glogl+= estpi.at(g)*exp(sump.at(g) - eps);
  }
  // code for taus can go here
  glogl = log(glogl) + eps;
 
  return(glogl);

}
