#include"negative_binomial_glm_cpp.h"

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

