#include"sam_cpp_pred_v1.h"
#include"sam_cpp2.h"

extern "C" { SEXP sam_cpp_pred(SEXP RX, SEXP RW, SEXP Roffset,
                               SEXP RG, SEXP RS, SEXP RnObs, SEXP Rpx, SEXP Rpw,
                               SEXP Rdisty, SEXP Rlinky, SEXP Rtype, SEXP Rpredtype,
                               SEXP Ralpha, SEXP Rbeta, SEXP Rgamma, SEXP Rtau,
                               SEXP Rbootalpha, SEXP Rbootbeta, SEXP Rbootgamma, SEXP Rboottau,
                               SEXP Rnboot){

  // define the sam pred classes
  sam_pred_classes pred_class;

  //  the data and the parameter sets
  pred_class.data.setVals(RX, RW, Roffset, RS, RG, Rpx, Rpw, RnObs, Rdisty, Rlinky, Rtype, Rpredtype, Rnboot);	//read in the data
  pred_class.params.setParams(pred_class.data, Ralpha, Rbeta, Rgamma, Rtau); //set the parameters for prediction.
  pred_class.bootparams.setParams(pred_class.data, Rbootalpha, Rbootbeta, Rbootgamma, Rboottau); // set the bootstrap parameters

  // vector to catch the point predictions
  int predCols = (pred_class.data.predtype==0) ? pred_class.data.nG : pred_class.data.nS;
  vector<double> ptPreds(pred_class.data.nObs*predCols, 0);

  // predict the point predictions
  if(pred_class.data.predtype==0)
    pt_predict_fun(pred_class.data, pred_class.params, ptPreds);
  else
    pt_predict_species_fun(pred_class.data, pred_class.params, ptPreds);

  // predict the bootstrap predictions
  vector<double> bootPreds(pred_class.data.nObs*predCols*max(pred_class.data.nboot,1),0);
  if(pred_class.data.nboot>0){
    if(pred_class.data.predtype==0)
      boot_predict_fun(pred_class.data, pred_class.bootparams, bootPreds);
    else
      boot_predict_species_fun(pred_class.data, pred_class.bootparams, bootPreds);
  }

  SEXP Rpreds = PROTECT(allocVector(REALSXP, pred_class.data.nObs*predCols));
  for( int ii=0; ii<(pred_class.data.nObs*predCols); ii++) REAL(Rpreds)[ii] = ptPreds[ii];

  SEXP RbootPredsOut = PROTECT(allocVector(REALSXP, pred_class.data.nObs*predCols*max(pred_class.data.nboot,1)));
  for( int ii=0; ii<(pred_class.data.nObs*predCols*max(pred_class.data.nboot,1)); ii++) REAL(RbootPredsOut)[ii] = bootPreds[ii];

  const char *names[] = {"preds", "bootPreds", ""};
  SEXP Rres = PROTECT(mkNamed(VECSXP, names));
  SET_VECTOR_ELT(Rres, 0, Rpreds);
  SET_VECTOR_ELT(Rres, 1, RbootPredsOut);
  UNPROTECT(3);

  return( Rres);
  }
} // end of the C external wrapper!

// shared linear-predictor / mean computation, used for both point and (per-draw) bootstrap prediction,
// and for both archetype- and species-level output (they differ only in how mus gets weighted by tau).
void calc_sam_pred_mus(vector<double> &mus, const sam_pred_data &dat,
                       const double *alpha, const double *beta, const double *gamma){

  double lp;

  for( int g=0; g<dat.nG; g++){
    for( int s=0; s<dat.nS; s++){
      for( int i=0; i<dat.nObs; i++){
        lp = alpha[s] + dat.offset[i];
        for( int j=0;j<dat.nPX; j++){
          lp += beta[MATREF2D(g,j,(dat.nG))] * dat.X[MATREF2D(i,j,dat.nObs)];
        }
        if(dat.nPW>0){
          for( int l=0;l<dat.nPW; l++){
            lp += gamma[MATREF2D(s,l,(dat.nS))] * dat.W[MATREF2D(i,l,dat.nObs)];
          }
        }

        //If the response is wanted transform the lp.
        if(dat.type==0){
          if(dat.disty==1){//bernoulli
            if(dat.linky==0)
              mus[MATREF3D(i,s,g,dat.nObs,dat.nS)] = inverse_logit(lp);
            if(dat.linky==1)
              mus[MATREF3D(i,s,g,dat.nObs,dat.nS)] = inverse_cloglog(lp);
          }
          if(dat.disty==2){ //poisson
            mus[MATREF3D(i,s,g,dat.nObs,dat.nS)] = inverse_log(lp);
          }
          if(dat.disty==3){ //negative binomial
            mus[MATREF3D(i,s,g,dat.nObs,dat.nS)] = inverse_log(lp);
          }
          if(dat.disty==4){ //tweedie
            mus[MATREF3D(i,s,g,dat.nObs,dat.nS)] = inverse_log(lp);
          }
          if(dat.disty==5){//normal
            mus[MATREF3D(i,s,g,dat.nObs,dat.nS)] = lp;
          }
          if(dat.disty==6){//binomial
            if(dat.linky==0)
              mus[MATREF3D(i,s,g,dat.nObs,dat.nS)] = inverse_logit(lp);
            if(dat.linky==1)
              mus[MATREF3D(i,s,g,dat.nObs,dat.nS)] = inverse_cloglog(lp);
          }
        }
        if(dat.type==1){
          mus[MATREF3D(i,s,g,dat.nObs,dat.nS)] = lp;
        }
      }
    }
  }
}

// archetype-level weighting: preds[i,g] = sum_s( mus[i,s,g]*tau[s,g] ) / sum_s( tau[s,g] )
void weight_archetype_preds(vector<double> &preds, const vector<double> &mus, const sam_pred_data &dat, const double *tau){

  vector<double> tau_g(dat.nG,0);

  for( int g=0; g<dat.nG; g++){
    for (int ss=0; ss < dat.nS; ss++)
      tau_g[g] += tau[MATREF2D(ss,g,dat.nS)];

    for (int ii=0; ii < dat.nObs; ii++ ){
      for(int ss=0; ss < dat.nS;ss++){
        preds[MATREF2D(ii,g,dat.nObs)] += mus[MATREF3D(ii,ss,g,dat.nObs,dat.nS)]*tau[MATREF2D(ss,g,dat.nS)];
      }
      preds[MATREF2D(ii,g,dat.nObs)] = preds[MATREF2D(ii,g,dat.nObs)]/tau_g[g];
    }
  }
}

// species-level weighting: preds[i,s] = sum_g( mus[i,s,g]*tau[s,g] )  -- tau rows already sum to one, no
// further normalisation needed.
void weight_species_preds(vector<double> &preds, const vector<double> &mus, const sam_pred_data &dat, const double *tau){

  for( int g=0; g<dat.nG; g++){
    for (int ii=0; ii < dat.nObs; ii++ ){
      for(int ss=0; ss < dat.nS;ss++){
        preds[MATREF2D(ii,ss,dat.nObs)] += mus[MATREF3D(ii,ss,g,dat.nObs,dat.nS)]*tau[MATREF2D(ss,g,dat.nS)];
      }
    }
  }
}

// point prediction for sams -- archetype level.
void pt_predict_fun(const sam_pred_data &dat, const sam_pred_params &params, vector<double> &preds){

  vector<double> mus(dat.nG*dat.nS*dat.nObs, 0);
  calc_sam_pred_mus(mus, dat, params.Alpha, params.Beta, params.Gamma);
  weight_archetype_preds(preds, mus, dat, params.Tau);
}

// point prediction for sams -- species level.
void pt_predict_species_fun(const sam_pred_data &dat, const sam_pred_params &params, vector<double> &preds){

  vector<double> mus(dat.nG*dat.nS*dat.nObs, 0);
  calc_sam_pred_mus(mus, dat, params.Alpha, params.Beta, params.Gamma);
  weight_species_preds(preds, mus, dat, params.Tau);
}

//bootstrap prediction for sams -- archetype level. tau is recomputed (in R) for every bootstrap draw,
//so a full (nboot x nS x nG) tau array is expected here -- unlike an earlier version of this function,
//the same tau is NOT recycled across draws.
void boot_predict_fun(const sam_pred_data &dat, const sam_pred_bootparams &bootparams, vector<double> &bootpreds){

  vector<double> mus(dat.nG*dat.nS*dat.nObs, 0);
  vector<double> alpha_b(dat.nS,0), beta_b(dat.nG*dat.nPX,0), gamma_b(dat.nS*dat.nPW,0), tau_b(dat.nS*dat.nG,0);
  vector<double> preds_b(dat.nObs*dat.nG,0);

  for(int b=0; b<dat.nboot; b++){

    for(int s=0; s<dat.nS; s++) alpha_b[s] = bootparams.bootAlpha[MATREF2D(b,s,dat.nboot)];
    for(int i=0; i<(dat.nG*dat.nPX); i++) beta_b[i] = bootparams.bootBeta[MATREF2D(b,i,dat.nboot)];
    if(dat.nPW>0)
      for(int i=0; i<(dat.nS*dat.nPW); i++) gamma_b[i] = bootparams.bootGamma[MATREF2D(b,i,dat.nboot)];
    for(int g=0; g<dat.nG; g++)
      for(int s=0; s<dat.nS; s++)
        tau_b[MATREF2D(s,g,dat.nS)] = bootparams.bootTau[MATREF3D(b,s,g,dat.nboot,dat.nS)];

    mus.assign(mus.size(),0.0);
    preds_b.assign(preds_b.size(),0.0);
    calc_sam_pred_mus(mus, dat, &alpha_b[0], &beta_b[0], &gamma_b[0]);
    weight_archetype_preds(preds_b, mus, dat, &tau_b[0]);

    for( int g=0; g<dat.nG; g++)
      for (int ii=0; ii < dat.nObs; ii++ )
        bootpreds[MATREF3D(ii,g,b,dat.nObs,dat.nG)] = preds_b[MATREF2D(ii,g,dat.nObs)];
  }
}

//bootstrap prediction for sams -- species level.
void boot_predict_species_fun(const sam_pred_data &dat, const sam_pred_bootparams &bootparams, vector<double> &bootpreds){

  vector<double> mus(dat.nG*dat.nS*dat.nObs, 0);
  vector<double> alpha_b(dat.nS,0), beta_b(dat.nG*dat.nPX,0), gamma_b(dat.nS*dat.nPW,0), tau_b(dat.nS*dat.nG,0);
  vector<double> preds_b(dat.nObs*dat.nS,0);

  for(int b=0; b<dat.nboot; b++){

    for(int s=0; s<dat.nS; s++) alpha_b[s] = bootparams.bootAlpha[MATREF2D(b,s,dat.nboot)];
    for(int i=0; i<(dat.nG*dat.nPX); i++) beta_b[i] = bootparams.bootBeta[MATREF2D(b,i,dat.nboot)];
    if(dat.nPW>0)
      for(int i=0; i<(dat.nS*dat.nPW); i++) gamma_b[i] = bootparams.bootGamma[MATREF2D(b,i,dat.nboot)];
    for(int g=0; g<dat.nG; g++)
      for(int s=0; s<dat.nS; s++)
        tau_b[MATREF2D(s,g,dat.nS)] = bootparams.bootTau[MATREF3D(b,s,g,dat.nboot,dat.nS)];

    mus.assign(mus.size(),0.0);
    preds_b.assign(preds_b.size(),0.0);
    calc_sam_pred_mus(mus, dat, &alpha_b[0], &beta_b[0], &gamma_b[0]);
    weight_species_preds(preds_b, mus, dat, &tau_b[0]);

    for( int s=0; s<dat.nS; s++)
      for (int ii=0; ii < dat.nObs; ii++ )
        bootpreds[MATREF3D(ii,s,b,dat.nObs,dat.nS)] = preds_b[MATREF2D(ii,s,dat.nObs)];
  }
}



// Classes to setup data for sam prediction.
// A wrapper class
sam_pred_classes::sam_pred_classes(){};
sam_pred_classes::~sam_pred_classes(){};

// A sam pred data class
sam_pred_data::sam_pred_data(){};
sam_pred_data::~sam_pred_data(){};

void sam_pred_data::setVals(SEXP &RX, SEXP &RW, SEXP &Roffset,
                            SEXP &RS, SEXP &RG, SEXP &Rpx, SEXP &Rpw,
                            SEXP &RnObs, SEXP &Rdisty, SEXP &Rlinky, SEXP &Rtype, SEXP &Rpredtype, SEXP &Rnboot){

  //integers
  nS = *(INTEGER( RS));
  nG = *(INTEGER( RG));
  nPX = *(INTEGER( Rpx));
  nPW = *(INTEGER( Rpw));
  nObs = *(INTEGER( RnObs));
  disty = *(INTEGER( Rdisty));
  linky = *(INTEGER( Rlinky));
  type = *(INTEGER( Rtype));
  predtype = *(INTEGER( Rpredtype));
  nboot = *(INTEGER( Rnboot));

  // input data - doubles
  X = REAL( RX);
  W = REAL( RW);
  offset = REAL( Roffset);

}

// A sam pred params class
sam_pred_params::sam_pred_params(){};
sam_pred_params::~sam_pred_params(){};


void sam_pred_params::setParams(const sam_pred_data &dat, SEXP &Ralpha, SEXP &Rbeta, SEXP &Rgamma, SEXP &Rtau){

  // Actual parameters
  Alpha = REAL( Ralpha);
  Beta = REAL( Rbeta);
  Gamma = REAL( Rgamma);
  Tau = REAL( Rtau);

  // how many of each indices
  nalpha = dat.nS;
  nbeta = (dat.nG*dat.nPX);
  ngamma = (dat.nS*dat.nPW);
  ntau = (dat.nG*dat.nS);

  nTot = nalpha + nbeta + ngamma + ntau;
}

// A sam pred params class
sam_pred_bootparams::sam_pred_bootparams(){};
sam_pred_bootparams::~sam_pred_bootparams(){};

void sam_pred_bootparams::setParams(const sam_pred_data &dat, SEXP &Rbootalpha, SEXP &Rbootbeta, SEXP &Rbootgamma, SEXP &Rboottau){

  // Actual parameters
  bootAlpha = REAL( Rbootalpha);
  bootBeta = REAL( Rbootbeta);
  bootGamma = REAL( Rbootgamma);
  bootTau = REAL( Rboottau); // (nboot x nS x nG), one tau matrix per bootstrap draw

  // How many of each indices
  nbootalpha = dat.nS;
  nbootbeta = (dat.nG*dat.nPX);
  nbootgamma = (dat.nS*dat.nPW);
  nboottau = (dat.nboot*dat.nG*dat.nS);

  nbootTot = nbootalpha + nbootbeta + nbootgamma + nboottau;
}
