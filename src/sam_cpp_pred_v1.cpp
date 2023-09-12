#include"sam_cpp_pred_v1.h"
#include"sam_cpp2.h"

extern "C" { SEXP sam_cpp_pred(SEXP RX, SEXP RW, SEXP RU, SEXP Roffset,
                               SEXP RG, SEXP RS, SEXP RnObs, SEXP Rpx, SEXP Rpw, SEXP Rpu,
                               SEXP Rdisty, SEXP Rlinky, SEXP Rtype,
                               SEXP Ralpha, SEXP Rbeta, SEXP Rgamma, SEXP Rdelta, SEXP Rtau,
                               SEXP Rbootalpha, SEXP Rbootbeta, SEXP Rbootgamma, SEXP Rbootdelta, SEXP Rboottau, // currently recycling the same tau for all bootstraps.
                               SEXP Rnboot){//}, SEXP RptPreds){//}, SEXP RbootPreds){

  // define the sam pred classes
  sam_pred_classes pred_class;

  //  the data and the parameter sets
  pred_class.data.setVals(RX, RW, RU, Roffset, RS, RG, Rpx, Rpw, Rpu, RnObs, Rdisty, Rlinky, Rtype, Rnboot);	//read in the data
  pred_class.params.setParams(pred_class.data, Ralpha, Rbeta, Rgamma, Rdelta, Rtau); //set the parameters for prediction.
  pred_class.bootparams.setParams(pred_class.data, Rbootalpha, Rbootbeta, Rbootgamma, Rbootdelta, Rboottau); // set the bootstrap parameters

  // vector to catch the point predictions
  vector<double> ptPreds(pred_class.data.nObs*pred_class.data.nG, 0);

  // predict the point predictions
  pt_predict_fun(pred_class.data, pred_class.params, ptPreds);

  // predict the bootstrap predictions
  // vector<double> bootPreds(pred_class.data.nObs*pred_class.data.nG*pred_class.data.nboot,0);
  if(pred_class.data.nboot>0){
    vector<double> bootPreds(pred_class.data.nObs*pred_class.data.nG*pred_class.data.nboot,0);
    boot_predict_fun(pred_class.data, pred_class.bootparams, bootPreds);
  }



  SEXP Rpreds = PROTECT(allocVector(REALSXP, pred_class.data.nObs*pred_class.data.nG));
  for( int ii=0; ii<(pred_class.data.nObs*pred_class.data.nG); ii++) REAL(Rpreds)[ii] = ptPreds[ii];
  UNPROTECT(1);

  const char *names[] = {"preds"};                   /* note the null string */
  SEXP Rres = PROTECT(mkNamed(VECSXP, names));  /* list of length 3 */
  SET_VECTOR_ELT(Rres, 0, Rpreds);        // predictions
  UNPROTECT(1);

  return( Rres);
  }
} // end of the C external wrapper!

// point prediction for sams - currently only archetypes - need to do conditional species prediction too.
void pt_predict_fun(const sam_pred_data &dat, const sam_pred_params &params, vector<double> &preds){

  vector<double> tau_g(dat.nG,0);
  vector<double> mus(dat.nG*dat.nS*dat.nObs, 0);
  vector<double> lps(dat.nG*dat.nS, 0);	//the nG x nS intercepts
  vector<double> etaAll(dat.nObs,0); //linear predictor for bias/all parameters.
  double lp;


  // If there is an all/bias parameter in the model what is overall contribution of this into the model.
if(dat.nPU>0){
  for(int i=0; i<dat.nObs; i++){
    for(int d=0; d<dat.nPU; d++){
      etaAll[i] += dat.U[MATREF2D(i,d,dat.nObs)]*params.Delta[d];
    }
  }
}

//predict for each archetype group.
for( int g=0; g<dat.nG; g++){
  for( int s=0; s<dat.nS; s++){
    lps.at(MATREF2D(g,s,dat.nG)) = params.Alpha[s];
    for( int i=0; i<dat.nObs; i++){
      lp = lps.at(MATREF2D(g,s,dat.nG)) + dat.offset[i] + etaAll[i];
      for( int j=0;j<dat.nPX; j++){
        lp += params.Beta[MATREF2D(g,j,(dat.nG))] * dat.X[MATREF2D(i,j,dat.nObs)];
      }
      if(dat.nPW>0){
        for( int l=0;l<dat.nPW; l++){
          lp += params.Gamma[MATREF2D(s,l,(dat.nS))] * dat.W[MATREF2D(i,l,dat.nObs)];
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

  // do the tau sums to get the sum tau per archetype
  for (int ss=0; ss < dat.nS; ss++){
    tau_g[g] += params.Tau[MATREF2D(ss,g,dat.nS)];

  }

  for (int ii=0; ii < dat.nObs; ii++ ){
    for(int ss=0; ss < dat.nS;ss++){
      mus[MATREF3D(ii,ss,g,dat.nObs,dat.nS)] = mus[MATREF3D(ii,ss,g,dat.nObs,dat.nS)]*params.Tau[MATREF2D(ss,g,dat.nS)];
      preds[MATREF2D(ii,g,dat.nObs)] += mus[MATREF3D(ii,ss,g,dat.nObs,dat.nS)];
    }
  }

  for (int ii=0; ii < dat.nObs; ii++ ){
      preds[MATREF2D(ii,g,dat.nObs)] = preds[MATREF2D(ii,g,dat.nObs)]/tau_g[g];
    }
  }

}


//predict for bootstrap objects
void boot_predict_fun(const sam_pred_data &dat, const sam_pred_bootparams &bootparams, vector<double> &bootpreds){

  vector<double> tau_g(dat.nG,0);
  vector<double> mus(dat.nG*dat.nS*dat.nObs, 0);
  vector<double> lps(dat.nG*dat.nS, 0);	//the nG x nS intercepts
  vector<double> etaAll(dat.nObs,0); //linear predictor for bias/all parameters.
  double lp;

  for(int b=0; b<dat.nboot; b++){

    // tau_g.assign(tau_g.size(),0);
    // vector<double> tau_g(dat.nG,0);
    // vector<double> mus(dat.nG*dat.nS*dat.nObs, 0);
    // vector<double> lps(dat.nG*dat.nS, 0);	//the nG x nS intercepts
    // vector<double> etaAll(dat.nObs,0); //linear predictor for bias/all parameters.
    // double lp;

  // If there is an all/bias parameter in the model what is overall contribution of this into the model.
  if(dat.nPU>0){
    for(int i=0; i<dat.nObs; i++){
      for(int d=0; d<dat.nPU; d++){
        etaAll[i] += dat.U[MATREF2D(i,d,dat.nObs)]*bootparams.bootDelta[MATREF2D(b,d,dat.nboot)];
      }
    }
  }

  //predict for each archetype group.
  for( int g=0; g<dat.nG; g++){
    for( int s=0; s<dat.nS; s++){
      lps.at(MATREF2D(g,s,dat.nG)) = bootparams.bootAlpha[MATREF2D(b,s,dat.nboot)];
      for( int i=0; i<dat.nObs; i++){
        lp = lps.at(MATREF2D(g,s,dat.nG)) + dat.offset[i] + etaAll[i];
        for( int j=0;j<dat.nPX; j++){
          lp += bootparams.bootBeta[MATREF3D(b,g,j,dat.nboot,dat.nG)] * dat.X[MATREF2D(i,j,dat.nObs)];
        }
        if(dat.nPW>0){
          for( int l=0;l<dat.nPW; l++){
            lp += bootparams.bootGamma[MATREF3D(b,s,l,dat.nboot,dat.nS)] * dat.W[MATREF2D(i,l,dat.nObs)];
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

    // do the tau sums to get the sum tau per archetype
    for (int ss=0; ss < dat.nS; ss++){
      tau_g[g] += bootparams.bootTau[MATREF2D(ss,g,dat.nS)];
    }

    for (int ii=0; ii < dat.nObs; ii++ ){
      for(int ss=0; ss < dat.nS;ss++){
        mus[MATREF3D(ii,ss,g,dat.nObs,dat.nS)] = mus[MATREF3D(ii,ss,g,dat.nObs,dat.nS)]*bootparams.bootTau[MATREF2D(ss,g,dat.nS)];
        bootpreds[MATREF3D(ii,g,b,dat.nObs,dat.nG)] += mus[MATREF3D(ii,ss,g,dat.nObs,dat.nS)];
      }
    }

    for (int ii=0; ii < dat.nObs; ii++ ){
      bootpreds[MATREF3D(ii,g,b,dat.nObs,dat.nG)] = bootpreds[MATREF3D(ii,g,b,dat.nObs,dat.nG)]/tau_g[g];
    }
  }
  }

}



// Classes to setup data for sam prediction.
// A wrapper class
sam_pred_classes::sam_pred_classes(){};
sam_pred_classes::~sam_pred_classes(){};

// A sam pred data class
sam_pred_data::sam_pred_data(){};
sam_pred_data::~sam_pred_data(){};

void sam_pred_data::setVals(SEXP &RX, SEXP &RW, SEXP &RU, SEXP &Roffset,
                            SEXP &RS, SEXP &RG, SEXP &Rpx, SEXP &Rpw, SEXP &Rpu,
                            SEXP &RnObs, SEXP &Rdisty, SEXP &Rlinky, SEXP &Rtype, SEXP &Rnboot){

  //integers
  nS = *(INTEGER( RS));
  nG = *(INTEGER( RG));
  nPX = *(INTEGER( Rpx));
  nPW = *(INTEGER( Rpw));
  nPU = *(INTEGER( Rpu));
  nObs = *(INTEGER( RnObs));
  disty = *(INTEGER( Rdisty));
  linky = *(INTEGER( Rlinky));
  type = *(INTEGER( Rtype));
  nboot = *(INTEGER( Rnboot));

  // input data - doubles
  X = REAL( RX);
  W = REAL( RW);
  U = REAL( RU);
  offset = REAL( Roffset);

}

// A sam pred params class
sam_pred_params::sam_pred_params(){};
sam_pred_params::~sam_pred_params(){};


void sam_pred_params::setParams(const sam_pred_data &dat, SEXP &Ralpha, SEXP &Rbeta, SEXP &Rgamma, SEXP &Rdelta, SEXP &Rtau){

  // Actual parameters
  Alpha = REAL( Ralpha);
  Beta = REAL( Rbeta);
  Gamma = REAL( Rgamma);
  Delta = REAL( Rdelta);
  Tau = REAL( Rtau);

  // how many of each indices
  nalpha = dat.nS;
  nbeta = (dat.nG*dat.nPX);
  ngamma = (dat.nS*dat.nPW);
  ndelta = dat.nPU;
  ntau = (dat.nG*dat.nS);

  nTot = nalpha + nbeta + ngamma + ndelta + ntau;
}

// A sam pred params class
sam_pred_bootparams::sam_pred_bootparams(){};
sam_pred_bootparams::~sam_pred_bootparams(){};

void sam_pred_bootparams::setParams(const sam_pred_data &dat, SEXP &Rbootalpha, SEXP &Rbootbeta, SEXP &Rbootgamma,SEXP &Rbootdelta, SEXP &Rboottau){

  // Actual parameters
  bootAlpha = REAL( Rbootalpha);
  bootBeta = REAL( Rbootbeta);
  bootGamma = REAL( Rbootgamma);
  bootDelta = REAL( Rbootdelta);
  bootTau = REAL( Rboottau);

  // How many of each indices
  nbootalpha = dat.nS;
  nbootbeta = (dat.nG*dat.nPX);
  nbootgamma = (dat.nS*dat.nPW);
  nbootdelta = dat.nPU;
  nboottau = (dat.nG*dat.nS);

  nbootTot = nbootalpha + nbootbeta + nbootgamma + nbootdelta + nboottau;
}
