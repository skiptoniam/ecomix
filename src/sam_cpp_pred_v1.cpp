#include"sam_cpp_pred_v1.h"
#include"sam_cpp2.h"

extern "C" { SEXP sam_cpp_pred(SEXP RX, SEXP RW, SEXP RU, SEXP Roffset,
                               SEXP RG, SEXP RS, SEXP RnObs, SEXP Rpx, SEXP Rpw, SEXP Rpu,
                               SEXP Rdisty, SEXP Rlinky, SEXP Rtype,
                               SEXP Ralpha, SEXP Rbeta, SEXP Rgamma, SEXP Rdelta, SEXP Rtau,
                               // SEXP Rbootalpha, SEXP Rbootbeta, SEXP Rbootgamma, SEXP Rbootdelta, SEXP Rboottau,
                               SEXP Rnboot){//}, SEXP RptPreds){//}, SEXP RbootPreds){

  sam_pred_classes pred_class;

  // vector to catch tau sum.
  vector<double> tau_g(pred_class.data.nG,0);

  // bootalpha = REAL(Rbootalpha);
  // bootbeta = REAL(Rbootbeta);
  // bootgamma = REAL(Rbootgamma);
  // bootdelta = REAL(Rbootdelta);
  // boottau = REAL(Rboottau);

  // int kount;

  //  the data and the parameter sets
  pred_class.data.setVals(RX, RW, RU, Roffset, RS, RG, Rpx, Rpw, Rpu, RnObs, Rdisty, Rlinky, Rtype, Rnboot);	//read in the data
  pred_class.params.setParams(pred_class.data, Ralpha, Rbeta, Rgamma, Rdelta, Rtau); //set the parameters for prediction.


  vector<double> mus(pred_class.data.nG*pred_class.data.nS*pred_class.data.nObs, 0);
  vector<double> lps(pred_class.data.nG*pred_class.data.nS, 0);	//the nG x nS intercepts
  vector<double> etaAll(pred_class.data.nObs,0); //linear predictor for bias/all parameters.
  vector<double> ptPreds(pred_class.data.nObs*pred_class.data.nG, 0);
  double lp;	//the lin pred for the gth group, sth species and ith site
  // double lp_sppEta=0.0;   // split the linear predictor into two components.

  // need a flag which will decide if the user need archetype or species level predictions.

  // If there is an all/bias parameter in the model what is overall contribution of this into the model.
  if(pred_class.data.nPU>0){
    for(int i=0; i<pred_class.data.nObs; i++){
      for(int d=0; d<pred_class.data.nPU; d++){
        etaAll[i] += pred_class.data.U[MATREF2D(i,d,pred_class.data.nObs)]*pred_class.params.Delta[d];
      }
    }
  }

  //predict for each archetype group.
  for( int g=0; g<pred_class.data.nG; g++){
    for( int s=0; s<pred_class.data.nS; s++){
      lps.at(MATREF2D(g,s,pred_class.data.nG)) = pred_class.params.Alpha[s];
      for( int i=0; i<pred_class.data.nObs; i++){
        lp = lps.at(MATREF2D(g,s,pred_class.data.nG)) + pred_class.data.offset[i] + etaAll[i];
        for( int j=0;j<pred_class.data.nPX; j++){
          lp += pred_class.params.Beta[MATREF2D(g,j,(pred_class.data.nG))] * pred_class.data.X[MATREF2D(i,j,pred_class.data.nObs)];
        }
        if(pred_class.data.nPW>0){
          for( int l=0;l<pred_class.data.nPW; l++){
            lp += pred_class.params.Gamma[MATREF2D(s,l,(pred_class.data.nS))] * pred_class.data.W[MATREF2D(i,l,pred_class.data.nObs)];
          }
        }

        //If the response is wanted transform the lp.
        if(pred_class.data.type==0){
        if(pred_class.data.disty==1){//bernoulli
          if(pred_class.data.linky==0)
            mus[MATREF3D(i,s,g,pred_class.data.nObs,pred_class.data.nS)] = inverse_logit(lp);
          if(pred_class.data.linky==1)
            mus[MATREF3D(i,s,g,pred_class.data.nObs,pred_class.data.nS)] = inverse_cloglog(lp);
        }
        if(pred_class.data.disty==2){ //poisson
          mus[MATREF3D(i,s,g,pred_class.data.nObs,pred_class.data.nS)] = inverse_log(lp);
        }
        if(pred_class.data.disty==3){ //negative binomial
          mus[MATREF3D(i,s,g,pred_class.data.nObs,pred_class.data.nS)] = inverse_log(lp);
        }
        if(pred_class.data.disty==4){ //tweedie
          mus[MATREF3D(i,s,g,pred_class.data.nObs,pred_class.data.nS)] = inverse_log(lp);
        }
        if(pred_class.data.disty==5){//normal
          mus[MATREF3D(i,s,g,pred_class.data.nObs,pred_class.data.nS)] = lp;
        }
        if(pred_class.data.disty==6){//binomial
          if(pred_class.data.linky==0)
            mus[MATREF3D(i,s,g,pred_class.data.nObs,pred_class.data.nS)] = inverse_logit(lp);
          if(pred_class.data.linky==1)
            mus[MATREF3D(i,s,g,pred_class.data.nObs,pred_class.data.nS)] = inverse_cloglog(lp);
        }
        }
        if(pred_class.data.type==1){
          mus[MATREF3D(i,s,g,pred_class.data.nObs,pred_class.data.nS)] = lp;
        }
      }
    }

     // do the tau sums to get the sum tau per archetype
     for (int ss=0; ss < pred_class.data.nS; ss++){
       tau_g[g] += pred_class.params.Tau[MATREF2D(ss,g,pred_class.data.nS)];

     }

     for (int ii=0; ii < pred_class.data.nObs; ii++ ){
       for(int ss=0; ss < pred_class.data.nS;ss++){
       mus[MATREF3D(ii,ss,g,pred_class.data.nObs,pred_class.data.nS)] = mus[MATREF3D(ii,ss,g,pred_class.data.nObs,pred_class.data.nS)]*pred_class.params.Tau[MATREF2D(ss,g,pred_class.data.nS)];
       ptPreds[MATREF2D(ii,g,pred_class.data.nObs)] += mus[MATREF3D(ii,ss,g,pred_class.data.nObs,pred_class.data.nS)];
      }
    }

     for (int ii=0; ii < pred_class.data.nObs; ii++ ){
         ptPreds[MATREF2D(ii,g,pred_class.data.nObs)] = ptPreds[MATREF2D(ii,g,pred_class.data.nObs)]/tau_g[g];
     }
  }

  SEXP Rpreds = PROTECT(allocVector(REALSXP, pred_class.data.nObs*pred_class.data.nG));
  for( int ii=0; ii<(pred_class.data.nObs*pred_class.data.nG); ii++) REAL(Rpreds)[ii] = ptPreds[ii];
  UNPROTECT(1);

  const char *names[] = {"preds",""};                   /* note the null string */
  SEXP Rres = PROTECT(mkNamed(VECSXP, names));  /* list of length 3 */
  SET_VECTOR_ELT(Rres, 0, Rpreds);        // predictions
  UNPROTECT(1);

  return( Rres);
  }
} // end of the C external wrapper!

// //write a quick function for doing the internal loop.
// void sam_cpp_archetype_pred(){
//
// }


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


void sam_pred_params::setParams(const sam_pred_data &dat, SEXP &Ralpha, SEXP &Rbeta,SEXP &Rgamma,SEXP &Rdelta,SEXP &Rtau){

  // Actual parameters
  Alpha = REAL( Ralpha);
  Beta = REAL( Rbeta);
  Gamma = REAL( Rgamma);
  Delta = REAL( Rdelta);
  Tau = REAL( Rtau);

  // how many of each indices
  nalpha = dat.nS*dat.nboot;
  nbeta = (dat.nG*dat.nPX)*dat.nboot;
  ngamma = (dat.nS*dat.nPW)*dat.nboot;
  ndelta = dat.nPU*dat.nboot;
  ntau = (dat.nG*dat.nS)*dat.nboot;

  nTot = nalpha + nbeta + ngamma + ndelta + ntau;
}
