#include"sam_cpp.h"
#include"sam_cpp_pred.h"

extern "C" { SEXP SAM_predict_C(SEXP Ry, SEXP RX, SEXP Roffset, SEXP Rspp_wts,
						   SEXP Rsite_spp_wts, SEXP Ry_not_na, SEXP Rtaus,
				           SEXP RnS, SEXP RnG, SEXP Rp, SEXP RnObs, SEXP Rdisty,
				           SEXP Ralpha, SEXP Rbeta, SEXP Reta, SEXP Rdisp,
				           SEXP RalphaBoot, SEXP RbetaBoot, SEXP RetaBoot, SEXP RdispBoot,
				           SEXP Rnboot, SEXP Rspp_pt_preds, SEXP Rgrp_pt_preds,
				           SEXP Rspp_boot_preds, SEXP Rgrp_boot_preds, SEXP RoptiDisp){
							   
	sam_cpp_all_classes all;
	int nboot = *(INTEGER(Rnboot));
	double *bootalpha, *bootbeta, *booteta, *bootParms;
	double *spp_pt_preds, *grp_pt_preds, *spp_boot_preds, *grp_boot_preds, *taus;
//	int type = *(INTEGER( Rtype));
//	int bootCount;

	//initialise the data structures -- they are mostly just pointers to REAL()s...
	all.data.setVals(Ry, RX, Roffset, Rspp_wts, Rsite_spp_wts, Ry_not_na, RnS, RnG, Rp, RnObs, Rdisty, RoptiDisp);	//read in the data
	all.params.setVals(all.data, Ralpha, Rbeta, Reta, Rdisp);	//read in the parameters

	//not creating a myFits object, as the data structure for the bootstrap fits is not present.
	vector<double> all_fits(all.data.nG*all.data.nObs*all.data.nS, all.data.NAnum);

	//ptPreds = REAL( RptPreds);
	//for( int i=0; i<all.data.nObs; i++)
		//calcLogPis( logPis, allPis.at(i), all.data, all.parms, i);
	//for( int g=0; g<all.data.nG; g++)
		//for( int i=0; i<all.data.nObs; i++)
			//ptPreds[MATREF2D(i,g,all.data.nObs)] = allPis.at(i).at(g);
			
	calc_mu_fits(all_fits, all.params, all.data); // This should give the g,s,i fits. 		
    spp_pt_preds = REAL( Rspp_pt_preds);
    grp_pt_preds = REAL( Rgrp_pt_preds);
    spp_boot_preds = REAL( Rspp_boot_preds);
    grp_boot_preds = REAL( Rgrp_boot_preds);
    taus = REAL( Rtaus);
	
	for(int s=0; s<all.data.nS; s++)
			for( int i=0; i<all.data.nObs; i++)
			    	for( int g=0; g<all.data.nG; g++)
				        spp_pt_preds[MATREF2D(i,s,all.data.nObs)] += all_fits[MATREF3D(i,g,s,all.data.nObs,all.data.nG)]*taus[MATREF2D(s,g,all.data.nS)];

    vector<double> tau_sum(all.data.nG,all.data.NAnum);
   	for( int g=0; g<all.data.nG; g++){
		for(int s=0; s<all.data.nS; s++){
		    tau_sum[g] +=  taus[MATREF2D(s,g,all.data.nS)];	
			for( int i=0; i<all.data.nObs; i++){
			        grp_pt_preds[MATREF2D(i,g,all.data.nObs)] += all_fits[MATREF3D(i,g,s,all.data.nObs,all.data.nG)]*taus[MATREF2D(s,g,all.data.nS)];			        
			}
	 	}
	}

	for( int g=0; g<all.data.nG; g++) 	
        for( int i=0; i<all.data.nObs; i++)
			        grp_pt_preds[MATREF2D(i,g,all.data.nObs)] = grp_pt_preds[MATREF2D(i,g,all.data.nObs)]/tau_sum[g];
    			

	//setting up the bootstrap values for alpha, tau, beta, disps
	bootalpha = REAL( RalphaBoot);
	bootbeta = REAL( RbetaBoot);
	booteta = REAL( RetaBoot);
	bootParms = (double *) R_alloc(all.params.nTot,sizeof(double));
	int kount;

	for( int b=0; b<nboot; b++){
		kount = 0;
		for( int i=0; i<all.params.nalpha; i++){
			bootParms[kount] = bootalpha[MATREF2D(b,i,nboot)];
			kount++;
		}
		for( int i=0; i<all.params.nbeta; i++){
			bootParms[kount] = bootbeta[MATREF2D(b,i,nboot)];
			kount++;
		}
		for( int i=0; i<all.params.neta; i++){
			bootParms[kount] = booteta[MATREF2D(b,i,nboot)];
			kount++;
		}
		all.params.update( bootParms, all.data);

	//generatte a new set of values for fitting.
	vector<double> all_fits(all.data.nG*all.data.nObs*all.data.nS, all.data.NAnum);
	
	calc_mu_fits(all_fits, all.params, all.data); // This should give the g,s,i fits. 		

    taus = REAL( Rtaus);
	for ( int j=0; j<nboot; j++)
		for( int s=0; s<all.data.nS; s++)
			for( int i=0; i<all.data.nObs; i++)
			    	for( int g=0; g<all.data.nG; g++)
			    	     spp_boot_preds[MATREF3D(i,s,j, all.data.nObs, nboot)] += all_fits[MATREF3D(i,g,s,all.data.nObs,all.data.nG)]*taus[MATREF2D(s,g,all.data.nS)];

    for ( int j=0; j<nboot; j++){
		for( int g=0; g<all.data.nG; g++){
			for(int s=0; s<all.data.nS; s++){
				for( int i=0; i<all.data.nObs; i++){
						grp_boot_preds[MATREF3D(i,g,j,all.data.nObs,nboot)] += all_fits[MATREF3D(i,g,s,all.data.nObs,all.data.nG)]*taus[MATREF2D(s,g,all.data.nS)];			        
				}
			}
		}
	}
	
	for ( int j=0; j<nboot; j++)
		for( int g=0; g<all.data.nG; g++) 	
			for( int i=0; i<all.data.nObs; i++)
			        grp_boot_preds[MATREF3D(i,g,j,all.data.nObs,nboot)] = grp_pt_preds[MATREF2D(i,g,all.data.nObs)]/tau_sum[g];
	}
	
	double foo = -999999;
	SEXP Rfoo = PROTECT(allocVector(REALSXP, 1));
    REAL(Rfoo)[0] = foo;
   	UNPROTECT(1);
	SEXP Rspp_pt_preds_out = PROTECT(allocVector(REALSXP, all.data.nS*all.data.nObs));
	for(int i=0; i<all.data.nS*all.data.nObs; i++) REAL(Rspp_pt_preds_out)[i] = spp_pt_preds[i];
	UNPROTECT(1);		
	SEXP Rgrp_pt_preds_out = PROTECT(allocVector(REALSXP, all.data.nG*all.data.nObs));
	for(int i=0; i<all.data.nG*all.data.nObs; i++) REAL(Rgrp_pt_preds_out)[i] = grp_pt_preds[i];
	UNPROTECT(1);	
	SEXP Rspp_boot_preds_out = PROTECT(allocVector(REALSXP, all.data.nS*all.data.nObs*nboot));
	for(int i=0; i<all.data.nS*all.data.nObs*nboot; i++) REAL(Rspp_boot_preds_out)[i] = spp_boot_preds[i];
	UNPROTECT(1);				
	SEXP Rgrp_boot_preds_out = PROTECT(allocVector(REALSXP, all.data.nG*all.data.nObs*nboot));
	for(int i=0; i<all.data.nG*all.data.nObs*nboot; i++)REAL(Rgrp_boot_preds_out)[i] = grp_boot_preds[i];
	UNPROTECT(1);			
	
	//SEXP Rspp_pt_preds_out = PROTECT(allocVector(REALSXP, all.data.nS*all.data.nObs));
	//for(int s=0; s<all.data.nS; s++)
		//for(int i=0; i<all.data.nObs; i++) 
			//REAL(Rspp_pt_preds_out)[MATREF2D(i,s,all.data.nObs)] = spp_pt_preds[MATREF2D(i,s,all.data.nObs)];
	//UNPROTECT(1);		
	//SEXP Rgrp_pt_preds_out = PROTECT(allocVector(REALSXP, all.data.nG*all.data.nObs));
	//for(int g=0; g<all.data.nG; g++)
		//for(int i=0; i<all.data.nObs; i++) 
			//REAL(Rgrp_pt_preds_out)[MATREF2D(i,g,all.data.nObs)] = grp_pt_preds[MATREF2D(i,g,all.data.nObs)];
	//UNPROTECT(1);	
	//SEXP Rspp_boot_preds_out = PROTECT(allocVector(REALSXP, all.data.nS*all.data.nObs*nboot));
	//for(int s=0; s<all.data.nS; s++)
		//for(int i=0; i<all.data.nObs; i++) 
			//for(int j=0; j<nboot; j++)
				//REAL(Rspp_boot_preds_out)[MATREF3D(i,s,j, all.data.nObs,all.data.nS)] = spp_boot_preds[MATREF3D(i, s, j, all.data.nObs,all.data.nS)];
	//UNPROTECT(1);				
	//SEXP Rgrp_boot_preds_out = PROTECT(allocVector(REALSXP, all.data.nG*all.data.nObs*nboot));
	//for(int g=0; g<all.data.nG; g++)
		//for(int i=0; i<all.data.nObs; i++) 
			//for(int j=0; j<nboot; j++)
				//REAL(Rgrp_boot_preds_out)[MATREF3D(i,g,j, all.data.nObs,all.data.nG)] = grp_boot_preds[MATREF3D(i, g, j, all.data.nObs,all.data.nG)];
	//UNPROTECT(1);				


    const char *names[] = {"foo","spp_pt_pred", "grp_pt_pred", "spp_boot_pred", "grp_boot_pred",""};                   /* note the null string */
	SEXP Rres = PROTECT(mkNamed(VECSXP, names));  /* list of length 4 */
	SET_VECTOR_ELT(Rres, 0, Rfoo);      // species point predictions
	SET_VECTOR_ELT(Rres, 1, Rspp_pt_preds_out);      // species point predictions
	SET_VECTOR_ELT(Rres, 2, Rgrp_pt_preds_out);   	// group point predictions
	SET_VECTOR_ELT(Rres, 3, Rspp_boot_preds_out);   // species boot predictions
	SET_VECTOR_ELT(Rres, 4, Rgrp_boot_preds_out);   // group boot predictions
	UNPROTECT(1);
	return( Rres);

}

}

