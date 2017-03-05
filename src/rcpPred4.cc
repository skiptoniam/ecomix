#include"rcp4.h"
#include"rcpPred4.h"

extern "C" { SEXP RCP_predict_C( SEXP Ry, SEXP RX, SEXP RW, SEXP Roffset, SEXP Rwts,
				SEXP RS, SEXP RG, SEXP Rp, SEXP Rpw, SEXP RnObs, SEXP Rdisty,
				SEXP Ralpha, SEXP Rtau, SEXP Rbeta, SEXP Rgamma, SEXP Rdisps, SEXP Rpowers,
				SEXP Rconc, SEXP Rsd, SEXP RsdGamma, SEXP RdispLocat, SEXP RdispScale,
				SEXP RalphaBoot, SEXP RtauBoot, SEXP RbetaBoot,	//only using parameters associated with pis. That is what we are predicting.
				SEXP Rnboot,
				SEXP RptPreds, SEXP RbootPreds, SEXP RoptiDisp)
{
	allClasses all;
	int nboot = *(INTEGER(Rnboot));
	double *bootalpha, *boottau, *bootbeta, *bootParms;
	double *ptPreds, *bootPreds;
//	int type = *(INTEGER( Rtype));
//	int bootCount;

	//initialise the data structures -- they are mostly just pointers to REAL()s...
	all.data.setVals( Ry, RX, RW, Roffset, RS, RG, Rp, Rpw, RnObs, Rdisty, RoptiDisp, Rwts);	//read in the data
	all.parms.setVals( all.data, Ralpha, Rbeta, Rtau, Rgamma, Rdisps, Rpowers, Rconc, Rsd, RsdGamma, RdispLocat, RdispScale);	//read in the parameters

	//not creating a myFits object, as the data structure for the bootstrap fits is not present.
	vector<double> logPis(all.data.nG, all.data.NAnum);//, pis( dat.nG, dat.NAnum);
	vector< vector<double> > allPis( all.data.nObs, vector<double> (all.data.nG, all.data.NAnum));	//(nObs x nG) matrix --2D

	ptPreds = REAL( RptPreds);
	for( int i=0; i<all.data.nObs; i++)
		calcLogPis( logPis, allPis.at(i), all.data, all.parms, i);
	for( int g=0; g<all.data.nG; g++)
		for( int i=0; i<all.data.nObs; i++)
			ptPreds[MATREF(i,g,all.data.nObs)] = allPis.at(i).at(g);

	//setting up the bootstrap values for alpha, tau, beta, disps
	bootalpha = REAL( RalphaBoot);
	boottau = REAL( RtauBoot);
	bootbeta = REAL( RbetaBoot);
	bootParms = (double *) R_alloc(all.parms.nTot,sizeof(double));
	int kount;

	bootPreds = REAL( RbootPreds);
	for( int b=0; b<nboot; b++){
		kount = 0;
		for( int i=0; i<all.parms.nalpha; i++){
			bootParms[kount] = bootalpha[MATREF(b,i,nboot)];
			kount++;
		}
		for( int i=0; i<all.parms.ntau; i++){
			bootParms[kount] = boottau[MATREF(b,i,nboot)];
			kount++;
		}
		for( int i=0; i<all.parms.nbeta; i++){
			bootParms[kount] = bootbeta[MATREF(b,i,nboot)];
			kount++;
		}
		all.parms.update( bootParms, all.data);

		//calculate pis and store them
		int place;
		for( int i=0; i<all.data.nObs; i++)
			calcLogPis( logPis, allPis.at(i), all.data, all.parms, i);
		for( int g=0; g<all.data.nG; g++)
			for( int i=0; i<all.data.nObs; i++){
				place = MATREF3D(i,g,b,all.data.nObs, all.data.nG);//MATREF( MATREF(i,g,all.data.nObs),b,(all.data.nObs*all.data.nG));
				bootPreds[place] = allPis.at(i).at(g);
			}
	}

	SEXP Rres;	//R object to return -- it is meaningless!
	Rres = allocVector(REALSXP,1);
	double *R_pointer = REAL(Rres);
	R_pointer[0] = -9999;

	return( Rres);


}

}

/*void calcMargFits( double *ptPreds, int bootCount, const vector<double> &allMus, const vector< vector<double> > &allPis, const myData &dat)
{
	for( int i=0; i<dat.nObs; i++)
		for( int s=0; s<dat.nS; s++){
			ptPreds[MATREF(MATREF(i,s,dat.nObs),bootCount,(dat.nObs*dat.nS))] = 0.0;
			for( int g=0; g<dat.nG; g++)
				ptPreds[MATREF(MATREF(i,s,dat.nObs),bootCount,(dat.nObs*dat.nS))] += allMus.at(MATREF(g,s,dat.nG))*allPis.at(i).at(g);
		}

}*/
