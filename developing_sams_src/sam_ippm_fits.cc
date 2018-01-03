#include"sam.h"

sam_ippm_fits::sam_ippm_fits(){};
sam_ippm_fits::~sam_ippm_fits(){};

void sam_ippm_fits::initialise( const int &nObs, const int &nG, const int &nS, const int &NAnum){
	
	allPis.resize(nObs);
	for( int i=0; i<nObs; i++)
		allPis.at(i).resize(nG, NAnum);
		
	allMus.resize(nObs*nS*nG, NAnum);
	
	allLogDens.resize(nObs);
	for( int i=0; i<nObs; i++)
		allLogDens.at(i).resize(nG, NAnum);
	
	allLogls.resize( nObs, NAnum);
};


void sam_ippm_fits::zero(const int &NAnum){
	for( int i=0; i<(int)allPis.size(); i++)
		allPis.at(i).assign(allPis.at(i).size(),NAnum); 
	allMus.assign( allMus.size(), NAnum);
	for( int i=0; i<(int)allLogDens.size(); i++)
		allLogDens.at(i).assign( allLogDens.at(i).size(), NAnum);
	allLogls.assign( allLogls.size(), NAnum);			
};
