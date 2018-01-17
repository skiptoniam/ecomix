#include"sam.h"

sam_ippm_fits::sam_ippm_fits(){};
sam_ippm_fits::~sam_ippm_fits(){};

void sam_ippm_fits::initialise( const int &nObs, const int &nG, const int &nS, const int &NAnum){
	
	// vector for storing pis
	allPis.resize(nG, NAnum);
	
	// array for catching Mu	
	allMus.resize(nObs*nS*nG, NAnum);
	
	// array for catching species loglike contribution to the model.
	allSum_f_species.resize(nG*nS, NAnum);
	
};

void sam_ippm_fits::zero(const int &NAnum){
	
	allPis.assign(allPis.size(),NAnum); 
	allMus.assign( allMus.size(), NAnum);
	allSum_f_species.assign(allSpecies_loglike_contrib.size(), NAnum);			
};
