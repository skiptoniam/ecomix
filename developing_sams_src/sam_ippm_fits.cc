#include"sam_ippm.h"

sam_ippm_fits::sam_ippm_fits(){};
sam_ippm_fits::~sam_ippm_fits(){};

void sam_ippm_fits::initialise( const int &nObs, const int &nG, const int &nS, const int &NAnum){
	
	// vector for storing pis
	allPis.resize(nG, NAnum);

	// vector for storing pis
	estpi.resize(nG-1, NAnum);
	
	// array for catching Mu	
	allMus.resize(nObs*nS*nG, NAnum);
	
	// array for catching species loglike contribution to the model.
	log_like_species_contrib.resize(nS, NAnum);

	// array for catching species/group loglike contribution to the model.
	log_like_species_group_contrib.resize(nS*nG, NAnum);
	
};

void sam_ippm_fits::zero(const int &NAnum){
	
	allPis.assign(allPis.size(),NAnum); 
	estpi.assign(estpi.size(),NAnum); 
	allMus.assign(allMus.size(), NAnum);
	log_like_species_contrib.assign(log_like_species_contrib.size(), NAnum);
	log_like_species_group_contrib.assign(log_like_species_group_contrib.size(), NAnum);
			
};
