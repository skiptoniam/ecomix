#include"sam_bernoulli_sp_ints.h"

sam_bernoulli_sp_ints_fits::sam_bernoulli_sp_ints_fits(){};
sam_bernoulli_sp_ints_fits::~sam_bernoulli_sp_ints_fits(){};

void sam_bernoulli_sp_ints_fits::initialise( const int &nObs, const int &nG, const int &nS, const int &nP, const int &NAnum){
	
	// array for catching Mu	
	allMus.resize(nObs*nS*nG, NAnum);
	
	// array for catching species loglike contribution to the model.
	log_like_species_contrib.resize(nS, NAnum);

	// array for catching species/group loglike contribution to the model.
	log_like_species_group_contrib.resize(nS*nG, NAnum);
	
	// array for catching for dlogdalpha
	dlogdalpha.resize(nS*nG, NAnum);
	
	// array for catching for dlogdbeta
	dlogdbeta.resize(nG*nP*nS, NAnum); // This is actually nG*nP*nSP
	
	// array for catching for dlogdpi
	dlogdpi.resize(nG, NAnum);
	
};

void sam_bernoulli_sp_ints_fits::zero(const int &NAnum){
	
	allMus.assign(allMus.size(), NAnum);
	log_like_species_contrib.assign(log_like_species_contrib.size(), NAnum);
	log_like_species_group_contrib.assign(log_like_species_group_contrib.size(), NAnum);
	dlogdalpha.assign(dlogdalpha.size(), NAnum);
	dlogdbeta.assign(dlogdbeta.size(), NAnum);
	dlogdpi.assign(dlogdpi.size(), NAnum);
			
};
