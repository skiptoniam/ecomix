#include"sam_cpp2.h"

sam_fits::sam_fits(){};
sam_fits::~sam_fits(){};

void sam_fits::initialise( const int &nObs, const int &nG, const int &nS, const int &nPX, const int &nPW, const int &nPU, const int &NAnum){
	
	// array for catching Mu	
	allMus.resize(nObs*nS*nG, NAnum);
	
	// array for catching species loglike contribution to the model.
	log_like_species_contrib.resize(nS, NAnum);

	// array for catching species/group loglike contribution to the model.
	log_like_species_group_contrib.resize(nS*nG, NAnum);
	
	// array for catching dMu	
	all_derivs_mu.resize(nObs*nS*nG, NAnum);
	
	// array for catching for dlogdalpha
	dlogdalpha.resize(nS*nG, NAnum);
	
	// array for catching for dlogdbeta
	dlogdbeta.resize(nG*nPX*nS, NAnum); // This is actually nG*nP*nSP

	// array for catching dlogdgamma
	dlogdgamma.resize(nG*nPW*nS, NAnum); // This is actually nG*nP*nSP
	
	// array for catching dlogdgamma
	dlogddelta.resize(nPU, NAnum); // This is actually nPU
	
	// array for catching for dlogdpi
	dlogdpi.resize(nG, NAnum);
	
	// array for catching for dlogdtheta
	dlogdtheta.resize(nS*nG, NAnum);
	
};

void sam_fits::zero(const int &NAnum){
	
	allMus.assign(allMus.size(), NAnum);
	log_like_species_contrib.assign(log_like_species_contrib.size(), NAnum);
	log_like_species_group_contrib.assign(log_like_species_group_contrib.size(), NAnum);
	all_derivs_mu.assign(all_derivs_mu.size(), NAnum);
	dlogdalpha.assign(dlogdalpha.size(), NAnum);
	dlogdbeta.assign(dlogdbeta.size(), NAnum);
	dlogdgamma.assign(dlogdgamma.size(), NAnum);
	dlogdpi.assign(dlogdpi.size(), NAnum);
	dlogdtheta.assign(dlogdtheta.size(), NAnum);
			
};

