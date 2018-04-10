#include"sam_ippm.h"

sam_ippm_fits::sam_ippm_fits(){};
sam_ippm_fits::~sam_ippm_fits(){};

void sam_ippm_fits::initialise( const int &nObs, const int &nG, const int &nS, const int &nP, const int &NAnum){
	
	// vector for of pis on natural scale
	pis_loglike.resize(nG-1, 0);
	
	// vector for of pis on addiative logitic scale
	//pis_grads.resize(nG-1, 0);
	
	// vector for of pis on natural scale
	//grad_pis_trans.resize(nG-1, 0);
	
	// vector for of pis on addiative logitic scale
	//grad_pis_non_trans.resize(nG-1, 0);
	
	// array for catching Mu	
	allMus.resize(nObs*nS*nG, 0);
	
	// array for catching species loglike contribution to the model.
	log_like_species_contrib.resize(nS,0);

	// array for catching species/group loglike contribution to the model.
	log_like_species_group_contrib.resize(nS*nG, 0);
	
	// array for catching for dlogdalpha
	dlogdalpha.resize(nS*nG, 0);
	
	// array for catching for dlogdbeta
	dlogdbeta.resize(nG*nP*nS, 0); // This is actually nG*nP*nSP
	
	// array for catching for dlogdpi
	dlogdpi.resize(nG, 0);
	
};

void sam_ippm_fits::zero(const int &NAnum){
	

	pis_loglike.assign(pis_loglike.size(),NAnum); 
	//pis_grads.assign(pis_grads.size(),NAnum); 
	//grad_pis_trans.assign(grad_pis_trans.size(),NAnum); 
	//grad_pis_non_trans.assign(grad_pis_non_trans.size(),NAnum); 
	allMus.assign(allMus.size(), NAnum);
	log_like_species_contrib.assign(log_like_species_contrib.size(), NAnum);
	log_like_species_group_contrib.assign(log_like_species_group_contrib.size(), NAnum);
	dlogdalpha.assign(dlogdalpha.size(), NAnum);
	dlogdbeta.assign(dlogdbeta.size(), NAnum);
	dlogdpi.assign(dlogdpi.size(), NAnum);
			
};

