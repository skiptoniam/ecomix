#include"sam_ippm.h"

sam_ippm_fits::sam_ippm_fits(){};
sam_ippm_fits::~sam_ippm_fits(){};

void sam_ippm_fits::initialise( const int &nObs, const int &nG, const int &nS, const int &nP, const int &NAnum){
	
	// vector for of pis on natural scale
	//pis_loglike.resize(nG-1, 0);
	
	// vector for of pis on natural scale
	par_pis.resize(nG-1, NAnum);
	
	// vector for of pis on addiative logitic scale
	par_etas.resize(nG-1, NAnum);
	
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

void sam_ippm_fits::zero(const int &NAnum){
	

	//pis_loglike.assign(pis_loglike.size(),NAnum); 
	par_pis.assign(par_pis.size(),NAnum); 
	par_etas.assign(par_etas.size(),NAnum); 
	//grad_pis_non_trans.assign(grad_pis_non_trans.size(),NAnum); 
	allMus.assign(allMus.size(), NAnum);
	log_like_species_contrib.assign(log_like_species_contrib.size(), NAnum);
	log_like_species_group_contrib.assign(log_like_species_group_contrib.size(), NAnum);
	dlogdalpha.assign(dlogdalpha.size(), NAnum);
	dlogdbeta.assign(dlogdbeta.size(), NAnum);
	dlogdpi.assign(dlogdpi.size(), NAnum);
			
};

