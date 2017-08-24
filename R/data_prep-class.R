#' @title data_prep objects
#' @rdname data_prep
#' @name table_to_species_data
#' @description prepares table of observational records into a sites x species matrix for use in \code{ecomix}.
#' @param observation_data data to be imported into R. A table of species occurrences and site records.
#' @param site_id Name of column that contains site names /IDs.
#' @param species_id Name of column that contains species information/ID.
#' @param measurement_id Name of column that contains species observational data. eg, 'occurrence', 'abundance' or 'biomass'.
#' @return A sites x species matrix to used as species_data in \code{ecomix} models.
#' @export

table_to_species_data <- function(observation_data, site_id = "site_id", species_id = "species_id", measurement_id = NULL) {

    nr <- length(levels(as.factor(observation_data[, site_id])))
    rn <- levels(as.factor(observation_data[, site_id]))

    cn <- levels(as.factor(observation_data[, species_id]))
    nc <- length(cn)
    nm <- matrix(0, nr, nc, dimnames = list(rn, cn))
    for (i in 1:length(observation_data[, 1])) {
        m <- as.character(observation_data[i, site_id])
        n <- as.character(observation_data[i, species_id])
        if (is.na(m) == TRUE | is.null(m) == TRUE | is.na(n) == TRUE | is.null(n) == TRUE){
            (next)(i)
        }
        if (m == "" | m == " " | n == "" | n == " "){
            (next)(i)
        }
        if (!is.null(measurement_id)){
            nm[m, n] <- nm[m, n] + observation_data[i, measurement_id]
            } else {
            nm[m, n] <- 1
            }
    }
    fm <- nm[rowSums(nm) > 0, ]
        return(as.matrix(fm))
}

#' @rdname data_prep
#' @name make_mixture_data
#' @param species_data A character string providing the name of the column in \verb{data} that contains the counts
#' @param covariate_data A character string providing the name of the column in \verb{data} that contains the binary occurence data
# #' @param weights A formula decribing the variables for which delta will be the coefficients.
# #' @param offset A formula giveng additional variables to be used for all groups.
#'
#' @return
#' A list containing the following elements:
#' \item{species_data}{A matrix of occurrence, counts, biomass, or...}
#' \item{covariate_data}{The design matrix for the global covariates}
# #' \item{weights}{An optional vector of 'prior weights' to be used in the fitting process}
# #' \item{offset}{An optinal vector that can be used to specify an a priori known component to be included in the linear predictor during fitting. }
#'
#' @author Skipton Woolley
#' @export

make_mixture_data <- function(species_data,covariate_data){

  #check species data
  species_data_check(species_data)

  #check covariate data
  covariate_data_check(covariate_data)

  if(!identical(dim(species_data)[1],dim(covariate_data)[1]))
    stop('dimensions of species matrix sites (rows) and covariate data at sites (rows) do not match, have another look')

  out <- list('species_data'=species_data,'covariate_data'=covariate_data)
  return(out)
}

species_data_check <- function(x){
    stopifnot(is.matrix(x)|is.data.frame(x))
    stopifnot(all(is.finite(x)))
}

covariate_data_check <- function(x){
  stopifnot(is.matrix(x)|is.data.frame(x))
  stopifnot(all(is.finite(x)))
}

# make_ppm_weights <- function()
