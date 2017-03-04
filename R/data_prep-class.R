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

