context('data_preperation-class')
library(ecomix)

test_that('data_preperation classes work', {

  ## setup the species data
  spp <- paste0("spp_",rep(letters[1:10],each=10),rep(letters[1:10],10))
  spp_abund <- rpois(100,5)
  site <- paste0("site_",rep(letters[1:5],each=20),rep(letters[1:5],20))
  sptable <- data.frame(spp=spp,site=site,spp_abund=spp_abund)
  mat_pa <- table_to_species_data(sptable,site_id = "site",species_id = "spp")
  mat_abund <- table_to_species_data(sptable,site_id = "site",species_id = "spp",measurement_id = "spp_abund")
  testthat::expect_is(mat_pa,"matrix")
  testthat::expect_is(mat_abund,"matrix")


  ## make a mixture model dataframe
  X <- data.frame(x1=runif(25),x2=rnorm(25))
  dat <- make_mixture_data(mat_pa,X)
  testthat::expect_is(dat,"data.frame")

  ## expect an error if the data is miss matched
  X <- data.frame(x1=runif(50),x2=rnorm(50))
  testthat::expect_error(make_mixture_data(mat_pa,X))

  ## expect error if na passed.
  testthat::expect_error(make_mixture_data(NA,NA))

})
