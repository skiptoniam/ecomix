context('data_preperation-class')
library(ecomix)

test_that('data_preperation classes work', {

  testthat::expect_error(make_mixture_data(NA,NA))

})
