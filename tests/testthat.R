Sys.setenv("R_TESTS" = "")
library(testthat)
library(ecomix)
library(raster)
library(scales)

test_check("ecomix")
