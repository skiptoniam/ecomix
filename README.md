ecomix: Finite mixture models for multiple species grouping of ecological data <img src="man/figures/logo.png" align="right" alt="" width="160" />
==================================================================================================================================================

[![Build
Status](https://travis-ci.org/skiptoniam/ecomix.svg?branch=master)](https://travis-ci.org/skiptoniam/ecomix.svg?branch=master)
[![Coverage
Status](https://img.shields.io/codecov/c/github/skiptoniam/ecomix/master.svg)](https://codecov.io/github/skiptoniam/ecomix?branch=master)
[![Project Status: Active ??? The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)

ecomix is a package to implement statistical models on multivariate
species data. Two main type of finite mixture models are available which
group multivariate data at the species level (Species Archetype Models;
using the `species_mix` functions) or site level (Region of Common
Profiles; using the `regional_mix` functions).

#### Installation

<!-- `ecomix` version 1.0.0 can be install from cran using `install.package("ecomix")`. -->

The development version of `ecomix` can be installed from GitHub using
the `devtools` package:

``` r
devtools::install_github('skiptoniam/ecomix')
library(ecomix)
```

An example vignette on how to run and interpret a Species Archetype
Models (species\_mix) is provided within the package or can be viewed at
<a href="https://skiptoniam.github.io/ecomix/articles/SAM-example.html" class="uri">https://skiptoniam.github.io/ecomix/articles/SAM-example.html</a>

An example vignette on how to run and interpret a Region of Common
Profiles model (regional\_mix) is provided with the package or can be
viewed at
<a href="https://skiptoniam.github.io/ecomix/articles/RCP-example.html" class="uri">https://skiptoniam.github.io/ecomix/articles/RCP-example.html</a>
