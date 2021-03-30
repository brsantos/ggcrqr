
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ggcrqr

<!-- badges: start -->

<!-- badges: end -->

The goal of ggcrqr is to …

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("brsantos/ggcrqr")
```

## Example

This is a basic example which shows you how to estimate the model
considering the data available in the package: `m_breast_cancer`.

``` r
library(ggcrqr)
#> Loading required package: LaplacesDemon
#> Loading required package: truncnorm
## basic example code
## one should check whether the values for burn and jump are adequate.
model <- bayesGG(time_to_d ~ age_group + stage_c, 
                 data = m_breast_cancer, 
                 q = 0.5, d = "cens", burn = 30000, jump = 40, 
                 guess = c(-0.1, 0.5, rep(0, 6)))
#> 
#> Laplace's Demon was called on Tue Mar 30 10:02:59 2021
#> 
#> Performing initial checks...
#> Algorithm: Adaptive Metropolis 
#> 
#> Laplace's Demon is beginning to update...
#> Iteration: 20000,   Proposal: Multivariate,   LP: -715
#> Iteration: 40000,   Proposal: Multivariate,   LP: -718
#> Iteration: 60000,   Proposal: Multivariate,   LP: -716.5
#> 
#> Assessing Stationarity
#> Assessing Thinning and ESS
#> Creating Summaries
#> Creating Output
#> 
#> Laplace's Demon has finished.
```
