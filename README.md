
<!-- README.md is generated from README.Rmd. Please edit that file -->

# CRV3J

<!-- badges: start -->

[![R-CMD-check](https://github.com/s3alfisc/clusterjack/workflows/R-CMD-check/badge.svg)](https://github.com/s3alfisc/clusterjack/actions)
<!-- badges: end -->

`CRV3J` implements the CRV3 Jackknife algorithm proposed in [MacKinnon,
Nielsen and
Webb](https://scholar.google.de/citations?view_op=view_citation&hl=de&user=PdkdfhMAAAAJ&sortby=pubdate&citation_for_view=PdkdfhMAAAAJ:8VtEwCQfWZkC).

## Installation

To install the package, run

``` r
devtools::install_github("s3alfisc/CRV3J")
```

## Example

``` r
library(CRV3J)
library(fwildclusterboot)
library(fixest)

set.seed(98765)
# few large clusters (around 10000 obs)
N <- 100000
N_G1 <-100
data <- fwildclusterboot:::create_data(
    N = N,
    N_G1 = N_G1,
    icc1 = 0.8,
    N_G2 = 10,
    icc2 = 0.8,
    numb_fe1 = 10,
    numb_fe2 = 10,
    seed = 12
  )

feols_fit <- feols(proposition_vote ~ treatment  + log_income |Q1_immigration + Q2_defense, cluster = ~group_id1 , data = data)
```

Calculate Jackknife CVR3 Variance-Covariance Matrix

``` r
vcov <- CRV3J::vcov_CR3J(
  obj = feols_fit, 
  cluster = data$group_id1
)
```

Results can now be collected via the `coeftest` package:

``` r
library(sandwich)
library(lmtest)

coeftest(feols_fit, vcov)
#> 
#> z test of coefficients:
#> 
#>              Estimate Std. Error  z value  Pr(>|z|)    
#> treatment   0.0033009  0.0012081   2.7324  0.006287 ** 
#> log_income -0.0120565  0.0005545 -21.7432 < 2.2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```
