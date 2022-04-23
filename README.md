
<!-- README.md is generated from README.Rmd. Please edit that file -->

# clusterjack

<!-- badges: start -->
<!-- badges: end -->

`clusterjack` implements the CRV3 Jackknife algorithm proposed in
[MacKinnon, Nielsen and
Webb](https://scholar.google.de/citations?view_op=view_citation&hl=de&user=PdkdfhMAAAAJ&sortby=pubdate&citation_for_view=PdkdfhMAAAAJ:8VtEwCQfWZkC).

## Installation

To install the package, run

``` r
devtools::install_github("s3alfisc/clusterjack")
```

## Example

``` r
library(clusterjack)
library(fwildclusterboot)
library(clubSandwich)
library(lmtest)
library(pracma)
library(modelsummary)

set.seed(98765)
# few large clusters (around 10000 obs)
N <- 10000
N_G1 <-10
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
lm_fit <- lm(proposition_vote ~ treatment  + log_income , data = data)

if(N == N_G1){
  data$group_id1 <- 1:N
}
```

Calculate Jackknife CVR3 Variance-Covariance Matrix

``` r
pracma::tic()
vcovJN <- clusterjack::vcovJN(
  model = lm_fit, 
  clustid = data$group_id1
)
pracma::toc()
#> elapsed time is 0.020000 seconds
```

Calculate

``` r
pracma::tic()
vcovCR3 <- 
clubSandwich::vcovCR(
  lm_fit, 
  cluster = data$group_id1, 
  type = "CR3")
pracma::toc()
#> elapsed time is 11.030000 seconds
```

``` r
library(sandwich)

vcovHC <- sandwich::vcovHC(
  lm_fit, 
  type = "HC3"
)
```

Confidence Intervals:

``` r
cat("vcovJN", "\n")
#> vcovJN
lmtest::coefci(lm_fit, parm = c("treatment", "log_income"), vcov = vcovJN)
#>                   2.5 %       97.5 %
#> treatment  -0.001325088  0.019158044
#> log_income -0.011346096 -0.005543023

cat("vcovCR3", "\n")
#> vcovCR3
lmtest::coefci(lm_fit, parm = c("treatment", "log_income"), vcov = vcovCR3)
#>                  2.5 %       97.5 %
#> treatment  -0.00187908  0.019712037
#> log_income -0.01150305 -0.005386071

cat("vcovHC3", "\n")
#> vcovHC3
lmtest::coefci(lm_fit, parm = c("treatment", "log_income"), vcov = vcovHC)
#>                    2.5 %       97.5 %
#> treatment   0.0009465414  0.016886415
#> log_income -0.0102328076 -0.006656311
```

Hypothesis test:

``` r
vc <- list(vcovHC = vcovHC, vcovCR3 = vcovCR3, vcovJN = vcovJN)
mlist <- lapply(vc, coeftest, x = lm_fit)
msummary(mlist)
```

|             |  vcovHC  | vcovCR3  |  vcovJN  |
|:------------|:--------:|:--------:|:--------:|
| (Intercept) |  1.055   |  1.055   |  1.055   |
|             | (0.011)  | (0.017)  | (0.016)  |
| treatment   |  0.009   |  0.009   |  0.009   |
|             | (0.004)  | (0.006)  | (0.005)  |
| log_income  |  -0.008  |  -0.008  |  -0.008  |
|             | (0.001)  | (0.002)  | (0.001)  |
| Num.Obs.    |  10000   |  10000   |  10000   |
| AIC         | -3480.2  | -3480.2  | -3480.2  |
| BIC         | -3451.4  | -3451.4  | -3451.4  |
| Log.Lik.    | 1744.119 | 1744.119 | 1744.119 |
