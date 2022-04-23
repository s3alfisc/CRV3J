
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
#> elapsed time is 0.320000 seconds
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
#> elapsed time is 156.190000 seconds
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
#> treatment  -0.000289118  0.004948919
#> log_income -0.013038535 -0.010844641

cat("vcovCR3", "\n")
#> vcovCR3
lmtest::coefci(lm_fit, parm = c("treatment", "log_income"), vcov = vcovCR3)
#>                    2.5 %       97.5 %
#> treatment  -0.0003023121  0.004962113
#> log_income -0.0130440609 -0.010839115

cat("vcovHC3", "\n")
#> vcovHC3
lmtest::coefci(lm_fit, parm = c("treatment", "log_income"), vcov = vcovHC)
#>                    2.5 %       97.5 %
#> treatment  -0.0003999807  0.005059782
#> log_income -0.0125171348 -0.011366041
```

Regression table:

``` r
vc <- list(vcovHC = vcovHC, vcovCR3 = vcovCR3, vcovJN = vcovJN)
mlist <- lapply(vc, coeftest, x = lm_fit)
msummary(mlist)
```

|             |  vcovHC  | vcovCR3  |  vcovJN  |
|:------------|:--------:|:--------:|:--------:|
| (Intercept) |  1.087   |  1.087   |  1.087   |
|             | (0.003)  | (0.006)  | (0.006)  |
| treatment   |  0.002   |  0.002   |  0.002   |
|             | (0.001)  | (0.001)  | (0.001)  |
| log_income  |  -0.012  |  -0.012  |  -0.012  |
|             | (0.000)  | (0.001)  | (0.001)  |
| Num.Obs.    |  100000  |  100000  |  100000  |
| AIC         | -18838.4 | -18838.4 | -18838.4 |
| BIC         | -18800.3 | -18800.3 | -18800.3 |
| Log.Lik.    | 9423.188 | 9423.188 | 9423.188 |
