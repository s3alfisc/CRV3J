test_that("CV3J equal to HC3 when all clusters are singletons", {
  
  set.seed(20190513)
  library(fixest)
  library(sandwich)
  
  m <- 8
  cluster <- factor(rep(LETTERS[1:m], 3 + rpois(m, 5)))
  n <- length(cluster)
  const <- rep("whuzzup.", n)
  X <- matrix(rnorm(3 * n), n, 3)
  nu <- rnorm(m)[cluster]
  e <- rnorm(n)
  w <- rgamma(n, shape = 3, scale = 3)
  y <- X %*% c(.4, .3, -.3) + nu + e
  
  dat <- data.frame(y, X, cluster, const, w, row = 1:n)
  
  lm_fit <- lm(y ~ X1 + X2 + X3 , data = dat)
  WLS_fit <- lm(y ~ X1 + X2 + X3 , data = dat, weights = w)
  
  CR3f <- vcov_CR3J(lm_fit, cluster = dat$row)
  expect_equal(sandwich::vcovHC(lm_fit, type = "HC3"), nobs(lm_fit) / (nobs(lm_fit) - 1) *as.matrix(CR3f), ignore_attr = TRUE)
  CR3f <- vcov_CR3J(WLS_fit, cluster = dat$row)
  expect_equal(vcovHC(WLS_fit, type = "HC3"), nobs(WLS_fit) / (nobs(WLS_fit) - 1) *as.matrix(CR3f), ignore_attr = TRUE)

  feols_fit <- feols(y ~ X1 + X2 + X3 , data = dat)
  feols_WLS_fit <- feols(y ~ X1 + X2 + X3 , data = dat, weights = w)
  
  CR3f <- vcov_CR3J(feols_fit, cluster = dat$row)
  expect_equal(vcovHC(feols_fit, type = "HC3"), nobs(feols_fit) / (nobs(feols_fit) - 1) *as.matrix(CR3f), ignore_attr = TRUE)
  
  # CR3f <- vcov_CR3J(feols_WLS_fit, cluster = dat$row)
  # expect_equal(vcovHC(feols_WLS_fit, type = "HC3"), nobs(feols_WLS_fit) / (nobs(feols_WLS_fit) - 1) *as.matrix(CR3f), ignore_attr = TRUE)
  
})
