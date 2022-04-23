vcovJN <- function(model, clustid, type = "CRV3", ...){

  #' Implements the fast algorithm for CRV3 Jackknife estimation in MacKinnon, Nielsen & Webb (2022)
  #' @param model An object of type lm
  #' @param clustid A clustering vector
  #' @param type 'CRV3' or 'CRVJ3'
  #' @param ... other function arguments passed to 'vcov'
  #' @export
  #' @importFrom stats coefficients model.frame model.matrix model.response

  if(!inherits(model, "lm")){
    stop("clusterjack currently only works for models of class lm.")
  }

  fml <- model$call$formula
  model_frame <- model.frame(fml, eval(model$call$data))

  X <- model.matrix(model, model_frame)
  y <- model.response(model_frame)
  beta_hat <- coefficients(model)

  unique_clusters <- unique(clustid)
  G <- length(unique_clusters)
  small_sample_correction <- (G-1)/G


  #calculate X_g'X_g
  tXgXg <- lapply(
    seq_along(unique_clusters),
    function(x) crossprod(X[clustid == x,])
  )
  tXX <- Reduce("+", tXgXg)
  #all.equal(tXX, crossprod(X))

  tXgyg <- lapply(
    seq_along(unique_clusters),
    function(x)
      t(X[clustid == x,]) %*% y[clustid == x]
  )
  tXy <- Reduce("+", tXgyg)
  #all.equal(tXy, t(X) %*% y)
  #mean(solve(tXX) %*% tXy - coef(model))
  # initiate jackknife

  beta_jack <-
    lapply(
      seq_along(unique_clusters),
      function(x){
        solve(tXX - tXgXg[[x]]) %*% (tXy - (t(X[clustid == x,]) %*% y[clustid == x]))
    })

  #solve(tXX - tXgXg[[x]]) %*% (tXy - (t(X[clustid == x,]) %*% y[clustid ==x]))
  #lm.fit(y = y[clustid != x], x = X[clustid != x,])$coef

  if(type == "CRV3J"){
    beta_bar <- beta_center <- Reduce("+", beta_jack) / G
  } else if(type == "CRV3"){
    beta_center <- beta_hat
  }

  V3 <- lapply(
    seq_along(unique_clusters),
    function(x)
       tcrossprod(beta_jack[[x]] - beta_center)
  )

  vcov <- Reduce("+", V3) * small_sample_correction

  vcov
}
