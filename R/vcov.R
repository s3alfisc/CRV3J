vcov_CR3J.lm <- function(obj, cluster, ...) {
  
  #' Compute CR3 Jackknive variance covariance matrices of objects of type lm
  #' @param object An object of type lm
  #' @param cluster A clustering vector
  #' @param ... other function arguments passed to 'vcov'
  #' @export

  cluster <- droplevels(as.factor(cluster))
  
  alias <- is.na(coef_CS(obj))
  X <- model_matrix(obj)
  if (any(alias)) {
    X <- X[, !alias, drop = FALSE]
  }  
  
  p <- NCOL(X)
  N <- NROW(X)
  
  cluster_length <- length(cluster)
  
  if (cluster_length != N) {
    
    cluster <- droplevels(handle_vectors(cluster, obj))
    
    if (length(cluster) != N) {
      stop("Clustering variable must have length equal to the number of rows in the data used to fit obj.")
    }
    
  } 
  
  if (any(is.na(cluster))) stop("Clustering variable cannot have missing values.")
  
  J <- nlevels(cluster)
  if (J < 2) stop("Cluster-robust variance estimation will not work when the data only includes a single cluster.")
  
  # get list of regression residuals, u_g in MNW notation
  resid <- residuals_CS(obj)
  res_list <- split(resid, cluster)
  
  y <- model.response(model.frame(obj))
  y_list <- split(y, cluster)
  
  # X_g's in MNW notation
  X_list <- matrix_list(X, cluster, "row")
  # W_g's in MNW notation
  W_list <- weightMatrix(obj, cluster)
  
  # list of weight adjusted X's
  XW_list <- Map(function(x, w) as.matrix(t(x) %*% sqrt(w)), x = X_list, w = W_list)
  yW_list <- Map(function(y, w) as.matrix(t(y) %*% sqrt(w)), y = y_list, w = W_list)
 
  # get small sample adjustments
  adjustments <- 1
  # multiply design matrices XW by sqrt(small_sample_adjustments)
  E_list <- lapply(XW_list, function(e) e * adjustments)
  
  XWg_XWg <- Map(function(e) tcrossprod(e), e = E_list)
  yWg_XWg <-  Map(function(e,y) e %*% y, e = E_list, y = yW_list)
  tXX <- Reduce("+", XWg_XWg)
  tXy <- Reduce("+", yWg_XWg)
  # beta_g - beta
  coef_list <- Map(function(a = tXX, b, c = tXy, d, beta = coef(obj)) MASS::ginv(a - b) %*% (c - d) - beta, b = XWg_XWg, d = yWg_XWg)
  vcov <- Reduce("+",Map(function(x) tcrossprod(x), x = coef_list) )
    
  vcov <- vcov * (J-1)/J
    
  v_scale <- v_scale(obj)
  w_scale <- attr(W_list, "w_scale")
  if (is.null(w_scale)) w_scale <- 1L
    
    
  # instead of sandwich::bread(obj), faster as tXX already computed
  # bread <- N * solve(tXX)
  bread <- sandwich::bread(obj)
    
  rownames(vcov) <- colnames(vcov) <- colnames(X)
  attr(vcov, "cluster") <- cluster
  attr(vcov, "bread") <- bread
  attr(vcov, "v_scale") <- v_scale
  attr(vcov, "est_mats") <- 
    # required because CR3f computes t(X) %*% sqrt(w)
           Map(function(x, w) as.matrix(t(x) %*% w), x = X_list, w = W_list)
  attr(vcov, "adjustments") <- adjustments
  class(vcov) <- c("vcovCJ","clubSandwich")
  return(vcov)
  
}


vcov_CR3J.fixest <- function(obj, cluster, ...) {
  
  #' Compute CR3 Jackknive variance covariance matrices of objects of type fixest
  #' @param object An object of type lm
  #' @param cluster A clustering vector
  #' @param ... other function arguments passed to 'vcov'
  #' @export
  #' @importFrom fixest demean
  
  cluster <- droplevels(as.factor(cluster))
  
  has_fe <- length(obj$fixef_vars) > 0
  
  X <- model.matrix(obj, type = "rhs")
  if(has_fe){
    fe <- model.matrix(obj, type = "fixef")
    X <- fixest::demean(X, fe)
  }

  p <- NCOL(X)
  N <- NROW(X)
  
  cluster_length <- length(cluster)
  
  if (cluster_length != N) {
    
    cluster <- droplevels(handle_vectors(cluster, obj))
    
    if (length(cluster) != N) {
      stop("Clustering variable must have length equal to the number of rows in the data used to fit obj.")
    }
    
  } 
  
  if (any(is.na(cluster))) stop("Clustering variable cannot have missing values.")
  
  J <- nlevels(cluster)
  if (J < 2) stop("Cluster-robust variance estimation will not work when the data only includes a single cluster.")
  
  # get list of regression residuals, u_g in MNW notation
  resid <- resid(obj)
  res_list <- split(resid, cluster)
  
  y <- model.matrix(obj, type = "lhs")
  y_list <- split(y, cluster)
  
  # X_g's in MNW notation
  X_list <- matrix_list(X, cluster, "row")
  # W_g's in MNW notation
  W_list <- weightMatrix(obj, cluster)
  
  # list of weight adjusted X's
  XW_list <- Map(function(x, w) as.matrix(t(x) %*% sqrt(w)), x = X_list, w = W_list)
  yW_list <- Map(function(y, w) as.matrix(t(y) %*% sqrt(w)), y = y_list, w = W_list)
  
  # get small sample adjustments
  adjustments <- 1
  # multiply design matrices XW by sqrt(small_sample_adjustments)
  E_list <- lapply(XW_list, function(e) e * adjustments)
  
  XWg_XWg <- Map(function(e) tcrossprod(e), e = E_list)
  yWg_XWg <-  Map(function(e,y) e %*% t(y), e = E_list, y = yW_list)
  tXX <- Reduce("+", XWg_XWg)
  tXy <- Reduce("+", yWg_XWg)
  # beta_g - beta
  coef_list <- Map(function(a = tXX, b, c = tXy, d, beta = coef(obj)) MASS::ginv(a - b) %*% (c - d) - beta, b = XWg_XWg, d = yWg_XWg)
  vcov <- Reduce("+",Map(function(x) tcrossprod(x), x = coef_list) )
  
  vcov <- vcov * (J-1)/J
  
  v_scale <- v_scale(obj)
  w_scale <- attr(W_list, "w_scale")
  if (is.null(w_scale)) w_scale <- 1L
  
  
  # instead of sandwich::bread(obj), faster as tXX already computed
  # bread <- N * solve(tXX)
  bread <- sandwich::bread(obj)
  
  rownames(vcov) <- colnames(vcov) <- colnames(X)
  attr(vcov, "cluster") <- cluster
  attr(vcov, "bread") <- bread
  attr(vcov, "v_scale") <- v_scale
  attr(vcov, "est_mats") <- 
    # required because CR3f computes t(X) %*% sqrt(w)
    Map(function(x, w) as.matrix(t(x) %*% w), x = X_list, w = W_list)
  attr(vcov, "adjustments") <- adjustments
  class(vcov) <- c("vcovCJ","clubSandwich")
  return(vcov)
  
}


vcov_CR3J <- function(object, ...) {
  
  #' Compute CR3 Jackknive variance covariance matrices of objects of type lm and fixest
  #' @param object An object of class `lm` or `fixest``
  #' @param ... Other arguments
  #' @export
  
  UseMethod("vcov_CR3J")
  
}


