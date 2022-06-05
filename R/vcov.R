vcov_CR3J.lm <- function(obj, cluster, ...) {
  
  #' Compute CR3 Jackknive variance covariance matrices of objects of type lm
  #' @param obj An object of type lm
  #' @param cluster A clustering vector
  #' @param ... other function arguments passed to 'vcov'
  #' @importFrom stats coef model.frame model.response model.matrix weights
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
    
  rownames(vcov) <- colnames(vcov) <- colnames(X)
  attr(vcov, "cluster") <- cluster
  # required because CR3f computes t(X) %*% sqrt(w)
  class(vcov) <- c("vcov_CV3J")
  
  return(vcov)
  
}


vcov_CR3J.fixest <- function(obj, cluster, ...) {
  
  #' Compute CR3 Jackknive variance covariance matrices of objects of type fixest
  #' @param obj An object of type lm
  #' @param cluster A clustering vector
  #' @param ... other function arguments passed to 'vcov'
  #' @export

  cluster <- droplevels(as.factor(cluster))
  
  has_fe <- length(obj$fixef_vars) > 0
  
  X <- model.matrix(obj, type = "rhs")
  y <- model.matrix(obj, type = "lhs")
  
  w <- weights(obj)
  if(!is.null(w)){
    X <- sqrt(w) * X
    y <- sqrt(w) * y 
    stop("Weighted least squares (WLS) is currently not supported for objects of type fixest.")
  }

  if(has_fe){
    fe <- model.matrix(obj, type = "fixef")
    X <- fixest::demean(X = X, f = fe)
    y <- fixest::demean(X = y, f = fe)
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
  
  y_list <- split(y, cluster)
  
  # X_g's in MNW notation
  X_list <- matrix_list(X, cluster, "row")

  # get small sample adjustments
  adjustments <- 1
  # multiply design matrices XW by sqrt(small_sample_adjustments)
  E_list <- lapply(X_list, function(e) e * adjustments)
  
  XWg_XWg <- Map(function(e) crossprod(e), e = E_list)
  yWg_XWg <-  Map(function(e,y) t(e) %*% as.matrix(y), e = E_list, y = y_list)
  tXX <- Reduce("+", XWg_XWg)
  tXy <- Reduce("+", yWg_XWg)
  # beta_g - beta
  coef_list <- Map(function(a = tXX, b, c = tXy, d, beta = coef(obj)) MASS::ginv(a - b) %*% (c - d) - beta, b = XWg_XWg, d = yWg_XWg)
  vcov <- Reduce("+",Map(function(x) tcrossprod(x), x = coef_list) )
  
  vcov <- vcov * (J-1)/J

  rownames(vcov) <- colnames(vcov) <- colnames(X)
  class(vcov) <- c("vcov_CV3J")
  return(vcov)
  
}


vcov_CR3J <- function(obj, ...) {
  
  #' Compute CR3 Jackknive variance covariance matrices of objects of type lm and fixest
  #' @param obj An object of class `lm` or `fixest``
  #' @param ... Other arguments
  #' @export
  
  UseMethod("vcov_CR3J")
  
}


