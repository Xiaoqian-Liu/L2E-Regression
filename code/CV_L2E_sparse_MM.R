#'  Cross Validation for L2E Sparse Regression with distance penalization
#' 
#' \code{CV_L2E_sparse_MM} performs k-fold cross-validation for robust sparse regression undere the L2 criterion with
#' distance penalty
#' 
#' @param y Response vector
#' @param X Design matrix
#' @param beta0 Initial vector of regression coefficients, can be missed
#' @param tau0 Initial precision estimate
#' @param kSeq  A sequence of tuning parameter k,the number of nonzero entries in the estimated coefficients
#' @param rhoSeq A sequence of tuning parameter rho, can be missed
#' @param nfolds The number of cross-validation folds. Default is 10.
#' @param seed Users can set the seed of the random number generator to obtain reproducible results.
#' @param method Median or mean to compute the objective
#' @param max_iter Maximum number of iterations
#' @param tol Relative tolerance
#' @param trace Whether to trace the progress of the cross-validation
#' @export
#' 
CV_L2E_sparse_MM <- function(y, X, beta0, tau0, kSeq, rhoSeq,  nfolds=5, seed=1234, method="median",
                             max_iter=1e2, tol=1e-4, trace=TRUE) {
  
  
  if(missing(rhoSeq)){
    rhoSeq <- 10^seq(0, 4, length.out = 20)  # set a sequence of rho
  }
  
  if(missing(beta0)){
    beta0 <- rnorm(ncol(X))  # initial beta
  }
  
  if(missing(tau0)){
    tau0 <- 1/mad(y)   # initial tau
  }
  
  # Set up folds
  if (!missing(seed)) set.seed(seed)
  n <- length(y)
  fold <- sample(1:n %% nfolds)
  fold[fold==0] <- nfolds
  
  
  # Do cross-validation
  
  cv.args <- list()
  cv.args$beta0 <- beta0
  cv.args$tau0 <- tau0
  cv.args$kSeq <- kSeq
  cv.args$rhoSeq <- rhoSeq
  cv.args$max_iter <- max_iter
  cv.args$tol <- tol
  cv.args$Show.Time <- FALSE
  
  
  Loss <- matrix(0, nrow = nfolds, ncol = length(kSeq))
  for (i in 1:nfolds) {
    if (trace) cat("Starting CV fold #", i, sep="","\n")
    res <- cv_fold_l2e_MM(i, y, X, fold, cv.args, method=method)
    Loss[i, ] <- res
  }
  
  # Return
  cve <- apply(Loss, 2, mean)  
  cvse <- apply(Loss, 2, sd)/sqrt(nfolds)
  min <- which.min(round(cve, 8))
  
  
  #find the lambda.1se
  for (i in min:1) {
    if(cve[i]>cve[min]+cvse[min])
      break
  }
  if(min==1){
    k.1se <- kSeq[1]
    min_1se <- 1
  }else{
    k.1se <- kSeq[i+1]
    min_1se <- i+1
  }

  
  return(list(cve=cve, cvse=cvse, min=min, k.min=kSeq[min], min_1se=min_1se, k.1se=k.1se,
              kSeq=kSeq, rhoSeq = rhoSeq, fold=fold))
  
}





cv_fold_l2e_MM <- function(i, y, X, fold, cv.args, method="median") {
  cv.args$y <- y[fold!=i]
  cv.args$X <- X[fold!=i, , drop=FALSE]
  fit.i <- do.call("L2E_sparse_MM", cv.args)
  
  # data in hold-out
  y_out <- y[fold==i]
  X_out <- X[fold==i, , drop=FALSE]
 
  
  L <- length(fit.i$kSeq)
  loss <- double(L)
  
  for (l in 1:L) {
    bhat <- fit.i$Beta[, l]  # get the estimated beta with the l-th k
    Xbeta <- X_out %*% bhat
    r <- y_out - Xbeta
    tauhat <- fit.i$Tau[l]
    
    loss[l] <- objective_tau(tau = tauhat, r = r, method=method) ### use median istead of mean to account for outliers
  }
  
  return(loss)
}

