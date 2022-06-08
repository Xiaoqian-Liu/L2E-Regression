#' Cross Validation for L2E Sparse Regression with existing penalization methods
#' 
#' \code{CV_L2E_sparse_ncv} performs k-fold cross-validation for robust sparse regression under the L2 criterion.
#'  Available penalties include lasso, MCP and SCAD.
#' 
#' @param y Response vector
#' @param X Design matrix
#' @param beta0 Initial vector of regression coefficients, can be omitted
#' @param tau0 Initial precision estimate, can be omitted
#' @param lambdaSeq A decreasing sequence of tuning parameter lambda, can be omitted
#' @param penalty Available penalties include lasso, MCP and SCAD.
#' @param nfolds The number of cross-validation folds. Default is 5.
#' @param seed Users can set the seed of the random number generator to obtain reproducible results.
#' @param method Median or mean to compute the objective
#' @param max_iter Maximum number of iterations
#' @param tol Relative tolerance
#' @param trace Whether to trace the progress of the cross-validation
#' @export
#' 
CV_L2E_sparse_ncv <- function(y, X, beta0, tau0, lambdaSeq,  penalty="MCP", nfolds=5, seed=1234, method="median", 
                          max_iter=1e2, tol=1e-4, trace=TRUE) {
  
  
  if(missing(lambdaSeq)){
    lambdaSeq <- 10^seq(1, -4, length.out = 20)  # set a sequence of lambda
  }
  
  if(missing(beta0)){
    beta0 <- double(ncol(X))  # initial beta
  }
  
  if(missing(tau0)){
    tau0 <- 1/mad(y)  # initial tau
  }
  
  if (tau0 <= 0) stop("Entered non-positive tau0")
  
  
  # Set up folds
  if (!missing(seed)) set.seed(seed)
  n <- length(y)
  fold <- sample(1:n %% nfolds)
  fold[fold==0] <- nfolds
  
  
  # Do cross-validation
  
  cv.args <- list()
  cv.args$b <- beta0
  cv.args$tau <- tau0
  cv.args$lambdaSeq <- lambdaSeq
  cv.args$penalty <- penalty
  cv.args$max_iter <- max_iter
  cv.args$tol <- tol
  cv.args$Show.Time <- FALSE
  
  
  Loss <- matrix(0, nrow = nfolds, ncol = length(lambdaSeq))
  for (i in 1:nfolds) {
    if (trace) cat("Starting CV fold #", i, sep="","\n")
    res <- cv_fold_l2e_ncv(i, y, X, fold, cv.args, method=method)
    Loss[i, ] <- res
  }
  
  # Return
  cve <- apply(Loss, 2, mean) 
  cvse <- apply(Loss, 2, sd)/sqrt(nfolds)
  min <- which.min(round(cve,8))
  
  
  # find the lambda.1se
  for (i in min:1) {
    if(cve[i]>cve[min]+cvse[min])
      break
  }
  
  if(min==1){
    lambda.1se <- lambdaSeq[1]
    min_1se <- 1
  }else{
    lambda.1se <- lambdaSeq[i+1]
    min_1se <- i+1
  }
  

  
  return(list(cve=cve, cvse=cvse, min=min, lambda.min=lambdaSeq[min], min_1se=min_1se, lambda.1se=lambda.1se,
              lambdaSeq=lambdaSeq,  fold=fold))
  
}





cv_fold_l2e_ncv <- function(i, y, X, fold, cv.args, method="median") {
  cv.args$y <- y[fold!=i]
  cv.args$X <- X[fold!=i, , drop=FALSE]
  fit.i <- do.call("L2E_sparse_ncv", cv.args)
  
  # data in hold-out
  y_out <- y[fold==i]
  X_out <- X[fold==i, , drop=FALSE]
  
  
  L <- length(fit.i$lambdaSeq)
  loss <- double(L)
  
  for (l in 1:L) {
    bhat <- fit.i$Beta[, l]  # get the estimated beta with the l-th k
    Xbeta <- X_out %*% bhat
    r <- y_out - Xbeta
    tauhat <- fit.i$Tau[l]
    
    loss[l] <- objective_tau(tau = tauhat, r = r, method=method) ### use median instead of mean to account for outliers
  }

  return(loss)
}

