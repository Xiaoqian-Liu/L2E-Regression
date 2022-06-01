#' Solution path of the L2E Trend Filtering Regression with Lasso
#' 
#' \code{L2E_TF_lasso} Computes the solution path of the robust trend filtering regression under the L2 criterion with Lasso penalty
#' 
#' @param y Response vector
#' @param X Design matrix
#' @param beta0 Initial vector of regression coefficients, can be omitted
#' @param tau0 Initial precision estimate, can be omitted
#' @param D The fusion matrix
#' @param lambdaSeq A decreasing sequence of values for the tuning parameter lambda, can be omitted
#' @param max_iter Maximum number of iterations
#' @param tol Relative tolerance
#' @param Show.Time Report the computing time
#' @export
#' 
L2E_TF_lasso <- function(y, X, beta0, tau0, D, lambdaSeq, max_iter=1e2, tol=1e-4, Show.Time=TRUE){
  
  if(missing(X)){
    X <- diag(nrow = length(y))  # initial X is identity matrix by default
  }
  
  if(missing(beta0)){
    beta0 <- rep(mean(y), ncol(X))  # initial beta
  }
  
  if(missing(tau0)){
    tau0 <- 1/mad(y)  # initial tau
  }
  
  if(missing(lambdaSeq)){
    lambdaSeq <- 10^seq(1, -4, length.out = 20)  # set a sequence of lambda
  }
  
  
  Nlambda <- length(lambdaSeq)
  
  # save outputs
  Beta <- matrix(0, nrow = ncol(X), ncol = Nlambda)
  Tau <- double(Nlambda)
  
  time <- proc.time()
  for (j in 1:Nlambda) {
    
    res <- l2e_regression_TF_lasso(y, X, beta=beta0, tau=tau0, D=D, lambda=lambdaSeq[j], 
                                 max_iter=max_iter, tol=tol, Show.Time = FALSE)
    Beta[, j] <- beta0 <- res$beta
    Tau[j] <- res$tau # warm start tau is not good
    
  }
  runtime <-  proc.time() - time
  if(Show.Time) print(runtime)
  
  return(list(Beta=Beta, Tau=Tau, runtime=runtime, lambdaSeq=lambdaSeq))
  
}
