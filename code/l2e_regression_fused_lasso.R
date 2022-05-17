#' L2E Fused Regression with Lasso
#' 
#' \code{l2e_regression_fused_lasso} Performs robust fused regression under the L2 criterion with Lasso penalty
#' 
#' @param y Response vector
#' @param X Design matrix
#' @param beta Initial vector of regression coefficients
#' @param tau Initial precision estimate
#' @param D The fusion matrix
#' @param lambda The tuning parameter 
#' @param max_iter Maximum number of iterations
#' @param tol Relative tolerance
#' @param Shou.Time Report the computing time
#' @export
#' 
l2e_regression_fused_lasso <- function(y, X, beta, tau, D, lambda=1, max_iter=1e2,
                                      tol=1e-4, Show.Time=TRUE) {
  
  time <- proc.time()
  Obj <- double(max_iter)
  for (i in 1:max_iter) {
    
    # update beta
    beta_last <- beta
    sol_beta <- update_beta_fused_lasso(y, X, beta_last, tau, D, lambda, max_iter=1e2,tol=1e-4)
    beta <- sol_beta$beta 
    # update tau
    r <- y - X%*%beta
    eta_last <- log(tau)  # get eta as lin line 9
    res_eta <- update_eta_bktk(r,eta_last, max_iter=10, tol=tol) # update eta as in line 10-12
    eta <- res_eta$eta
    tau <- exp(eta) # update tau as in line 13
    
    
    Obj[i] <- objective(eta, r) + lambda*sum(abs(D%*%beta))
    # Check for convergence
    A <- norm(as.matrix(beta_last-beta),'f') < tol*(1 + norm(as.matrix(beta_last),'f'))
    B <- abs(eta_last-eta) < tol*(1 + abs(eta_last))
    if (A & B) break
  }
  
  if(Show.Time) print(proc.time() - time)
  
  return(list(beta=beta,tau=tau, iter=i, Obj=Obj[1:i]))
  
}
