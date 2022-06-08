#' L2E Trend Filtering Regression with distance penalty
#' 
#' \code{l2e_regression_TF_dist} Performs robust trend filtering regression under the L2 criterion with distance penalty
#' 
#' @param y Response vector
#' @param X Design matrix
#' @param beta Initial vector of regression coefficients
#' @param tau Initial precision estimate
#' @param D The fusion matrix
#' @param k The number of nonzero entries in `D %*% beta`
#' @param rho The parameter in the proximal distance algorithm
#' @param max_iter Maximum number of iterations
#' @param tol Relative tolerance
#' @param Show.Time Report the computing time
#' 
l2e_regression_TF_dist <- function(y, X, beta, tau, D, k, rho=1,  max_iter=1e2,
                                     tol=1e-4,  Show.Time=TRUE) {
  
  if (tau <= 0) stop("Entered non-positive initial tau")
  
  time <- proc.time()
  Obj <- double(max_iter)
  for (i in 1:max_iter) {
    
    # update beta
    beta_last <- beta
    sol_beta <- update_beta_TF_dist(y, X, beta_last, tau, D, k, rho, max_iter=1e2,tol=1e-4)
    beta <- sol_beta$beta 
    # update tau
    r <- y - X%*%beta
    eta_last <- log(tau)  # get eta as in line 9
    res_eta <- update_eta_bktk(r,eta_last, tol=tol) # update eta as in line 10-12
    eta <- res_eta$eta
    tau <- exp(eta) # update tau as in line 13
 

    Dbeta <- as.vector(D%*%beta)
    Dbeta_proj <- Proj_sparse(Dbeta, k)
    Obj[i] <- objective(eta, r) + rho*norm(Dbeta - Dbeta_proj, "2")^2/2
    # Check for convergence
    A <- norm(as.matrix(beta_last-beta),'f') < tol*(1 + norm(as.matrix(beta_last),'f'))
    B <- abs(eta_last-eta) < tol*(1 + abs(eta_last))
    if (A & B) break
  }
  
  if(Show.Time) print(proc.time() - time)
  
  return(list(beta=beta,tau=tau, iter=i, Obj=Obj[1:i]))
  
}
