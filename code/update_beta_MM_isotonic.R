#' Update beta with isotonic regularization using MM
#' 
#' \code{update_beta_MM_isotonic} updates the regression coefficients for l2e isotonic regression using MM
#' 
#' @param y response
#' @param beta initial vector of regression coefficients
#' @param tau precision estimate
#' @param max_iter maximum number of iterations
#' @param tol relative tolerance
#' @export

update_beta_MM_isotonic <- function(y,beta,tau,max_iter=1e2,tol=1e-4) {  
  
  n <- length(y)
  for (i in 1:max_iter) {
    beta_last <- beta
    # Compute the weights
    r <- y - beta_last
    w <- exp(-0.5* (tau*r)**2 )
    w <- ifelse(w==0, 1e-20, w) # 0 weights will cause NA in beta
    
    # Now solve a weighted isotonic regression,
    # Iso package has a function pava, solving \sum_{i=1}^n w_ii (y_i - \beta_i)^2
    beta <- gpava(z=1:n, y =y, weights = sqrt(w))$x
    if (norm(as.matrix(beta_last-beta),'f') < tol*(1 + norm(as.matrix(beta_last),'f'))) break
  }
  return(list(beta=beta,iter=i))
}
