#' Update beta with convex regularization using MM
#' 
#' \code{update_beta_MM_convex} updates the regression coefficients for l2e convex regression using MM
#' 
#' @param y response
#' @param beta initial vector of regression coefficients
#' @param tau precision estimate
#' @param max_iter maximum number of iterations
#' @param tol relative tolerance
#' @export

update_beta_MM_convex <- function(y,beta,tau,max_iter=1e2,tol=1e-4) {  
  
  n <- length(y)
  for (i in 1:max_iter) {
    beta_last <- beta
    ## Compuete the weightes
    r <- y - beta_last
    w <- exp(-0.5* (tau*r)**2 )
    
    ## Now solve a weighted convex regression,
    ## cobs pacakge has a function conreg, solving \sum_{i=1}^n w_ii (y_i - \beta_i)^2
    beta <- conreg(x=1:n, y, sqrt(w), convex = TRUE, method = "Duembgen06_R")$yf
    # beta <- fitted(shapereg(ytilde ~ conv(1:n), weights = sqrt(w)))
    # beta <- beta/sqrt(w)
    #print(i)
    if (norm(as.matrix(beta_last-beta),'f') < tol*(1 + norm(as.matrix(beta_last),'f'))) break
  }
  return(list(beta=beta,iter=i))
}
