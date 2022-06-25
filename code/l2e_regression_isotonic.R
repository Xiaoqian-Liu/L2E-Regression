#' LTE Isotonic Regression
#' 
#' \code{l2e_regression_isotonic} Performs L2E isotonic regression via inexact coordinate descent.
#' 
#' @param y Response vector
#' @param b Initial vector of regression coefficients
#' @param tau Initial precision estimate
#' @param max_iter Maximum number of iterations
#' @param tol Relative tolerance
#' @param Show.Time Report the computing time
#' @export
#' @examples
#' set.seed(12345)
#' n <- 200
#' tau <- 1
#' x <- seq(-2.5, 2.5, length.out=n)
#' f <- x^3
#' y <- f + (1/tau)*rnorm(n)
#' 
#' # Clean Data
#' plot(x, y, pch=16, cex.lab=1.5, cex.axis=1.5, cex.sub=1.5, col='gray')
#' lines(x, f, lwd=3)
#' 
#' tau <- 1
#' b <- y
#' sol <- l2e_regression_isotonic(y, b, tau)
#' 
#' plot(x, y, pch=16, cex.lab=1.5, cex.axis=1.5, cex.sub=1.5, col='gray')
#' lines(x, f, lwd=3)
#' iso <- isotone::gpava(1:n, y)$x
#' lines(x, iso, col='blue', lwd=3)
#' lines(x, sol$beta, col='dark green', lwd=3)
#' 
#' # Contaminated Data
#' ix <- 0:9
#' y[45 + ix] <- 14 + rnorm(10)
#' 
#' plot(x, y, pch=16, cex.lab=1.5, cex.axis=1.5, cex.sub=1.5, col='gray')
#' lines(x, f, lwd=3)
#' 
#' tau <- 1
#' b <- y
#' sol <- l2e_regression_isotonic(y, b, tau)
#' 
#' plot(x, y, pch=16, cex.lab=1.5, cex.axis=1.5, cex.sub=1.5, col='gray')
#' lines(x, f, lwd=3)
#' iso <- isotone::gpava(1:n, y)$x
#' lines(x, iso, col='blue', lwd=3)
#' lines(x, sol$beta, col='dark green', lwd=3)
#' 
l2e_regression_isotonic <- function(y,b,tau,max_iter=1e2,tol=1e-4,Show.Time=TRUE) {
  
  time <- proc.time()
  sd_y <- sd(y)
  
  # save the inner loop iters
  iter_beta <- rep(0, max_iter)
  iter_tau <- rep(0, max_iter)
  
  for (i in 1:max_iter) {
    b_last <- b
    tau_last <- tau
    
    res_b <- update_beta_isotonic(y,b,tau,max_iter=10,tol=1e-4)
    b <- res_b$beta 
    iter_beta[i] <- res_b$iter
    
    r <- y - b
    res_tau <- update_tau_R(r, tau_last, sd_y)
    tau <- res_tau$tau
    iter_tau[i] <- res_tau$iter
    
    A <- norm(as.matrix(b_last-b),'f') < tol*(1 + norm(as.matrix(b_last),'f'))
    B <- abs(tau_last-tau) < tol*(1 + tau_last)
    if (A & B) break
  }
  
  if(Show.Time) print(proc.time() - time)
  
  return(list(beta=b, tau=tau, iter=i,
              iter_beta=iter_beta[1:i], iter_tau=iter_tau[1:i]))
}