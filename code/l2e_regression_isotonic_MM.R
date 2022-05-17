#' L2E Isotonic Regression - MM
#' 
#' \code{l2e_regression_isotonic} performs L2E isotonic regression via block coordinate descent 
#' with MM for updating beta andtau
#' 
#' @param y response
#' @param beta initial vector of regression coefficients
#' @param tau initial precision estimate
#' @param max_iter maximum number of iterations
#' @param tol relative tolerance
#' @param Show.Time Report the computing time
#' @export
#' @examples
#' set.seed(12345)
#' n <- 200
#' tau <- 1
#' x <- seq(-2.5,2.5,length.out=n)
#' f <- x^3
#' y <- f + (1/tau)*rnorm(n)
#' 
#' ## Clean
#' pdf(file='isotonic_clean0.pdf')
#' plot(x,y,pch=16,cex.lab=1.5, cex.axis=1.5,cex.sub=1.5)
#' lines(x,f,col='blue',lwd=3)
#' dev.off()
#' 
#' tau <- 1
#' b <- y
#' sol <- l2e_regression_isotonic_MM(y,b,tau)
#' 
#' pdf(file='isotonic_clean1.pdf')
#' plot(x,y,pch=16,cex.lab=1.5, cex.axis=1.5,cex.sub=1.5)
#' lines(x,f,col='blue',lwd=3)
#' iso <- gpava(1:n,y)$x
#' lines(x,iso,col='red',lwd=3)
#' lines(x,sol$beta,col='green',lwd=3)
#' dev.off()
#' 
#' ## Contaminated
#' ix <- 0:9
#' y[45 + ix] <- 14 + rnorm(10)
#' 
#' pdf(file='isotonic_contaminated0.pdf')
#' plot(x,y,pch=16,cex.lab=1.5, cex.axis=1.5,cex.sub=1.5)
#' lines(x,f,col='blue',lwd=3)
#' dev.off()
#' 
#' tau <- 1
#' b <- y
#' sol <- l2e_regression_isotonic_MM(y,b,tau)
#' 
#' pdf(file='isotonic_contaminated1.pdf')
#' plot(x,y,pch=16,cex.lab=1.5, cex.axis=1.5,cex.sub=1.5)
#' lines(x,f,col='blue',lwd=3)
#' iso <- gpava(1:n,y)$x
#' lines(x,iso,col='red',lwd=3)
#' lines(x,sol$beta,col='green',lwd=3)
#' dev.off()
l2e_regression_isotonic_MM <- function(y,beta,tau,max_iter=1e2,tol=1e-4, Show.Time=TRUE) {
  
  time <- proc.time()
  # save the inner loop iters
  iter_beta <- rep(0, max_iter)
  iter_eta <- rep(0, max_iter)
  
  for (i in 1:max_iter) {
    # update beta
    beta_last <- beta
    res_beta <- update_beta_MM_isotonic(y,beta_last,tau,max_iter=1e2,tol=tol)
    beta <- res_beta$beta
    iter_beta[i] <- res_beta$iter
    
    # update tau
    r <- y - beta
    eta_last <- log(tau)  # get eta as lin line 9
    res_eta <- update_eta_bktk(r,eta_last, tol=tol) # update eta as in line 10-12
    eta <- res_eta$eta
    tau <- exp(eta) # update tau as in line 13
    iter_eta[i] <- res_eta$iter
    
    A <- norm(as.matrix(beta_last-beta),'f') < tol*(1 + norm(as.matrix(beta_last),'f'))
    B <- abs(eta_last-eta) < tol*(1 + abs(eta_last))
    if (A & B) break
  }
  
  if(Show.Time) print(proc.time() - time)
  return(list(beta=beta,tau=tau, iter=i,
              iter_beta = iter_beta[1:i], iter_eta = iter_eta[1:i]))
}
