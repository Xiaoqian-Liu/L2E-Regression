#' L2E Sparse Regression with distance penalization
#' 
#' \code{L2E_sparse_dist} Computes the solution path of the robust sparse regression under the L2 criterion with distance penalty
#' 
#' @param y Response vector
#' @param X Design matrix
#' @param beta0 Initial vector of regression coefficients, can be omitted
#' @param tau0 Initial precision estimate, can be omitted
#' @param kSeq A sequence of tuning parameter k, the number of nonzero entries in the estimated coefficients
#' @param rhoSeq An increasing sequence of tuning parameter rho, can be omitted
#' @param stepsize The stepsize parameter for the MM algorithm (0, 1)
#' @param sigma The halving parameter sigma (0, 1)
#' @param max_iter Maximum number of iterations
#' @param tol Relative tolerance
#' @param Show.Time Report the computing time
#' @export
#' @examples
#' set.seed(12345)
#' n <- 200
#' tau <- 1
#' f <- matrix(c(rep(2,5), rep(0,45)), ncol = 1)
#' X <- X0 <- matrix(rnorm(n*50), nrow = n)
#' y <- y0 <- X0 %*% f + (1/tau)*rnorm(n)
#' x <- 1:length(f)
#' 
#' ## Clean Data 
#' k <- c(4,5,6)
#' sol <- L2E_sparse_dist(y=y, X=X, kSeq=k)
#' 
#' plot(x, f, type='b', pch=1, ylim=c(-2.5,3))
#' points(x, sol$Beta[,1], col='blue', type='b', pch=0) ## k=4
#' points(x, sol$Beta[,2], col='red', type='b', pch=2) ## k=5
#' points(x, sol$Beta[,3], col='dark green', type='b', pch=8) ## k=6
#' 
#' ## Contaminated Data
#' ix <- 1:20
#' y[ix] <- 2 + y0[ix] 
#' X[ix,] <- 2 + X0[ix,]
#' 
#' sol <- L2E_sparse_dist(y=y, X=X, kSeq=k)
#' 
#' plot(x, f, type='b', pch=1, ylim=c(-2.5,3))
#' points(x, sol$Beta[,1], col='blue', type='b', pch=0) ## k=4
#' points(x, sol$Beta[,2], col='red', type='b', pch=2) ## k=5
#' points(x, sol$Beta[,3], col='dark green', type='b', pch=8) ## k=6
#' 
L2E_sparse_dist <- function(y, X, beta0, tau0, kSeq, rhoSeq, stepsize = 0.9, sigma=0.5, max_iter=1e2,
                                     tol=1e-4, Show.Time=TRUE) {
  
  if(missing(beta0)){
    beta0 <- rnorm(ncol(X))  # initial beta
  }
  
  if(missing(tau0)){
    tau0 <- 1/mad(y)  # initial tau
  }
  
  if (tau0 <= 0) stop("Entered non-positive tau0")
  
  if(missing(rhoSeq)){
    rhoSeq <- 10^seq(0, 4, length.out = 20)  # set a sequence of rho
  }
  

  Nk <- length(kSeq)
  Nrho <- length(rhoSeq)
  
  # save outputs
  Beta <- matrix(0, nrow = ncol(X), ncol = Nk)
  Beta_path <- list()
  Tau <- double(Nk)
  Tau_path <- matrix(0, Nk, Nrho)
  time <- double(Nk)
  
  for (i in 1:Nk) {
    
    k <- kSeq[i]
    beta <- matrix(0, nrow = ncol(X), ncol = Nrho)
    tau <- double(Nrho)
    
    time <- proc.time()
    for (j in 1:Nrho) {
      
      res <- l2e_regression_sparse_MM(y, X, beta = beta0, tau = tau0, k=k, rho=rhoSeq[j], stepsize = stepsize, sigma=sigma, 
                                      max_iter=max_iter, tol=tol, Show.Time = FALSE)
      beta[, j] <- beta0 <- res$beta
      tau[j] <- tau0 <- res$tau
      
    }
    runtime <-  proc.time() - time
    if(Show.Time) print(runtime)
    
    Beta_path[[i]] <- beta
    Beta[, i] <- beta0
    Tau_path[i, ] <- tau
    Tau[i] <- tau0
    time[i] <- runtime[[3]]
  }
  
 
  return(list(Beta=Beta, Beta_path=Beta_path, Tau=Tau, Tau_path=Tau_path, time = time,
              rhoSeq = rhoSeq, kSeq = kSeq))
  
}
