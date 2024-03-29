#' Solution path of the L2E Fused Regression with distance penalization
#' 
#' \code{L2E_fused_dist} computes the solution path of the robust sparse regression under the L2 criterion with distance penalty
#' 
#' @param y Response vector
#' @param X Design matrix
#' @param beta0 Initial vector of regression coefficients
#' @param tau0 Initial precision estimate
#' @param D The fusion matrix
#' @param kSeq A sequence of tuning parameter k, the number of nonzero entries in the estimated coefficients
#' @param rhoSeq A increasing sequence of tuning parameter rho
#' @param max_iter Maximum number of iterations
#' @param tol Relative tolerance
#' @param Shou.Time Report the computing time
#' @export
#' 
L2E_fused_dist <- function(y, X, beta0, tau0, D, kSeq, rhoSeq, max_iter=1e2,tol=1e-4, Show.Time=TRUE) {
  
  if(missing(beta0)){
    beta0 <- filter(MedianFilter(9), y)  # initial beta, random initial is very bad
  }
  
  if(missing(tau0)){
    tau0 <- 1/mad(y)  # initial tau
  }
  
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
      
      res <- l2e_regression_fused_dist(y, X, beta = beta0, D=D, tau = tau0, k=k, rho=rhoSeq[j], 
                                      max_iter=max_iter, tol=tol, Show.Time = FALSE)
      beta[, j] <- beta0 <- res$beta
      tau[j] <- tau0 <- res$tau # but warm start tau here helps reduce running time
      
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
