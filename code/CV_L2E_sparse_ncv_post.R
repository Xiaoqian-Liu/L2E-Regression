CV_L2E_sparse_ncv_post <- function(y, X, beta0, tau0, lambda, penalty, beta_star){
  
  sol_final <- l2e_regression_sparse_ncv(y, X, beta0, tau0, lambda = lambda, 
                                         penalty=penalty, Show.Time = FALSE)
  
  ## Get active set to refine the raw L2E betahat
  true_active_set <- which(beta_star != 0)
  activesetL2E <- which(sol_final$beta != 0)
  activeset_truepos_L2E <- length(which(activesetL2E %in% true_active_set))
  activeset_falsepos_L2E <- length(which(!(activesetL2E %in% true_active_set)))
  X1 <- X[,activesetL2E,drop=FALSE] # make new covariates with just active set
  
  ## Use optim::bfgs to get unbiased estimate
  ## Objective function for optim
  fxx <- function(xx) {
    n <- length(y)
    beta <- matrix(xx[-1], ncol=1)
    tau <- xx[1]
    r <- y - X1%*%beta
    h <- 0.5*tau/sqrt(pi) - (tau/n)*sqrt(2/pi)*sum(exp(-0.5*(tau*r)^2))
    return(h)
  }
  
  ## Gradient function for optim
  gxx <- function(xx) {
    n <- length(y)
    tau <- xx[1]
    beta <- matrix(xx[-1], ncol=1)
    r <- y - X1%*%beta
    w <- exp(-0.5*(tau*r)^2)
    dh_beta <- -((tau^3)/n)*sqrt(2/pi)*t(X1)%*%matrix(w*r, ncol=1)
    dh_tau <- 0.5/sqrt(pi) - (1/n)*sqrt(2/pi)*sum(w*(1 - (tau*r)^2))
    return(c(dh_tau, dh_beta))
  }
  
  tau1 <- 1/mad(y)
  ## multi-try
  best <- Inf
  for (i in 1:1e3) {
    sol1 <- optim(c(rexp(1, rate=1/tau1), mad(y)*rnorm(length(activesetL2E))),
                  fn=fxx, gr=gxx, method="L-BFGS-B",lower = c(0, rep(Inf, length(activesetL2E))))
    if (sol1$value < best) {
      best <- sol1$value
      sol_best <- sol1
    }
  }
  betaL2E_v2 <- numeric(length=length(beta0))
  betaL2E_v2[activesetL2E] <- sol_best$par[-1]
  #relmse_L2E <- sqrt( sum( (betaL2E_v2 - beta0) ** 2)/sum(beta0**2))
  #residuals_L2E <- y - X1 %*% sol_best$par[-1] # residuals
  return(list(beta=betaL2E_v2, tau=sol_best$par[1]))
}
