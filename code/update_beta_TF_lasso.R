update_beta_TF_lasso <- function(y,X,beta,tau,D,lambda,max_iter=1e2,tol=1e-4) {  
  
  n <- nrow(X)
  for (i in 1:max_iter) {
    
    beta_last <- beta
    Xbeta <- X %*% beta
    r <- y - Xbeta
    w <- as.vector(exp(-0.5* (tau*r)**2 ))
    
    beta <- sol_TF_lasso(y, X, w, D, lambda)$beta
                                                                                                                                                                                                                                             
    
    if (norm(as.matrix(beta_last-beta),'f') < tol*(1 + norm(as.matrix(beta_last),'f'))) break
  }
  
  return(list(beta=beta,iter=i))
}









## QP to solve trend filtering lasso problem
sol_TF_lasso <- function(y, X, w, D, lambda) {
  
  n <- length(y)
  p <- ncol(X)
  r <- nrow(D)
  
  W <- Diagonal(n=n, x = w)
  XtW <- t(X)%*%W
  XtWy <- XtW%*%y
  XtWX <- XtW%*%X
  
  P <- Matrix(0, nrow = p+2*r, ncol = p+2*r)
  P[1:p, 1:p] <- XtWX/2
  
  q <- c(as.vector(-XtWy), rep(lambda, 2*r))

  A <- cbind(D, Diagonal(n=r, x=-1), Diagonal(n=r, x=1))
  
  b <- rep(0,r)
  
  
  model <- list()
  
  model$A          <- A
  model$Q          <- P
  model$obj        <- as.vector(q)
  model$modelsense <- 'min'
  model$rhs        <- b
  model$sense      <- rep('=', r)
  model$lb        <- c(rep(-Inf, p), rep(0, 2*r))
  
  
  params <- list(OutputFlag=0)
  
  result <- gurobi(model, params)
  
  beta <- result$x[1:p]
  
  v <- result$x[(p+1):(p+r)]-result$x[(p+r+1):(p+2*r)]
  return(list(beta=beta, v=v))
  
}
