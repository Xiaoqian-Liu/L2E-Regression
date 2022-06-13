update_beta_sparse_MM <- function(y,X,beta,tau,k,rho,stepsize=0.9,sigma=0.5,max_iter=1e2,tol=1e-4) {  
  
  n <- nrow(X)

  for (i in 1:max_iter) {
    
    beta_last <- beta
    Xbeta <- X %*% beta
    r <- y - Xbeta
    w <- as.vector(exp(-0.5* (tau*r)**2 ))
    
    
    beta <- ProxDist(beta0=beta, rho=rho, X=X, y=y, w=w, k=k, stepsize=stepsize, sigma=sigma,
                     max_iter=max_iter, tol=tol)$beta
    
    if (norm(as.matrix(beta_last-beta),'f') < tol*(1 + norm(as.matrix(beta_last),'f'))) break
  }

  return(list(beta=beta,iter=i))
}






## Proximal distance algorithm to update beta
ProxDist <- function(beta0, rho, X, y, w, k, stepsize=0.9, sigma=0.5, max_iter=1e2, tol=1e-4){
  n <- nrow(X)
  p <- ncol(X)
  
  # prepare the inverse matrix
  Ip <- Matrix::Diagonal(n=p, x=1)
  W <- Matrix::Diagonal(n=n, x = w)
  XtW <- t(X)%*%W
  XtWy <- XtW%*%y
  XtWX <- XtW%*%X
  
  # compute S=H_inv
  if(n>=p){
    S <- solve(XtWX + rho*Ip)
  }else{
    
    XXt <- X%*%t(X)
    H <- solve(Matrix::Diagonal(n=n, x=1/w) + XXt/rho)
    S <- Ip/rho - t(X)%*%H%*%X/(rho^2)
  }
  
  # now start the interation
  for(j in 1:max_iter){
    z0 <- Proj_sparse(beta0, k)
    gradf <- as.matrix(rho*(beta0 - z0) + XtWX%*%beta0 - XtWy)
    v <- as.vector(-S%*%gradf)
    eta <- 1
    
    f1 <- f(X, y, beta=beta0+eta*v, k, rho)
    f0 <- f(X, y, beta=beta0, k, rho)
    diff <- as.numeric(stepsize*eta*t(gradf)%*%beta0)
    
    while(f1> f0+diff){
      eta <- sigma*eta
      f1 <- f(X, y, beta=beta0+eta*v, k, rho)
      f0 <- f(X, y, beta=beta0, k, rho)
      diff <- as.numeric(stepsize*eta*t(gradf)%*%beta0)
      
    }
    
    beta1 <- beta0+eta*v
    
    if (norm(as.matrix(beta1-beta0),'f') < tol*(1 + norm(as.matrix(beta0),'f'))) {
      break
    }else{
      beta0 <- beta1
    }
  }
  
  return(list(beta = beta1, iter = j))
}






## Project x onto the set C={x| x has at most k zero entries}
Proj_sparse <- function(x, k){
  
  res <- sort(abs(x), method="quick", index.return=TRUE)
  ind <- tail(res$ix, k)
  xnew <- rep(0, length(x))
  xnew[ind] <- x[ind]
  return(xnew)
}

f <- function(X, y, beta, k, rho){
  pjbeta <- Proj_sparse(beta, k)
  s1 <- rho*norm(beta - pjbeta, "2")^2/2
  Xbeta <- X%*%beta
  s2 <- norm(y-Xbeta, "2")^2/2
  return(s1+s2)
}
