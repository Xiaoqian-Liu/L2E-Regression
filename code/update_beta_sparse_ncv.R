update_beta_sparse_ncv <- function(y,X,beta,tau,lambda, penalty, max_iter=1e2,tol=1e-4) {  
  
  n <- nrow(X)
 
  for (i in 1:max_iter) {
    
    beta_last <- beta
    Xbeta <- X %*% beta
    r <- y - Xbeta
    w <- as.vector(exp(-0.5* (tau*r)**2 ))
    
    
    W <- Matrix::Diagonal(n=n, x = sqrt(as.vector(w)))
    Xtilde <- as.matrix(W%*%X)
    ytilde <- as.vector(W%*%y)
    beta <- as.vector(ncvreg::ncvfit(Xtilde, ytilde, init = beta_last, penalty=penalty,lambda = lambda,
                             max.iter = 100, warn = FALSE)$beta)
    
    if (norm(as.matrix(beta_last-beta),'f') < tol*(1 + norm(as.matrix(beta_last),'f'))) break
  }

  return(list(beta=beta,iter=i))
}





