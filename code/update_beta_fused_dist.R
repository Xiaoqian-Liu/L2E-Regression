update_beta_fused_dist <- function(y,X,beta,tau,D,k,rho,max_iter=1e2,tol=1e-4) {  
  
  n <- nrow(X)
  
  for (i in 1:max_iter) {
    beta_last <- beta
    Xbeta <- X %*% beta
    r <- y - Xbeta
    w <- as.vector(exp(-0.5* (tau*r)**2 ))
    
    
    beta <- gradient_descent(y, X, w, beta_last, D, rho, k)$beta

    if (norm(as.matrix(beta_last-beta),'f') < tol*(1 + norm(as.matrix(beta_last),'f'))) break
  }
  
  return(list(beta=beta,iter=i))
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



#' Conjugate gradient algorithm for solving high-dimenisonal linear equation system Ax=b
#' 
#' @param A  the coefficient matrix
#' @param b  the vector on the right-side hand
#' @param x0  the initial vector for the solution
#' @export

cg <- function(A, b, x0, tol=1e-3, maxiter=1e3){
  A=Matrix(A, sparse = TRUE)
  #first check A
  if(isSymmetric(A)=="FALSE"){
    print("Error: The A matrix is not symmetric!")
  }
  
  x_k <- x0
  r_k <- A%*%x_k-b
  p_k <- -r_k
  for (k in 1:maxiter) {
    y <- A%*%p_k
    a_k <- as.numeric(t(r_k)%*%(r_k)/t(p_k)%*%y)
    x_1 <- x_k+a_k*p_k
    r_1 <- r_k+a_k*y
    
    dif <- norm(r_1,"2")
    if(dif>tol){
      beta_1 <- as.numeric(t(r_1)%*%r_1/t(r_k)%*%(r_k))
      p_1 <- -r_1+beta_1*p_k
      
      p_k <- p_1
      r_k <- r_1
      x_k <- x_1
      
    }else{
      break()
    }
  }
  
  return(list(x=x_1, iters = k, diffnorm = dif))
  
}



gradient_descent <- function(y, X, w, beta, D, rho, k, max_iter=1e2, tol=1e-5){
  n <- length(w)
  W <- Diagonal(n=n, x = w)
  XtW <- t(X)%*%W
  XtWy <- XtW%*%y
  XtWX <- XtW%*%X
  
  for (i in 1:max_iter) {
    #print(i)
    beta_last <- beta
    a <- XtWX%*%beta_last
    Dbeta <- as.vector(D%*%beta_last)
    Dbeta_proj <- Proj_sparse(Dbeta, k)
    
    gradient <- as.vector(a - XtWy+rho*t(D)%*%(Dbeta- Dbeta_proj))
    
    Av <- XtWX%*%gradient
    vAv <- t(gradient)%*%Av
    Dv <- as.vector(D%*%gradient)

    stepsize <- as.numeric(norm(gradient, "2")^2/(vAv + rho*norm(Dv, "2")^2))
    #print(stepsize)
    beta <- beta_last - stepsize*gradient
    
    if (norm(as.matrix(beta_last-beta),'f') < tol*(1 + norm(as.matrix(beta_last),'f'))) break
  }
  
  return(list(beta=beta))
}
