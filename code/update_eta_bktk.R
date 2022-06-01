#' Update eta using Newton's method with backtracking
#' 
#' \code{update_eta_bktk} Updates the precision parameter tau = e^eta for L2E regression using Newton's method
#' 
#' @param y Response vector
#' @param r Vector of residual
#' @param eta Initial estimate of eta
#' @param max_iter Maximum number of iterations
#' @param tol Relative tolerance
#' @export
#' 
update_eta_bktk <- function(r, eta, max_iter=1e2, tol=1e-10) {
  
  
  n <- length(r)
  r_sq <- r^2
  stepsize <- 1
  stepsize_shrinkage <- 0.9
  
  
  first_derivative_seq <- double(max_iter)
  second_derivative_seq <- double(max_iter)
  Eta <- double(max_iter)
  Obj <- double(max_iter)
  
  if ( exp(eta)<=0 ) print("non-positive initial tau in update_eta_bktk")
  
  for (i in 1:max_iter) {
    
    eta_last <- eta
    tau_last <- exp(eta_last)  # avoid computing tau_last in the following
    
    # some elements for computing the derivatives
    v1 <- exp(-0.5*tau_last^2*r_sq) # the w_i's
    v2 <- r_sq*v1
    
    if(all(round(v1,2)==1)) break # Check if v1 is close to 1
    
    first_derivative <- tau_last/(2*sqrt(pi)) - tau_last*sqrt(2/pi)*mean(v1)+
      tau_last^3*sqrt(2/pi)*mean(v2)
    first_derivative_seq[i] <- first_derivative
    
    second_derivative <- tau_last/(2*sqrt(pi))+ 4*tau_last^3*sqrt(2/pi)*mean(v2)
    second_derivative_seq[i] <- second_derivative
    
    
    lam <-  first_derivative^2/second_derivative
    if(lam < tol) break
    
    ### backtracking
    dd <- -first_derivative/second_derivative
    f1 <- objective(eta_last + stepsize*dd, r)
    f0 <- objective(eta_last, r)
    
    
    while (f1>f0-0.5*stepsize*lam) {
      stepsize <- stepsize_shrinkage*stepsize
      f1 <- objective(eta_last + stepsize*dd, r)
    }
    
    eta <- eta_last - stepsize*first_derivative
    stepsize <- 1
    
    Eta[i] <- eta
    Obj[i] <- objective(eta, r)
  }
  
  if(i>max_iter) i=i-1  
  
  return(list(eta=as.numeric(eta),iter=i, Eta=Eta[1:i], Obj=Obj[1:i],
              first_derivative_seq = first_derivative_seq[1:i],
              second_derivative_seq = second_derivative_seq[1:i]))
  
}



