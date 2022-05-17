#' Objective function of the L2E regression
#' 
#' \code{objective} compute the obejctive of the L2E regression 
#' 
#' @param eta the current estimate of tau
#' @param r vector of residual
#' @param method mean or median 
#' @export
objective <- function(eta, r, method="mean"){
  
  v1 <- exp(-0.5*exp(2*eta)*r^2)
  
  s1 <- exp(eta)/(2*sqrt(pi))
  
  if(method=="mean"){
    s2 <- exp(eta)* sqrt(2/pi)*mean(v1)
  }else{
    s2 <- exp(eta)* sqrt(2/pi)*median(v1)
  }
  
  return(s1-s2)
}


#' Objective function of the L2E regression
#' 
#' \code{objective_tau} compute the obejctive of the L2E regression 
#' 
#' @param tau the current estimate of tau
#' @param r vector of residual
#' @param method mean or median 
#' @export
objective_tau <- function(tau, r, method="mean"){
  
  v1 <- exp(-0.5*tau^2*r^2)
  
  s1 <- tau/(2*sqrt(pi))
  
  if(method=="mean"){
    s2 <- tau* sqrt(2/pi)*mean(v1)
  }else{
    s2 <- tau* sqrt(2/pi)*median(v1)
  }
  
  return(s1-s2)
}


