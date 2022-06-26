#' Objective function of the L2E regression
#' 
#' \code{objective} compute the objective of the L2E regression in terms of eta
#' 
#' @param eta The current estimate of eta
#' @param r Vector of residuals
#' @param method Mean or median 
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
#' \code{objective_tau} compute the objective of the L2E regression in terms of tau
#' 
#' @param tau The current estimate of tau
#' @param r Vector of residuals
#' @param method Mean or median 
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


