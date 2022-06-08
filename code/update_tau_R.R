#' Tau update function
#' 
#' \code{update_tau_R} Function for updating tau
#' 
#' @param r Residual vector
#' @param tau Current estimate for tau
#' @param sd_y Standard deviation of y
#' @param max_iter Maximum number of iterations
#' @param tol Relative tolerance
#' @return Returns a list object containing the new estimate for tau (scalar) and the number of iterations (scalar) the update step utilized
#' @export
update_tau_R <- function(r, tau, sd_y, max_iter=1e2, tol=1e-10) {
  tau_min <- as.numeric(1/sd_y)
  tau_max <- Inf
  tau <- as.numeric(tau)
  tau_sq <- as.numeric(tau^2)
  r_sq <- r^2
  r_2norm_sq <- sum(r_sq)
  if (r_2norm_sq == 0) { 
    tau_optimal = 0 
    } else {
    tau_optimal <- 1/min(abs(r[abs(r) > 0]))
  }
  
  n <- length(r)
  L <- (3/n) * sqrt(2/pi) * tau_optimal * r_2norm_sq * exp(-0.5) 
  
  for (i in 1:max_iter) {
    tau_last <- tau
    tau_sq <- tau_last^2
    
    e <- exp(-0.5 * tau_sq * r_sq)
    df <- 1/(2*sqrt(pi)) - (sqrt(2/pi)/n)*(sum(e) - tau_sq*sum(r_sq*e))
    tau <- tau - (1/L)*df
    
    if (tau < tau_min) {tau <- tau_min}
    if (tau > tau_max) {tau <- tau_max}
    if (abs(df) < tol) {break}
  }
  
  return(list(tau=as.numeric(tau),iter=i))
}

