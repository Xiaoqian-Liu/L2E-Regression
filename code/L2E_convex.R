#' Wrapper function for L2E convex regression
#' 
#' \code{L2E_convex} Performs convex regression under the L2 criterion. Available methods include projected gradient descent (PG) and majorization-minimization (MM).
#' 
#' @param y Response vector
#' @param beta Initial vector of regression coefficients
#' @param tau Initial precision estimate
#' @param method Available methods include PG and MM. MM by default.
#' @param max_iter Maximum number of iterations
#' @param tol Relative tolerance
#' @param Show.Time Report the computing time
#' @export
#' @examples
#' set.seed(12345)
#' n <- 200
#' tau <- 1
#' x <- seq(-2, 2, length.out=n)
#' f <- x^4 + x
#' y <- f + (1/tau) * rnorm(n)
#' 
#' ## Clean Data
#' plot(x, y, pch=16, cex.lab=1.5, cex.axis=1.5, cex.sub=1.5, col='gray')
#' lines(x, f, lwd=3)
#' 
#' tau <- 1
#' b <- y
#' ## Least Squares method
#' cvx <- fitted(cobs::conreg(y, convex=TRUE)) 
#' ## MM method
#' sol_mm <- L2E_convex(y, b, tau)
#' ## PG method
#' sol_pg <- L2E_convex(y, b, tau, method='PG') 
#' 
#' plot(x, y, pch=16, cex.lab=1.5, cex.axis=1.5, cex.sub=1.5, col='gray')
#' lines(x, f, lwd=3)
#' lines(x, cvx, col='blue', lwd=3) ## LS
#' lines(x, sol_mm$beta, col='red', lwd=3) ## MM
#' lines(x, sol_pg$beta, col='green', lwd=3) ## PG
#' 
#' ## Contaminated Data
#' ix <- 0:9
#' y[45 + ix] <- 14 + rnorm(10)
#' 
#' plot(x, y, pch=16, cex.lab=1.5, cex.axis=1.5, cex.sub=1.5, col='gray')
#' lines(x, f, lwd=3)
#' 
#' tau <- 1
#' b <- y
#' cvx <- fitted(cobs::conreg(y, convex=TRUE))
#' sol_mm <- L2E_convex(y, b, tau)
#' sol_pg <- L2E_convex(y, b, tau, method='PG')
#' 
#' plot(x, y, pch=16, cex.lab=1.5, cex.axis=1.5, cex.sub=1.5, col='gray')
#' lines(x, f, lwd=3)
#' lines(x, cvx, col='blue', lwd=3) ## LS
#' lines(x, sol_mm$beta, col='red', lwd=3) ## MM
#' lines(x, sol_pg$beta, col='green', lwd=3) ## PG 
L2E_convex <- function(y, beta, tau, method = "MM", max_iter = 1e2, tol = 1e-4, Show.Time = TRUE) {
  
  if (tau <= 0) stop("Entered non-positive initial tau")
  
  if (method == "MM") {
    l2e_regression_convex_MM(y=y, beta=beta, tau=tau, max_iter=max_iter, tol=tol, Show.Time=Show.Time)
  }
  else if (method == "PG") {
    l2e_regression_convex(y=y, b=beta, tau=tau, max_iter=max_iter, tol=tol, Show.Time=Show.Time)
  }
  else {
    stop("Method unavailable. Available methods are PG and MM.")
  }
}