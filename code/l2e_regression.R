#' LTE multivariate regression
#' 
#' \code{l2e_regression} Performs L2E regression via inexact coordinate descent
#' 
#' @param y Response vector
#' @param X Design matrix
#' @param tau Initial precision estimate
#' @param b Initial vector of regression coefficients
#' @param max_iter Maximum number of iterations
#' @param tol Relative tolerance
#' @return Returns a list object containing the estimates for beta (vector) and tau (scalar)
#' @importFrom stats sd
#' @export
#' @examples
#' 
#' # Bank data example
#' y <- bank$y
#' X <- as.matrix(bank[,1:13])
#' X0 <- as.matrix(cbind(rep(1,length(y)), X))
#' tauinit <- 1/mad(y)
#' binit <- matrix(0, 14, 1)#' 
#' sol <- l2e_regression(y, X0, tauinit, binit)
#' r <- y - X0 %*% sol$beta
#' ix <- which(abs(r) > 3/sol$tau)
#' l2e_fit <- X0 %*% sol$beta
#' plot(y, l2e_fit, ylab='Predicted values', pch=16, cex=0.8)
#' points(y[ix], l2e_fit[ix], pch=16, col='blue', cex=0.8)
#' 
l2e_regression <- function(y,X,tau,b,max_iter=1e2,tol=1e-4) {
  n <- nrow(X); p <- ncol(X)
  if (p >= n) stop('Current implementation can only handle n > p')
  QRF <- qr(X)
  sd_y <- as.numeric(sd(y))
  for (i in 1:max_iter) {
    b_last <- b
    tau_last <- tau
    b <- update_beta_qr(y,X,QRF,tau_last,b_last)$beta
    r <- y - X%*%b
    tau <- update_tau_R(r,tau_last, sd_y)$tau
    A <- norm(as.matrix(b_last-b),'f') < tol*(1 + norm(as.matrix(b_last),'f'))
    B <- abs(tau_last-tau) < tol*(1 + tau_last)
    if (A & B) break
  }
  return(list(beta=b,tau=tau))
}