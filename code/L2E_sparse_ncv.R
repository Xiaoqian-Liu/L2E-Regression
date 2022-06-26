#' Solution path of L2E Sparse regression with existing penalization methods
#' 
#' \code{L2E_sparse_ncv} computes the solution path of robust sparse regression under the L2 criterion. Available penalties include lasso, MCP and SCAD.
#' 
#' @param y Response vector
#' @param X Design matrix
#' @param b Initial vector of regression coefficients, can be omitted
#' @param tau Initial precision estimate, can be omitted
#' @param lambdaSeq A decreasing sequence of values for the tuning parameter lambda, can be omitted
#' @param penalty Available penalties include lasso, MCP and SCAD.
#' @param max_iter Maximum number of iterations
#' @param tol Relative tolerance
#' @param Show.Time Report the computing time
#' @export
#' @examples
#' set.seed(12345)
#' n <- 200
#' tau <- 1
#' f <- matrix(c(rep(2,5), rep(0,45)), ncol = 1)
#' X <- X0 <- matrix(rnorm(n*50), nrow = n)
#' y <- y0 <- X0 %*% f + (1/tau)*rnorm(n)
#' x <- 1:length(f)
#' 
#' ## Clean Data 
#' lambda <- 10^-1
#' sol_mcp <- L2E_sparse_ncv(y=y, X=X, lambdaSeq=lambda)
#' sol_lasso <- L2E_sparse_ncv(y=y, X=X, lambdaSeq=lambda, penalty="lasso")
#' sol_scad <- L2E_sparse_ncv(y=y, X=X, lambdaSeq=lambda, penalty="SCAD")
#' 
#' plot(x, f, type='b', pch=1, ylim=c(-1,3))
#' points(x, sol_mcp$Beta, col='blue', type='b',pch=0) ## MCP
#' points(x, sol_lasso$Beta, col='red', type='b',pch=2) ## LASSO
#' points(x, sol_scad$Beta, col='dark green', type='b',pch=8) ## SCAD
#' 
#' ## Contaminated Data
#' ix <- 1:20
#' y[ix] <- 2 + y0[ix] 
#' X[ix,] <- 2 + X0[ix,]
#' 
#' sol_mcp <- L2E_sparse_ncv(y=y, X=X, lambdaSeq=lambda)
#' sol_lasso <- L2E_sparse_ncv(y=y, X=X, lambdaSeq=lambda, penalty="lasso")
#' sol_scad <- L2E_sparse_ncv(y=y, X=X, lambdaSeq=lambda, penalty="SCAD")
#' 
#' plot(x, f, type='b', pch=1, ylim=c(-1,3))
#' points(x, sol_mcp$Beta, col='blue', type='b', pch=0) ## MCP
#' points(x, sol_lasso$Beta, col='red', type='b', pch=2) ## LASSO
#' points(x, sol_scad$Beta, col='dark green', type='b', pch=8) ## SCAD
#' 
L2E_sparse_ncv <- function(y,X,b,tau,lambdaSeq,penalty="MCP",max_iter=1e2,tol=1e-4,Show.Time=TRUE) {
  
  if(missing(b)){
    b <- double(ncol(X))  # initial beta
  }
  
  if(missing(tau)){
    tau <- 1/mad(y)   # initial tau
  }
  
  if (tau <= 0) stop("Entered non-positive initial tau")
  
  if(missing(lambdaSeq)){
    lambdaSeq <- 10^seq(1, -4, length.out = 20)  # set a sequence of rho
  }
  
  
  Nlambda <- length(lambdaSeq)
  
  # save outputs
  Beta <- matrix(0, nrow = ncol(X), ncol = Nlambda)
  Tau <- double(Nlambda)
  
  time <- proc.time()
  for (j in 1:Nlambda) {
    
    res <- l2e_regression_sparse_ncv(y, X, b, tau, lambda=lambdaSeq[j], penalty=penalty,
                                 max_iter=max_iter, tol=tol, Show.Time = FALSE)
    Beta[, j] <- b <- res$beta
    Tau[j] <- tau <- res$tau
    
  }
  runtime <-  proc.time() - time
  if(Show.Time) print(runtime)
  
  return(list(Beta=Beta, Tau=Tau, runtime=runtime, lambdaSeq=lambdaSeq))
  
}
