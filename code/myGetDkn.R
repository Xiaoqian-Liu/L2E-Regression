#' Compute kth order differencing matrix
#' 
#' @param k order of the differencing matrix
#' @param n Number of time points
#' @export
myGetDkn <- function(k, n) {
  
  D1n<-function(n){         #a function generates 1-st order differnce matrix with n cloumns
    D1n=matrix(0, nrow = (n-1), ncol = n)
    for (i in 1:(n-1)) {
      D1n[i,i]=-1
      D1n[i,i+1]=1
    }
    return(D1n)
  }
  if(k==1){
    return(D1n(n))
  } else{
    Dkn=D1n(n)
    for (i in 1:(k-1)) {
      m=n-i
      Dkn=D1n(m)%*%Dkn
    }
    A=Matrix(Dkn, sparse=TRUE)
    return(A)
  }
}