#' This is a wrapper for the positive GOLAZO problem.
#'
#' This is simply running golazo() with L_ij=0 and U_ij=rho.
#' @param S the sample covariance matrix
#' @param rho the penalty on the positive entries of K (can be Inf)
#' @param tol the convergence tolerance (default tol=1e-7)
#' @param diagonal.pen if FALSE (default) the diagonal of K is not penalized.
#' @param verbose if TRUE (default) the output will be printed.
#' @return the optimal value of the concentration matrix
#' @return the number of iterations the algorithm needed to converge
#' @return the corresponding value of the log-likelihood
#' @keywords coordinate descent, concentration matrix.
#' @export
#' @examples
#' data(ability.cov)
#' S <- ability.cov$cov
#' R <- stats::cov2cor(S)
#' d <- nrow(R)
#' res <- positive.golazo(R,rho=0.1)
#' Khat <- res$K
#' print(Khat)
positive.golazo <- function(S,rho,tol=1e-7,diagonal.pen=FALSE,verbose=TRUE){
  d <- nrow(S)
  L <- matrix(0,d,d)
  U <- matrix(rho,d,d)
  if (diagonal.pen==FALSE){
    diag(U) <- 0
  }
  res <- golazo(S,L,U,tol=tol,verbose=verbose)
  return(res)
}

