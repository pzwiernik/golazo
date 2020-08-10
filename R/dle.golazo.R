#' This is a wrapper for the dual maximum likelihood GOLAZO problem.
#'
#' The function simply running golazo() with L_ij=-Inf and U_ij=0 (for certain ij).
#' @param K Positive semidefinite matrix. This will be typically the inverse of the sample covariance matrix or the output of a previous GOLAZO run.
#' @param edges It takes three possible values: "all", "input" (default), or it is equal to an adjacency  matrix. If it is equal to input teh algorithm first finds the non-zero entries of the input matrix K and onnly the corresponding entries of Sigma are assumed nonnegative. If edges is equal to "all" all entries of Sigma are bound to be nnonnegative. This corresponds to computing the MLE under (positively) associated distributions. Finally, edges equal to any adjacency matrix (1 for edges, 0 otherwise) allows to connstraint arbitrary entries of Sigma to be nonnegative.
#' @param tol The convergence tolerance (default tol=1e-7). The algorithm termininnates when teh dual gap (guaranteed to be nonnegative) is less than tol.
#' @param verbose if TRUE (default) the output will be printed.
#' @return Sig the optimal value of the covariannce matrix
#' @return K the optimal value of the concentration matrix
#' @return it the number of iterations
#' @keywords coordinate descent, dual MLE, local association.
#' @export
#' @examples
#' data(ability.cov)
#' S <- ability.cov$cov
#' R <- stats::cov2cor(S)
#' d <- nrow(R)
#' res <- positive.golazo(R,rho=0.1)
#' res2 <- dle.golazo(res$K)
#' print(res2$Sig)
dle.golazo <- function(K,edges="input",tol=1e-7,verbose=TRUE){
  d <- nrow(K)
  U <- matrix(0,d,d)
  if (edges[1]=="input"){
    Rhat <- stats::cov2cor(K)
    diag(Rhat) <- 0
    Rhat <- abs(Rhat)
    edges <- which(Rhat > 1e-6)
    L <- matrix(0,d,d)
    L[edges] <- -Inf
  } else if (edges[1]=="all"){
    L <- matrix(-Inf,d,d)
    diag(U) <- 0
  } else{
    L <- ifelse(edges==1,-Inf,edges)
  }
  res <- golazo(K,L,U,tol=tol,verbose=verbose)
  return(list(Sig=res$K,K=res$Sig,it=res$it))
}

