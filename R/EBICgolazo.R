#' Finds  the optimal penalty parameter for the GOLAZO approach with respect to the EBIC criterion.
#'
#' This function computes the EBIC criterion for a seried of rho parameters. The penalty matrices are \eqn{\rho L}{rho L}, \eqn{\rho U}{rho U},
#' where L and U are fixed in  advance. If L, U are not specified, they will be set to the default value for the positive GOLAZO problem,
#' that is, \eqn{L_{ij}=0}{L_ij=0} and \eqn{U_{ij}=1}{U_ij=1} for all \eqn{i\neq j}{i!=j}.
#' @param S a positive semidefinite  matrix
#' @param n the sample size
#' @param L matrix of lower penalties (may contain -Inf)
#' @param U matrix of upper penalties (may contain Inf)
#' @param tol the convergence tolerance (default tol=1e-8) for the dual gap
#' @param edge.tol the threshold for the normalized entries of K to be treated treated as zero
#' @param gamma the EBIC parameter between 0 and 1 (default gamma=0.5)
#' @param rhomin minimal rho to check
#' @param rhomax maximal rho to check
#' @param nrhos number of rhos in the interval (rhomin,rhomax)
#' @param verbose if TRUE (default) the output will be printed.
#' @return the  list of rhos for which EBIC was computed
#' @return the list of teh corresponding EBIC values
#' @return the optimal rho
#' @keywords EBIC, optimal penalty
#' @export
#' @examples
#' EBICgolazo(diag(10),n=10)
#'
##### Algorithm 3
EBICgolazo <- function(S,n=NULL,L= NULL,U=NULL,tol=1e-6,edge.tol=1e-6,gamma=0.5,rhomin=0.01,rhomax=1,nrhos=10,verbose=TRUE){
  # rhomax=1 makes a lot of sense if S is a correlation matrix
  d <- nrow(S)
  if (is.null(L)||is.null(U)){
      if (verbose==TRUE){
        cat("The value of L and U set to the default value as for the pGOLAZO: L_ij=0, U_ij=1.\n")
      }
     L <- matrix(0,d,d)
     U <- matrix(1,d,d)-diag(d)
  }
  rhos <- seq(from=rhomin,to=rhomax,length.out=nrhos)
  ebic <- rhos
  ebic.gamma <- gamma
  for (rho in rhos){
    LL <- rho*L
    UU <- rho*U
    res <- golazo(S,L=LL,U=UU,tol=tol,verbose=FALSE)
    # compute EBIC
    K <- res$K
    KR <- stats::cov2cor(K) #to make edge count independend of scalings
    nedg <- length(which(abs(KR[upper.tri(abs(KR),diag=FALSE)])> edge.tol))
    ebic[which(rhos==rho)] <- -(n)*(log(det(K))-sum(S*K))+nedg * (log(n) + 4 * ebic.gamma * log(d))
  }
  return(list(rhos=rhos,ebic=ebic,opt.rho=rhos[min(which(ebic==min(ebic)))]))
}

