#' Performs Graphical Oriented LAZy Optimization by optimizing the dual problem.
#'
#' This function implements a simple block-coordinate descent algorithm to find the maximum of the regularized
#' Gaussiann log-likelihood  with  a an assymetric penalty of lasso type.
#' @param S the sample covariance matrix
#' @param L matrix of lower penalties (can be -Inf)
#' @param U matrix of upper penalties (can be Inf)
#' @param tol the convergence tolerance (default tol=1e-7)
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
#' L <- matrix(0,d,d)
#' U <- matrix(0.2,d,d)
#' diag(U) <- 0
#' res <- golazo(R,L=L,U=U)
#' Khat <- res$K
#' print(Khat)
golazo <- function(S,L,U,tol=1e-7,verbose=TRUE){
  d <- nrow(S)
  if (verbose==TRUE){
      cat("** The function maximizes the log-likelihood function penalized with the general GOLAZO penalty.\n")
  }
  # as explained in the paper, without loss of generality big entries in L and U can be thresholded
  if (Inf %in% diag(U)){
    cat("Error: U contains Inf on the diagonal. The optimum does not exist.\n")
    return()
  } else{
    aux <- sqrt(outer(diag(S+U),diag(S+U)))
    # save the diagonals for later
    dU <- diag(U)
    dL <- diag(L)
    U <- pmin(U,-S+aux)
    V <- pmax(L,-S-aux)
    diag(U) <- dU
    diag(L) <- dL
  }
  ### compute the starting point
  S <- S+diag(diag(U))
  # the easy case is when S is bounded away from the boundary of the PSD cone
  if (min(eigen(S)$values)>1e-4){
    if (verbose==TRUE){
      cat("The algorithm is initiated at the input matrix.\n")
    }
    Sig <- S
  } else{
    if (verbose==TRUE){
      cat("The input matrix is not positive definite. Computing the starting point..")
    }
    # this is the GLASSO case
    if ((min(U+diag(d))>0 && max(L-diag(d))<0)){
      Z <- diag(diag(S))
    } else {
      if (min(U+diag(d))==0){
        cat("\n \n **Warning: This combination of L,U (with U having zero off-diagonal entries) is not supported unless S is PD..\n")
        return()
      }
      Z <- Zmatrix(S)
    }
    t <- 1
    # perform simple backtracing to find a feasible point
    while(!(min(L<=round(t*(Z-S),8)) && min(round(t*(Z-S),8)<=U))){
      t <- t/2
    }
    Sig <- (1-t)*S+t*Z
    Sig <- (1-t)*S+t*Z # this is the starting point
    if (verbose==TRUE){
      cat(" DONE\n")
    }
  }
  K <- solve(Sig)
  it <- 0
  if (verbose==TRUE){
    cat("\n The algorithm will stop when the dual gap is below: ",tol,"\b.\n\n")
    cat("Iteration | Dual Gap\n")
  }
  dualgap <- Inf
  while(dualgap > tol){
    for (j in 1:d){
      A <- rbind(diag(d-1),-diag(d-1))
      b <- c(S[j,-j]+L[j,-j],-(S[j,-j]+U[j,-j]))
      active <- which((b>-Inf & b<Inf))
      A <- A[active,]
      b <- b[active]
      y <- quadprog::solve.QP(Dmat=solve(Sig[-j,-j]),dvec=rep(0,d-1),Amat=t(A),bvec=b)$solution
      Sig[j,-j] <- Sig[-j,j] <- y
    }
    it <- it+1
    K <- solve(Sig)
    #roundK <- K * (abs(K)>tol)
    dualgap <- sum(S*K)-d+sum(pmax(L*(abs(K)>tol)*K,U*(abs(K)>tol)*K))
    if (verbose==TRUE){
      cat(it,"\t  | ",dualgap,"\n")
    }
  }
  return(list(K=(K+t(K))/2,Sig=(Sig+t(Sig))/2,it=it))
}

