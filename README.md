# golazo
Flexible regularised likelihood estimation using the GOLAZO approach. 


GOLAZO is a generalized approach for penalized log-likelihood estimation where the l1-penalty used in GLASSO is replaced with a general assymetric type of penalty. More precisely for a non-positive matrix L and a nonnegative matrix U we maximize `logdet(K)-tr(SK)+sum_{i,j} max{L_{ij}K_{ij},U_{ij}K_{ij}}`. The algorithm optimizes the dual problem using a block coordinate descent updates. 


```{r message = FALSE, warning = FALSE}
library(devtools)
install_github("pzwiernik/golazo", build_vignettes=TRUE)
library(golazo)
```

# Basic positive GOLAZO problem

If the dataset that we analyse has some positive dependence structure, it may be good to penalize only negative partial correlations in order not to loose any positive dependence information. Many datasets encountered in psychology have this type of strong positive dependence. 

To have a simple example we load the following "Ability and inteligence" datasets where six tests were given to 112 individuals. The covariance matrix is given in this object.

```{r}
data(ability.cov)
n <- ability.cov$n.obs
S <- ability.cov$cov
R <- stats::cov2cor(S)
d <- nrow(R)
```

The simplest thing to do is to learn the underlying inverse covariance matrix K for a fixed value of the penalty parameter. We use the function `positive.golazo()` with some penalty value.

```{r}
res <- positive.golazo(R,rho=0.01)
Khat <- res$K
print(round(Khat,3))
```


# About the `golazo()` function

The function `positive.golazo()` is only a wrapper that calls function `golazo()`. The function `golazo()` is the main function of this package. Apart from S, its main parameters are the penalty matrices L and U. L has non-positive entries and encodes penalties on negative entries of the innverse covariance matrix. U has non-negative entries and encodes penalties on negative entries of the innverse covariance matrix. 

## Positive GOLAZO again

The previous computation using `positive.golazo()` can be done directly with `golazo()`. By setting L to zero so that negative entries of K (positive partial correlations) are not penalized.

```{r }
rho <- 0.01
L <- matrix(0,d,d)
U <- matrix(rho,d,d)
diag(U) <- 0
res <- golazo(R,L=L,U=U)
print(round(res$K,3))
```

## Graphical LASSO

Graphical LASSO is a special case of GOLAZO. Here we optimize `logdet(K) +tr(SK)+rho sum_{i,j}|K_{ij}|`.

```{r }
L <- matrix(-rho,d,d)
U <- matrix(rho,d,d)
res <- golazo(R,L=L,U=U)
print(round(res$K,3))
```

This version puts the same penalty on all entries of the inverse covariance matrix. This introduces an unnecessary bias. It is better to run a version with no penalty on the diagonal (both L and U are zero on the diagonal). Here we optimize $\log\det K +{\rm tr}(SK)+\rho\sum_{i\neq j}|K_{ij}|$.

```{r }
L <- matrix(-rho,d,d)
U <- matrix(rho,d,d)
diag(U) <- diag(L) <- 0
res <- golazo(R,L=L,U=U)
print(round(res$K,3))
```

## Gaussian graphical models

It is straightforward to set certain entries of the inverse covariance to zero in advance. Here for example we will do this with entries (1,2) and (3,4).

```{r }
L <- matrix(-rho,d,d)
U <- matrix(rho,d,d)
diag(U) <- diag(L) <- 0
U[1,2] <- U[2,1] <- U[3,4] <- U[4,3] <- Inf
L[1,2] <- L[2,1] <- L[3,4] <- L[4,3] <- -Inf
res <- golazo(R,L=L,U=U)
print(round(res$K,3))
```

Our code also computes the MLE under the corresponding graphical model.

```{r }
L <- matrix(0,d,d)
U <- matrix(0,d,d)
U[1,2] <- U[2,1] <- U[3,4] <- U[4,3] <- Inf
L[1,2] <- L[2,1] <- L[3,4] <- L[4,3] <- -Inf
res <- golazo(R,L=L,U=U)
print(round(res$K,3))
```

## Inequality  constraints

Suppose that the goal is to estimate under inequality constraints. For our dataset, a natural assumption is that `K_{ij}<= 0` (K is an M-matrix). Incorporating these inequalities in our case is straightforward and the resulting estimate is the maximum likelihood estimator under the M-matrix constraint. 

```{r }
L <- matrix(0,d,d)
U <- matrix(Inf,d,d)
diag(U) <- 0
res <- golazo(R,L=L,U=U)
print(round(res$K,3))
```

# Positive GOLAZO with EBIC score

It is convenient to have an automatic way of choosing in some sense optimal penalty.  We first compute the positive GOLAZO estimate for a number of possible penalty parameters. Our proposal is to fix $L,U$ and for each $\rho$ we compute the EBIC score for penalties $\rho L$ and $\rho U$. 

Warning: Although the GOLAZO algorithm runs for $>500$ nodes, EBIC reruns this algorithm many times and so it may take a long time to terminate.

We first fix L,U. Here specified for the positive GOLAZO problem.

```{r message=FALSE,results='hide'}
L <- matrix(0,d,d)
U <- matrix(1,d,d)
diag(U) <- 0
```

Now we run the wrapper that calls `golazo` multiple times to find the best  EBIC penalty.

```{r}
ebic <- EBICgolazo(R,n=n,L=L,U=U,rhomin=0.001,rhomax=1,nrhos=10,gamma=0.5)
plot(ebic$rhos,ebic$ebic,type="l",xlab="rho",ylab="EBIC",main="EBIC for gamma=0.5")
rho <- ebic$opt.rho
```

We quickly rerun the algorithm for this optimal $\rho$.

```{r}
res <- positive.golazo(R,rho=rho)
Khat <- res$K
print(round(Khat,3))
```

The estimated matrix is already an inverse M-matrix (all partial correlations are nonnegative) and it is equal to the MLE under the M-matrix constraint as computed above.  

We run the graphical lasso for the same data. The optimal penalty is chosen again with EBIC. 

```{r}
L <- matrix(-1,d,d)
U <- matrix(1,d,d)
diag(U) <- diag(L) <- 0
gebic <- EBICgolazo(R,n=n,L=L,U=U,rhomin=0.001,rhomax=1,nrhos=10,gamma=0.5)
plot(gebic$rhos,gebic$ebic,type="l",xlab="rho",ylab="EBIC",main="EBIC for gamma=0.5")
grho <- gebic$opt.rho
```

Again, rerun the calculations for the optimal rho.

```{r}
L <- matrix(-grho,d,d)
U <- matrix(grho,d,d)
diag(U) <- diag(L) <- 0
res <- golazo(R,L=L,U=U)
Kglasso <- res$K
print(round(Kglasso,3))
```


It turns out that here the estimator is also an M-matrix. The estimated graphs look similar but the positive GOLAZO estimate gives the higher likelihood value.

```{r message=FALSE}
c(log(det(Khat))-sum(diag(R%*%Khat)),log(det(Kglasso))-sum(diag(R%*%Kglasso)))
```

Penalizing for size of the graph the EBIC criterion prefers the positive glasso estimate.

```{r}
# number of edges
Rhat <- cov2cor(Khat)
Rglasso <- cov2cor(Kglasso)
edgesK = length(which(abs(Rhat[upper.tri(abs(Rhat),diag=FALSE)])> 1e-6))
edgesKglasso = length(which(abs(Rglasso[upper.tri(abs(Rglasso),diag=FALSE)])> 1e-6))
#EBIC parameter
ebic.gamma <- 0.5
# EBIC
ebicK <- -(n)*(log(det(Khat))-sum(diag(R%*%Khat)))+edgesK * (log(n) + 4 * ebic.gamma * log(d))
ebicKglasso <- -(n)*(log(det(Kglasso))-sum(diag(R%*%Kglasso)))+edgesKglasso * (log(n) + 4 * ebic.gamma * log(d))
cat("The positive glasso EBIC: ",ebicK,"\n")
cat("The glasso EBIC: ",ebicKglasso,"\n")
```

# Dual MLE estimate

To learn locally associated models, Lauritzen, Zwiernik supplement positive GOLAZO with the second step where the   optimal point $\hat K^\rho$ obtained from positive GOLAZO is used as the input data to compute the dual maximum likelihood estimator under the  constraints $\Sigma_{ij}\geq 0$ for all $ij\in E(\hat G)$ with $\hat G$ being the graph of $\hat K^\rho$. The dual likelihood takes the form $\log\det\Sigma-{\rm tr}(\hat K^\rho \Sigma)$. 

By replacing the  role of $K$ and $\Sigma$, this again can be optimized using the GOLAZO approach. We first show the longer way to obtain this estimator in order to explain what  is under the hood. Then we will present a wrapper that makes these computations easier.

First identify the  graph of $\hat  K^\rho$.

```{r}
Rhat <- stats::cov2cor(Khat)
diag(Rhat) <- 0
Rhat <- abs(Rhat)
edges <- which(Rhat > 1e-6)
```

Then  we create penalty matrix L that  takes the edge constraints into account and run `golazo()`. Since `\hat K^\rho` is already locally associated, this step of the algorithm does nothing.

```{r}
L <- matrix(0,d,d)
L[edges] <- -Inf
U <- matrix(0,d,d)
resDLE <- golazo(Khat,L=L,U=U)
print(resDLE$K)
```

Here, we need to remember that the parameter `resDLE$K` is equal to the  optimal `\check \Sigma` and `resDLE$Sig` is the optimal `\check K`. Since this may  be  confusing, we get another good reason to use the wrapper  `la.golazo()` instead (`la` stands for Locally  Associated). 

One extra parameter of `la.golazo()` is `edges`. If  `edges="input"` then the edges of $\hat G$ are computed as above and the algorithm constrains only  the corresponding entries of $\Sigma$ to be nonnnegative. If  `edges="all"` then  all entries of $\Sigma$ are assumed to be nonnnnegative. In this case the algorithm computes the MLE for (positively) associated Gaussian distributions. Finally, if `edges` is an adjacency  matrix (1 for edge entries and 0 otherwise) an arbitrary  subset of the entries of $\Sigma$ can be constrained to be nonnegative.

```{r}
resDLE <- dle.golazo(Khat)
print(resDLE$Sig)
```


The positive graphical lasso procedure already gave an edge positive distribution so  the second step does not do anything.

# Locally associated distributions in one line
Using the wrappers `positive.golazo()` and `dle.golazo()` we can now easily  compute the DLE for local association. Given a parameter rho (which can be found using EBIC) we get:

```{r}
res <- positive.golazo(R,rho=rho)
resDLE <- dle.golazo(res$K)
```

In our case the dataset had a very strong positive dependence information so each part of the algorithm converged in one step. 


