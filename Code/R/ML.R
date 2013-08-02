## Copyright 2013 Nick Polson, James Scott, and Jesse Windle.

## This file is part of BayesLogit, distributed under the GNU General Public
## License version 3 or later and without ANY warranty, implied or otherwise.


norm.llh <- function(beta, m, V)
{
  ## beta : P x J-1; coefficients, column assumed to be beta_J = 0.
  ## m : P*(J-1); vectorized mean.
  ## V : P*(J-1) x P*(J-1); Variance matrix for vectorized mean.
  ## P : Precision, in case you do not want to use variance.
  beta = as.numeric(beta)
  ev   = eigen(V);

  P = ev$vectors %*% ( t(ev$vectors) / ev$values );  
  llh = -0.5 * sum( log(abs(ev$values)) ) - 0.5 * t(beta - m) %*% P %*% (beta - m);

  llh
}

norm.llh.2 <- function(beta, m, evalues, P)
{
  ## beta : P x J-1; coefficients, column assumed to be beta_J = 0.
  ## m : P*(J-1); vectorized mean.
  ## V : P*(J-1) x P*(J-1); Variance matrix for vectorized mean.
  ## P : Precision, in case you do not want to use variance.
  beta = as.numeric(beta)
  
  llh = -0.5 * sum( log(abs(evalues)) ) - 0.5 * t(beta - m) %*% P %*% (beta - m);

  llh
}

mlogit.llh <- function(beta, y, X, n,
                       m.0=array(0, dim=c(ncol(X), ncol(y))),
                       P.0=array(0, dim=c(ncol(X), ncol(X), ncol(y))))
{
  ## beta : P x J-1; coefficients, column assumed to be beta_J = 0.
  ## y : N x J-1; fraction of responses.
  ## X : N x P; design matrix
  ## P.0 : P x P x J-1; array of matrices for independent prior. 
  beta = matrix(beta, ncol(X), ncol(y));
  XB = X %*% beta;
  ll = sum( n * (rowSums(y * XB) - log( rowSums(exp(XB))+1 )) )
  lr = 0; for (j in 1:ncol(y))
    lr = lr - 0.5 * t(beta[,j] - m.0[,j]) %*% P.0[,,j] %*% (beta[,j] - m.0[,j]);
  ll + lr
}

grad.mlogit.llh <- function(beta, y, X, n,
                       m.0=array(0, dim=c(ncol(X), ncol(y))),
                       P.0=array(0, dim=c(ncol(X), ncol(X), ncol(y))))
{
  beta  = matrix(beta, ncol(X), ncol(y));
  expXB = exp(X %*% beta);
  grad  = t(X) %*% ( y - expXB * n / (rowSums(expXB)+1) );
  grad.lr = matrix(0, ncol(X), ncol(y)); for(j in 1:ncol(y))
    grad.lr[,j] = -P.0[,,j] %*% (beta[,j] - m.0[,j]);
  grad = grad + grad.lr
    
  as.numeric(grad)
}

hessian.mlogit.llh <- function(beta, y, X, n, P.0=array(0, dim=c(ncol(X), ncol(X), ncol(y))))
{
  beta  = matrix(beta, ncol(X), ncol(y));
  P = ncol(X);
  J = ncol(y) + 1;
  hess = array(0, dim=c(P, J-1, P, J-1));

  XB = X %*% beta
  expXB = exp(XB)
  sum.expXB = rowSums(expXB) + 1;

  for (p in 1:(J-1)) {
    for (j in 1:(J-1)) {

      hess[,j,,p] = t(X) %*% (X * expXB[,p] * expXB[,j] * n / sum.expXB^2)

      if (p == j) {
        hess[,j,,p] = hess[,j,,p] - t(X) %*% (X * expXB[,j] * n / sum.expXB) - P.0[,,j]
      }
            
    }
  }

  hess.vec = array(hess, dim=c(P*(J-1), P*(J-1)));
}

ind.metroplis <- function(y, X, n, m0, P0, m, V, samp=1000, burn=100, tune = 0.25, verbose=1000)
{
  evd = eigen(V);
  evalues = evd$values;
  Prec = evd$vectors %*% ( t(evd$vectors) / evd$values );  
  L  = chol(V);

  N = nrow(X)
  P = ncol(X)
  J = ncol(y) + 1
  PJ1 = P * (J - 1);

  output = list(
    beta = array(0, dim=c(samp, P*(J-1))),
    alpha = rep(0, samp)
    )

  beta = output$beta[1,]
  
  for (i in 1:(samp+burn)) {

    ppsl = m + tune * L %*% rnorm(PJ1);
    ## log.ratio = mlogit.llh(ppsl, y, X, n, m0, P0) + norm.llh.2(beta, m, evalues, Prec) -
    ##   mlogit.llh(beta, y, X, n, m0, P0) - norm.llh.2(ppsl, m, evalues, Prec);
    log.ratio = mlogit.llh(ppsl, y, X, n, m0, P0) + norm.llh(beta, m, V) -
      mlogit.llh(beta, y, X, n, m0, P0) - norm.llh(ppsl, m, V);
    alpha = min(exp(log.ratio), 1);
    if (runif(1) < alpha) {
      beta = ppsl
    }
    
    if (i > burn) {
      output$beta[i-burn,] = beta;
      output$alpha[i-burn] = alpha
    }

    if (i %% verbose == 0) {
      if (i > burn) cat("Ave alpha:", mean(output$alpha[1:(i-burn)]), ", ");
      cat("Iteration:", i, "\n");
    }
  }

  output
}

################################################################################
                                   ## MAIN ##
################################################################################

y = as.matrix(y)
X = as.matrix(X)

P = ncol(X);
J = ncol(y) + 1;

beta.0 = matrix(0, P, J-1);

m.0 = array(0, dim=c(P, J-1));
P.0 = array(diag(0.0001, P), dim=c(P, P, J-1));

optim.out = optim(beta.0, mlogit.llh, gr=grad.mlogit.llh, method="BFGS", hessian=TRUE,
                  y=y, X=X, n=n, P.0=P.0, control=list(fnscale=-1));

beta.ml = matrix(optim.out$par, P, J-1);

hess.num = optim.out$hessian
hess.pen = hessian.mlogit.llh(beta.ml, y, X, n, P.0)

V.num = solve(-1*optim.out$hessian);
V.pen = solve(-1*optim.out$hessian);
Chol.ml = chol(V.pen);

m = as.numeric(beta.ml);
V = V.pen;

mh = ind.metroplis(y, X, n, m.0, P.0, m, V, samp=10000, burn=1000, tune=0.1)
beta.met.ave = colMeans(mh$beta);
beta.met.ave = array(beta.met.ave, dim=c(P, J-1));
