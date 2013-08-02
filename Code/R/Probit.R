## Copyright 2013 Nick Polson, James Scott, and Jesse Windle.

## This file is part of BayesLogit, distributed under the GNU General Public
## License version 3 or later and without ANY warranty, implied or otherwise.


## library("msm");
## dyn.unload("~/RPackage/BayesBridge/Code/C/Bridge.so");
## dyn.load("~/RPackage/BayesBridge/Code/C/Bridge.so");
## source("~/RPackage/BayesBridge/Code/C/BridgeWrapper.R");

################################################################################
                                  ## PROBIT ##
################################################################################

probit <- function(y, X, samp=1000, burn=100, verbose=100, m.0=NULL, C.0=NULL)
{
  X = as.matrix(X);
  
  N = nrow(X);
  P = ncol(X);
  M = samp;
  
  out = list(
    beta = array(dim=c(M, P)),
    z    = array(dim=c(M, N))
    )

  left  = rep(0, N);  left [y==0]=-Inf;
  right = rep(0, N);  right[y==1]= Inf;

  ## PREPROCESS ##
  if (is.null(C.0)) {
    Prc = t(X) %*% X;
    b.0 = 0.0
  }
  else {
    C.0Inv = solve(C.0);
    b.0    = C.0Inv %*% m.0;
    Prc = t(X) %*% X + C.0Inv;
  }
  
  beta.V = solve(Prc);
  beta.L = t(chol(beta.V));
  
  beta = rep(0.0, P);
  
  for (i in 1:(samp+burn)) {

    psi = X %*% beta;

    ## Draw z.
    z = rtnorm(N, psi, 1.0, left, right);

    ## Draw beta.
    beta.b = t(X) %*% z + b.0;
    beta.m = beta.V %*% beta.b;
    beta = beta.m + beta.L %*% rnorm(P);

    if (i > burn) {
      out$beta[i-burn, ] = beta;
      out$z   [i-burn, ] = z;
    }

    if (i%%verbose==0) {
      cat("Probit: Iteration", i, "\n");
    }
    
  }

  out
}
