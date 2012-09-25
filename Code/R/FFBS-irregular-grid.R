
## This FFBS routine works when there are irregularly spaced observations.

FFBS <- function(y.u, X, tpred, mu, phi, Omega, W, m0, C0, tm)
{
  ## y: vector: y^u
  ## X: matrix: Design Matrix
  ## mu: vector: mean, if it exists.
  ## phi: vector:  autoregression parameter
  ## r: matrix: mixture index.
  ## W: vector: variance of omega_t, i.e. omega_t \sim N(0, diag(W))
  ## m0: vector: mean of beta_0.
  ## C0: matrix: var  of beta_0.
  ## tm: maps time to index.

  K = length(tm);
  T = tpred[tm[K]];
  P = ncol(X);
  N = nrow(X);

  m = array(m0, dim=c(P, T+1));
  C = array(C0, dim=c(P, P, T+1));
  R = array(0., dim=c(P, P, T+1));
  a = array(0., dim=c(P, T+1));
  UR = array(0., dim=c(P, P, T+1));

  beta = array(0, dim=c(P, T+1));

  Phi = diag(phi, P);
  idx = 1;
  
  ## FEED FORWARD ##
  for (i in 1:T) {
    i.h = i+1 # Index of beta (and all other variables that start at 0)
    a[,i.h]  = (1 - phi) * mu + phi * m[,i.h-1];
    R[,,i.h] = Phi %*% (C[,,i.h-1] * phi) + W;
    UR[,,i.h] = chol(R[,,i.h]);

    ## Check if we have an observation.
    if (tpred[tm[idx]]==i) {
      ## Indicies and pick correct parts of y, X, and r.
      n.i = ifelse(idx==K, N+1 - tm[idx], tm[idx+1] - tm[idx]);
      idc = tm[idx]:(tm[idx] + n.i - 1);
      y.i = y.u[idc];
      X.i = as.matrix(X[idc,], nrow=n.i);
      
      ## Mean and Variance of nu -- CAN CHANGE DEPENDING PG OR FS.
      omega.i = Omega[idc];
      b.i = rep(0, length(omega.i));
      V.i = omega.i
      
      ## Forecast
      f.i = X.i %*% a[,i.h] + b.i;
      rho.i = X.i %*% R[,,i.h];
      e.i = y.i - f.i;
      
      ## Joint to Conditional
      ## if (n.i > 1.5 * P) {
        ## tXdivV = t(X.i) / V.i;
        ## M  = chol2inv(UR[,,i.h]) + tXdivV %*% X.i;
        ## L  = t(chol(M));
        ## xi = backsolve(t(L), tXdivV)
        ## QInv = diag(1/V.i, n.i) - t(xi) %*% xi;
        ## A.i  = t(rho.i) %*% QInv;
        ## m[,i.h]  = a[,i.h] + A.i %*% e.i;
        ## C[,,i.h] = R[,,i.h] - A.i %*% rho.i;
      ## }
      ## else {
        Q.i = rho.i %*% t(X.i) + diag(V.i, n.i);
        L  = t(chol(Q.i));
        xi = forwardsolve(L, rho.i);
        tA.i  = backsolve(t(L), xi);
        m[,i.h]  = a[,i.h] + t(tA.i) %*% e.i;
        C[,,i.h] = R[,,i.h] - t(xi) %*% xi;
      ## }

      idx = idx + 1;
    }
    else {

      cat("Warning: missing data at time", i, ".\n");
      m[,i.h]  = a[,i.h]
      C[,,i.h] = R[,,i.h]
      
    }

  }

  ## BACKWARDS SAMPLE ##
  i.h = T+1;
  L = t(chol(C[,,i.h]));
  beta[,i.h] = m[,i.h] + L %*% rnorm(P);

  for (i in (T+1):2) {
    i.y = i-1;
    rho.i = Phi %*% C[,,i-1];
    xi = forwardsolve(t(UR[,,i]), rho.i);
    tB.i = backsolve(UR[,,i], xi);
    e.i  = beta[,i] - a[,i];
    ell  = m[,i-1] + t(tB.i) %*% e.i;
    U    = C[,,i-1] - t(xi) %*% xi;

    L = t(chol(U));
    beta[,i-1] = ell + L %*% rnorm(P);
    
  }

  beta
}
