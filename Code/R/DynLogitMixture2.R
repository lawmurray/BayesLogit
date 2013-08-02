## Copyright 2013 Nick Polson, James Scott, and Jesse Windle.

## This file is part of BayesLogit, distributed under the GNU General Public
## License version 3 or later and without ANY warranty, implied or otherwise.

## This is a variant upon FS&F (2007) where we allow for multiple observations
## of a binary response at each time step.

## Independent AR(1)'s.  Maybe should change this.
source("Stationary.R");
source("LogitByMixture.R")

## Define normal mixture.  FS&F (2007) p. 3511.
normal.mixture = list(
  w = c(0.00397, 0.0396, 0.168, 0.147, 0.125, 0.101, 0.104, 0.116, 0.107, 0.088),
  m = c(5.09, 3.29, 1.82, 1.24, 0.764, 0.391, 0.0431, -0.306, -0.673, -1.06),
  v = c(4.50, 2.02, 1.10, 0.422, 0.198, 0.107, 0.0778, 0.0766, 0.0947, 0.146)
  )
  normal.mixture$s = sqrt(normal.mixture$v)
  normal.mixture$N = length(normal.mixture$w)

ev.samp = -log(rexp(10000));
normal.mixture$mar.mean = mean(ev.samp);
normal.mixture$mar.var  = var(ev.samp);

NM = normal.mixture;

FFBS <- function(y.u, X, tpred, mu, phi, r, W, m0, C0, tm)
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
      X.i = matrix(X[idc,], nrow=n.i);
      
      ## Mean and Variance of nu -- CAN CHANGE DEPENDING PG OR FS.
      r.i = r[idc];
      b.i = NM$m[r.i];
      V.i = NM$v[r.i];
      
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

      ## cat("Warning: missing data at time", i, ".\n");
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

  out = list("beta"=beta, "m"=m);
}

dyn.logit.mix <- function(y, X, tpred, samp=1000, burn=100, verbose=10000,
                          m.0=NULL, C.0=NULL,
                          mu.m0=NULL, mu.V0=NULL,
                          phi.m0=NULL, phi.V0=NULL,
                          W.a0=NULL, W.b0=NULL,
                          beta.true=NULL, y.u.true=NULL, r.true=NULL,
                          mu.true=NULL, phi.true=NULL, W.true=NULL)
{
  ## X: design matrix
  ## tpred: time

  ## Assumptions:
  ## tpred is ordered from smallest to largest with tpred[1] = 1 and tpred[N] = T.
  ## X must match the ordering of tpred.
  
  N = nrow(X)
  P = ncol(X)
  M = samp
  T = tpred[N]

  ## Set up tm, which maps from time back to the first instance of that time
  ## found in tpred.
  tm = rep(NA, N);
  prev = tpred[1]-1;
  for (i in 1:N) {
    if (tpred[i] != prev) tm[i] = i;
    prev = tpred[i]
  }
  tm = tm[!is.na(tm)];
  
  ## Default prior parameters.
  if (is.null(m.0))    m.0    = rep(0.0, P);
  if (is.null(C.0))    C.0    = diag(1.0, P);
  if (is.null(mu.m0))  mu.m0  = rep(0.0 , P);
  if (is.null(mu.V0))  mu.V0  = rep(0.01, P);
  if (is.null(phi.m0)) phi.m0 = rep(0.99, P);
  if (is.null(phi.V0)) phi.V0 = rep(0.01, P);
  if (is.null(W.a0))   W.a0   = rep(1.0, P);
  if (is.null(W.b0))   W.b0   = rep(1.0, P);

  ## Output data structure.
  out = list(
    beta = array(0, dim=c(M, P, T+1)),
    y.u  = array(0, dim=c(M, N)),
    r    = array(0, dim=c(M, N)),
    mu   = array(0, dim=c(M, P)),
    phi  = array(0, dim=c(M, P)),
    W    = array(0, dim=c(M, P)),
    tm   = tm
    )

  beta = matrix(0.0, P, T+1);  # even odds.
  y.u  = matrix(0.0, N);
  r    = matrix(1.0, N);
  mu   = matrix(0.0, P);
  phi  = matrix(0.99, P);
  W    = W.b0 / W.a0;
  
  ## Always be careful with indicators.  You don't want to start with a very
  ## unlikely value (given the data).  It doesn't matter here because we can
  ## sample y.u without respect to r to start.  Nonetheless, seed randomly out
  ## of principle.
  r = sample.int(NM$N, N, prob=NM$w, replace=TRUE);

  ## In case we are doing testing or we want to constrain to local level model.
  if (!is.null(beta.true)) beta = beta.true;
  if (!is.null(y.u.true))  y.u  = y.u.true;
  if (!is.null(r.true))    r    = r.true;
  if (!is.null(mu.true))   mu   = mu.true;
  if (!is.null(phi.true))  phi  = phi.true;
  if (!is.null(W.true))    W    = W.true;

  start.time = proc.time();
  
  for (i in 1:(samp+burn)) {

    ## Use tpred + 1 because we keep beta[0];
    f = apply(X * t(matrix(beta[,tpred+1], nrow=P)), 1, sum);
    lambda = exp(f);

    ## (y.u, r | beta)  = (y.u | beta) (r | y.u, beta)
    ## (y.u | beta, r) != (y.u | beta)
    y.u  = draw.utility(y, lambda)
    r    = draw.indicators(y.u, lambda)

    ## Draw beta;
    ffbs = FFBS(y.u, X, tpred, mu, phi, r, diag(W, P), m.0, C.0, tm);
    beta = ffbs$beta;
    
    ## AR(1) - phi, W assumed to be diagonal !!!
    ## mu  = draw.mu.R(beta, phi, W, mu.m0, mu.V0) 
    ## phi = draw.phi.R(beta, mu, W, phi.m0, phi.V0, phi)
    W   = draw.W.R  (beta, mu, phi, W.a0, W.b0)

    ## cat("phi:", phi, "W:", W, "\n");
    
    if (i > burn) {
      out$beta[i-burn,,]  = beta;
      out$y.u[i-burn,]    = y.u;
      out$r[i-burn,]      = r;
      out$mu[i-burn, ]    = mu;
      out$phi[i-burn, ]   = phi;
      out$W[i-burn, ]     = W;
    }

    if (i %% verbose == 0) cat("Iteration", i, "\n");    
    
  }

  end.time = proc.time();
  diff.time = end.time - start.time;
  out$diff.time = diff.time
  
  out
}

################################################################################
                                  ## TEST 1 ##
################################################################################

## SYNTHETIC DATA ##
if (FALSE) {

  T = 100;
  N = 1;

  beta = array(0, dim=c(N, T+1));
  X = matrix(1, T, N);
  tpred = 1:T;

  ## Parameters
  W   = 0.1;
  mu  = 0.0;
  phi = 0.95

  ## Prior
  b.m0 = 0.0;
  b.C0 = 2.0;
  phi.m0 = 0.9
  phi.V0 = 0.1;
  W.a0   = 100;
  W.b0   = W.a0 * W;

  ## Synthetic
  beta[,1] = 0.0;
  for (i in 2:(T+1)) {
    beta[,i] = phi * beta[,i-1] + sqrt(W) * rnorm(1);
  }

  f = apply(t(X) * beta[,tpred+1], 2, sum);
  lambda = exp(f);

  r    = sample.int(NM$N, T, prob=NM$w, replace=TRUE);
  ep   = rnorm(T, NM$m[r], NM$s[r]);
  y.u  = log(lambda) + ep;
  y.u0 = -1 * log(rexp(T));
  y    = as.numeric(y.u > y.u0);

  ## Simulate
  out = dyn.logit.mix(y, X, tpred, samp=600, burn=0,
                      m.0=b.m0, C.0=b.C0,
                      mu.m0=NULL, mu.V0=NULL,
                      phi.m0=NULL, phi.V0=NULL,
                      W.a0=W.a0, W.b0=W.b0,
                      verbose=100,
                      beta.true=NULL,
                      y.u.true=NULL, r.true=NULL, # Both must be set or NULL.
                      mu.true=0.0, phi.true=1.0, W.true=NULL)
  
}

## PLOT ##
if (FALSE) {

  beta.mean = apply(out$beta, c(2,3), mean);
  adj.y.u   = y.u - NM$m[r];

  y.u.mean  = apply(out$y.u, 2, mean);
  y.u.upper = apply(out$y.u, 2, function(x){quantile(x, 0.95)});
  y.u.lower = apply(out$y.u, 2, function(x){quantile(x, 0.05)});
  
  ymin = min(beta.mean, beta, adj.y.u);
  ymax = max(beta.mean, beta, adj.y.u);
  
  plot(0:T, beta[1,], type="l", ylim=c(ymin,ymax));
  points(0:T, beta[1,])
  
  points(y*(ymax-ymin) + ymin, cex=0.1);

  lines(y.u.upper - NM$m[r], col="gray");
  lines(y.u.lower - NM$m[r], col="gray");
  lines(y.u.mean - NM$m[r], col=5);
  
  abline(h=mean(adj.y.u), col="gray");
  lines(0:T, beta.mean[1,], col=2);
  abline(h=mean(beta.mean[1,]), col=2, lty=c(2,2));
  
  lines(0:T, out$beta[50,1,], col=3);
  
}

if (FALSE) {

  beta.mean = apply(out$beta, c(2,3), mean);
  beta.95   = apply(out$beta, c(2,3), function(x){quantile(x, 0.95)});
  beta.05   = apply(out$beta, c(2,3), function(x){quantile(x, 0.05)});

  y.u.mean  = apply(out$y.u, 2, mean);
  y.u.upper = apply(out$y.u, 2, function(x){quantile(x, 0.95)});
  y.u.lower = apply(out$y.u, 2, function(x){quantile(x, 0.05)});
  
  ymin = min(beta.mean, out$beta, y.u.lower);
  ymax = max(beta.mean, out$beta, y.u.upper);

  plot(0:T, beta.mean[1,], col=2, type="l", ylim=c(ymin,ymax));
  lines(0:T, beta.95[1,], col="pink")
  lines(0:T, beta.05[1,], col="pink")
  abline(h=mean(beta.mean[1,]), col=2, lty=c(2,2));

  lines(1:T, y.u.mean, col=4);
  lines(1:T, y.u.lower, col="gray");
  lines(1:T, y.u.upper, col="gray");

  points(1:T, y*(ymax-ymin) + ymin, cex=0.1);

  ## lines(y.u, col="gray")

  lines(0:T, out$beta[300,1,], col=3);
  lines(1:T, out$y.u[300,], col=5);

}

################################################################################
                                  ## TEST 2 ##
################################################################################

## SYNTHETIC DATA ##
if (FALSE) {

  T = 100;
  P = 1;
  N = 2*T;

  beta = array(0, dim=c(P, T+1));
  X = matrix(1, 2*T, P);
  tpred = ceiling(seq(0.5,T,0.5));
  
  ## Parameters
  W   = 0.1;
  mu  = 0.0;
  phi = 0.95

  ## Prior
  b.m0 = 0.0;
  b.C0 = 2.0;
  phi.m0 = 0.9
  phi.V0 = 0.1;
  W.a0   = 10;
  W.b0   = W.a0 * W;

  ## Synthetic
  beta[,1] = 0.0;
  for (i in 2:(T+1)) {
    beta[,i] = phi * beta[,i-1] + sqrt(W) * rnorm(1);
  }

  f = apply(X * t(matrix(beta[,tpred+1], nrow=P)), 1, sum);
  lambda = exp(f);

  r    = sample.int(NM$N, N, prob=NM$w, replace=TRUE);
  ep   = rnorm(N, NM$m[r], NM$s[r]);
  y.u  = log(lambda) + ep;
  y.u0 = -1 * log(rexp(N));
  y    = as.numeric(y.u > y.u0);

  ## Simulate
  out = dyn.logit.mix(y, X, tpred, samp=500, burn=0,
                      m.0=b.m0, C.0=b.C0,
                      mu.m0=NULL, mu.V0=NULL,
                      phi.m0=NULL, phi.V0=NULL,
                      W.a0=W.a0, W.b0=W.b0,
                      verbose=100,
                      beta.true=NULL,
                      y.u.true=NULL, r.true=NULL, # Both must be set or NULL.
                      mu.true=0.0, phi.true=1.0, W.true=NULL)
  
}

## PLOT ##
if (FALSE) {

  beta.mean = apply(out$beta, c(2,3), mean);
  adj.y.u   = y.u - NM$m[r];

  y.u.mean  = apply(out$y.u, 2, mean);
  y.u.upper = apply(out$y.u, 2, function(x){quantile(x, 0.95)});
  y.u.lower = apply(out$y.u, 2, function(x){quantile(x, 0.05)});
  
  ymin = min(beta.mean, beta, adj.y.u);
  ymax = max(beta.mean, beta, adj.y.u);
  
  plot(0:T, beta[1,], type="l", ylim=c(ymin,ymax));
  points(0:T, beta[1,])
  
  points(seq(0.5,T,0.5), y*(ymax-ymin) + ymin, cex=0.1);

  lines(seq(0.5,T,0.5), y.u.upper - NM$m[r], col="gray");
  lines(seq(0.5,T,0.5), y.u.lower - NM$m[r], col="gray");
  lines(seq(0.5,T,0.5), y.u.mean - NM$m[r], col=5);
  
  abline(h=mean(adj.y.u), col=4);
  lines(0:T, beta.mean[1,], col=2);
  points(0:T, beta.mean[1,], col=2)
  abline(h=mean(beta.mean[1,]), col=2, lty=c(2,2));
  
  lines(0:T, out$beta[50,1,], col=3);
  
}

if (FALSE) {

  beta.mean = apply(out$beta, c(2,3), mean);
  beta.95   = apply(out$beta, c(2,3), function(x){quantile(x, 0.95)});
  beta.05   = apply(out$beta, c(2,3), function(x){quantile(x, 0.05)});

  y.u.mean  = apply(out$y.u, 2, mean);
  y.u.upper = apply(out$y.u, 2, function(x){quantile(x, 0.95)});
  y.u.lower = apply(out$y.u, 2, function(x){quantile(x, 0.05)});
  
  ymin = min(beta.mean, out$beta, y.u.lower);
  ymax = max(beta.mean, out$beta, y.u.upper);

  plot(0:T,   beta.mean[1,], col=2, type="l", ylim=c(ymin,ymax));
  lines(0:T,  beta.95[1,], col="pink")
  lines(0:T,  beta.05[1,], col="pink")
  abline(h=mean(beta.mean[1,]), col=2, lty=c(2,2));

  lines(seq(0.5,T,0.5), y.u.mean, col=4);
  lines(seq(0.5,T,0.5), y.u.lower, col="gray");
  lines(seq(0.5,T,0.5), y.u.upper, col="gray");

  points(seq(0.5,T,0.5), y*(ymax-ymin) + ymin, cex=0.1);

  ## lines(y.u, col="gray")

  lines(0:T, out$beta[300,1,], col=3);
  lines(seq(0.5,T,0.5), out$y.u[300,], col=5);
  
}

################################################################################
                              ## TOKYO RAINFALL ##
################################################################################

## TOKYO RAINFALL ##
if (FALSE) {
  tkrain = as.numeric(na.omit(as.numeric(as.matrix(tokyo))))
  
  ## tkrain = tkrain[1:100]
  T = length(tkrain);
  df = data.frame("Time"=c(1:T,1:T), "y"=c(as.numeric(tkrain>=1), as.numeric(tkrain>=2)));
  df = df[order(df$Time),]

  N = nrow(df)
  tpred = df$Time;
  y = df$y
  X = matrix(1, N, 1);

  ## Prior
  b.m0 = -1.0;
  b.C0 = 1.0;
  W      = 0.1;
  W.a0   = 300;
  W.b0   = W.a0 * W;

  out = dyn.logit.mix(y, X, tpred, samp=100, burn=100,
                      m.0=b.m0, C.0=b.C0,
                      mu.m0=NULL, mu.V0=NULL,
                      phi.m0=NULL, phi.V0=NULL,
                      W.a0=W.a0, W.b0=W.b0,
                      verbose=100,
                      beta.true=NULL,
                      y.u.true=NULL, r.true=NULL, # Both must be set or NULL.
                      mu.true=0.0, phi.true=1.0, W.true=NULL)
  
}

if (FALSE) {

  beta.mean = apply(out$beta, c(2,3), mean);
  beta.95   = apply(out$beta, c(2,3), function(x){quantile(x, 0.95)});
  beta.05   = apply(out$beta, c(2,3), function(x){quantile(x, 0.05)});

  y.u.mean  = apply(out$y.u, 2, mean);
  y.u.upper = apply(out$y.u, 2, function(x){quantile(x, 0.95)});
  y.u.lower = apply(out$y.u, 2, function(x){quantile(x, 0.05)});
  
  ymin = min(beta.mean, out$beta, y.u.lower);
  ymax = max(beta.mean, out$beta, y.u.upper);

  plot(0:T, beta.mean[1,], col=2, type="l", ylim=c(ymin,ymax));
  points(0:T, beta.mean[1,], col=2);
  lines(0:T, beta.95[1,], col="pink")
  lines(0:T, beta.05[1,], col="pink")
  abline(h=mean(beta.mean[1,]), col=2, lty=c(2,2));

  lines(seq(0.5,T,0.5), y.u.mean, col=4);
  lines(seq(0.5,T,0.5), y.u.lower, col="gray");
  lines(seq(0.5,T,0.5), y.u.upper, col="gray");

  points(seq(0.5,T,0.5), y*(ymax-ymin) + ymin, cex=0.1);

  ## lines(y.u, col="gray")

  lines(0:T, out$beta[300,1,], col=3);
  ## lines(seq(0.5,T,0.5), out$y.u[300,], col=5);
  
}

if (FALSE) {

  lambda = exp(out$beta[,1,]);
  prob = lambda / (1 + lambda);
  
  prob.mean = apply(prob, 2, mean);
  prob.95   = apply(prob, 2, function(x){quantile(x, 0.95)});
  prob.05   = apply(prob, 2, function(x){quantile(x, 0.05)});
  
  ymin = min(prob.mean, prob)
  ymax = max(prob.mean, prob)

  plot(0:T, prob.mean, col=2, type="l", ylim=c(ymin,ymax));
  lines(0:T, prob.95, col="pink")
  lines(0:T, prob.05, col="pink")
  abline(h=mean(prob.mean), col=2, lty=c(2,2));

  points(seq(0.5,T,0.5), y*(ymax-ymin) + ymin, cex=0.1);
  
  ## lines(0:T, prob[50,], col=3);
  
}

################################################################################
                          ## TRIAL TOKAY RAINFALL 2 ##
################################################################################

## TOKYO RAINFALL - Just do if it rained or not ##
if (FALSE) {
  tkrain = as.numeric(na.omit(as.numeric(as.matrix(tokyo))))
  
  tkrain = tkrain[1:100]
  T = length(tkrain);
  df = data.frame("Time"=c(1:T), "y"=c(as.numeric(tkrain>=1)));
  df = df[order(df$Time),]

  N = nrow(df)
  tpred = df$Time;
  y = df$y
  X = matrix(1, N, 1);

  ## Prior
  b.m0 = 0.0;
  b.C0 = 1.0;
  W      = 1.0;
  W.a0   = 10;
  W.b0   = W.a0 * W;

  out = dyn.logit.mix(y, X, tpred, samp=300, burn=0,
                      m.0=b.m0, C.0=b.C0,
                      mu.m0=NULL, mu.V0=NULL,
                      phi.m0=NULL, phi.V0=NULL,
                      W.a0=W.a0, W.b0=W.b0,
                      verbose=100,
                      beta.true=NULL,
                      y.u.true=NULL, r.true=NULL, # Both must be set or NULL.
                      mu.true=0.0, phi.true=1.0, W.true=NULL)
  
}

if (FALSE) {

  beta.mean = apply(out$beta, c(2,3), mean);
  beta.95   = apply(out$beta, c(2,3), function(x){quantile(x, 0.95)});
  beta.05   = apply(out$beta, c(2,3), function(x){quantile(x, 0.05)});

  y.u.mean  = apply(out$y.u, 2, mean);
  y.u.upper = apply(out$y.u, 2, function(x){quantile(x, 0.95)});
  y.u.lower = apply(out$y.u, 2, function(x){quantile(x, 0.05)});
  
  ymin = min(beta.mean, out$beta, y.u.lower);
  ymax = max(beta.mean, out$beta, y.u.upper);

  plot(0:T, beta.mean[1,], col=2, type="l", ylim=c(ymin,ymax));
  lines(0:T, beta.95[1,], col="pink")
  lines(0:T, beta.05[1,], col="pink")
  abline(h=mean(beta.mean[1,]), col=2, lty=c(2,2));

  lines(1:T, y.u.mean, col=4);
  lines(1:T, y.u.lower, col="gray");
  lines(1:T, y.u.upper, col="gray");

  points(1:T, y*(ymax-ymin) + ymin, cex=0.1);

  ## lines(y.u, col="gray")

  lines(0:T, out$beta[300,1,], col=3);
  lines(1:T, out$y.u[300,], col=5);
  
}


if (FALSE) {

  lambda = exp(out$beta[,1,]);
  prob = lambda / (1 + lambda);
  
  prob.mean = apply(prob, 2, mean);
  prob.95   = apply(prob, 2, function(x){quantile(x, 0.95)});
  prob.05   = apply(prob, 2, function(x){quantile(x, 0.05)});
  
  ymin = min(prob.mean, prob)
  ymax = max(prob.mean, prob)

  plot(0:T, prob.mean, col=2, type="l", ylim=c(ymin,ymax));
  lines(0:T, prob.95, col="pink")
  lines(0:T, prob.05, col="pink")
  abline(h=mean(prob.mean), col=2, lty=c(2,2));

  points(y*(ymax-ymin) + ymin ~ Time, data=df);
  
  ## lines(y.u, col="gray")

  lines(0:T, prob[50,], col=3);
  
}

################################################################################
                                ## Marketing ##
################################################################################

if (FALSE) {

  bank = read.csv("~/Downloads/bank.csv", sep=";")
  bank = bank[order(bank$duration),];

  y = as.numeric(bank$y=="yes")
  X = model.matrix(y ~ loan, data=bank);
  tpred = bank$duration-3; ## Simple fix for now.

  N = nrow(X);
  P = ncol(X);
  T = tpred[N]
  
  ## Prior
  b.m0 = rep(-1.0, P);
  b.C0 = diag(1.0, P);
  W      = 0.01;
  W.a0   = 300;
  W.b0   = W.a0 * W;

  out = dyn.logit.mix(y, X, tpred, samp=100, burn=0, verbose=10,
                      m.0=b.m0, C.0=b.C0,
                      mu.m0=NULL, mu.V0=NULL,
                      phi.m0=NULL, phi.V0=NULL,
                      W.a0=W.a0, W.b0=W.b0,
                      beta.true=NULL,
                      y.u.true=NULL, r.true=NULL, # Both must be set or NULL.
                      mu.true=0.0, phi.true=1.0, W.true=NULL)
  
}

if (FALSE) {

  beta.mean = apply(out$beta, c(2,3), mean);
  beta.95   = apply(out$beta, c(2,3), function(x){quantile(x, 0.95)});
  beta.05   = apply(out$beta, c(2,3), function(x){quantile(x, 0.05)});
  
  ymin = min(beta.mean, out$beta);
  ymax = max(beta.mean, out$beta);
  
  par(mfrow=c(2,1))

  for (j in 1:2) {
  
    plot(0:T, beta.mean[j,], col=2, type="l", ylim=c(ymin,ymax));
    lines(0:T, beta.95[j,], col="pink")
    lines(0:T, beta.05[j,], col="pink")
    abline(h=mean(beta.mean[j,]), col=2, lty=c(2,2));
    lines(0:T, out$beta[100,j,], col=3);
    
    points(jitter(tpred), y*(ymax-ymin) + ymin, cex=0.1);
    
  }
  
}

if (FALSE) {

  lambda = exp(out$beta[,1,]);
  prob = lambda / (1 + lambda);
  
  prob.mean = apply(prob, 2, mean);
  prob.95   = apply(prob, 2, function(x){quantile(x, 0.95)});
  prob.05   = apply(prob, 2, function(x){quantile(x, 0.05)});
  
  ymin = min(prob.mean, prob)
  ymax = max(prob.mean, prob)

  par(mfrow=c(1,1));
  
  plot(0:T, prob.mean, col=2, type="l", ylim=c(ymin,ymax));
  lines(0:T, prob.95, col="pink")
  lines(0:T, prob.05, col="pink")
  abline(h=mean(prob.mean), col=2, lty=c(2,2));

  points(jitter(tpred), y*(ymax-ymin) + ymin, cex=0.1);
  
}
