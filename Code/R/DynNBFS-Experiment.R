## Copyright 2013 Nick Polson, James Scott, and Jesse Windle.

## This file is part of BayesLogit, distributed under the GNU General Public
## License version 3 or later and without ANY warranty, implied or otherwise.


## This is the first Script I wrote to dynamic negative binomial regression.

## Independent AR(1)'s.  Maybe should change this.
source("~/RV-Project/Code/C_Examples/MyLib/Gibbs/Ind/AR1/Stationary/Stationary.R");
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

## FFBS.dyn
FFBS.dyn <- function(y.u, X, tpred, mu, phi, W, m0, C0, tm, Omega, r=NULL)
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
      ## omega.i = Omega[idc];
      ## b.i = rep(0, length(omega.i));
      ## V.i = 1.0 / omega.i
      r.i = r[idc];
      b.i = -1 * NM$m[r.i];
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

  out = list("beta"=beta, "m"=m)
  
  out
} ## FFBS.dyn

## From notes.  lambda is mean of Pois.  psi is log-odds.
draw.indicators <- function(lambda, psi)
{
  ## y.u - N x 1 - latent variable y^u in paper.
  ## lambda = X beta

  res = log(lambda) - psi
  log.wds = log(NM$w) - log(NM$s);

  ## Safer to work on log scale.  Columns correspond to outcome index i!
  ## NOTICE THAT WE ADD m.  This is due to using -1 * EV distribution.
  log.post = -0.5 * outer(1*NM$m, res, "+")^2 / NM$v + log.wds;
  unnrm.post = exp(log.post);

  ## Now sample.
  r = apply(unnrm.post, 2, function(x){sample.int(n=NM$N, size=1, prob=x)})
}  ## draw.indicators

## y is response.  psi is log-odds.  Maybe just include this in MCMC.
draw.lambda <- function(y, d, psi)
{
  alpha = exp(psi)
  lambda = rgamma(y+d, rate=1+1/alpha);
} ## draw.lambda

##------------------------------------------------------------------------------

dyn.NB.mix <- function(y, X, tpred, 
                         samp=1000, burn=100, verbose=100000,
                         m.0=NULL, C.0=NULL,
                         mu.m0=NULL, mu.V0=NULL,
                         phi.m0=NULL, phi.V0=NULL,
                         W.a0=NULL, W.b0=NULL,
                         beta.true=NULL, r.true=NULL, lambda.true=NULL,
                         mu.true=NULL, phi.true=NULL, W.true=NULL, d=1)
{
  ## X: N by P matrix
  ## y: N by 1 vector, avg response
  ## tpred: vector: the vector of observation times.
  ##                We assume that tpred starts at 1 and ends at T.

  ## PREPROCESS ##
  
  X = as.matrix(X)
  N = nrow(X);
  P = ncol(X);
  ## I'm not sure why I had this here.  Maybe something with BayesLogit?
  ## Or so that I could combine times.  I'd say combining things isn't that useful.
  ## X = cbind(X, tpred);

  y.prior = 0.0;
  x.prior = rep(0, P+1);
  n.prior = 0.0;
  
  ## ## Combine data.
  ## ## new.data = logit.combine(y, X, n, y.prior, x.prior, n.prior);
  ## new.data = list("y"=y, "X"=X, "n"=n)
  ## y = new.data$y;
  ## X = as.matrix(new.data$X[,1:P]);
  ## n = new.data$n;
  ## n.prior = 0.0;
  ## tpred = as.matrix(new.data$X[,P+1]);

  ## X = as.matrix(X);

  ## SETUP ##

  N = nrow(X)
  p = ncol(X)
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
    lambda = array(0, dim=c(M, T)),
    r    = array(0, dim=c(M, T)),
    mu   = array(0, dim=c(M, P)),
    phi  = array(0, dim=c(M, P)),
    W    = array(0, dim=c(M, P)),
    tm   = tm
    )

  beta = matrix(0.0, P, T+1);  # even odds.
  r    = matrix(0.0, T);
  mu   = matrix(0.0, P);
  phi  = matrix(0.99, P);
  W    = W.b0 / W.a0;
  lambda = rep(0, N);
  
  ## In case we are doing testing or we want to constrain to local level model.
  if (!is.null(beta.true)) beta = beta.true;
  if (!is.null(mu.true))   mu   = mu.true;
  if (!is.null(phi.true))  phi  = phi.true;
  if (!is.null(W.true))    W    = W.true;
  if (!is.null(r.true))    r    = r.true;
  if (!is.null(lambda.true)) lambda = lambda.true;
  
  start.time = proc.time();
  
  ## SAMPLE ##
  
  for ( j in 1:(samp+burn) )
  {
    ## draw lambda, r
    psi = apply(X * t(matrix(beta[,tpred+1], nrow=P)), 1, sum);
    alpha = exp(psi)

    lambda = rgamma(N, y+d, rate=1+1/alpha);
    r = draw.indicators(lambda, psi);

    ## Draw beta;
    z = log(lambda);
    ffbs = FFBS.dyn(z, X, tpred, mu, phi, diag(W, P), m.0, C.0, tm, NULL, r);
    beta = ffbs$beta;
    
    ## AR(1) - phi, W assumed to be diagonal !!!
    ## mu  = draw.mu.R(beta, phi, W, mu.m0, mu.V0) 
    ## phi = draw.phi.R(beta, mu, W, phi.m0, phi.V0, phi)
    W   = draw.W.R  (beta, mu, phi, W.a0, W.b0)
    
    # Record if we are past burn-in.
    if (j > burn) {
      out$lambda[j-burn,] = lambda
      out$r[j-burn,]      = r;
      out$beta[j-burn,,]  = beta;
      out$mu[j-burn, ]    = mu;
      out$phi[j-burn, ]   = phi;
      out$W[j-burn, ]     = W;
    }

    if (j %% verbose == 0) { cat("Dyn NB Mix: Iteration", j, "\n"); }
  }

  end.time = proc.time();
  diff.time = end.time - start.time;
  out$diff.time = diff.time
  
  out
} ## dyn.NB.mix

################################################################################
                                   ## MAIN ##
################################################################################

if (FALSE) {

  T = 500;
  N = 1;

  beta = array(0, dim=c(N, T+1));
  X = matrix(1, T, N);

  ## Parameters
  W = 0.1;
  phi = 0.95;
  mu = 4.0

  ## Prior
  m0 = mu;
  C0 = 4.0;
  mu.m0 = 0.0;
  mu.V0 = 4.0;
  phi.m0 = 0.9;
  phi.V0 = 0.1;
  W.a0   = 10;
  W.b0   = 1.0;

  ## Synthetic
  beta[,1] = m0;
  for (i in 2:(T+1)) {
    beta[,i] = mu + phi* (beta[,i-1] - mu) + sqrt(W) * rnorm(1);
  }

  psi = apply(X * t(beta)[-1,], 1, sum)
  alpha = exp(psi)

  y = rep(0, T);
  r = rep(0, T);
  lambda = rep(0, T);
  ep = rep(0, T);

  for (i in 1:T) {
    r[i] = sample.int(NM$N, 1);
    ep[i] = rnorm(1, -1*NM$m[r[i]], NM$s[r[i]]);
    lambda[i] = exp(psi[i] + ep[i]);
    y[i] = rpois(1, lambda[i]);
  }

  tpred = 1:T;

  ## Simulate
  source("DynNB-Mixture.R")
  out <- dyn.NB.mix(y, X, tpred, 
                    samp=100, burn=0, verbose=50,
                    m.0=m0, C.0=C0,
                    mu.m0=mu.m0, mu.V0=mu.V0,
                    phi.m0=phi.m0, phi.V0=phi.V0,
                    W.a0=W.a0, W.b0=W.b0,
                    beta.true=NULL, r.true=NULL, lambda.true=NULL,
                    mu.true=mu, phi.true=phi, W.true=NULL)
}

if (FALSE) {

  beta.m = array(0, dim=c(N, T+1));
  for (i in 1:N) {
    beta.m[i,] = apply(out$beta[,i,], 2, mean);
  }

  ymin = min(beta.m, beta);
  ymax = max(beta.m, beta);

  plot(beta[1,], type="l", ylim=c(0, ymax))
  lines(beta.m[1,], col=2)
  points(y, col="gray")
  
}

################################################################################
                                 ## APPENDIX ##
################################################################################

## A note on the FFBS used herein.  This was adapted from <DynLogitPG.R>.  It is
## essentially identical.  I changed it so that it could accomodate the normal
## mixture, but the PG case is just commented out.  It is easy to go between the
## two.

## The original appiclation was binary logistic, but it turns out the FFBS is
## essentially the same.  It runs slow because I have adapted things so that you
## can have missing data.  I did this to try and accomodate the advertising
## dataset, which is irregularly spaced in the ``time'' variable, which I
## believe was duration of the phone call or soemthing like that.

## You can use a local level model by simple setting mu to 0 and phi to 1.
## Otherwise you can use a stationary AR(1) process.
