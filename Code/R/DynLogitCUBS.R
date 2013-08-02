## Copyright 2013 Nick Polson, James Scott, and Jesse Windle.

## This file is part of BayesLogit, distributed under the GNU General Public
## License version 3 or later and without ANY warranty, implied or otherwise.

binom.obs.dens <- function(y, n, psi, log.sc=TRUE)
{
  ## log.dens =  y * psi - n * log(1+exp(psi));
  ## out = ifelse(log.sc, log.dens, exp(log.dens)) # WRONG USE OF ifelse
  p = 1 / (1 + exp(-1 * psi))  
  dbinom(y, n, p, log=log.sc);
}

target.dens <- function(y, X, n, beta, mu, phi, W, m0, C0, log.sc=FALSE, alpha=NULL, obs="binom")
{
  ## The first rows of X = X.stc, X.dyn
 
  T   = length(y)
  N.b = nrow(as.matrix(W))
  N   = length(m0)
  N.a = N - N.b
  b.idc = 1:N.b + N.a
  a.idc = 1:N.a
  
  X.dyn = matrix(X[,b.idc], ncol=N.b);
  if (N.a > 0) X.stc = matrix(X[,a.idc], ncol=N.a);
  
  psi = apply(X.dyn * t(beta)[-1,], 1, sum);
  if (N.a > 0) psi = psi + X.stc %*% alpha;

  ## ar1.llh = ar1.dens(beta, mu, phi, W, m0, C0, log.sc=TRUE, alpha);
  ## obs.llh = sum(binom.obs.dens(y, n, psi, log.sc=TRUE))
  ## obs.llh = sum(gauss.obs.dens(y, n, psi, log.sc=TRUE))
  ## obs.llh = sum(poisson.obs.dens(y, psi, log.sc=TRUE))
  ## obs.llh = sum(neg.binom.obs.dens(y, n, psi, log.sc=TRUE))
  
  ar1.llh <- ar1.llh.C(beta, mu, phi, W, m0, C0, alpha);
  if (is.na(ar1.llh)) cat("ar1.llh is na.\n");
  ## Using a switch statement slowed things down by 3.5%.
  obs.llh <- switch(obs,
                    binom  = sum(binom.obs.dens(y, n, psi, log.sc=TRUE)),
                    nbinom = sum(nbinom.obs.dens(y, n, psi, log.sc=TRUE)),
                    norm   = sum(gauss.obs.dens(y, n, psi, log.sc=TRUE)))
  if (is.na(obs.llh)) cat("obs.llh is na.\n");

  ## obs.llh = do.call("binom.obs.dens", list(y, n, psi, log.sc=TRUE))

  llh = ar1.llh + obs.llh;
  out = ifelse(log.sc, llh, exp(llh))
  out
}

##------------------------------------------------------------------------------

dyn.logit.CUBS <- function(y, X.dyn, n, m0, C0,
                           samp=1000, burn=100, verbose=100,
                           mu.m0=NULL, mu.P0=NULL,
                           phi.m0=NULL, phi.P0=NULL,
                           W.a0=NULL, W.b0=NULL, X.stc=NULL,
                           mu.true = NULL, phi.true=NULL, W.true=NULL)
{
  ## Dim and restructure.
  X   = cbind(X.stc, X.dyn)
  C0  = as.matrix(C0);
  T   = length(y)
  N.b = ncol(X.dyn)
  N   = ncol(X)
  N.a = N - N.b
  M   = samp

  ## Default prior parameters -- almost a random walk for beta ##
  if (is.null(m0)     || is.null(C0))     { m0     = rep(0.0, N)   ; C0     = diag(1.0, N  ); }
  if (is.null(mu.m0)  || is.null(mu.P0))  { mu.m0  = rep(0.0 ,N.b) ; mu.P0  = rep(100 , N.b); }
  if (is.null(phi.m0) || is.null(phi.P0)) { phi.m0 = rep(0.99,N.b) ; phi.P0 = rep(100 , N.b); }
  if (is.null(W.a0)   || is.null(W.b0))   { W.a0   = rep(1.0, N.b) ; W.b0   = rep(1.0,  N.b);  }
  
  ## Output
  out <- list("beta" = array(0, dim=c(M, N.b, T+1)),
              "mu"   = array(0, dim=c(M, N.b)),
              "phi"  = array(0, dim=c(M, N.b)),
              "W"    = array(0, dim=c(M,N.b)),
              ## "ppsl.b"=array(0, dim=c(M, N.b, T+1)),
              "l.fdivq" = rep(0, M),
              "l.ratio" = rep(0, M),
              "lf.ppsl" = rep(0, M),
              "lq.ppsl" = rep(0, M),
              "ac.rate" = rep(0, M))
  if (N.a > 0) out$alpha=array(0, dim=c(M, N.a))

  ## Initialize
  beta    = array(0, dim=c(N.b, T+1));
  l.fdivq = -Inf;
  naccept = 0
  if (N.a > 0) ppsl.a = NULL else ppsl.a = rep(0, N.a);
  mu   = mu.m0
  phi  = phi.m0
  W    = W.b0 / W.a0;
  om   = rep(0, T);

  ## Check if known.
  know.phi <- know.mu <- know.W <- FALSE
  if (!is.null(phi.true))  { phi  = phi.true;  know.phi  = TRUE;
                             if (any(phi==1)) {  mu.true   = rep(0, N.b); } }
  if (!is.null(mu.true))   { mu   = mu.true ;  know.mu   = TRUE; }
  if (!is.null(W.true))    { W    = W.true  ;  know.W    = TRUE; }
  ## if (!is.null(iota.true)) { iota = iota.true; know.iota = TRUE; }

  start.time = proc.time()
  
  ## MCMC
  for(i in 1:(samp+burn)) {
    if (i==burn+1) {
      start.ess = proc.time();
      naccept = 0;
    }
    
    ## Draw beta
    W.mat = diag(W, N.b)
    ## draw = CUBS.R(y, X, n, mu, phi, W.mat, m0, C0);
    draw = CUBS.C(y, X, n, mu, phi, W.mat, m0, C0, obs="binom");
    ppsl.b  = draw$beta
    if (N.a > 0) ppsl.a = draw$alpha
    lf.ppsl = target.dens(y, X, n, ppsl.b, mu, phi, W.mat, m0, C0, log.sc=TRUE, alpha=ppsl.a, obs="binom")
    lq.ppsl = draw$log.dens;
    l.fdivq.ppsl = lf.ppsl - lq.ppsl
    l.ratio = l.fdivq.ppsl - l.fdivq
    ## l.fdivq.ppsl = 0; l.ratio = 0; lf.ppsl = 0; lq.ppsl = 0

    if (is.na(l.ratio) || is.nan(l.ratio)) {
      ## PROBLEM: YOU CAN HAVE SINGULAR MATRICES IN CUBS.C.
      ## THIS LEADS TO lq.ppsl=NaN.
      cat("l.ratio", l.ratio, "l.fdivq", l.fdivq, "lf.ppsl", lf.ppsl, "lq.ppsl", lq.ppsl, "\n");
      l.ratio = -Inf; ## Force rejection.
    }
    
    a.prob = min(1, exp(l.ratio))

    if (runif(1) < a.prob) {
      beta  = ppsl.b
      alpha = ppsl.a
      l.fdivq = l.fdivq.ppsl
      naccept = naccept + 1
    }
    ## End draw beta

    ## AR(1) - phi, W assumed to be diagonal !!!
    ## mu  = draw.mu.R(beta, phi, W, mu.m0, mu.V0) 
    ## phi = draw.phi.R(beta, mu, W, phi.m0, phi.V0, phi)
    ## W   = draw.W.R  (beta, mu, phi, W.a0, W.b0)
    if (!know.mu)  mu  = draw.mu.ar1.ind (beta, phi, W, mu.m0, mu.P0)
    if (!know.phi) phi = draw.phi.ar1.ind(beta, mu, W, phi.m0, phi.P0, phi)
    if (!know.W)   W   = draw.W.ar1.ind  (beta, mu, phi, W.a0, W.b0)
    
    if (i > burn) {
      ## out$ppsl.b[i-burn,,] = ppsl.b
      out$beta[i-burn,,]   = beta
      if (N.a > 0)
        out$alpha[i-burn,] = alpha
      out$mu[i-burn,]      = mu
      out$phi[i-burn,]     = phi
      out$W[i-burn,]       = W;
      out$l.fdivq[i-burn]  = l.fdivq
      out$l.ratio[i-burn]  = l.ratio
      out$lf.ppsl[i-burn]  = lf.ppsl
      out$lq.ppsl[i-burn]  = lq.ppsl
      out$ac.rate[i-burn]  = naccept / (i-burn)
    }

    if (i %% verbose == 0) {
      cat("DynLogitCUBS: ");
      if (i > burn) cat("Accept rate:", naccept / (i-burn), ", ");
      cat("Iteration:", i, "\n");
    }
    
  }

  end.time = proc.time()
  out$total.time = end.time - start.time
  out$ess.time   = end.time - start.ess
  out$last       = draw
  out$a.rate     = naccept / (i - burn);
  
  out
  
}

################################################################################
                                ## TEST BINOM ##
################################################################################

if (FALSE) {

  dyn.unload("BayesLogit.so")
  dyn.load("BayesLogit.so")
  source("LogitWrapper.R")
  
  ## source("DynLogitCUBS.R")
  source("DynLogitPG.R")
  
}

if (FALSE) {

  T = 300;
  P = 1;

  beta = array(0, dim=c(P, T+1));
  X = matrix(1, nrow=T, ncol=P);
  if (P != 1) X = matrix(rnorm(T*P), nrow=T);
  
  N = nrow(X);

  ## Parameters
  iota = rep(0, P);
  W   = rep(0.1, P);
  mu  = rep(0, P);
  phi = rep(0.95, P)

  ## Prior
  b.m0 = rep(0.0, P);
  b.C0 = diag(2.0, P);
  mu.m0 = mu
  mu.V0 = rep(1.0, P);
  phi.m0 = rep(0.9, P)
  phi.V0 = rep(0.1, P)
  W.a0   = rep(10, P);
  W.b0   = W.a0 * W;

  ## Synthetic
  beta[,1] = mu;
  for (i in 2:(T+1)) {
    beta[,i] = mu + phi * (beta[,i-1]-mu) + sqrt(W) * rnorm(P);
  }

  n = rep(1, T)
  psi = iota + apply(X * t(beta)[-1,], 1, sum); ## log-odds
  p   = exp(psi) / (1 + exp(psi));
  y = rbinom(T, n, p);
  w = rpg.devroye(T, n, psi);

  m0 = mu
  C0 = b.C0
  samp = 2000
  burn = 200
  verbose = 100

  ## source("DynLogitCUBS.R")
  one.R <- CUBS.R(y, X, n, mu, phi, diag(W, P), m0, C0);
  one.C <- CUBS.C(y, X, n, mu, phi, diag(W, P), m0, C0);

  ## source("DynLogitCUBS.R")
  out.cubs <- dyn.logit.CUBS(y, X.dyn=X, n, m0=m0, C0=C0,
                             samp=samp, burn=burn, verbose=verbose,
                             mu.m0=mu.m0, mu.P0=1/mu.V0,
                             phi.m0=phi.m0, phi.P0=1/phi.V0,
                             W.a0=W.a0, W.b0=W.b0,
                             mu.true=mu, phi.true=phi, W.true=NULL)

  out = list("cubs"=list("sstat"=list("beta"=list())));
  out$pg = out$CUBS;

  out$cubs$ess.time = out.cubs$ess.time[3];
  out$cubs$sstat$beta[[1]] = sum.stat.dyn(out.cubs$beta, out$cubs$ess.time);
  out$cubs$sstat$beta = simplify2array(out$cubs$sstat$beta);
  
}

if (FALSE) {

  source("DynLogitPG.R")
  ## samp = 10000; burn=2000; verbose = 1000;
  out.pg <- dyn.logit.PG(y, X, n, samp=samp, burn=burn, verbose=verbose,
                         m.0=m0, C.0=C0,
                         mu.m0=mu.m0, mu.P0=1/mu.V0,
                         phi.m0=phi.m0, phi.P0=1/phi.V0,
                         W.a0=W.a0, W.b0=W.b0,
                         beta.true=NULL, iota.true=NULL, w.true=NULL,
                         mu.true=mu, phi.true=phi, W.true=NULL)

  out$pg$ess.time = out.pg$ess.time[3];
  out$pg$sstat$beta[[1]] = sum.stat.dyn(out.pg$beta, out$pg$ess.time)
  out$pg$sstat$beta = simplify2array(out$pg$sstat$beta);
  
  out.table = setup.table.dyn(out, "beta")
}

if (FALSE) {

  ## par(mfrow=c(1,1))
  beta.mean = apply(plot.out$beta, c(2,3), mean);
  beta.95   = apply(plot.out$beta, c(2,3), function(x){quantile(x, 0.95)});
  beta.05   = apply(plot.out$beta, c(2,3), function(x){quantile(x, 0.05)});

  ymin = min(beta.mean, plot.out$beta);
  ymax = max(beta.mean, plot.out$beta);

  plot(0:T, beta.mean[1,], col=2, type="l", ylim=c(ymin,ymax));
  lines(0:T, beta.95[1,], col="pink")
  lines(0:T, beta.05[1,], col="pink")
  abline(h=mean(beta.mean[1,]), col=2, lty=c(2,2));

  lines(0:T, beta[1,]);
  
  if (n[1] > 3) { points(1:T, log(y / (n-y)), cex=0.1) } else
  { points(1:T, (y-min(y)) / (max(y)-min(y)) * (ymax-ymin) + ymin, cex=0.1) }
  
  lines(0:T, plot.out$beta[100,1,], col=3);
  
}

if (FALSE) {

  par(mfrow=c(1,3))
  hist(plot.out$mu[,1])
  hist(plot.out$phi[,1])
  hist(plot.out$W[,1])

  apply(plot.out$mu,  2, mean)
  apply(plot.out$phi, 2, mean)
  apply(plot.out$W,   2, mean)
  
}
