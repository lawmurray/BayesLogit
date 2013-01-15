gaussian.dens <- function(x, m, LorV, log.sc=FALSE, is.V=FALSE)
{
  if (is.V) LorV = t(chol(LorV));
  ## L = chol(Var)
  log.dens = 0.0
  e = (x - m);
  e = backsolve(LorV, e, upper.tri=FALSE);
  ## e = forwardsolve(LorV, e)
  ## cat(e, "\n")
  log.dens = log.dens - 0.5 * (t(e) %*% e) - sum(log(diag(LorV)));
  out = ifelse(log.sc, log.dens, exp(log.dens));
  out
}

binom.obs.dens <- function(y, n, psi, log.sc=TRUE)
{
  ## log.dens =  y * psi - n * log(1+exp(psi));
  ## out = ifelse(log.sc, log.dens, exp(log.dens)) # WRONG USE OF ifelse
  p = 1 / (1 + exp(-1 * psi))  
  dbinom(y, n, p, log=log.sc);
}

gauss.obs.dens <- function(y, sig2, psi, log.sc=TRUE)
{
  dnorm(y, psi, sqrt(sig2), log=log.sc)
}

poisson.obs.dens <- function(y, psi, log.sc=TRUE)
{
  lambda = exp(psi)
  dpois(y, lambda, log=log.sc)
}

neg.binom.obs.dens <- function(y, n, zeta, log.sc=TRUE)
{
  ## n = shape
  psi = zeta - log(n)
  p = 1 / (1 + exp(-1 * psi))
  dnbinom(y, n, 1-p, log=TRUE)
}

ar1.dens <- function(beta, mu, phi, W, m0, C0, log.sc=FALSE, alpha=NULL)
{
  T = ncol(beta) - 1;
  L = t(chol(W))
  log.dens = 0;
  C1 = 

  ## llh
  for (i in 2:(T+1)) {
    mean.i = mu + phi * (beta[,i-1] - mu);
    ## log.dens = log.dens + gaussian.dens(beta[,i], mean.i, L, log.sc=TRUE, is.V=FALSE);
    e = beta[,i] - mean.i
    e = backsolve(L, e, upper.tri=FALSE)
    log.dens = log.dens - 0.5 * (t(e) %*% e) - sum(log(diag(L)));
  }
  ## ## log.dens = log.dens - T * sum(log(diag(L)));

  ## prior
  L = t(chol(C0))
  ## log.dens = log.dens + gaussian.dens(c(beta[,1], alpha), m0, L, log.sc=TRUE, is.V=FALSE);
  e = c(beta[,1], alpha) - m0;
  e = backsolve(L, e, upper.tri=FALSE)
  log.dens = log.dens - 0.5 * (t(e) %*% e) - sum(log(diag(L)));
  
  out = ifelse(log.sc, log.dens, exp(log.dens));
  out
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
  ## Using a switch statement slowed things down by 3.5%.
  obs.llh <- switch(obs,
                    binom  = sum(binom.obs.dens(y, n, psi, log.sc=TRUE)),
                    nbinom = sum(nbinom.obs.dens(y, n, psi, log.sc=TRUE)),
                    norm   = sum(gauss.obs.dens(y, n, psi, log.sc=TRUE)))

  ## obs.llh = do.call("binom.obs.dens", list(y, n, psi, log.sc=TRUE))

  llh = ar1.llh + obs.llh;
  out = ifelse(log.sc, llh, exp(llh))
  out
}

##------------------------------------------------------------------------------

cubs.mh <- function(y, X.dyn, n, m0, C0,
                    samp=1000, burn=100, verbose=100,
                    mu.m0=NULL, mu.P0=NULL,
                    phi.m0=NULL, phi.P0=NULL,
                    W.a0=NULL, W.b0=NULL, X.stc=NULL,
                    mu.true = NULL, phi.true=NULL, W.true=NULL,
                    obs=c("binom", "nbinom", "norm"))
{
  ## Dim and restructure.
  X   = cbind(X.stc, X.dyn)
  C0  = as.matrix(C0);
  T   = length(y)
  N.b = ncol(X.dyn)
  N   = ncol(X)
  N.a = N - N.b
  M   = samp
  obs = obs[1];

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
              "a.rate"  = rep(0, M))
  if (N.a > 0) out$alpha=array(0, dim=c(M, N.a))

  ## Initialize
  beta    = array(0, dim=c(N.b, T+1));
  l.fdivq = -Inf;
  naccept = 0
  if (N.a > 0) ppsl.a = NULL else ppsl.a = rep(0, N.a);

  ## Check if known.
  know.phi <- know.mu <- know.W <- FALSE
  if (!is.null(phi.true))  { phi  = phi.true;  know.phi  = TRUE;
                             if (phi[1]==1) {  mu.true   = rep(0, N.b); } }
  if (!is.null(mu.true))   { mu   = mu.true ;  know.mu   = TRUE; }
  if (!is.null(W.true))    { W    = W.true  ;  know.W    = TRUE; }
  ## if (!is.null(iota.true)) { iota = iota.true; know.iota = TRUE; }

  start.time = proc.time()
  
  ## MCMC
  for(i in 1:(samp+burn)) {
    if (i==burn+1) start.ess = proc.time();
    
    ## Draw beta
    W.mat = diag(W, N.b)
    ## draw = CUBS.R(y, X, n, mu, phi, W.mat, m0, C0);
    draw = CUBS.C(y, X, n, mu, phi, W.mat, m0, C0, obs=obs);
    ppsl.b  = draw$beta
    if (N.a > 0) ppsl.a = draw$alpha
    lf.ppsl = target.dens(y, X, n, ppsl.b, mu, phi, W.mat, m0, C0, log.sc=TRUE, alpha=ppsl.a, obs=obs)
    lq.ppsl = draw$log.dens;
    l.fdivq.ppsl = lf.ppsl - lq.ppsl
    l.ratio = l.fdivq.ppsl - l.fdivq
    ## l.fdivq.ppsl = 0; l.ratio = 0; lf.ppsl = 0; lq.ppsl = 0
    
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
    if (!know.phi) phi = draw.phi.ar1.ind(beta, phi, W, phi.m0, phi.P0, phi)
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
      out$a.rate[i-burn]   = naccept / i
    }

    if (i %% verbose == 0) cat("CUBS: iteration", i, "a.rate:", naccept / i, "\n");
    
  }

  end.time = proc.time()
  out$total.time = end.time - start.time
  out$ess.time   = end.time - start.ess
  out$last          = draw
  
  out
  
}

                 

################################################################################
                                ## TEST BINOM ##
################################################################################

if (FALSE) {

  dyn.unload("BayesLogit.so")
  dyn.load("BayesLogit.so")
  source("LogitWrapper.R")
  
  ## source("CUBS-MH.R")
  source("FFBS.R")
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

  ## source("CUBS-MH.R")
  one.R <- CUBS.R(y, X, n, mu, phi, diag(W, P), m0, C0);
  one.C <- CUBS.C(y, X, n, mu, phi, diag(W, P), m0, C0);

  ## source("CUBS-MH.R")
  out.cubs <- cubs.mh(y, X.dyn=X, n, m0=m0, C0=C0,
                      samp=samp, burn=burn, verbose=verbose,
                      mu.m0=mu.m0, mu.P0=1/mu.V0,
                      phi.m0=phi.m0, phi.P0=1/phi.V0,
                      W.a0=W.a0, W.b0=W.b0,
                      mu.true=mu, phi.true=phi, W.true=NULL,
                      obs="binom")

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

################################################################################
                                ## TEST NORM ##
################################################################################

if (FALSE) {

  T = 100;
  P = 1;

  beta = array(0, dim=c(P, T+1));
  X = matrix(1, nrow=T, ncol=P);

  N = nrow(X);

  ## Parameters
  iota = 0;
  W   = 0.1;
  mu  = 2.0;
  phi = 0.95

  ## Prior
  b.m0 = 0.0;
  b.C0 = 2.0;
  phi.m0 = 0.9
  phi.V0 = 0.1;
  W.a0   = 10;
  W.b0   = W.a0 * W;

  ## Synthetic
  beta[,1] = mu;
  for (i in 2:(T+1)) {
    beta[,i] = mu + phi * (beta[,i-1]-mu) + sqrt(W) * rnorm(P);
  }

  sig2 = rep(0.1^2, T);
  psi = iota + apply(X * t(beta)[-1,], 1, sum);
  y = rnorm(T, psi, sqrt(sig2));

  m0 = mu
  C0 = diag(b.C0, P);
  M = 1000
  burn = 100
  verbose = 100
  
  ## source("CUBS-MH.R")
  one.R <- CUBS.R(y, X, sig2, mu, phi, W, m0, C0);
  one.C <- CUBS.C(y, X, sig2, mu, phi, diag(W, P), m0, C0, obs="norm");

  ## source("CUBS-MH.R")
  out.cubs <- cubs.mh(y, X.dyn=X, sig2, m0, C0,
                      samp=M, burn=burn, verbose=verbose,
                      mu.true=mu, phi.true=phi, W.true=W,
                      obs="norm")
  
}

if (FALSE) {

  beta.mean = apply(out$beta, c(2,3), mean);
  beta.95   = apply(out$beta, c(2,3), function(x){quantile(x, 0.95)});
  beta.05   = apply(out$beta, c(2,3), function(x){quantile(x, 0.05)});

  ymin = min(beta.mean, out$beta);
  ymax = max(beta.mean, out$beta);

  plot(0:T, beta.mean[1,], col=2, type="l", ylim=c(ymin,ymax));
  lines(0:T, beta.95[1,], col="pink")
  lines(0:T, beta.05[1,], col="pink")
  abline(h=mean(beta.mean[1,]), col=2, lty=c(2,2));

  lines(0:T, beta[1,]);
  
  points(1:T, y, cex=0.1);
  
  lines(0:T, out$beta[100,1,], col=3);
  
}

################################################################################
                               ## TEST POISSON ##
################################################################################

if (FALSE) {

  T = 100;
  P = 1;

  beta = array(0, dim=c(P, T+1));
  X = matrix(1, nrow=T, ncol=P);

  N = nrow(X);

  ## Parameters
  iota = 0;
  W   = 0.1;
  mu  = 10;
  phi = 0.95

  ## Prior
  b.m0 = mu;
  b.C0 = 2.0;
  phi.m0 = 0.9
  phi.V0 = 0.1;
  W.a0   = 10;
  W.b0   = W.a0 * W;

  ## Synthetic
  beta[,1] = mu;
  for (i in 2:(T+1)) {
    beta[,i] = mu + phi * (beta[,i-1]-mu) + sqrt(W) * rnorm(P);
  }

  psi = iota + apply(X * t(beta)[-1,], 1, sum);
  y = rpois(T, exp(psi));

  m0 = mu
  C0 = diag(b.C0, P);
  M = 1000
  verbose = 100
  
  ## source("CUBS-MH.R")
  one.draw <- CUBS.R(y, X, rep(0, length(y)), mu, phi, W, m0, C0);

  ## source("CUBS-MH.R")
  out.cubs <- dyn.logit.cubs(y, X.dyn=X, sig2, m0, C0,
                             samp=M, verbose=verbose,
                             mu.true=mu, phi.true=phi, W.true=W)
  
}

if (FALSE) {

  beta.mean = apply(out$beta, c(2,3), mean);
  beta.95   = apply(out$beta, c(2,3), function(x){quantile(x, 0.95)});
  beta.05   = apply(out$beta, c(2,3), function(x){quantile(x, 0.05)});

  ymin = min(beta.mean, out$beta);
  ymax = max(beta.mean, out$beta);

  plot(0:T, beta.mean[1,], col=2, type="l", ylim=c(ymin,ymax));
  lines(0:T, beta.95[1,], col="pink")
  lines(0:T, beta.05[1,], col="pink")
  abline(h=mean(beta.mean[1,]), col=2, lty=c(2,2));

  lines(0:T, beta[1,]);
  
  points(1:T, log(y+1e-6), cex=0.1);
  
  lines(0:T, out$beta[100,1,], col=3);
  
}

################################################################################
                              ## TEST NEG BINOM ##
################################################################################

if (FALSE) {

  T = 100;
  P = 1;

  beta = array(0, dim=c(P, T+1));
  X = matrix(1, nrow=T, ncol=P);

  N = nrow(X);

  ## Parameters
  iota = 0;
  W   = 0.1;
  mu  = 0.5;
  phi = 0.95
  d = 4

  ## Prior
  b.m0 = 0.0;
  b.C0 = 2.0;
  phi.m0 = 0.9
  phi.V0 = 0.1;
  W.a0   = 10;
  W.b0   = W.a0 * W;

  ## Synthetic
  beta[,1] = mu;
  for (i in 2:(T+1)) {
    beta[,i] = mu + phi * (beta[,i-1]-mu) + sqrt(W) * rnorm(P);
  }

  n = rep(d, T)
  zeta = iota + apply(X * t(beta)[-1,], 1, sum); ## log-mean
  psi = zeta - log(d)
  p   = exp(psi) / (1 + exp(psi));
  y = rnbinom(T, n, 1-p);
  w = rpg.devroye(T, y+n, psi);

  m0 = mu
  C0 = diag(b.C0, P);
  M = 1000
  verbose = 100

  ## source("CUBS-MH.R")
  one.draw <- CUBS.R(y, X, n, mu, phi, W, m0, C0);

  ## source("CUBS-MH.R")
  out.cubs <- cubs.mh(y, X.dyn=X, n, m0, C0,
                      samp=M, verbose=verbose,
                      mu.true=mu, phi.true=phi, W.true=W,
                      obs="nbinom")
  
}

if (FALSE) {

  source("DynNBPG.R")
  samp = 1000
  burn = 0
  out.pg <- dyn.NB.PG(y, X.dyn=X, X.stc=NULL,
                      samp=samp, burn=burn, verbose=50,
                      m.0=m0, C.0=C0,
                      mu.m0=NULL, mu.V0=NULL,
                      phi.m0=NULL, phi.V0=NULL,
                      W.a0=NULL, W.b0=NULL,
                      d.true=n[1], w.true=NULL,
                      beta.true=NULL, iota.true=NULL,
                      mu.true=mu, phi.true=phi, W.true=W)
  
}

if (FALSE) {

  beta.mean = apply(out$beta, c(2,3), mean);
  beta.95   = apply(out$beta, c(2,3), function(x){quantile(x, 0.95)});
  beta.05   = apply(out$beta, c(2,3), function(x){quantile(x, 0.05)});

  ymin = min(beta.mean, out$beta);
  ymax = max(beta.mean, out$beta);

  plot(0:T, beta.mean[1,], col=2, type="l", ylim=c(ymin,ymax));
  lines(0:T, beta.95[1,], col="pink")
  lines(0:T, beta.05[1,], col="pink")
  abline(h=mean(beta.mean[1,]), col=2, lty=c(2,2));

  lines(0:T, beta[1,]);

  ly = log(y+1e-1)
  ## points(1:T, (ly-min(ly))/(max(ly)-min(ly))*(ymax-ymin) + ymin, cex=0.1);
  points(1:T, ly, cex=0.1);
  
  lines(0:T, out$beta[100,1,], col=3);
  
}

################################################################################
                                ## UNIT TEST ##
################################################################################

## AR1 DENS
if (FALSE) {

  T = 1000;
  P = 1;

  beta = array(0, dim=c(P, T+1));

  N = nrow(X);

  ## Parameters
  W   = 0.1
  mu  = 2.0
  phi = 0.9

  ## Prior
  b.m0 = mu;
  b.C0 = 2.0;

  ## Synthetic
  beta[,1] = b.m0;
  for (i in 2:(T+1)) {
    beta[,i] = mu + phi * (beta[,i-1]-mu) + sqrt(W) * rnorm(P);
  }  

  ar1.llh = function(theta, beta, m0, C0) {
    N = length(m0);
    mu = theta[1:N];
    phi = theta[1:N + N];
    W   = matrix(theta[1:(N*N) + 2 * N], ncol=N, nrow=N);
    ar1.dens(beta, mu, phi, W, m0, C0, log.sc=TRUE);
  }

  theta.0 = c(1.0, 0.5, 0.2)
  optim.out = optim(theta.0, ar1.llh, method="L-BFGS-B", hessian=TRUE,
    lower=c(-Inf, 0, 0.0001), upper=c(Inf, 1, Inf),
    beta=beta, m0=b.m0, C0=b.C0, control=list(fnscale=-1));
  
}

if (FALSE) {

  M = 10000
  r = 2
  s = 1

  samp.true = rbeta(1000, r, s);
  samp.true = log(samp.true) - log(1-samp.true)

  fq = binom.solve(c(r,s), c(0,0))

  m.ppsl = fq[1]
  v.ppsl = fq[2]

  out = list(beta=rep(0, M), l.fdviq=rep(0, M));

  beta = 0
  naccept = 0
  l.fdivq = -Inf
  
  for (i in 1:M) {
    ppsl = rnorm(1, m.ppsl, sqrt(v.ppsl));
    l.fdivq.ppsl = r * ppsl - (r+s) * log(1+exp(ppsl)) - dnorm(ppsl, m.ppsl, sqrt(v.ppsl), log=TRUE);
    l.ratio = l.fdivq.ppsl - l.fdivq
    a.prob = min(1, exp(l.ratio))
    if (runif(1) < a.prob) {
      beta = ppsl
      l.fdivq = l.fdivq.ppsl
      naccept = naccept + 1
    }
    out$beta[i] = beta
    out$l.fdivq[i] = l.fdivq
  }

  cat("a.rate:", naccept / M, "\n")

  hist(samp.true)
  hist(out$beta)

  mean(samp.true); var(samp.true)
  mean(out$beta); var(out$beta)
  
}
