library("rootSolve")

## CUBS: conjugate updating backward sampling
## Use a MH step to to correct.

binom.solve <- function(rs, fq)
{
  r = rs[1]
  s = rs[2]
  E = digamma(r) - digamma(s)
  V = trigamma(r) + trigamma(s)
  F1 = E - fq[1]
  F2 = V - fq[2]
  out = c("F1"=F1, "F2"=F2)
  out
}

poisson.solve <- function(rs, fq)
{
  F1 = digamma(rs[1]) - log(abs(rs[2])) - fq[1]
  F2 = trigamma(rs[1]) - fq[2]
  out = c("F1"=F1, "F2"=F2)
  out
}

cubs.draw <- function(y, X, n, mu, phi, W, m0, C0)
{
  ## When tracking (beta_t, alpha).  It may be the case that there is no alpha.
  ## z_t = x_t (beta_t, alpha_t) + ep_t, ep_t \sim N(0, V_t).
  ## beta_t = mu + phi * (beta_t - mu) + omega_t, omega_t \sim N(0,W).
  ## alpha_t = alpha_{t-1}
  
  ## y : vector of observations (T)
  ## X : design matrix (T x N)
  ## mu : mu (K)
  ## phi : vector (K)
  ## W : covariance MATRIX of innovations of beta (K x K)
  ## m0 : prior mean on (beta_0, alpha_0) (N)
  ## C0 : prior var on (beta_0, alpha_0) (N x N).

  W = as.matrix(W)
  
  T = length(y);
  N.b = ncol(W);
  N   = ncol(X);
  N.a = N - N.b;
  b.idc = 1:N.b;
  a.idc = 1:N.a+N.b

  with.alpha = N.a > 0;
  
  m = array(m0, dim=c(N, T+1));
  C = array(C0, dim=c(N, N, T+1));
  R = array(0., dim=c(N, N, T+1));
  a = array(0., dim=c(N, T+1));
  rs = array(0., dim=c(2, T));

  beta = array(0, dim=c(N.b, T+1));
  
  d = c( phi, rep(1, N.a) );
  D = diag(d, N);
  mu = c( mu, rep(0, N.a) );
  big.W = matrix(0, N, N); big.W[b.idc, b.idc] = W;
  if (length(n)==1) n = array(n, dim=T);
  
  ## Filter Forward
  for (i in 2:(T+1)) {
    i.l = i-1;

    a[,i]  = d * m[,i-1] + (1-d) * mu;
    R[,,i] = D %*% C[,,i-1] %*% D  + big.W;

    x.i = t(X[i.l,])
    f.i = x.i %*% a[,i];
    q.i = ( x.i %*% R[,,i] %*% X[i.l,] )[1];
    
    rho.i = R[,,i] %*% X[i.l,];
    A.i = rho.i / q.i;

    ## Conjugate update
    ## ## Binomial
    rs.i = multiroot(binom.solve, start=c(1,1), fq=c(f.i, q.i));
    rs[,i.l] = rs.i$root
    rstar.i = rs.i$root[1] + y[i.l];
    sstar.i = n[i.l] - y[i.l] + rs.i$root[2];
    fqstar.i = binom.solve(c(rstar.i, sstar.i), c(0, 0));
    ## ## cat("f,q:", f.i, q.i, "root:", rs.i$root, "f.root:", rs.i$f.root, "\n");
    ## ## Gaussian
    ## qstar.i = (1 / (1 / q.i + 1 / n[i.l]))
    ## fstar.i = (f.i / q.i + y[i.l] / n[i.l]) * qstar.i
    ## fqstar.i = c(fstar.i, qstar.i)
    ## ## cat(c(f.i, q.i), fqstar.i, "\n")
    ## ## Poisson
    ## r.i = multiroot(function(r,q){trigamma(r)-q}, start=1, q=q.i);
    ## s.i = exp(digamma(r.i$root[1]) - f.i)
    ## rs.i = list(root=c(r.i$root[1], s.i));
    ## rs[,i.l] = rs.i$root
    ## fstar.i = digamma(rs.i$root[1] + y[i.l]) - log(rs.i$root[2] + 1)
    ## qstar.i = trigamma(rs.i$root[1] + y[i.l])
    ## fqstar.i = c(fstar.i, qstar.i)
    ## Neg. Binomial
    ## fhat.i = f.i - log(n[i.l])
    ## rs.i = multiroot(binom.solve, start=c(1,1), fq=c(fhat.i, q.i));
    ## rs[,i.l] = rs.i$root
    ## rstar.i = rs.i$root[1] + y[i.l];
    ## sstar.i = n[i.l] + rs.i$root[2];
    ## fqstar.i = binom.solve(c(rstar.i, sstar.i), c(0, 0));
    ## fqstar.i[1] = fqstar.i[1] + log(n[i.l])
    
    m[,i]  = a[,i] + A.i * (fqstar.i[1] - f.i);
    C[,,i] = R[,,i] + rho.i %*% t(rho.i) * ( (fqstar.i[2] / q.i - 1) / q.i );

  }

  ## Keep track of log density.
  log.dens = 0
  
  ## Backwards sample
  L = t( chol(C[,,T+1]) );
  ## evd = eigen(C[,,T+1]);
  ## Rt = evd$vectors %*% diag(sqrt(evd$values), N) %*% t(evd$vectors);
  ep    = rnorm(N);
  theta = m[,T+1] + L %*% ep;
  ## alpha = ifelse(with.alpha, theta[a.idc], 0);
  if (with.alpha) alpha = theta[a.idc] else alpha = 0
  beta[,T+1] = theta[b.idc];

  log.dens = log.dens - 0.5 * (t(ep) %*% ep) - sum(log(diag(L)));

  for (i in (T+1):2) {

    A.bs = C[b.idc, b.idc, i-1] %*% (solve(R[b.idc, b.idc, i]) * phi);
    V.bs = C[b.idc, b.idc, i-1] - A.bs %*% R[b.idc, b.idc, i] %*% t(A.bs);
    m.bs = m[b.idc, i-1] + A.bs %*% (beta[,i] - a[b.idc,i]);

    L  = t(chol(V.bs));
    ep = rnorm(N.b)
    beta[,i-1] = m.bs + L %*% ep;

    log.dens = log.dens - 0.5 * (t(ep) %*% ep) - sum(log(diag(L)));
    
  }

  ## Need to return log density as well or r,s to calculate log-density.
  out = list("alpha"=alpha, "beta"=beta, "log.dens"=log.dens, "m"=m, "C"=C, "rs"=rs);
  out
}

##------------------------------------------------------------------------------

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

  ## llh
  for (i in 2:(T+1)) {
    mean.i = mu + phi * (beta[,i-1] - mu);
    log.dens = log.dens + gaussian.dens(beta[,i], mean.i, L, log.sc=TRUE, is.V=FALSE);
    ## e = (beta[,i] - mu) - phi * (beta[,i-1] - mu);
    ## e = backsolve(L, e, upper.tri=FALSE)
    ## log.dens = log.dens - 0.5 * (t(e) %*% e);
  }
  ## log.dens = log.dens - T * sum(log(diag(L)));
  
  ## prior
  L = t(chol(C0))
  ## e = c(beta[,1], alpha) - m0;
  ## e = backsolve(L, e, upper.tri=FALSE)
  ## log.dens = log.dens - 0.5 * (t(e) %*% e) - sum(log(diag(L)));
  log.dens = log.dens + gaussian.dens(c(beta[,1], alpha), m0, L, log.sc=TRUE, is.V=FALSE);
  
  out = ifelse(log.sc, log.dens, exp(log.dens));
  out
}

target.dens <- function(y, X, n, beta, mu, phi, W, m0, C0, log.sc=FALSE, alpha=NULL)
{
  ## The first rows of X correspond to beta.
 
  T = length(y)
  N.b = nrow(as.matrix(W))
  N = ncol(X)
  N.a = N - N.b
  b.idc = 1:N.b
  a.idc = N.b + 1:N.a
  
  X.dyn = matrix(X[,b.idc], ncol=N.b);
  if (N.a > 0) matrix(X.stc[,a.idc], ncol=N.a);
  
  psi = apply(X.dyn * t(beta)[-1,], 1, sum);
  if (N.a > 0) psi = psi + X.stc %*% alpha;
  
  ar1.llh = ar1.dens(beta, mu, phi, W, m0, C0, log.sc=TRUE, alpha);
  obs.llh = sum(binom.obs.dens(y, n, psi, log.sc=TRUE))
  ## obs.llh = sum(gauss.obs.dens(y, n, psi, log.sc=TRUE))
  ## obs.llh = sum(poisson.obs.dens(y, psi, log.sc=TRUE))
  ## obs.llh = sum(neg.binom.obs.dens(y, n, psi, log.sc=TRUE))

  llh = ar1.llh + obs.llh;
  out = ifelse(log.sc, llh, exp(llh))
  out
}

##------------------------------------------------------------------------------

dyn.logit.cubs <- function(y, X.dyn, n, m0, C0,
                           samp=1000, verbose=100,
                           mu.m0=NULL, mu.V0=NULL,
                           phi.m0=NULL, phi.V0=NULL,
                           W.a0=NULL, W.b0=NULL, X.stc=NULL,
                           mu.true = NULL, phi.true=NULL, W.true=NULL)
{

  X = cbind(X.dyn, X.stc)
  T = length(y)
  N.b = ncol(X.dyn)
  N = ncol(X)
  N.a = N - N.b
  M = samp

  ## Output
  out <- list("beta"=array(0, dim=c(M, N.b, T+1)),
              "ppsl.b"=array(0, dim=c(M, N.b, T+1)),
              "l.fdivq"=rep(0, M),
              "l.ratio"=rep(0, M),
              "lf.ppsl"=rep(0, M),
              "lq.ppsl"=rep(0, M),
              "a.rate"=rep(0, M))
  if (N.a > 0) out$alpha=array(0, dim=c(M, N.a))

  ## Initialize
  beta = array(0, dim=c(P, T+1));
  l.fdivq = -Inf;
  naccept = 0
  ppsl.a = NULL
  
  ## Check if known.
  if (!is.null(phi.true))  { phi  = phi.true;  know.phi  = TRUE;
                             if (phi==1) {mu.true = rep(0, N.b);}}
  if (!is.null(mu.true))   { mu   = mu.true;   know.mu   = TRUE; }
  if (!is.null(W.true))    { W    = W.true;    know.W    = TRUE; }
  ## if (!is.null(iota.true)) { iota = iota.true; know.iota = TRUE; }

  ## MCMC
  for(i in 1:M) {

    ## Draw beta
    draw = cubs.draw(y, X, n, mu, phi, W, m0, C0);
    ppsl.b  = draw$beta
    if (N.a > 0) ppsl.a = draw$alpha
    lf.ppsl = target.dens(y, X, n, ppsl.b, mu=mu, phi=phi, W=W, m0=m0, C0=C0, log.sc=TRUE, alpha=ppsl.a)
    lq.ppsl = draw$log.dens;
    l.fdivq.ppsl = lf.ppsl - lq.ppsl
    l.ratio = l.fdivq.ppsl - l.fdivq
    
    a.prob = min(1, exp(l.ratio))

    if (runif(1) < a.prob) {
      beta  = ppsl.b
      alpha = ppsl.a
      l.fdivq = l.fdivq.ppsl
      naccept = naccept + 1
    }
    ## End draw beta

    out$ppsl.b[i,,] = ppsl.b
    out$beta[i,,] = beta
    if (N.a > 0) out$alpha[i,] = alpha
    out$l.fdivq[i] = l.fdivq
    out$l.ratio[i] = l.ratio
    out$lf.ppsl[i] = lf.ppsl
    out$lq.ppsl[i] = lq.ppsl
    out$a.rate[i]  = naccept / i
    out$last       = draw

    if (i %% verbose == 0) cat("CUBS: iteration", i, "a.rate:", naccept / i, "\n");
    
  }

  out
  
}

                 

################################################################################
                                ## TEST BINOM ##
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
  mu  = 0;
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

  n = rep(50, T)
  psi = iota + apply(X * t(beta)[-1,], 1, sum); ## log-odds
  p   = exp(psi) / (1 + exp(psi));
  y = rbinom(T, n, p);
  w = rpg.devroye(T, n, psi);

  m0 = mu
  C0 = diag(b.C0, P);
  M = 1000
  verbose = 100

  ## source("DynLogitCUBS.R")
  one.draw <- cubs.draw(y, X, n, mu, phi, W, m0, C0);

  ## source("DynLogitCUBS.R")
  out.cubs <- dyn.logit.cubs(y, X.dyn=X, n, m0, C0,
                             samp=M, verbose=verbose,
                             mu.true=mu, phi.true=phi, W.true=W)
  
}

if (FALSE) {

  source("DynLogitPG.R")
  samp = 1000
  burn = 0
  out.pg <- dyn.logit.PG(y, X, n, samp=samp, burn=burn, verbose=100,
                         m.0=m0, C.0=C0,
                         mu.m0=NULL, mu.V0=NULL,
                         phi.m0=NULL, phi.V0=NULL,
                         W.a0=NULL, W.b0=NULL,
                         beta.true=NULL, iota.true=NULL, w.true=NULL,
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
  
  if (n[1] > 3) { points(1:T, log(y / (n-y)), cex=0.1) } else
  { points(1:T, (y-min(y)) / (max(y)-min(y)) * (ymax-ymin) + ymin, cex=0.1) }
  
  lines(0:T, out$beta[100,1,], col=3);
  
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
  M = 100
  verbose = 100
  
  ## source("DynLogitCUBS.R")
  one.draw <- cubs.draw(y, X, sig2, mu, phi, W, m0, C0);
  
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
  
  ## source("DynLogitCUBS.R")
  one.draw <- cubs.draw(y, X, rep(0, length(y)), mu, phi, W, m0, C0);

  ## source("DynLogitCUBS.R")
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

  ## source("DynLogitCUBS.R")
  one.draw <- cubs.draw(y, X, n, mu, phi, W, m0, C0);

  ## source("DynLogitCUBS.R")
  out.cubs <- dyn.logit.cubs(y, X.dyn=X, n, m0, C0,
                             samp=M, verbose=verbose,
                             mu.true=mu, phi.true=phi, W.true=W)
  
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
