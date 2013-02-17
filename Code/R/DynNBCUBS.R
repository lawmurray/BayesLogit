nbinom.obs.dens <- function(y, n, zeta, log.sc=TRUE)
{
  ## n = shape
  psi = zeta - log(n)
  p = 1 / (1 + exp(-1 * psi))
  dnbinom(y, n, 1-p, log=TRUE)
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

dyn.NB.CUBS <- function(y, X.dyn, m0, C0,
                        samp=1000, burn=100, verbose=100,
                        mu.m0=NULL, mu.P0=NULL,
                        phi.m0=NULL, phi.P0=NULL,
                        W.a0=NULL, W.b0=NULL, X.stc=NULL,
                        mu.true = NULL, phi.true=NULL, W.true=NULL, d.true=NULL)
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
              "d"    = array(0, dim=c(M)),
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
  d       = 1
  mu   = mu.m0;
  phi  = phi.m0;
  W    = W.b0 / W.a0;
  if (N.a > 0) ppsl.a = rep(0, N.a) else ppsl.a = NULL ;
  alpha = ppsl.a

  ## Check if known.
  know.phi <- know.mu <- know.W <- know.d <- FALSE
  if (!is.null(d.true))    { d    = d.true;    know.d    = TRUE; }
  if (!is.null(phi.true))  { phi  = phi.true;  know.phi  = TRUE;
                             if (phi[1]==1) {  mu.true   = rep(0, N.b); } }
  if (!is.null(mu.true))   { mu   = mu.true ;  know.mu   = TRUE; }
  if (!is.null(W.true))    { W    = W.true  ;  know.W    = TRUE; }
  ## if (!is.null(iota.true)) { iota = iota.true; know.iota = TRUE; }

  ## Preprocess ## 
  ymax = max(y);
  F = cumsum(hist(y, breaks=0:(ymax+1)-0.5, plot=FALSE)$counts)
  G = T - F;
  
  start.time = proc.time()
  
  ## MCMC
  for(i in 1:(samp+burn)) {
    if (i==burn+1) start.ess = proc.time();

    ## draw (d | beta)
    if (!know.d) {
      log.mean = apply(X.dyn * t(beta)[-1,], 1, sum);
      if (N.a > 0) log.mean = log.mean + X.stc %*% alpha;
      mu.lambda  = exp(log.mean)
      d = draw.df(y, d, mu.lambda, G, ymax);
    }
    n = rep(d, T);
    
    ## Draw beta
    W.mat = diag(W, N.b)
    ## draw = CUBS.R(y, X, n, mu, phi, W.mat, m0, C0, obs="nbinom");
    draw = CUBS.C(y, X, n, mu, phi, W.mat, m0, C0, obs="nbinom");
    ppsl.b  = draw$beta
    if (N.a > 0) ppsl.a = draw$alpha
    ## cat("N:", N, "N.a", N.a, "N.b", N.b, "\n");
    ## cat("dim beta:", dim(draw$beta), "\n")
    ## plot(draw$beta[1,])
    ## readline("")
    lf.ppsl = target.dens(y, X, n, ppsl.b, mu, phi, W.mat, m0, C0, log.sc=TRUE, alpha=ppsl.a, obs="nbinom")
    ## cat("lf.ppsl:", lf.ppsl, "\n")
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
      out$d[i-burn]        = d;
      out$l.fdivq[i-burn]  = l.fdivq
      out$l.ratio[i-burn]  = l.ratio
      out$lf.ppsl[i-burn]  = lf.ppsl
      out$lq.ppsl[i-burn]  = lq.ppsl
      out$a.rate[i-burn]   = naccept / i
    }

    if (i %% verbose == 0) cat("CUBS: iteration", i, "a.rate:", naccept / i, "\n");
    
  }

  if(!is.null(draw$rs)) out$rs = draw$rs

  end.time = proc.time()
  out$total.time = end.time - start.time
  out$ess.time   = end.time - start.ess
  out$last          = draw
  
  out
  
}

################################################################################
                              ## TEST NEG BINOM ##
################################################################################

if (FALSE) {

  T = 1000;
  P = 3;

  beta = array(0, dim=c(P, T+1));
  X = matrix(rnorm(T*P), nrow=T, ncol=P);

  N = nrow(X);

  ## Parameters
  iota = 0
  W   = rep(0.1, P);
  mu  = rep(0.5, P);
  phi = rep(0.95, P)
  d = 4

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

  n = rep(d, T)
  zeta = iota + apply(X * t(beta)[-1,], 1, sum); ## log-mean
  psi = zeta - log(d)
  p   = exp(psi) / (1 + exp(psi));
  y = rnbinom(T, n, 1-p);
  w = rpg.devroye(T, y+n, psi);

  m0 = mu
  C0 = b.C0
  scale.up = 1
  samp = 1000 * scale.up
  burn = 200 * scale.up
  verbose = 100 * scale.up

  ## source("DynNBCUBS.R")
  one.draw <- CUBS.R(y, X, n, mu, phi, diag(W, P), m0, C0, obs="nbinom");

  ## source("DynNBCUBS.R")
  out.cubs <- dyn.NB.CUBS(y, X.dyn=X, m0, C0,
                          samp=samp, verbose=verbose,
                          mu.true=NULL, phi.true=NULL, W.true=NULL,
                          d.true=NULL)
  
}

if (FALSE) {

  source("DynNBPG.R")
  ## samp = 1000; burn = 0
  out.pg <- dyn.NB.PG(y, X.dyn=X, X.stc=NULL,
                      samp=samp, burn=burn, verbose=verbose,
                      m.0=m0, C.0=C0,
                      mu.m0=NULL, mu.P0=NULL,
                      phi.m0=NULL, phi.P0=NULL,
                      W.a0=NULL, W.b0=NULL,
                      d.true=NULL, w.true=NULL,
                      beta.true=NULL, iota.true=NULL,
                      mu.true=NULL, phi.true=NULL, W.true=NULL)

  beta.cubs = apply(out.cubs$beta, c(2,3), mean)
  beta.pg   = apply(out.pg$beta, c(2,3), mean)

  plot(beta.cubs[1,], type="l");
  lines(beta.pg[1,], col=2)
  
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
