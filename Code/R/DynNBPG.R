## Negative binomial regression using PG augmentation.
## We model the log-mean here.

source("NB-Shape.R")
source("AR1.R"); ## Independent AR(1)'s.  Maybe should change this.

################################################################################
                             ## Dynamc NB by PG ##
################################################################################

dyn.NB.PG <- function(y, X.dyn, X.stc=NULL,
                      samp=1000, burn=100, verbose=100000,
                      m.0=NULL, C.0=NULL,
                      mu.m0=NULL, mu.P0=NULL,
                      phi.m0=NULL, phi.P0=NULL,
                      W.a0=NULL, W.b0=NULL,
                      d.true=NULL, w.true=NULL,
                      beta.true=NULL, iota.true=NULL,
                      mu.true=NULL, phi.true=NULL, W.true=NULL)
{
  ## m.0 = prior mean for (iota,beta_0) or (beta_0).
  ## C.0 = prior var  for (iota,beta_0) or (beta_0).

  ## y: counts (T)
  ## X: the design matrix (including covariates for non-dynamic coef.) (T x P)
  ## n: the number of trials (T)

  ## m.0: prior mean for (iota,beta_0) or (beta_0).  (P)
  ## C.0: prior var  for (iota,beta_0) or (beta_0).  (P)

  ## mu: mean of (beta_t) (P.b)
  ## phi: ar coef. of (beta_t) (P.b)
  ## W: innovation variance of (beta_t) (P.b) -- ASSUMED DIAGONAL.

  ## A few things to keep in mind.
  ## 1) phi = 1 ==> mu = 0.
  ## 2) phi = 1 && X==1 ==> iota = 0.
  ## 3) iota = unknown ==> beta = unknown.

  ## Process for (alpha, beta_t):
  ## alpha_t = alpha_{t-1}
  ## beta_t \sim AR(1).
  
  ## Structure.
  y = as.matrix(y)
  X = cbind(X.stc, X.dyn)
  C.0 = as.matrix(C.0);
  
  ## Dimension ##
  T = nrow(X);
  P = ncol(X);
  P.b = ncol(X.dyn);
  P.a = P - P.b
  M = samp;

  ## Default prior parameters -- almost a random walk for beta ##
  if (is.null(m.0)    || is.null(C.0))    { m.0    = rep(0.0, P)   ; C.0    = diag(1.0, P  ); }
  if (is.null(mu.m0)  || is.null(mu.P0))  { mu.m0  = rep(0.0 ,P.b) ; mu.P0  = rep(100., P.b); }
  if (is.null(phi.m0) || is.null(phi.P0)) { phi.m0 = rep(0.99,P.b) ; phi.P0 = rep(100., P.b); }
  if (is.null(W.a0)   || is.null(W.b0))   { W.a0   = rep(1.0, P.b) ; W.b0   = rep(1.0,  P.b);  }

  ## Output data structure.
  out = list(
    d    = array(0, dim=c(M)),
    w    = array(0, dim=c(M, T)),
    beta = array(0, dim=c(M, P.b, T+1)),
    alpha= array(0, dim=c(M, 1)),
    mu   = array(0, dim=c(M, P.b)),
    phi  = array(0, dim=c(M, P.b)),
    W    = array(0, dim=c(M, P.b))
    )
  if (P.a > 0) out$alpha = array(0, dim=c(M, P.a))
  
  ## Initialize ##
  d    = 1
  beta = matrix(0.00, P.b, T+1);  # even odds.
  iota = rep(0.00, P.a);
  mu   = mu.m0;
  phi  = phi.m0;
  W    = W.b0 / W.a0;
  om   = rep(0, T);
  
  ## In case we are doing testing or we want to constrain to local level model.
  know.d    = FALSE;
  know.beta = FALSE;
  know.w    = FALSE;
  know.phi  = FALSE;
  know.mu   = FALSE;
  know.W    = FALSE;
  know.iota = FALSE;

  if (!is.null(d.true))    { d    = d.true;    know.d    = TRUE; }
  if (!is.null(beta.true)) { beta = beta.true; know.beta = TRUE; }
  if (!is.null(w.true))    { om   = w.true;    know.w    = TRUE; }
  if (!is.null(phi.true))  { phi  = phi.true;  know.phi  = TRUE;
                             if (any(phi==1)) {mu.true = rep(0, P.b);}}
  if (!is.null(mu.true))   { mu   = mu.true;   know.mu   = TRUE; }
  if (!is.null(W.true))    { W    = W.true;    know.W    = TRUE; }
  if (!is.null(iota.true)) { iota = iota.true; know.iota = TRUE; }

  ## Check that we are okay.
  if (know.beta && P.a > 0 && !know.iota) {printf("Know beta, not iota, X.stc!=NULL."); return(0);}
    
  ## Preprocess ## 
  ymax = max(y);
  F = cumsum(hist(y, breaks=0:(ymax+1)-0.5, plot=FALSE)$counts)
  G = T - F;
  
  start.time = proc.time();

  ## Sample
  for ( j in 1:(samp+burn) )
  {
    if (j==burn+1) start.ess = proc.time();
    
    ## draw (d, w | beta) --- WARNING: JOINT DRAW.
    log.mean = apply(X.dyn * t(beta)[-1,], 1, sum);
    if (P.a > 0) log.mean = log.mean + X.stc %*% iota;
    mu.lambda  = exp(log.mean)
    
    ## draw (d | beta) (w | d, beta)
    if (!know.d) d = draw.df(y, d, mu.lambda, G, ymax);
    psi = log.mean - log(d);
    w = rpg.devroye(T, y+d, psi);

    ## draw beta
    kappa = 0.5 * (y-d)
    z = kappa / w + log(d);
    ## ffbs = FFBS.C(z, X, 1/w, mu, phi, diag(W, P.b), m.0, C.0)
    ffbs = CUBS.C(z, X, 1/w, mu, phi, diag(W, P.b), m.0, C.0, obs="norm");
    iota = ffbs$alpha;
    beta = ffbs$beta;
    

    ## AR(1) - phi, W assumed to be diagonal !!!
    if (!know.mu)  mu  = draw.mu.ar1.ind (beta, phi, W, mu.m0, mu.P0)
    if (!know.phi) phi = draw.phi.ar1.ind(beta, mu, W, phi.m0, phi.P0, phi)
    if (!know.W)   W   = draw.W.ar1.ind  (beta, mu, phi, W.a0, W.b0)
    
    # Record if we are past burn-in.
    if (j>burn) {
      out$d[j-burn]      = d
      out$w[j-burn,]     = w
      out$alpha[j-burn,] = iota
      out$beta[j-burn,,] = beta
      out$mu[j-burn, ]   = mu;
      out$phi[j-burn, ]  = phi;
      out$W[j-burn, ]    = W;
    }

    if (j %% verbose == 0) { print(paste("Dyn NB PG: Iteration", j)); }
  }

  end.time  = proc.time();
  
  out$total.time = end.time - start.time;
  out$ess.time   = end.time - start.ess;
  
  out
} ## dyn.NB.PG

################################################################################
                                   ## Test ##
################################################################################

if (FALSE) {

  ## library("BayesLogit")
  ## dyn.unload("BayesLogit.so")
  ## dyn.load("BayesLogit.so")
  ## dyn.load("BayesLogit.so"); source("LogitWrapper.R"); source("FFBS.R");
  
  T = 500;
  P = 1;

  beta = array(0, dim=c(P, T+1));
  X.dyn = matrix(1, T, P);
  ## X.dyn = matrix(rnorm(T*P), nrow=T, ncol=P)
  ## X.stc = matrix(1, nrow=T, ncol=1)

  ## Parameters
  W = rep(0.05, P)
  phi = rep(0.95, P)
  mu = rep(0.5, P)
  d = 4

  ## Prior
  m.0    = mu;
  C.0    = diag(4.0, P);
  mu.m0  = mu;
  mu.P0  = rep(1/4, P);
  phi.m0 = rep(0.9, P);
  phi.P0 = rep(100, P);
  W.a0   = rep(10 , P);
  W.b0   = rep(10 , P);

  ## Synthetic
  beta[,1] = m.0;
  for (i in 2:(T+1)) {
    beta[,i] = mu + phi* (beta[,i-1] - mu) + sqrt(W) * rnorm(P);
  }

  log.mean  = apply(X.dyn * t(beta)[-1,], 1, sum)
  psi       = log.mean - log(d)
  mu.lambda = exp(log.mean)

  lambda = rgamma(T, d, scale=exp(psi));
  y = rpois(T, lambda);
  
  ## Simulate ##
  ## source("DynNBPG.R")
  samp = 500; burn=000; verbose=100;
  out <- dyn.NB.PG(y, X.dyn=X.dyn, X.stc=NULL,
                   samp=samp, burn=burn, verbose=verbose,
                   m.0=m.0, C.0=C.0,
                   mu.m0=NULL, mu.P0=NULL,
                   phi.m0=NULL, phi.P0=NULL,
                   W.a0=W.a0, W.b0=W.b0,
                   d.true=NULL, w.true=NULL,
                   beta.true=NULL, iota.true=NULL,
                   mu.true=mu, phi.true=rep(1,P), W.true=NULL)

  ## source("DynNBPG.R")
  samp = 500; burn=000; verbose=100;
  out <- dyn.NB.PG(y, X.dyn=X.dyn, X.stc=NULL,
                   samp=samp, burn=burn, verbose=verbose,
                   m.0=m.0, C.0=C.0,
                   mu.m0=mu.m0, mu.P0=mu.P0,
                   phi.m0=phi.m0, phi.P0=phi.P0,
                   W.a0=W.a0, W.b0=W.b0,
                   d.true=NULL, w.true=NULL,
                   beta.true=NULL, iota.true=NULL,
                   mu.true=NULL, phi.true=NULL, W.true=NULL)
}

if (FALSE) {

  beta.m = array(0, dim=c(P, T+1));
  for (i in 1:P) {
    beta.m[i,] = apply(out$beta[,i,], 2, mean);
  }

  ymin = min(beta.m, beta);
  ymax = max(beta.m, beta);

  plot(beta[1,], type="l", ylim=c(ymin, ymax))
  lines(beta.m[1,], col=2)
  points(log(y), col="gray")
  
}

################################################################################
                                 ## GEN DATA ##
################################################################################

if (FALSE)
{

  T = 500
  P = 4
  corr.type = "low"
  nb.mean = 1000  

  ## for (P in c(2,4)) {
  ##   for (corr.type in c("low", "high")) {
  ##     for (nb.mean in c(10, 100)) {
  
  c  = 0.5
  d.true  = 4
  marg.V   = 5 / sqrt(P) * c
  phi.true = rep(0.95, P)

  W.true = marg.V * (1 - phi.true^2)

  beta = matrix(0, nrow=P, ncol=T+1)
  beta[,1] = 0
  for (i in 2:(T+1))
    beta[,i] = phi.true * (beta[,i-1]) + rnorm(P, 0, sqrt(W.true))
  
  xgrid = seq(-1, 1, length.out=T)
  
  tX = matrix(0, nrow=P, ncol=T);
  if (corr.type=="low")  freq = c(1, 2, 3, 4)
  if (corr.type=="high") freq = c(1, 1.1, 1.2, 1.3)
  for (i in 1:P)
    tX[i,] = cos(freq[i] * pi * xgrid);

  tX = tX / sqrt(P) * (1-c)
  X  = t(tX)
  
  log.mean = log(nb.mean) + colSums(beta[,-1] * tX)
  psi = log.mean - log(d.true)
  p.success  = 1 / (1 + exp(-psi))
  y   = rnbinom(T, d.true, 1-p.success) ## p.success is prob of registering a single count.

  filename = paste("DynNB-synth-", corr.type, "-", P, "-mu-", nb.mean, ".RData", sep="")
  
  if (FALSE) {
    save(d.true, nb.mean, marg.V, phi.true, W.true, beta, tX, X, log.mean, psi, p.success, y, freq, xgrid,
         file=filename, compress=TRUE)
  }
  
}
