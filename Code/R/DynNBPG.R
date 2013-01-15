## Negative binomial regression using PG augmentation.
## We model the log-mean here.

source("NB-Shape.R")
source("Stationary.R"); ## Independent AR(1)'s.  Maybe should change this.

################################################################################
                             ## Dynamc NB by PG ##
################################################################################

dyn.NB.PG <- function(y, X.dyn, X.stc=NULL,
                      samp=1000, burn=100, verbose=100000,
                      m.0=NULL, C.0=NULL,
                      mu.m0=NULL, mu.V0=NULL,
                      phi.m0=NULL, phi.V0=NULL,
                      W.a0=NULL, W.b0=NULL,
                      d.true=NULL, w.true=NULL,
                      beta.true=NULL, iota.true=NULL,
                      mu.true=NULL, phi.true=NULL, W.true=NULL)
{
  ## m.0 = prior mean for (iota,beta_0) or (beta_0).
  ## C.0 = prior var  for (iota,beta_0) or (beta_0).

  ## y: the average response (T)
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
  W = as.matrix(W);
  C.0 = as.matrix(C.0);
  
  ## Dimension ##
  T = nrow(X);
  P = ncol(X);
  P.b = ncol(X.dyn);
  P.a = P - P.b
  M = samp;

  ## Default prior parameters -- almost a random walk for beta ##
  if (is.null(m.0)    || is.null(C.0))    { m.0    = rep(0.0, P)   ; C.0    = diag(1.0, P  ); }
  if (is.null(mu.m0)  || is.null(mu.V0))  { mu.m0  = rep(0.0 ,P.b) ; mu.V0  = rep(0.01, P.b); }
  if (is.null(phi.m0) || is.null(phi.V0)) { phi.m0 = rep(0.99,P.b) ; phi.V0 = rep(0.01, P.b); }
  if (is.null(W.a0)   || is.null(W.b0))   { W.a0   = rep(1.0, P.b) ; W.b0   = rep(1.0,  P.b);  }

  ## Output data structure.
  out = list(
    d    = array(0, dim=c(M)),
    w    = array(0, dim=c(M, T)),
    iota = array(0, dim=c(M, max(P.a, 1))),
    beta = array(0, dim=c(M, P.b, T+1)),
    mu   = array(0, dim=c(M, P.b)),
    phi  = array(0, dim=c(M, P.b)),
    W    = array(0, dim=c(M, P.b))
    )

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
                             if (phi[1]==1) {mu.true = rep(0, P.b);}}
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
    if (P.a > 0) psi = log.mean + X.stc %*% iota;
    mu.lambda  = exp(log.mean)
    ## draw (d | beta) (w | d, beta)
    d = draw.df(d, mu.lambda, G, ymax);
    psi = log.mean - log(d);
    w = rpg.devroye(T, y+d, psi);

    ## draw beta
    kappa = 0.5 * (y-d)
    z = kappa / w + log(d);
    ffbs = FFBS.C(z, X, 1/w, mu, phi, diag(W, P.b), m.0, C.0)
    iota = ffbs$alpha
    beta = ffbs$beta
    
    # Record if we are past burn-in.
    if (j>burn) {
      out$d[j-burn]      = d
      out$w[j-burn,]     = w
      out$iota[j-burn]   = iota
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
  P = 2;

  beta = array(0, dim=c(P, T+1));
  X = matrix(1, T, P);

  ## Parameters
  W = rep(0.1, P)
  phi = rep(0.95, P)
  mu = rep(3.0, P)
  d = 4

  ## Prior
  m.0    = mu;
  C.0    = diag(4.0, P);
  mu.m0  = rep(mu , P);
  mu.V0  = rep(4.0, P);
  phi.m0 = rep(0.9, P);
  phi.V0 = rep(0.1, P);
  W.a0   = rep(10 , P);
  W.b0   = rep(1.0, P);

  ## Synthetic
  beta[,1] = m.0;
  for (i in 2:(T+1)) {
    beta[,i] = mu + phi* (beta[,i-1] - mu) + sqrt(W) * rnorm(P);
  }

  log.mean  = apply(X * t(beta)[-1,], 1, sum)
  psi       = log.mean - log(d)
  mu.lambda = exp(log.mean)

  lambda = rgamma(T, d, scale=exp(psi));
  y = rpois(T, lambda);
  
  ## Simulate ##
  ## source("DynNBPG.R")
  samp = 500; burn=000; verbose=100;
  out <- dyn.NB.PG(y, X.dyn=X, X.stc=NULL,
                   samp=samp, burn=burn, verbose=verbose,
                   m.0=m.0, C.0=C.0,
                   mu.m0=NULL, mu.V0=NULL,
                   phi.m0=NULL, phi.V0=NULL,
                   W.a0=W.a0, W.b0=W.b0,
                   d.true=NULL, w.true=NULL,
                   beta.true=NULL, iota.true=NULL,
                   mu.true=mu, phi.true=rep(1,P), W.true=NULL)
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
