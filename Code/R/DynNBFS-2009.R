## Negative binomial regression using FSF's log-gamma parameterization.
## We model the log-mean here.

## if (exists("TESTING")) {
if (!is.loaded("BayesLogit.so")) dyn.load("../C/BayesLogit.so");
source("ComputeMixture.R")
source("FFBS.R")
source("NBDF.R") ## Routine for sampling shape.
source("Stationary.R"); ## Independent AR(1)'s.  Maybe should change this.
## } ## TESTING

################################################################################

##------------------------------------------------------------------------------

draw.indicators.R <- function(res, nmix)
{
  ## y.u - N x 1 - latent variable y^u in paper.
  ## lambda = X beta

  nmix$s = sqrt(nmix$v)
  
  log.wds = log(nmix$p) - log(nmix$s);
  
  ## Safer to work on log scale.  Columns correspond to outcome index i!  Watch
  ## out for sd in outer product when adapting to other routines.
  log.post = -0.5 * outer(-nmix$m, res, "+")^2 / nmix$v + log.wds;
  unnrm.post = exp(log.post);

  ## Now sample. 
  r = apply(unnrm.post, 2, function(x){sample.int(n=nmix$nc, size=1, prob=x)})
}  ## draw.indicators

draw.indicators.C <- function(res, nmix)
{
  n = length(res);
  r = rep(0, n);
  
  OUT <- .C("draw_indicators_generic",
            as.integer(r), as.double(res), as.integer(n),
            as.double(nmix$p), as.double(nmix$m), as.double(sqrt(nmix$v)), as.integer(nmix$nc))

  OUT[[1]]
} ## draw.indicators.C

################################################################################

##------------------------------------------------------------------------------

dyn.NB.FS <- function(y, X.dyn, X.stc=NULL,
                      samp=1000, burn=100, verbose=100000,
                      m.0=NULL, C.0=NULL,
                      mu.m0=NULL, mu.V0=NULL,
                      phi.m0=NULL, phi.V0=NULL,
                      W.a0=NULL, W.b0=NULL,
                      d.true=NULL, lambda.true=NULL, r.true=NULL,
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
  if (is.null(W.a0)   || is.null(W.b0))   { W.a0   = rep(1.0, P.b) ; W.b0   = rep(1.0,  P.b); }

  ## Output data structure.
  out = list(
    d    = array(0, dim=c(M)),
    iota   = array(0, dim=c(M, max(P.a, 1))),
    beta   = array(0, dim=c(M, P.b, T+1)),
    lambda = array(0, dim=c(M, T)),
    r      = array(0, dim=c(M, T)),
    mu     = array(0, dim=c(M, P.b)),
    phi    = array(0, dim=c(M, P.b)),
    W      = array(0, dim=c(M, P.b))
    )

  ## Initialize ##
  d    = 1
  beta = matrix(0.00, P.b, T+1);  # mean = 1
  iota = rep(0.00, P.a);
  mu   = matrix(mu.m0, P.b);
  phi  = matrix(phi.m0, P.b);
  W    = rep(W.b0 / W.a0, P.b);
  
  ## In case we are doing testing or we want to constrain to local level model.
  know.d    = FALSE;
  know.beta = FALSE;
  know.lambda = FALSE;
  know.r    = FALSE;
  know.phi  = FALSE;
  know.mu   = FALSE;
  know.W    = FALSE;
  know.iota = FALSE;

  if (!is.null(d.true))    { d    = d.true;    know.d    = TRUE; }
  if (!is.null(beta.true)) { beta = beta.true; know.beta = TRUE; }
  if (!is.null(lambda.true)) { lambda = lambda.true; know.lambda = TRUE; }
  if (!is.null(r.true))    { r    = r.true;    know.r    = TRUE; }
  if (!is.null(phi.true))  { phi  = phi.true;  know.phi  = TRUE;
                             if (phi==1) {mu.true = rep(0, P.b);}}
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

  ## SAMPLE ##
   
  for ( j in 1:(samp+burn) )
  {
    if (j==burn+1) start.ess = proc.time();

    ## draw (d, lambda, r | beta) --- WARNING: JOINT DRAW.
    log.mean = apply(X.dyn * t(beta)[-1,], 1, sum);
    if (P.a > 0) psi = log.mean + X.stc %*% iota;
    mu.lambda = exp(log.mean)
    
    ## draw (d | beta)
    d = draw.df(d, mu.lambda, G, ymax);

    ## draw (lambda | d, beta)
    psi = log.mean - log(d);
    p = 1 / (1 + exp(-psi))
    lambda = rgamma(T, y+d, scale=p)

    ## draw (r | d, lambda, beta)
    nmix = compute.mixture(d);
    res  = psi - log(lambda)
    r    = draw.indicators.C(res, nmix);

    ## draw beta
    z = log(lambda) + log(d) + nmix$m[r]; ## So that we model the log mean.
    ffbs = FFBS.C(z, X, mu, phi, diag(W,P.b), nmix$v[r], m.0, C.0)
    iota = ffbs$alpha
    beta = ffbs$beta

    ## AR(1) - phi, W assumed to be diagonal !!!
    ## mu  = draw.mu.R(beta, phi, W, mu.m0, mu.V0) 
    ## phi = draw.phi.R(beta, mu, W, phi.m0, phi.V0, phi)
    W   = draw.W.R  (beta, mu, phi, W.a0, W.b0)
    
    # Record if we are past burn-in.
    if (j>burn) {
      out$d[j-burn]       = d
      out$lambda[j-burn,] = lambda
      out$r[j-burn,]      = r
      out$beta[j-burn,,]  = beta
      out$mu[j-burn, ]    = mu;
      out$phi[j-burn, ]   = phi;
      out$W[j-burn, ]     = W;
    }

    if (j %% verbose == 0) cat("Dyn NB FS: Iteration", j, "\n");

  }

  end.time  = proc.time();
  out$total.time = end.time - start.time;
  out$ess.time   = end.time - start.ess;
  
  out
} ## dyn.NB.mix

################################################################################
                                   ## MAIN ##
################################################################################

if (FALSE) {

  dyn.load("FSF_nmix.so")
  
  T = 500;
  P = 1;

  beta = array(0, dim=c(P, T+1));
  X = matrix(1, T, P);

  ## Parameters
  W = 0.1;
  phi = 0.95;
  mu = 0.0
  d = 4
  nmix = compute.mixture(d)

  ## Prior
  m.0     = 1.0;
  C.0     = 4.0;
  mu.m0  = 0.0;
  mu.V0  = 4.0;
  phi.m0 = 0.9;
  phi.V0 = 0.1;
  W.a0   = 10;
  W.b0   = 1.0;

  ## Synthetic
  beta[,1] = m.0;
  for (i in 2:(T+1)) {
    beta[,i] = mu + phi* (beta[,i-1] - mu) + sqrt(W) * rnorm(1);
  }

  log.mean  = apply(X * t(beta)[-1,], 1, sum)
  psi       = log.mean - log(d)
  mu.lambda = exp(log.mean)

  r = sample.int(nmix$nc, T, replace = TRUE, nmix$p);
  ep = rnorm(T, 1*nmix$m[r], sqrt(nmix$v[r]));
  lambda = (exp(psi - ep));
  y = rpois(T, lambda);
  
  ## Simulate ##
  source("DynNBFS-2009.R")
  samp = 500; burn=0;
  out <- dyn.NB.FS(y, X.dyn=X, X.stc=NULL,
                   samp=samp, burn=burn, verbose=50,
                   m.0=m.0, C.0=C.0,
                   mu.m0=mu.m0, mu.V0=mu.V0,
                   phi.m0=phi.m0, phi.V0=phi.V0,
                   W.a0=W.a0, W.b0=W.b0,
                   d.true=NULL, lambda.true=NULL, r.true=NULL, 
                   beta.true=NULL, iota.true=NULL,
                   mu.true=mu, phi.true=1.0, W.true=NULL)
  
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
