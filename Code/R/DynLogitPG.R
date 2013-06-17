## Binomial Logistic Regression using PolyaGamma augmentation.

## Depends:
## - FFBS.R or CUBS.R
## - Stationary.R

## if (exists("TESTING")) {
source("Stationary.R"); ## Independent AR(1)'s.  Maybe should change this.
## } ## TESTING

##------------------------------------------------------------------------------
## Binomial Logistic Regression.

dyn.logit.PG <- function(y, X.dyn, n=rep(1, length(y)), X.stc=NULL,
                         samp=1000, burn=100, verbose=100000,
                         m.0=NULL, C.0=NULL,
                         mu.m0=NULL, mu.P0=NULL,
                         phi.m0=NULL, phi.P0=NULL,
                         W.a0=NULL, W.b0=NULL,
                         beta.true=NULL, iota.true=NULL, w.true=NULL,
                         mu.true=NULL, phi.true=NULL, W.true=NULL)
{
  ## y: the counts (T)
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
  
  ## NOTE: We do not combine data. ##

  ## Structure.
  y = as.matrix(y)
  X = cbind(X.stc, X.dyn)
  ## W = as.matrix(W);
  C.0 = as.matrix(C.0);
  
  ## Dimension ##
  T = nrow(X);
  P = ncol(X);
  P.b = ncol(X.dyn);
  P.a = P - P.b
  M = samp;

  ## Default prior parameters -- almost a random walk for beta ##
  if (is.null(m.0)    || is.null(C.0))    { m.0    = rep(0.0, P)   ; C.0    = diag(1.0, P  ); }
  if (is.null(mu.m0)  || is.null(mu.P0))  { mu.m0  = rep(0.0 ,P.b) ; mu.P0  = rep(100 , P.b); }
  if (is.null(phi.m0) || is.null(phi.P0)) { phi.m0 = rep(0.99,P.b) ; phi.P0 = rep(100 , P.b); }
  if (is.null(W.a0)   || is.null(W.b0))   { W.a0   = rep(1.0, P.b) ; W.b0   = rep(1.0,  P.b);  }

  ## Output data structure ##
  out = list(
    w    = array(0, dim=c(M, T)),
    iota = array(0, dim=c(M, max(P.a, 1))),
    beta = array(0, dim=c(M, P.b, T+1)),
    mu   = array(0, dim=c(M, P.b)),
    phi  = array(0, dim=c(M, P.b)),
    W    = array(0, dim=c(M, P.b)),
    psi  = array(0, dim=c(M, T))
    )

  ## Initialize ##
  beta = matrix(0.00, P.b, T+1);  # even odds.
  iota = rep(0.00, P.a);
  mu   = matrix(0.00, P.b);
  phi  = matrix(0.99, P.b);
  W    = W.b0 / W.a0;
  om   = rep(0, T);
  
  ## In case we are doing testing or we want to constrain to local level model.
  know.beta = FALSE;
  know.w    = FALSE;
  know.phi  = FALSE;
  know.mu   = FALSE;
  know.W    = FALSE;
  know.iota = FALSE;
  
  if (!is.null(beta.true)) { beta = beta.true; know.beta = TRUE; }
  if (!is.null(w.true))    { om   = w.true;    know.w    = TRUE; }
  if (!is.null(phi.true))  { phi  = phi.true;  know.phi  = TRUE;
                             if (any(phi==1)) {mu.true = rep(0, P.b);}}
  if (!is.null(mu.true))   { mu   = mu.true;   know.mu   = TRUE; }
  if (!is.null(W.true))    { W    = W.true;    know.W    = TRUE; }
  if (!is.null(iota.true)) { iota = iota.true; know.iota = TRUE; }

  ## Check that we are okay.
  if (know.beta && P.a > 0 && !know.iota) {printf("Know beta, not iota, X.stc!=NULL."); return(0);}
  
  kappa = (y-n/2)
  
  ## SAMPLE ##

  start.time = proc.time();
 
  for ( j in 1:(samp+burn) )
  {
    if (j==burn+1) start.ess = proc.time();

    psi = apply(X.dyn * t(beta)[-1,], 1, sum);
    if (P.a > 0) { psi = psi + X.stc %*% iota; }
    
    ## draw om
    om = rpg.devroye(T, n, psi);

    ## Draw beta;
    z = kappa / om;
    ffbs = CUBS.C(z, X, 1/om, mu, phi, diag(W, P.b), m.0, C.0, obs="norm");
    iota = ffbs$alpha;
    beta = ffbs$beta;
    
    ## AR(1) - phi, W assumed to be diagonal !!!
    ## mu  = draw.mu.R(beta, phi, W, mu.m0, mu.V0) 
    ## phi = draw.phi.R(beta, mu, W, phi.m0, phi.V0, phi)
    ## W   = draw.W.R  (beta, mu, phi, W.a0, W.b0)
    if (!know.mu)  mu  = draw.mu.ar1.ind (beta, phi, W, mu.m0, mu.P0)
    if (!know.phi) phi = draw.phi.ar1.ind(beta, mu, W, phi.m0, phi.P0, phi)
    if (!know.W)   W   = draw.W.ar1.ind  (beta, mu, phi, W.a0, W.b0)

    # Record if we are past burn-in.
    if (j > burn) {
      out$w[j-burn,]      = om;
      out$iota[j-burn,]   = iota;
      out$beta[j-burn,,]  = beta;
      out$mu[j-burn, ]    = mu;
      out$phi[j-burn, ]   = phi;
      out$W[j-burn, ]     = W;
      out$psi[j-burn,]    = psi;
    }

    if (j %% verbose == 0) { cat("Dyn Logit PG: Iteration", j, "\n"); }
  }

  end.time = proc.time();
  
  out$total.time = end.time - start.time;
  out$ess.time   = end.time - start.ess;
  out$alpha = out$iota

  out
} ## logit.gibbs.R

################################################################################
                                 ## TEST 01 ##
################################################################################

if (FALSE) {

  T = 400;
  P = 4;

  beta = array(0, dim=c(P, T+1));
  X = matrix(0.1+rnorm(T*P), nrow=T, ncol=P);

  N = nrow(X);

  ## Parameters
  iota = 0;
  W   = rep(0.1, P);
  mu  = rep(1.0, P);
  phi = rep(0.9, P)

  ## Prior
  b.m0 = rep(0.0, P);
  b.C0 = diag(2.0, P);
  mu.m0 = mu
  mu.P0 = rep(0.1, P)
  phi.m0 = rep(0.9, P)
  phi.P0 = rep(100, P);
  W.a0   = rep(10, P);
  W.b0   = W.a0 * W;

  ## Synthetic
  beta[,1] = mu;
  for (i in 2:(T+1)) {
    beta[,i] = mu + phi * (beta[,i-1]-mu) + sqrt(W) * rnorm(P);
  }

  psi = iota + apply(X * t(beta)[-1,], 1, sum); ## log-odds
  p   = exp(psi) / (1 + exp(psi));
  y = rbinom(N, 1, p);
  w = rpg.devroye(N, 1, psi);
  
  ## Simulate
  source("DynLogitPG.R")
  samp = 1000
  burn = 0
  out <- dyn.logit.PG(y, X, samp=samp, burn=burn, verbose=100,
                      m.0=b.m0, C.0=b.C0,
                      mu.m0=mu.m0, mu.P0=mu.P0,
                      phi.m0=phi.m0, phi.P0=phi.P0,
                      W.a0=W.a0, W.b0=W.b0,
                      beta.true=NULL, iota.true=NULL, w.true=NULL,
                      mu.true=mu, phi.true=phi, W.true=NULL)

  source("DynLogitPG.R")
  samp = 1000
  burn = 0
  out <- dyn.logit.PG(y, X, samp=samp, burn=burn, verbose=100,
                      m.0=b.m0, C.0=b.C0,
                      mu.m0=mu.m0, mu.P0=mu.P0,
                      phi.m0=phi.m0, phi.P0=phi.P0,
                      W.a0=W.a0, W.b0=W.b0,
                      beta.true=NULL, iota.true=NULL, w.true=NULL,
                      mu.true=NULL, phi.true=NULL, W.true=NULL)
  
  ess = apply(out$beta[, 1, ], 2, ESS);
  mean(ess)
  
}

if (FALSE) {

  beta.mean = apply(out$beta, c(2,3), mean);
  beta.95   = apply(out$beta, c(2,3), function(x){quantile(x, 0.95)});
  beta.05   = apply(out$beta, c(2,3), function(x){quantile(x, 0.05)});

  ymin = min(beta.mean, out$beta);
  ymax = max(beta.mean, out$beta);

  i = 4
  plot(0:T, beta.mean[i,], col=2, type="l", ylim=c(ymin,ymax));
  lines(0:T, beta.95[i,], col="pink")
  lines(0:T, beta.05[i,], col="pink")
  abline(h=mean(beta.mean[1,]), col=2, lty=c(2,2));

  lines(0:T, beta[i,]);
  
  points(1:T, y*(ymax-ymin) + ymin, cex=0.1);
  
  lines(0:T, out$beta[300,i,], col=3);
  
}

################################################################################
                                 ## TEST 02 ##
################################################################################

if (FALSE) {

  T = 500;
  P = 1;

  beta = array(0, dim=c(P, T+1));
  X = matrix(1, nrow=T, ncol=P);
  n = rep(2, T)

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

  psi = iota + apply(X * t(beta)[-1,], 1, sum); ## log-odds
  p   = exp(psi) / (1 + exp(psi));
  y = rbinom(T, n, p);
  w = rpg.devroye(T, n, psi);
  
  ## Simulate
  source("DynLogitPG.R")
  samp = 500
  burn = 0
  out <- dyn.logit.PG(y, X, n, samp=samp, burn=burn, verbose=100,
                      m.0=b.m0, C.0=b.C0,
                      mu.m0=NULL, mu.P0=NULL,
                      phi.m0=NULL, phi.P0=NULL,
                      W.a0=W.a0, W.b0=W.b0,
                      beta.true=NULL, iota.true=iota, w.true=NULL,
                      mu.true=0.0, phi.true=1.0, W.true=NULL)
  
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
  
  points(1:T, (y/n)*(ymax-ymin) + ymin, cex=0.1);
  
  lines(0:T, out$beta[300,1,], col=3);
  
}


################################################################################
                              ## TOKYO RAINFALL ##
################################################################################

## TOKYO RAINFALL ##
if (FALSE) {
  tokyo = read.csv("~/RPackage/BayesLogit/Code/R/DataSets/tokyo-rain.csv");
  tkrain = as.numeric(na.omit(as.numeric(as.matrix(tokyo))))

  T = length(tkrain)
  y = tkrain
  X = matrix(1, nrow=T, ncol=1)
  n = rep(2, T);
  
  ## Prior
  iota = 0.0
  b.m0 = -1.0;
  b.C0 = 1.0;
  W      = 0.1;
  W.a0   = 300;
  W.b0   = W.a0 * W;

  ## Simulate
  source("DynLogitPG.R")
  samp = 500
  burn = 0
  out <- dyn.logit.PG(y, X, n,
                      samp=samp, burn=burn, verbose=100,
                      m.0=b.m0, C.0=b.C0,
                      mu.m0=NULL, mu.P0=NULL,
                      phi.m0=NULL, phi.P0=NULL,
                      W.a0=W.a0, W.b0=W.b0,
                      beta.true=NULL, iota.true=0, w.true=NULL,
                      mu.true=0.0, phi.true=1.0, W.true=NULL)

  ess = apply(out$beta[, 1, ], 2, ESS);
  mean(ess)
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

  points(1:T, (y/n)*(ymax-ymin) + ymin, cex=0.1);
  
  lines(0:T, out$beta[300,1,], col=3);
  
}

################################################################################
                                ## GEN SYNTH ##
################################################################################

if (FALSE)
{

  T = 500
  P = 2
  corr.type = "low"
  ntrial = 1

  ## for (P in c(2,4)) {
  ##   for (corr.type in c("low", "high")) {
  ##     for (ntrial in c(1, 10, 100)) {

  c  = 0.5
  psi.mean = 1.4
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

  n = rep(ntrial, T)
  
  psi = psi.mean + colSums(beta[,-1] * tX)
  p   = 1 / (1 + exp(-psi))
  y   = rbinom(T, n, p)

  filename = paste("DynLogit-synth-", corr.type, "-", P, "-n-", ntrial, ".RData", sep="")
  
  if (FALSE) {
    save(T, P, psi.mean, marg.V, phi.true, W.true, beta, xgrid, X, n, psi, p, y, freq, c,
         file=filename, compress=TRUE)
  }
  
}
