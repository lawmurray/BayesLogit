## Dynamic Binomial logistic regression using FS and Fussl 2012.

source("compmix.R")

################################################################################

## Taken from binomlogit package.
draw.yStar <- function(lambda, y, n){
  N = length(y)
  U = rgamma(N, shape=n, rate=1+lambda)
  V = rgamma(N, shape=y, rate=1)
  W = rgamma(N, shape=n-y, rate=lambda)
  -log((U+(y < n)*W)/(U+(y > 0)*V))
}

################################################################################

dyn.logit.dRUM <- function(y, X.dyn, n=1, X.stc=NULL,
                           samp=1000, burn=100, verbose=100,
                           m.0=NULL, C.0=NULL,
                           mu.m0=NULL, mu.P0=NULL,
                           phi.m0=NULL, phi.P0=NULL,
                           W.a0=NULL, W.b0=NULL,
                           z.true=NULL, r.true=NULL,
                           beta.true=NULL, iota.true=NULL,
                           mu.true=NULL, phi.true=NULL, W.true=NULL)
{
  ## m.0 = prior mean for (iota,beta_0) or (beta_0).
  ## C.0 = prior var  for (iota,beta_0) or (beta_0).

  ## y: binomial repsonse in {0:n} (T)
  ## X: the design matrix (including covariates for non-dynamic coef.) (T x P)
  ## n: the number of trials (1)

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
    z    = array(0, dim=c(M, T)),
    r    = array(0, dim=c(M, T)),
    iota = array(0, dim=c(M, max(P.a, 1))),
    beta = array(0, dim=c(M, P.b, T+1)),
    mu   = array(0, dim=c(M, P.b)),
    phi  = array(0, dim=c(M, P.b)),
    W    = array(0, dim=c(M, P.b))
    )

  ## Initialize ##
  beta = matrix(0.00, P.b, T+1);  # even odds.
  iota = rep(0.00, P.a);
  mu   = matrix(mu.m0, P.b);
  phi  = matrix(phi.m0, P.b);
  W    = W.b0 / W.a0;
  with.iota = TRUE
  n.rep = rep(n, T)

  ## Using compmix nmix.n ~ Type III logis.
  nmix.n    = as.list(compmix(n));
  nmix.n$nc = length(nmix.n$p)
  nmix.n$m  = rep(0, nmix.n$nc);
  
  ## In case we are doing testing or we want to constrain to local level model.
  know.beta = FALSE;
  know.z    = FALSE;
  know.r    = FALSE;
  know.phi  = FALSE;
  know.mu   = FALSE;
  know.W    = FALSE;
  know.iota = FALSE;

  if (!is.null(beta.true)) { beta = beta.true; know.beta = TRUE; }
  if (!is.null(z.true))    { z    = z.true;    know.z    = TRUE; }
  if (!is.null(r.true))    { r    = r.true;    know.r    = TRUE; }
  if (!is.null(phi.true))  { phi  = phi.true;  know.phi  = TRUE;
                             if (any(phi==1)) {mu.true = rep(0, P.b);}}
  if (!is.null(mu.true))   { mu   = mu.true;   know.mu   = TRUE; }
  if (!is.null(W.true))    { W    = W.true;    know.W    = TRUE; }
  if (!is.null(iota.true)) { iota = iota.true; know.iota = TRUE; }

  ## Check that we are okay.
  if (know.beta && P.a > 0 && !know.iota) {printf("Know beta, not iota, X.stc!=NULL."); return(0);}
  
  ## SAMPLE ##
  
  start.time = proc.time();
  
  for (j in 1:(samp+burn)) {
    if (j==burn+1) start.ess = proc.time();
    
    ## JOINT DRAW: (z, r | beta, y) = (z | beta, y) (r | z, beta, y)
    psi = apply(X.dyn * t(beta)[-1,], 1, sum);
    if (P.a > 0) psi = psi + X.stc %*% iota;
    lambda = exp(psi);

    ## USING compmix
    z = draw.yStar(lambda, y, n.rep)
    r = draw.indicators.C(z-psi, nmix.n);
    V.hat = nmix.n$v[r];

    ## (iota, beta | r, z, y)
    ffbs = CUBS.C(z, X, V.hat, mu, phi, diag(W, P.b), m.0, C.0, obs="norm");
    iota = ffbs$alpha
    beta = ffbs$beta

    ## AR(1) - phi, W assumed to be diagonal, i.e. vectors in R!!!
    if (!know.mu)  mu  = draw.mu.ar1.ind (beta, phi, W, mu.m0, mu.P0)
    if (!know.phi) phi = draw.phi.ar1.ind(beta, mu, W, phi.m0, phi.P0, phi)
    if (!know.W)   W   = draw.W.ar1.ind  (beta, mu, phi, W.a0, W.b0)
    
    if (j > burn) {
      out$z   [j-burn,]  = z;
      out$r   [j-burn,]  = r;
      out$iota[j-burn]   = iota;
      out$beta[j-burn,,] = beta;
      out$mu  [j-burn, ] = mu;
      out$phi [j-burn, ] = phi;
      out$W   [j-burn, ] = W;
    }

    if (j %% verbose == 0) cat("Dyn Logit dRUM: Iteration", j, "\n");
    
  }

  end.time = proc.time()

  out$total.time = end.time-start.time
  out$ess.time   = end.time-start.ess
  out$alpha = out$iota
  
  out
} ## dyn.logit.dRUM

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
  source("DynLogitdRUM.R")
  samp = 500
  burn = 0
  out <- dyn.logit.FS(y, X.dyn=X, n=rep(2, T), X.stc=NULL,
                      samp=samp, burn=burn, verbose=100,
                      m.0=b.m0, C.0=b.C0,
                      mu.m0=NULL, mu.P0=NULL,
                      phi.m0=NULL, phi.P0=NULL,
                      W.a0=W.a0, W.b0=W.b0,
                      z.true=NULL, r.true=NULL,
                      beta.true=NULL, iota.true=NULL,
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
