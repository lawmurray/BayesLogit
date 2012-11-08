## Dynamic BINARY logistic regression using FS's 2010 discrete mixture of
## normals.

if (!is.loaded("FSF_nmix.so")) dyn.load("../C/FSF_nmix.so");

source("FFBS.R")
source("Stationary.R"); ## Independent AR(1)'s.  Maybe should change this.

################################################################################
                              ## NORMAL MIXTURE ##
################################################################################

## This R script uses Fruhwirth-Schnatter and Fruhwirth's normal-mixture
## approximation to logistic regression (2010).

## Based on Monahan & Stefanski means based on their method.  They do not use
## the numbers in their article.  I'm not sure what it is.  They use the same
## probabilities but different variances?

## ## Define normal mixture -- 6 comp.  FS&F p. 119, based on Monahan & Stefanski (1992).
## normal.mixture = list(
##   w = c(1.8446, 17.268, 37.393, 31.697, 10.89, 0.90745) / 100,
##   m = rep(0, 6),
##   v = c(0.68159, 1.2419, 2.2388, 4.0724, 7.4371, 13.772)
##   )
##   normal.mixture$s = sqrt(normal.mixture$v)
##   normal.mixture$N = length(normal.mixture$w)
## c(1.21126, 0.89734, 0.66833, 0.49553, 0.36668, 0.26946)

## ## Define normal mixture -- 3 comp.  FS&F p. 119, based on Monahan & Stefanski (1992).
## normal.mixture = list(
##   w = c(25.22, 58.523, 16.257) / 100,
##   m = rep(0, 3),
##   v = c(1.2131, 2.9955, 7.5458)
##   )
##   normal.mixture$s = sqrt(normal.mixture$v)
##   normal.mixture$N = length(normal.mixture$w)

## Define normal mixture -- 6 comp.  FS&F p. 119, based on K-L distance.
normal.mixture = list(
  w = c(5.8726, 28.74, 36.756, 22.427, 5.8701, 0.33466) / 100,
  m = rep(0, 6),
  v = c(0.84678, 1.61, 2.8904, 5.0772, 8.9109, 15.923)
  )
  normal.mixture$s = sqrt(normal.mixture$v)
  normal.mixture$N = length(normal.mixture$w)
  normal.mixture$p = normal.mixture$w

## lgs.samp = -log(rexp(10000)) + log(rexp(10000));

## ## Define normal mixture -- 8 comp.  Monahan & Stefanski (1992).  They use
## ## inverse-scale.  NOTE! Initially, I only used the first 5-digits.  That was
## ## WAY off (in terms of the estimated variance).  You need a lot of precision.
## normal.mixture = list(
##   w = c(0.0032463432, 0.0515174770, 0.1950779126, 0.3155698236,
##         0.2741495761, 0.1310768806, 0.0279241871, 0.0014495678),
##   m  = rep(0, 8),
##   is = c(1.3653408062, 1.0595239710, 0.8307913137, 0.6507321666,
##          0.5081354253, 0.3963122451, 0.3089042522, 0.2382126164)
##   )
##   normal.mixture$s = 1 / normal.mixture$is
##   normal.mixture$v = normal.mixture$s^2
##   normal.mixture$N = length(normal.mixture$w)

lgs.samp = rlogis(10000)
normal.mixture$mar.mean = 0
normal.mixture$mar.var  = pi^2 / 3;

## Make a copy.
NM = normal.mixture

################################################################################

##------------------------------------------------------------------------------

draw.indicators <- function(z, lambda)
{
  ## y.u - N x 1 - latent variable y^u in paper.
  ## lambda = X beta

  res = z - log(lambda)
  log.wds = log(NM$w) - log(NM$s);

  ## Safer to work on log scale.  Columns correspond to outcome index i!
  log.post = -0.5 * outer(1/NM$s, res, "*")^2 + log.wds;
  unnrm.post = exp(log.post);

  ## Now sample. 
  r = apply(unnrm.post, 2, function(x){sample.int(n=NM$N, size=1, prob=x)})
}  ## draw.indicators

draw.indicators.C <- function(z, lambda, nmix)
{
  n = length(z);
  r = rep(0, n);
  
  OUT <- .C("draw_indicators_logistic",
            as.integer(r), as.double(z), as.double(lambda), as.integer(n),
            as.double(nmix$w), as.double(nmix$s), as.integer(nmix$N))

  OUT[[1]]
} ## draw.indicators.C

##------------------------------------------------------------------------------

draw.z <- function(lambda, y){
  n = length(lambda)
  u = runif(n)
  z = log(lambda * u + y) - log(1 - u + lambda * (1-y));
  z
} ## draw.z

################################################################################

dyn.logit.FS <- function(y, X.dyn, n=1, X.stc=NULL,
                         samp=1000, burn=100, verbose=100,
                         m.0=NULL, C.0=NULL,
                         mu.m0=NULL, mu.V0=NULL,
                         phi.m0=NULL, phi.V0=NULL,
                         W.a0=NULL, W.b0=NULL,
                         z.true=NULL, r.true=NULL,
                         beta.true=NULL, iota.true=NULL,
                         mu.true=NULL, phi.true=NULL, W.true=NULL)
{
  ## m.0 = prior mean for (iota,beta_0) or (beta_0).
  ## C.0 = prior var  for (iota,beta_0) or (beta_0).

  ## y: binomial repsonse in {0:n} (T)
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

  ## Output data structure ##
  out = list(
    z    = array(0, dim=c(M, T*n)),
    r    = array(0, dim=c(M, T*n)),
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
  nmix = NM
  
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
                             if (phi==1) {mu.true = rep(0, P.b);}}
  if (!is.null(mu.true))   { mu   = mu.true;   know.mu   = TRUE; }
  if (!is.null(W.true))    { W    = W.true;    know.W    = TRUE; }
  if (!is.null(iota.true)) { iota = iota.true; know.iota = TRUE; }

  ## Check that we are okay.
  if (know.beta && P.a > 0 && !know.iota) {printf("Know beta, not iota, X.stc!=NULL."); return(0);}

  ## For binomial.
  expand = as.numeric( matrix(1:T, nrow=n, ncol=T, byrow=TRUE) );
  samp.y = function(x){ sample(c(rep(1,x),rep(0,n-x)), n) };
  
  ## SAMPLE ##
  
  start.time = proc.time();
  
  for (j in 1:(samp+burn)) {
    if (j==burn+1) start.ess = proc.time();
    
    ## JOINT DRAW: (z, r | beta, y) = (z | beta, y) (r | z, beta, y)
    psi = apply(X.dyn * t(beta)[-1,], 1, sum);
    if (P.a > 0) psi = psi + X.stc %*% iota;
    lambda = exp(psi);

    lambda.hat = lambda[expand]; ## Need for binomail.  Remove hats for pure binary.
    y.hat      = as.numeric(apply(y, 1, samp.y));
    
    z    = draw.z(lambda.hat, y.hat);
    r    = draw.indicators.C(z, lambda.hat, nmix)

    z.hat = apply(matrix(z        , nrow=n), 2, sum) / n;
    V.hat = apply(matrix(nmix$v[r], nrow=n), 2, sum) / n^2; 
    ## V.hat = nmix$v[r] ## when n = 1.
    
    ## (iota, beta | r, z, y)
    ffbs = FFBS.C(z.hat, X, mu, phi, diag(W,P), V.hat, m.0, C.0)
    iota = ffbs$alpha
    beta = ffbs$beta

    ## AR(1) - phi, W assumed to be diagonal !!!
    ## mu  = draw.mu.R(beta, phi, W, mu.m0, mu.V0) 
    ## phi = draw.phi.R(beta, mu, W, phi.m0, phi.V0, phi)
    W   = draw.W.R  (beta, mu, phi, W.a0, W.b0)
    
    if (j > burn) {
      out$z [j-burn,]    = z;
      out$r   [j-burn,]  = r;
      out$iota[j-burn]   = iota;
      out$beta[j-burn,,] = beta;
      out$mu[j-burn, ]   = mu;
      out$phi[j-burn, ]  = phi;
      out$W[j-burn, ]    = W;
    }

    if (j %% verbose == 0) cat("Dyn Binary Logit FS: Iteration", j, "\n");
    
  }

  end.time = proc.time()

  out$total.time = end.time-start.time
  out$ess.time   = end.time-start.ess
  
  out
} ## dyn.logit.FS

################################################################################
                                 ## TEST 01 ##
################################################################################

if (FALSE) {

  dyn.load("FSF_nmix.so")
  
  T = 500;
  P = 1;

  beta = array(0, dim=c(P, T+1));
  X = matrix(1, nrow=T, ncol=P);

  N = nrow(X);

  ## Parameters
  W   = 0.1;
  mu  = 2.0;
  phi = 0.95

  ## Prior
  m0     = mu;
  C0     = 2.0;
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

  nmix = NM
  psi = apply(X * t(beta)[-1,], 1, sum); ## log-odds
  p   = exp(psi) / (1 + exp(psi));
  r   = sample.int(nmix$N, N, replace=TRUE, prob=nmix$w);
  ep  = rnorm(N, 0, nmix$s[r]);
  z   = psi + ep;
  y   = as.numeric(z>0);
  
  ## Simulate
  source("DynLogitFS-2010.R")
  samp = 300
  burn = 0
  out <- dyn.logit.FS(y, X.dyn=X, n=1, X.stc=NULL,
                      samp=samp, burn=burn, verbose=100,
                      m.0=m0, C.0=C0,
                      mu.m0=NULL, mu.V0=NULL,
                      phi.m0=NULL, phi.V0=NULL,
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

  lines(0:T, beta[1,]);
  
  points(1:T, y*(ymax-ymin) + ymin, cex=0.1);
  
  lines(0:T, out$beta[300,1,], col=3);
  
}

################################################################################
                                 ## TEST 02 ##
################################################################################

if (FALSE) {

  dyn.load("FSF_nmix.so")
  
  T = 500;
  P = 1;
  n = 2;

  beta = array(0, dim=c(P, T+1));
  X = matrix(1, nrow=T, ncol=P);

  N = nrow(X);

  ## Parameters
  W   = 0.1;
  mu  = 2.0;
  phi = 0.95

  ## Prior
  m0     = mu;
  C0     = 2.0;
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

  epnd = as.numeric( matrix(1:T, nrow=n, ncol=T, byrow=TRUE) );
  
  nmix = NM
  psi = apply(X * t(beta)[-1,], 1, sum); ## log-odds
  p   = exp(psi) / (1 + exp(psi));
  r   = sample.int(nmix$N, N*n, replace=TRUE, prob=nmix$w);
  ep  = rnorm(N*n, 0, nmix$s[r]);
  z   = psi[epnd] + ep;
  y   = apply(matrix(as.numeric(z>0), nrow=n), 2, sum);
  
  ## Simulate
  source("DynLogitFS-2010.R")
  samp = 300
  burn = 0
  out <- dyn.logit.FS(y, X.dyn=X, n=n, X.stc=NULL,
                      samp=samp, burn=burn, verbose=100,
                      m.0=m0, C.0=C0,
                      mu.m0=NULL, mu.V0=NULL,
                      phi.m0=NULL, phi.V0=NULL,
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
  source("DynLogitFS-2010.R")
  samp = 500
  burn = 0
  out <- dyn.logit.FS(y, X.dyn=X, n=2, X.stc=NULL,
                      samp=samp, burn=burn, verbose=100,
                      m.0=b.m0, C.0=b.C0,
                      mu.m0=NULL, mu.V0=NULL,
                      phi.m0=NULL, phi.V0=NULL,
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
