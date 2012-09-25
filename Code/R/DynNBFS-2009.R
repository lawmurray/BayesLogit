
## Negative binomial regression using FSF's log-gamma parameterization.
## We model the log-mean here.

## Independent AR(1)'s.  Maybe should change this.
source("~/RV-Project/Code/C_Examples/MyLib/Gibbs/Ind/AR1/Stationary/Stationary.R");
source("~/RPackage/BayesLogit/Code/R/ComputeMixture.R")
if (!is.loaded("FSF_nmix.so")) dyn.load("~/RPackage/BayesLogit/Code/R/FSF_nmix.so");

################################################################################

## When tracking (alpha, beta_t) or (\beta_t)
FFBS <- function(z, X, mu, phi, W, r, m0, C0, NM, with.alpha=FALSE)
{
  ## When tracking (alpha, beta_t) or (\beta_t)
  ## z_t = alpha + x_t beta_t + ep_t, ep_t \sim N(m_{r_i}, v_{r_i})
  ## beta_t = mu + phi * (beta_t - mu) + omega_t, omega_t \sim N(0,W).
  
  ## z : vector of observations.
  ## X : design matrix
  ## phi : vector
  ## r : mixture indicators : list
  ## W : covariance MATRIX of innovations of beta.
  ## m0 : prior mean on (alpha, beta)
  ## C0 : prior var on (alpha, beta)

  T = length(z);
  N = ncol(X);
  N.b  = length(phi);
  N.a = N - N.b;
  if (with.alpha) a.idc = 1:N.a;
  b.idc = 1:N.b+N.a;
  
  m = array(m0, dim=c(N, T+1));
  C = array(C0, dim=c(N, N, T+1));
  R = array(0., dim=c(N, N, T+1));
  a = array(0., dim=c(N, T+1));

  beta = array(0, dim=c(N.b, T+1));
  
  d = c( rep(1, N.a), phi );
  D = diag(d, N);
  big.W = matrix(0, N, N); big.W[b.idc, b.idc] = W;
  
  ## Feed Forward
  for (i in 2:(T+1)) {
    i.l = i-1;

    a[,i]  = d * m[,i-1] + (1-d) * mu;
    R[,,i] = D %*% C[,,i-1] %*% D  + big.W;

    tF.i = t(X[i.l,])
    f.i  = tF.i %*% a[,i] + NM$m[r[i.l]];

    Q    = tF.i %*% R[,,i] %*% t(tF.i) + diag(NM$v[r[i.l]], 1);
    QI = solve(Q);

    ## nm.p = 1 / NM$v[r[[i.l]]];
    Rtx = R[,,i] %*% X[i.l,];
    ## xRtx = X[i.l,] %*% Rtx;
    ## QI  = diag(nm.p, n.i) - (nm.p / (1/xRtx + sum(nm.p))) %*% t(nm.p);

    e.i = z[i.l] - f.i;

    ## A.i = R[,,i] %*% t(tF.i) %*% QI;
    A.i = Rtx %*% QI;

    ## We could simplify further.
    m[,i] = a[,i] + A.i %*% e.i;
    ## C[,,i] = R[,,i] - A.i %*% Q %*% t(A.i);
    C[,,i] = R[,,i] - Rtx %*% QI %*% t(Rtx);
    
  }

  ## Backward Sample
  ## L = t( chol(C[,,T+1]) );
  evd = eigen(C[,,T+1]);
  Rt = evd$vectors %*% diag(sqrt(evd$values), N) %*% t(evd$vectors);
  theta = m[,T+1] + Rt %*% rnorm(N);
  alpha = ifelse(with.alpha, theta[a.idc], 0);
  beta[,T+1] = theta[b.idc];
  
  for (i in (T+1):2) {

    B = C[,,i-1] %*% (solve(R[,,i]) * d);
    theta.V = C[,,i-1] - B %*% R[,,i] %*% t(B);
    L = t( chol(theta.V[b.idc, b.idc]) );
    
    e = beta[,i] - a[b.idc,i];
    beta.m = m[b.idc,i-1] + B[b.idc, b.idc] %*% e;

    beta[,i-1] = beta.m + L %*% rnorm(N.b);
  }
  
  list("alpha"=alpha, "beta"=beta);
} ## FFBS

##------------------------------------------------------------------------------

df.llh <- function(d, mu, G, ymax)
{
  p =  1 / (1 + d / mu)
  sum( log(d+0:(ymax-1)) * G[1:ymax] ) + d * sum(log(1-p)) + sum(y * log(p)) ;
}

draw.df <- function(d.prev, mu, G, ymax)
{
  ## optim.out <- optim(d.prev, fn=df.llh, gr = NULL,
  ##                    mu, G, ymax,                              ## for minimization
  ##                    method="L-BFGS-B", lower=1, hessian=TRUE, control=list(fnscale=-1));
  
  ## mle = optim.out$par
  ## fim = -1 / optim.out$hessian

  ## cat("d.prev:", d.prev, "mle:", mle, "fim:", fim, "\n");

  d.new = d.prev
  
  ## Mixture of MH kernels.
  w = c(0.5, 0.5);
  ## k = sample.int(2, 1, prob=w);
  k = 1
  
  if (k==1) {
    ## Kernel 1: RW MH.
    rw.lower = max(d.prev - 1, 1);
    rw.upper = d.prev + 1;
    rw.grid  = rw.lower:rw.upper;
    rw.n     = length(rw.grid)
    rw.p     = rep(1/rw.n, rw.n);
    
    d.prop = sample(rw.grid, 1, prob=rw.p);
    
    ltarget = df.llh(d.prop, mu, G, ymax) - df.llh(d.prev, mu, G, ymax)
    lppsl = log(ifelse(d.prop==1, 1/2, 1/3)) - log(ifelse(d.prev==1, 1/2, 1/3));
    lratio = ltarget + lppsl
    
    if (runif(1) < exp(lratio)) d.new = d.prop
  }

  if (k==2) {
    ## Kernel 2: Ind MH.
    d.m = max(1, mle);
    d.prop = rpois(1, d.m);

    if (d.prop==0) return(d.prev);
    
    p.prop = 1 / (1 + d.prop / mu)
    p.prev = 1 / (1 + d.prev / mu)
    ltarget = df.llh(d.prop, mu, G, ymax) - df.llh(d.prev, mu, G, ymax)
    lppsl   = dpois(d.prev, d.m, log=TRUE) - dpois(d.prop, d.m, log=TRUE)
    lratio  = ltarget + lppsl

    if (runif(1) < exp(lratio)) d.new = d.prop
  }
    
  d.new
}

##------------------------------------------------------------------------------

draw.indicators <- function(res, nmix)
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

##------------------------------------------------------------------------------

dyn.NB.FS <- function(y, X,
                      samp=1000, burn=100, verbose=100000,
                      m.0=NULL, C.0=NULL,
                      mu.m0=NULL, mu.V0=NULL,
                      phi.m0=NULL, phi.V0=NULL,
                      W.a0=NULL, W.b0=NULL,
                      d.true=NULL, lambda.true=NULL, r.true=NULL,
                      beta.true=NULL, iota.true=NULL,
                      mu.true=NULL, phi.true=NULL, W.true=NULL)
{
  ## y: N by 1 vector, avg response
  ## X: N by P matrix

  ## m.0 = prior mean for (iota,beta_0) or (beta_0).
  ## C.0 = prior var  for (iota,beta_0) or (beta_0).
  
  ## Dimension ##

  y = as.matrix(y)
  X = as.matrix(X)
  N = nrow(X);
  P = ncol(X);
  M = samp
  with.iota = is.null(iota.true) ## Maybe should adjust prior.
  
  ## Default prior parameters.
  if (is.null(m.0))    m.0    = rep(0.0, P);
  if (is.null(C.0))    C.0    = diag(1.0, P);
  if (is.null(mu.m0))  mu.m0  = rep(0.0 , P);
  if (is.null(mu.V0))  mu.V0  = rep(0.01, P);
  if (is.null(phi.m0)) phi.m0 = rep(0.99, P);
  if (is.null(phi.V0)) phi.V0 = rep(0.01, P);
  if (is.null(W.a0))   W.a0   = rep(1.0, P);
  if (is.null(W.b0))   W.b0   = rep(1.0, P);

  ## Output data structure.
  out = list(
    iota = array(0, dim=c(M)),
    beta = array(0, dim=c(M, P, N+1)),
    lambda = array(0, dim=c(M, N)),
    r    = array(0, dim=c(M, N)),
    mu   = array(0, dim=c(M, P)),
    phi  = array(0, dim=c(M, P)),
    W    = array(0, dim=c(M, P))
    )

  ## Initialize ##
  beta = matrix(0.0, P, N+1);  ## mean = 1.
  iota = 0.0
  d    = 1
  mu   = matrix(mu.m0, P);
  phi  = matrix(phi.m0, P);
  W    = rep(W.b0 / W.a0, P)
  
  ## In case we are doing testing or we want to constrain to local level model.
  if (!is.null(d.true))    d    = d.true;
  if (!is.null(r.true))    r    = r.true;
  if (!is.null(lambda.true)) lambda = lambda.true;
  if (!is.null(iota.true)) {iota = iota.true; with.iota = FALSE}
  if (!is.null(beta.true)) beta = beta.true;
  if (!is.null(mu.true))   mu   = mu.true;
  if (!is.null(phi.true))  phi  = phi.true;
  if (!is.null(W.true))    W    = W.true;

  ## Preprocess ## 
  ymax = max(y);
  F = cumsum(hist(y, breaks=0:(ymax+1)-0.5, plot=FALSE)$counts)
  G = N - F;
  
  start.time = proc.time();

  ## SAMPLE ##
  for ( j in 1:(samp+burn) )
  {
    if (j==burn+1) start.ess = proc.time();

    ## draw (d, lambda, r | beta) --- WARNING: JOINT DRAW.
    log.mean  = iota + apply(X * t(beta)[-1,], 1, sum)
    mu.lambda = exp(log.mean)
    
    ## draw (d | beta)
    d = draw.df(d, mu.lambda, G, ymax);

    ## draw (lambda | d, beta)
    psi = log.mean - log(d);
    p = 1 / (1 + exp(-psi))
    lambda = rgamma(N, y+d, scale=p)

    ## draw (r | d, lambda, beta)
    nmix = compute.mixture(d);
    res  = psi - log(lambda)
    r    = draw.indicators.C(res, nmix);

    ## draw beta
    nmix$m = -1 * nmix$m;     ## We want log Gamma not - log Gamma.
    z = log(lambda) + log(d); ## So that we model the log mean.
    ffbs = FFBS(z, X, mu, phi, diag(W,P), r, m.0, C.0, nmix, with.alpha=with.iota)
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
  N = 1;

  beta = array(0, dim=c(N, T+1));
  X = matrix(1, T, N);

  ## Parameters
  W = 0.1;
  phi = 0.95;
  mu = 1.0
  d = 4
  nmix = compute.mixture(d)
  iota = 0

  ## Prior
  m.0     = mu;
  C.0     = 4.0;
  mu.m0  = 0.0;
  mu.V0  = 4.0;
  phi.m0 = 0.9;
  phi.V0 = 0.1;
  W.a0   = 10;
  W.b0   = 1.0;

  ## Synthetic
  beta[,1] = m0;
  for (i in 2:(T+1)) {
    beta[,i] = mu + phi* (beta[,i-1] - mu) + sqrt(W) * rnorm(1);
  }

  log.mean  = iota + apply(X * t(beta)[-1,], 1, sum)
  psi       = log.mean - log(d)
  mu.lambda = exp(log.mean)

  r = sample.int(nmix$nc, T, replace = TRUE, nmix$p);
  ep = rnorm(T, 1*nmix$m[r], sqrt(nmix$v[r]));
  lambda = (exp(psi - ep));
  y = rpois(T, lambda);
  
  ## Simulate ##
  source("DynNBFS-2009.R")
  samp = 500; burn=0; 
  out <- dyn.NB.FS(y, X,
                   samp=samp, burn=burn, verbose=50,
                   m.0=m.0, C.0=C.0,
                   mu.m0=mu.m0, mu.V0=mu.V0,
                   phi.m0=phi.m0, phi.V0=phi.V0,
                   W.a0=W.a0, W.b0=W.b0,
                   d.true=NULL, lambda.true=NULL, r.true=NULL, 
                   beta.true=NULL, iota.true=iota,
                   mu.true=mu, phi.true=1.0, W.true=NULL)
}

if (FALSE) {

  beta.m = array(0, dim=c(N, T+1));
  for (i in 1:N) {
    beta.m[i,] = apply(out$beta[,i,], 2, mean);
  }

  ymin = min(beta.m, beta);
  ymax = max(beta.m, beta);

  plot(beta[1,], type="l", ylim=c(0, ymax))
  lines(beta.m[1,], col=2)
  points(log(y), col="gray")
  
}

################################################################################
                                 ## APPENDIX ##
################################################################################

## Regarding the FFBS below.  This was adapted from <DynLogitPG.R>.  It is
## essentially identical.  I changed it so that it could accomodate the normal
## mixture, but the PG case is just commented out.  It is easy to go between the
## two.

## The original appiclation was binary logistic, but it turns out the FFBS is
## essentially the same.  It runs slow because I have adapted things so that you
## can have missing data.  I did this to try and accomodate the advertising
## dataset, which is irregularly spaced in the ``time'' variable, which I
## believe was duration of the phone call or soemthing like that.

## You can use a local level model by simple setting mu to 0 and phi to 1.
## Otherwise you can use a stationary AR(1) process.

## Old FFBS
## ## FFBS.dyn
## FFBS.dyn <- function(y.u, X, tpred, mu, phi, W, m0, C0, tm, Omega, r=NULL)
## {
##   ## y: vector: y^u
##   ## X: matrix: Design Matrix
##   ## mu: vector: mean, if it exists.
##   ## phi: vector:  autoregression parameter
##   ## r: matrix: mixture index.
##   ## W: vector: variance of omega_t, i.e. omega_t \sim N(0, diag(W))
##   ## m0: vector: mean of beta_0.
##   ## C0: matrix: var  of beta_0.
##   ## tm: maps time to index.

##   K = length(tm);
##   T = tpred[tm[K]];
##   P = ncol(X);
##   N = nrow(X);

##   m = array(m0, dim=c(P, T+1));
##   C = array(C0, dim=c(P, P, T+1));
##   R = array(0., dim=c(P, P, T+1));
##   a = array(0., dim=c(P, T+1));
##   UR = array(0., dim=c(P, P, T+1));

##   beta = array(0, dim=c(P, T+1));

##   Phi = diag(phi, P);
##   idx = 1;
  
##   ## FEED FORWARD ##
##   for (i in 1:T) {
##     i.h = i+1 # Index of beta (and all other variables that start at 0)
##     a[,i.h]  = (1 - phi) * mu + phi * m[,i.h-1];
##     R[,,i.h] = Phi %*% (C[,,i.h-1] * phi) + W;
##     UR[,,i.h] = chol(R[,,i.h]);

##     ## Check if we have an observation.
##     if (tpred[tm[idx]]==i) {
##       ## Indicies and pick correct parts of y, X, and r.
##       n.i = ifelse(idx==K, N+1 - tm[idx], tm[idx+1] - tm[idx]);
##       idc = tm[idx]:(tm[idx] + n.i - 1);
##       y.i = y.u[idc];
##       X.i = matrix(X[idc,], nrow=n.i);
      
##       ## Mean and Variance of nu -- CAN CHANGE DEPENDING PG OR FS.
##       ## omega.i = Omega[idc];
##       ## b.i = rep(0, length(omega.i));
##       ## V.i = 1.0 / omega.i
##       r.i = r[idc];
##       b.i = -1 * NM$m[r.i];
##       V.i = NM$v[r.i];

##       ## Forecast
##       f.i = X.i %*% a[,i.h] + b.i;
##       rho.i = X.i %*% R[,,i.h];
##       e.i = y.i - f.i;
      
##       ## Joint to Conditional
##       ## if (n.i > 1.5 * P) {
##         ## tXdivV = t(X.i) / V.i;
##         ## M  = chol2inv(UR[,,i.h]) + tXdivV %*% X.i;
##         ## L  = t(chol(M));
##         ## xi = backsolve(t(L), tXdivV)
##         ## QInv = diag(1/V.i, n.i) - t(xi) %*% xi;
##         ## A.i  = t(rho.i) %*% QInv;
##         ## m[,i.h]  = a[,i.h] + A.i %*% e.i;
##         ## C[,,i.h] = R[,,i.h] - A.i %*% rho.i;
##       ## }
##       ## else {
##         Q.i = rho.i %*% t(X.i) + diag(V.i, n.i);
##         L  = t(chol(Q.i));
##         xi = forwardsolve(L, rho.i);
##         tA.i  = backsolve(t(L), xi);
##         m[,i.h]  = a[,i.h] + t(tA.i) %*% e.i;
##         C[,,i.h] = R[,,i.h] - t(xi) %*% xi;
##       ## }

##       idx = idx + 1;
##     }
##     else {

##       ## cat("Warning: missing data at time", i, ".\n");
##       m[,i.h]  = a[,i.h]
##       C[,,i.h] = R[,,i.h]
      
##     }

##   }

##   ## BACKWARDS SAMPLE ##
##   i.h = T+1;
##   L = t(chol(C[,,i.h]));
##   beta[,i.h] = m[,i.h] + L %*% rnorm(P);

##   for (i in (T+1):2) {
##     i.y = i-1;
##     rho.i = Phi %*% C[,,i-1];
##     xi = forwardsolve(t(UR[,,i]), rho.i);
##     tB.i = backsolve(UR[,,i], xi);
##     e.i  = beta[,i] - a[,i];
##     ell  = m[,i-1] + t(tB.i) %*% e.i;
##     U    = C[,,i-1] - t(xi) %*% xi;

##     L = t(chol(U));
##     beta[,i-1] = ell + L %*% rnorm(P);
    
##   }

##   out = list("beta"=beta, "m"=m)
  
##   out
## } ## FFBS.dyn
