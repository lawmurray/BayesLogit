
## Negative binomial regression using PG augmentation.
## We model the log-mean here.

## Independent AR(1)'s.  Maybe should change this.
source("~/RV-Project/Code/C_Examples/MyLib/Gibbs/Ind/AR1/Stationary/Stationary.R");

################################################################################

## When tracking (alpha, beta_t) or (\beta_t)
FFBS <- function(z, X, mu, phi, W, w, m0, C0, with.alpha=FALSE)
{
  ## When tracking (alpha, beta_t) or (\beta_t)
  ## z_t = alpha + x_t beta_t + ep_t, ep_t \sim N(0, 1/w_t)
  ## beta_t = mu + phi * (beta_t - mu) + omega_t, omega_t \sim N(0,W).
  
  ## z : vector of observations.
  ## X : design matrix
  ## phi : vector
  ## r : mixture indicators : list
  ## W : covariance matrix of innovations of beta.
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
    f.i  = tF.i %*% a[,i];

    Q    = tF.i %*% R[,,i] %*% t(tF.i) + diag(1.0 / w[i.l], 1);
    QI = solve(Q);

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

dyn.NB.PG <- function(y, X,
                      samp=1000, burn=100, verbose=100000,
                      m.0=NULL, C.0=NULL,
                      mu.m0=NULL, mu.V0=NULL,
                      phi.m0=NULL, phi.V0=NULL,
                      W.a0=NULL, W.b0=NULL,
                      d.true=NULL, w.true=NULL,
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
  with.iota = is.null(iota.true)  ## Maybe should adjust prior.
  
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
    d    = array(0, dim=c(M)),
    w    = array(0, dim=c(M, N)),
    iota = array(0, dim=c(M)),
    beta = array(0, dim=c(M, P, N+1)),
    mu   = array(0, dim=c(M, P)),
    phi  = array(0, dim=c(M, P)),
    W    = array(0, dim=c(M, P))
    )

  ## Initialize ##
  beta = matrix(0.0, P, N+1);  ## mean = 1.
  iota = 0.0
  d    = 1
  w    = rep(1, M);
  mu   = matrix(mu.m0, P);
  phi  = matrix(phi.m0, P);
  W    = rep(W.b0 / W.a0, P)
  
  ## In case we are doing testing or we want to constrain to local level model.
  if (!is.null(d.true))    d = d.true;
  if (!is.null(w.true))    w = w.true;
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

  ## Sample
  for ( j in 1:(samp+burn) )
  {
    if (j==burn+1) start.ess = proc.time();
    
    ## draw (d, w | beta) --- WARNING: JOINT DRAW.
    log.mean  = iota + apply(X * t(beta)[-1,], 1, sum)
    mu.lambda  = exp(log.mean)

    ## draw (d | beta)
    d = draw.df(d, mu.lambda, G, ymax);
    
    ## draw (w | d, beta)
    psi = log.mean - log(d);
    w = rpg.devroye(N, y+d, psi);

    ## draw beta
    kappa = 0.5 * (y-d)
    z = kappa / w + log(d);
    ffbs = FFBS(z, X, mu, phi, diag(W, P), w, m.0, C.0, with.alpha=with.iota)
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

  library("BayesLogit")
  
  T = 500;
  P = 1;

  beta = array(0, dim=c(P, T+1));
  X = matrix(1, T, P);

  ## Parameters
  W = 0.1;
  phi = 0.95;
  mu = 1.0
  d = 4
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
  beta[,1] = m.0;
  for (i in 2:(T+1)) {
    beta[,i] = mu + phi* (beta[,i-1] - mu) + sqrt(W) * rnorm(1);
  }

  log.mean  = iota + apply(X * t(beta)[-1,], 1, sum)
  psi       = log.mean - log(d)
  mu.lambda = exp(log.mean)

  lambda = rgamma(T, d, exp(psi));
  y = rpois(T, lambda);
  
  ## Simulate ##
  source("DynNBPG.R")
  samp = 500; burn=0; 
  out <- dyn.NB.PG(y, X,
                   samp=samp, burn=burn, verbose=50,
                   m.0=m.0, C.0=C.0,
                   mu.m0=NULL, mu.V0=NULL,
                   phi.m0=NULL, phi.V0=NULL,
                   W.a0=W.a0, W.b0=W.b0,
                   d.true=NULL, w.true=NULL,
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
