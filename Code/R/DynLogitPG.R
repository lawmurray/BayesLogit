
## Binomial Logistic Regression using PolyaGamma augmentation.

## Independent AR(1)'s.  Maybe should change this.
source("~/RV-Project/Code/C_Examples/MyLib/Gibbs/Ind/AR1/Stationary/Stationary.R");

FFBS <- function(z, X, mu, phi, W, w, m0, C0, with.alpha=FALSE)
{
  ## When tracking (alpha, beta_t) or (\beta_t)
  ## z_t = alpha + x_t beta_t + ep_t, ep_t \sim N(0, 1/w_t)
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
## Binomial Logistic Regression.

dyn.logit.PG <- function(y, X, n=rep(1, length(y)),
                         samp=1000, burn=100, verbose=100000,
                         m.0=NULL, C.0=NULL,
                         mu.m0=NULL, mu.V0=NULL,
                         phi.m0=NULL, phi.V0=NULL,
                         W.a0=NULL, W.b0=NULL,
                         beta.true=NULL, iota.true=NULL, w.true=NULL,
                         mu.true=NULL, phi.true=NULL, W.true=NULL)
{
  ## y: the counts
  ## X: the design matrix.
  ## n: the number of trials.

  ## m.0 = prior mean for (iota,beta_0) or (beta_0).
  ## C.0 = prior var  for (iota,beta_0) or (beta_0).
  
  ## NOTE: We do note combine data. ##
  
  ## Dimension ##
  y = as.matrix(y)
  X = as.matrix(X)
  N = nrow(X);
  P = ncol(X);
  M = samp;
  
  ## Default prior parameters ##
  if (is.null(m.0))    m.0    = rep(0.0, P);
  if (is.null(C.0))    C.0    = diag(1.0, P);
  if (is.null(mu.m0))  mu.m0  = rep(0.0 , P);
  if (is.null(mu.V0))  mu.V0  = rep(0.01, P);
  if (is.null(phi.m0)) phi.m0 = rep(0.99, P);
  if (is.null(phi.V0)) phi.V0 = rep(0.01, P);
  if (is.null(W.a0))   W.a0   = rep(1.0, P);
  if (is.null(W.b0))   W.b0   = rep(1.0, P);

  ## Output data structure ##
  out = list(
    w    = array(0, dim=c(M, N)),
    iota = array(0, dim=c(M)),
    beta = array(0, dim=c(M, P, N+1)),
    mu   = array(0, dim=c(M, P)),
    phi  = array(0, dim=c(M, P)),
    W    = array(0, dim=c(M, P))
    )

  ## Initialize ##
  beta = matrix(0.0, P, N+1);  # even odds.
  mu   = matrix(0.0, P);
  phi  = matrix(0.99, P);
  W    = W.b0 / W.a0;
  om   = rep(0, N);
  with.iota = TRUE
  
  ## In case we are doing testing or we want to constrain to local level model.
  if (!is.null(beta.true)) beta = beta.true;
  if (!is.null(mu.true))   mu   = mu.true;
  if (!is.null(phi.true))  phi  = phi.true;
  if (!is.null(W.true))    W    = W.true;
  if (!is.null(w.true))    om   = w.true;
  if (!is.null(iota.true)) {iota = iota.true; with.iota=FALSE}
  
  kappa = (y-n/2)

  start.time = proc.time();
  
  ## SAMPLE ##
  
  for ( j in 1:(samp+burn) )
  {
    if (j==burn+1) start.ess = proc.time();
    
    ## draw om
    psi = iota + apply(X * t(beta)[-1,], 1, sum);
    om = rpg.devroye(N, n, psi);

    ## Draw beta;
    z = kappa / om;
    ffbs = FFBS(z, X, mu, phi, diag(W, P), om,  m.0, C.0, with.alpha=with.iota);
    iota = ffbs$alpha;
    beta = ffbs$beta;
    
    ## AR(1) - phi, W assumed to be diagonal !!!
    ## mu  = draw.mu.R(beta, phi, W, mu.m0, mu.V0) 
    ## phi = draw.phi.R(beta, mu, W, phi.m0, phi.V0, phi)
    W   = draw.W.R  (beta, mu, phi, W.a0, W.b0)

    # Record if we are past burn-in.
    if (j > burn) {
      out$w[j-burn,]      = om;
      out$iota[j-burn]    = iota;
      out$beta[j-burn,,]  = beta;
      out$mu[j-burn, ]    = mu;
      out$phi[j-burn, ]   = phi;
      out$W[j-burn, ]     = W;
    }

    if (j %% verbose == 0) { cat("Dyn Logit PG: Iteration", j, "\n"); }
  }

  end.time = proc.time();
  out$total.time = end.time - start.time;
  out$ess.time   = end.time - start.ess;
  
  out
} ## logit.gibbs.R

################################################################################
                                 ## TEST 01 ##
################################################################################

if (FALSE) {

  T = 500;
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

  psi = iota + apply(X * t(beta)[-1,], 1, sum); ## log-odds
  p   = exp(psi) / (1 + exp(psi));
  y = rbinom(N, 1, p);
  w = rpg.devroye(N, 1, psi);
  
  ## Simulate
  source("DynLogitPG.R")
  samp = 500
  burn = 0
  out <- dyn.logit.PG(y, X, samp=samp, burn=burn, verbose=100,
                      m.0=b.m0, C.0=b.C0,
                      mu.m0=NULL, mu.V0=NULL,
                      phi.m0=NULL, phi.V0=NULL,
                      W.a0=W.a0, W.b0=W.b0,
                      beta.true=NULL, iota.true=iota, w.true=NULL,
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
                      mu.m0=NULL, mu.V0=NULL,
                      phi.m0=NULL, phi.V0=NULL,
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
  out = dyn.logit.PG(y, X, n, samp=samp, burn=burn, verbose=100,
                     m.0=b.m0, C.0=b.C0,
                     mu.m0=NULL, mu.V0=NULL,
                     phi.m0=NULL, phi.V0=NULL,
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
