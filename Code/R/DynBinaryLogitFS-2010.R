## Dynamic BINARY logistic regression.

## Independent AR(1)'s.  Maybe should change this.
source("~/RV-Project/Code/C_Examples/MyLib/Gibbs/Ind/AR1/Stationary/Stationary.R");
if (!is.loaded("FSF_nmix.so")) dyn.load("~/RPackage/BayesLogit/Code/R/FSF_nmix.so");

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

dyn.logit.FS <- function(y, X,
                         samp=1000, burn=100, verbose=100,
                         m.0=NULL, C.0=NULL,
                         mu.m0=NULL, mu.V0=NULL,
                         phi.m0=NULL, phi.V0=NULL,
                         W.a0=NULL, W.b0=NULL,
                         z.true=NULL, r.true=NULL,
                         beta.true=NULL, iota.true=NULL,
                         mu.true=NULL, phi.true=NULL, W.true=NULL)
{
  ## y: the counts
  ## X: the design matrix.

  ## m.0 = prior mean for (iota,beta_0) or (beta_0).
  ## C.0 = prior var  for (iota,beta_0) or (beta_0).
  
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
    z    = array(0, dim=c(M, N)),
    r    = array(0, dim=c(M, N)),
    iota = array(0, dim=c(M)),
    beta = array(0, dim=c(M, P, N+1)),
    mu   = array(0, dim=c(M, P)),
    phi  = array(0, dim=c(M, P)),
    W    = array(0, dim=c(M, P))
    )

  ## Initialize ##
  beta = matrix(0.0, P, N+1);  # even odds
  iota = 0
  mu   = matrix(0.0, P);
  phi  = matrix(0.99, P);
  W    = W.b0 / W.a0;
  with.iota = TRUE
  nmix = NM

  ## In case we are doing testing.  May remove later.
  if (!is.null(z.true))    z = z.true;
  if (!is.null(r.true))    r = r.true;
  if (!is.null(beta.true)) beta = beta.true;
  if (!is.null(mu.true))   mu   = mu.true;
  if (!is.null(phi.true))  phi  = phi.true;
  if (!is.null(W.true))    W    = W.true;
  if (!is.null(iota.true)) {iota = iota.true; with.iota=FALSE}

  start.time = proc.time();
  
  for (j in 1:(samp+burn)) {
    if (j==burn+1) start.ess = proc.time();
    
    ## JOINT DRAW: (z, r | beta, y) = (z | beta, y) (r | z, beta, y)
    psi  = iota + apply(X * t(beta)[-1,], 1, sum)
    lambda = exp(psi);
    z    = draw.z(lambda, y);
    r    = draw.indicators.C(z, lambda, nmix)
    
    ## (iota, beta | r, z, y)
    ffbs = FFBS(z, X, mu, phi, diag(W,P), r, m.0, C.0, nmix, with.alpha=with.iota)
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
  iota = 0;
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
  psi = iota + apply(X * t(beta)[-1,], 1, sum); ## log-odds
  p   = exp(psi) / (1 + exp(psi));
  r   = sample.int(nmix$N, N, replace=TRUE, prob=nmix$w);
  ep  = rnorm(N, 0, nmix$s[r]);
  z   = psi + ep;
  y   = as.numeric(z>0);
  
  ## Simulate
  source("DynBinaryLogitFS-2010.R")
  samp = 300
  burn = 0
  out <- dyn.logit.FS(y, X,
                      samp=samp, burn=burn, verbose=100,
                      m.0=m0, C.0=C0,
                      mu.m0=NULL, mu.V0=NULL,
                      phi.m0=NULL, phi.V0=NULL,
                      W.a0=W.a0, W.b0=W.b0,
                      z.true=NULL, r.true=NULL,
                      beta.true=NULL, iota.true=0,
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

################################################################################
