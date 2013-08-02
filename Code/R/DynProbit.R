## Copyright 2013 Nick Polson, James Scott, and Jesse Windle.

## This file is part of BayesLogit, distributed under the GNU General Public
## License version 3 or later and without ANY warranty, implied or otherwise.


## if (exists("TESTING")) {
if (!is.loaded("BayesLogit.so")) dyn.load("../C/BayesLogit.so");
source("DynLogitPG.R")
## } ## TESTING

################################################################################

dyn.probit <- function(y, X, tpred, 
                       samp=1000, burn=100, verbose=100000,
                       m.0=NULL, C.0=NULL,
                       mu.m0=NULL, mu.V0=NULL,
                       phi.m0=NULL, phi.V0=NULL,
                       W.a0=NULL, W.b0=NULL,
                       beta.true=NULL, z.true=NULL,
                       mu.true=NULL, phi.true=NULL, W.true=NULL)
{
  ## y: vector: binary response.
  ## X: matrix: design matrix
  ## tpred: vector: the vector of observation times.
  ##                We assume that tpred starts at 1 and ends at T.
  ##                tpred and X must be ``aligned.''
  
  X = as.matrix(X);

  ## SETUP ##
  
  N = nrow(X);
  P = ncol(X);
  M = samp;
  T = tpred[N];

  ## Set up tm, which maps from time back to the first instance of that time
  ## found in tpred.
  tm = rep(NA, N);
  prev = tpred[1]-1;
  for (i in 1:N) {
    if (tpred[i] != prev) tm[i] = i;
    prev = tpred[i]
  }
  tm = tm[!is.na(tm)];
  
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
    z    = array(0, dim=c(M, N)),
    beta = array(0, dim=c(M, P, T+1)),
    mu   = array(0, dim=c(M, P)),
    phi  = array(0, dim=c(M, P)),
    W    = array(0, dim=c(M, P)),
    tm   = tm
    )

  left  = rep(0, N);  left [y==0]=-Inf;
  right = rep(0, N);  right[y==1]= Inf;
  
  beta = matrix(0.0, P, T+1);  # even odds.
  mu   = matrix(0.0, P);
  phi  = matrix(0.99, P);
  W    = W.b0 / W.a0;
  z    = rep(0, N);
  
  ## In case we are doing testing or we want to constrain to local level model.
  if (!is.null(beta.true)) beta = beta.true;
  if (!is.null(mu.true))   mu   = mu.true;
  if (!is.null(phi.true))  phi  = phi.true;
  if (!is.null(W.true))    W    = W.true;
  if (!is.null(z.true))    z    = z.true;

  start.time = proc.time();
  
  ## SAMPLE ##
  
  for (i in 1:(samp+burn)) {

    psi = apply(X * t(matrix(beta[,tpred+1], nrow=P)), 1, sum);
    
    ## Draw z.
    z = rtnorm(N, psi, 1.0, left, right);

    ## Draw beta;
    ffbs = FFBS.PG(z, X, tpred, mu, phi, rep(1.0, N), diag(W, P), m.0, C.0, tm);
    beta = ffbs$beta;
    
    ## AR(1) - phi, W assumed to be diagonal !!!
    ## mu  = draw.mu.R(beta, phi, W, mu.m0, mu.V0) 
    ## phi = draw.phi.R(beta, mu, W, phi.m0, phi.V0, phi)
    W   = draw.W.R  (beta, mu, phi, W.a0, W.b0)
    
    # Record if we are past burn-in.
    if (i > burn) {
      out$z[i-burn,]      = z
      out$beta[i-burn,,]  = beta;
      out$mu[i-burn, ]    = mu;
      out$phi[i-burn, ]   = phi;
      out$W[i-burn, ]     = W;
    }

    if (i%%verbose==0) {
      cat("Dyn Probit: Iteration", i, "\n");
    }
    
  }

  end.time = proc.time();
  diff.time = end.time - start.time;
  out$diff.time = diff.time

  out
}

################################################################################
                                  ## TEST 1 ##
################################################################################

if (FALSE) {

  T = 100;
  P = 1;
  N = 2*T;

  beta = array(0, dim=c(P, T+1));
  X = matrix(1, 2*T, P);
  tpred = ceiling(seq(0.5,T,0.5));
  
  ## Parameters
  W   = 0.1;
  mu  = 0.0;
  phi = 0.95

  ## Prior
  b.m0 = 0.0;
  b.C0 = 2.0;
  phi.m0 = 0.9
  phi.V0 = 0.1;
  W.a0   = 10;
  W.b0   = W.a0 * W;

  ## Synthetic
  beta[,1] = 0.0;
  for (i in 2:(T+1)) {
    beta[,i] = phi * beta[,i-1] + sqrt(W) * rnorm(1);
  }

  psi = apply(X * t(matrix(beta[,tpred+1], nrow=P)), 1, sum);
  z   = rnorm(N, psi, 1.0);
  y   = as.numeric(z >= 0);

  out <- dyn.probit(y, X, tpred, 
                    samp=500, burn=100, verbose=100,
                    m.0=b.m0, C.0=b.C0,
                    mu.m0=NULL, mu.V0=NULL,
                    phi.m0=NULL, phi.V0=NULL,
                    W.a0=W.a0, W.b0=W.b0,
                    beta.true=NULL, z.true=NULL,
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
  
  points(seq(0.5,T,0.5), y*(ymax-ymin) + ymin, cex=0.1);
  
  lines(0:T, out$beta[300,1,], col=3);
  
}

################################################################################
                              ## TOKYO RAINFALL ##
################################################################################

## TOKYO RAINFALL ##
if (FALSE) {
  tokyo = read.csv("DataSets/tokyo-rain.csv")
  tkrain = as.numeric(na.omit(as.numeric(as.matrix(tokyo))))
  
  ## tkrain = tkrain[1:100]
  T = length(tkrain);
  df = data.frame("Time"=c(1:T,1:T), "y"=c(as.numeric(tkrain>=1), as.numeric(tkrain>=2)));
  df = df[order(df$Time),]

  N = nrow(df)
  tpred = df$Time;
  y = df$y
  X = matrix(1, N, 1);

  ## Prior
  b.m0 = -1.0;
  b.C0 = 1.0;
  W      = 0.1;
  W.a0   = 300;
  W.b0   = W.a0 * W;

  ## Simulate
  out = dyn.probit(y, X, tpred, samp=500, burn=100, verbose=100,
                     m.0=b.m0, C.0=b.C0,
                     mu.m0=NULL, mu.V0=NULL,
                     phi.m0=NULL, phi.V0=NULL,
                     W.a0=W.a0, W.b0=W.b0,
                     beta.true=NULL, z.true=NULL,
                     mu.true=0.0, phi.true=1.0, W.true=NULL)

  ess = apply(out$beta[seq(1, 500, 1), 1, ], 2, ESS);
  
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

  points(seq(0.5,T,0.5), y*(ymax-ymin) + ymin, cex=0.1);
  
  lines(0:T, out$beta[300,1,], col=3);
  
}

if (FALSE) {

  psi = out$beta[,1,];
  prob = pnorm(psi, 0.0, 1.0);
  
  prob.mean = apply(prob, 2, mean);
  prob.95   = apply(prob, 2, function(x){quantile(x, 0.95)});
  prob.05   = apply(prob, 2, function(x){quantile(x, 0.05)});
  
  ymin = min(prob.mean, prob)
  ymax = max(prob.mean, prob)

  plot(0:T, prob.mean, col=2, type="l", ylim=c(ymin,ymax));
  lines(0:T, prob.95, col="pink")
  lines(0:T, prob.05, col="pink")
  abline(h=mean(prob.mean), col=2, lty=c(2,2));

  points(seq(0.5,T,0.5), y*(ymax-ymin) + ymin, cex=0.1);
  
  ## lines(0:T, prob[50,], col=3);
  
}

################################################################################
                                ## Marketing ##
################################################################################

if (FALSE) {

  bank = read.csv("DataSets/bank.csv", sep=";")
  bank = bank[order(bank$duration),];

  y = as.numeric(bank$y=="yes")
  X = model.matrix(y ~ loan, data=bank);
  tpred = bank$duration-3; ## Simple fix for now.

  N = nrow(X);
  P = ncol(X);
  T = tpred[N]
  
  ## Prior
  b.m0 = rep(-1.0, P);
  b.C0 = diag(1.0, P);
  W      = 0.01;
  W.a0   = 300;
  W.b0   = W.a0 * W;

  ## Simulate
  out <- dyn.probit(y, X, tpred, samp=500, burn=100, verbose=100,
                    m.0=b.m0, C.0=b.C0,
                    mu.m0=NULL, mu.V0=NULL,
                    phi.m0=NULL, phi.V0=NULL,
                    W.a0=W.a0, W.b0=W.b0,
                    beta.true=NULL, z.true=NULL,
                    mu.true=0.0, phi.true=1.0, W.true=NULL)
  
}

if (FALSE) {

  beta.mean = apply(out$beta, c(2,3), mean);
  beta.95   = apply(out$beta, c(2,3), function(x){quantile(x, 0.95)});
  beta.05   = apply(out$beta, c(2,3), function(x){quantile(x, 0.05)});

  ymin = min(beta.mean, out$beta);
  ymax = max(beta.mean, out$beta);

  par(mfrow=c(2,1))

  for (j in 1:2) {
  
    plot(0:T, beta.mean[j,], col=2, type="l", ylim=c(ymin,ymax));
    lines(0:T, beta.95[j,], col="pink")
    lines(0:T, beta.05[j,], col="pink")
    abline(h=mean(beta.mean[j,]), col=2, lty=c(2,2));
    lines(0:T, out$beta[100,j,], col=3);
    
    points(jitter(tpred), y*(ymax-ymin) + ymin, cex=0.1);
    
  }
  
}
