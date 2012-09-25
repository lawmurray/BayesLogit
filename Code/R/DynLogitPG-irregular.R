
## Independent AR(1)'s.  Maybe should change this.
source("~/RV-Project/Code/C_Examples/MyLib/Gibbs/Ind/AR1/Stationary/Stationary.R");

FFBS.PG <- function(y.u, X, tpred, mu, phi, Omega, W, m0, C0, tm)
{
  ## y: vector: y^u
  ## X: matrix: Design Matrix
  ## mu: vector: mean, if it exists.
  ## phi: vector:  autoregression parameter
  ## r: matrix: mixture index.
  ## W: vector: variance of omega_t, i.e. omega_t \sim N(0, diag(W))
  ## m0: vector: mean of beta_0.
  ## C0: matrix: var  of beta_0.
  ## tm: maps time to index.

  K = length(tm);
  T = tpred[tm[K]];
  P = ncol(X);
  N = nrow(X);

  m = array(m0, dim=c(P, T+1));
  C = array(C0, dim=c(P, P, T+1));
  R = array(0., dim=c(P, P, T+1));
  a = array(0., dim=c(P, T+1));
  UR = array(0., dim=c(P, P, T+1));

  beta = array(0, dim=c(P, T+1));

  Phi = diag(phi, P);
  idx = 1;
  
  ## FEED FORWARD ##
  for (i in 1:T) {
    i.h = i+1 # Index of beta (and all other variables that start at 0)
    a[,i.h]  = (1 - phi) * mu + phi * m[,i.h-1];
    R[,,i.h] = Phi %*% (C[,,i.h-1] * phi) + W;
    UR[,,i.h] = chol(R[,,i.h]);

    ## Check if we have an observation.
    if (tpred[tm[idx]]==i) {
      ## Indicies and pick correct parts of y, X, and r.
      n.i = ifelse(idx==K, N+1 - tm[idx], tm[idx+1] - tm[idx]);
      idc = tm[idx]:(tm[idx] + n.i - 1);
      y.i = y.u[idc];
      X.i = matrix(X[idc,], nrow=n.i);
      
      ## Mean and Variance of nu -- CAN CHANGE DEPENDING PG OR FS.
      omega.i = Omega[idc];
      b.i = rep(0, length(omega.i));
      V.i = 1.0 / omega.i

      ## Forecast
      f.i = X.i %*% a[,i.h] + b.i;
      rho.i = X.i %*% R[,,i.h];
      e.i = y.i - f.i;
      
      ## Joint to Conditional
      ## if (n.i > 1.5 * P) {
        ## tXdivV = t(X.i) / V.i;
        ## M  = chol2inv(UR[,,i.h]) + tXdivV %*% X.i;
        ## L  = t(chol(M));
        ## xi = backsolve(t(L), tXdivV)
        ## QInv = diag(1/V.i, n.i) - t(xi) %*% xi;
        ## A.i  = t(rho.i) %*% QInv;
        ## m[,i.h]  = a[,i.h] + A.i %*% e.i;
        ## C[,,i.h] = R[,,i.h] - A.i %*% rho.i;
      ## }
      ## else {
        Q.i = rho.i %*% t(X.i) + diag(V.i, n.i);
        L  = t(chol(Q.i));
        xi = forwardsolve(L, rho.i);
        tA.i  = backsolve(t(L), xi);
        m[,i.h]  = a[,i.h] + t(tA.i) %*% e.i;
        C[,,i.h] = R[,,i.h] - t(xi) %*% xi;
      ## }

      idx = idx + 1;
    }
    else {

      ## cat("Warning: missing data at time", i, ".\n");
      m[,i.h]  = a[,i.h]
      C[,,i.h] = R[,,i.h]
      
    }

  }

  ## BACKWARDS SAMPLE ##
  i.h = T+1;
  L = t(chol(C[,,i.h]));
  beta[,i.h] = m[,i.h] + L %*% rnorm(P);

  for (i in (T+1):2) {
    i.y = i-1;
    rho.i = Phi %*% C[,,i-1];
    xi = forwardsolve(t(UR[,,i]), rho.i);
    tB.i = backsolve(UR[,,i], xi);
    e.i  = beta[,i] - a[,i];
    ell  = m[,i-1] + t(tB.i) %*% e.i;
    U    = C[,,i-1] - t(xi) %*% xi;

    L = t(chol(U));
    beta[,i-1] = ell + L %*% rnorm(P);
    
  }

  out = list("beta"=beta, "m"=m)
  
  out
}

dyn.logit.PG <- function(y, X, tpred, n=rep(1, length(y)),
                         samp=1000, burn=100, verbose=100000,
                         m.0=NULL, C.0=NULL,
                         mu.m0=NULL, mu.V0=NULL,
                         phi.m0=NULL, phi.V0=NULL,
                         W.a0=NULL, W.b0=NULL,
                         beta.true=NULL,
                         mu.true=NULL, phi.true=NULL, W.true=NULL)
{
  ## X: N by P matrix
  ## y: N by 1 vector, avg response
  ## n: N by 1 vector, # of obs at distinct x
  ## tpred: vector: the vector of observation times.
  ##                We assume that tpred starts at 1 and ends at T.

  ## PREPROCESS ##
  
  X = as.matrix(X)
  N = nrow(X);
  P = ncol(X);
  X = cbind(X, tpred);

  y.prior = 0.0;
  x.prior = rep(0, P+1);
  n.prior = 0.0;
  
  ## Combine data.
  ## new.data = logit.combine(y, X, n, y.prior, x.prior, n.prior);
  new.data = list("y"=y, "X"=X, "n"=n)
  y = new.data$y;
  X = as.matrix(new.data$X[,1:P]);
  n = new.data$n;
  n.prior = 0.0;
  tpred = as.matrix(new.data$X[,P+1]);

  ## X = as.matrix(X);

  ## SETUP ##

  N = nrow(X)
  p = ncol(X)
  M = samp
  T = tpred[N]

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
    w    = array(0, dim=c(M, N)),
    beta = array(0, dim=c(M, P, T+1)),
    mu   = array(0, dim=c(M, P)),
    phi  = array(0, dim=c(M, P)),
    W    = array(0, dim=c(M, P)),
    tm   = tm
    )

  beta = matrix(0.0, P, T+1);  # even odds.
  mu   = matrix(0.0, P);
  phi  = matrix(0.99, P);
  W    = W.b0 / W.a0;
  om   = rep(0, N);
  
  ## In case we are doing testing or we want to constrain to local level model.
  if (!is.null(beta.true)) beta = beta.true;
  if (!is.null(mu.true))   mu   = mu.true;
  if (!is.null(phi.true))  phi  = phi.true;
  if (!is.null(W.true))    W    = W.true;
  
  alpha = (y-1/2)*n
  kappa = colSums(X*alpha)

  start.time = proc.time();
  
  ## SAMPLE ##
  
  for ( j in 1:(samp+burn) )
  {

    ## draw om
    psi = apply(X * t(matrix(beta[,tpred+1], nrow=P)), 1, sum);
    om = rpg.devroye(N, n, psi);

    ## Draw beta;
    z = alpha / om;
    ffbs = FFBS.PG(z, X, tpred, mu, phi, om, diag(W, P), m.0, C.0, tm);
    beta = ffbs$beta;
    
    ## AR(1) - phi, W assumed to be diagonal !!!
    ## mu  = draw.mu.R(beta, phi, W, mu.m0, mu.V0) 
    ## phi = draw.phi.R(beta, mu, W, phi.m0, phi.V0, phi)
    W   = draw.W.R  (beta, mu, phi, W.a0, W.b0)
    
    # Record if we are past burn-in.
    if (j > burn) {
      out$w[j-burn,]   = om
      out$beta[j-burn,,]  = beta;
      out$mu[j-burn, ]    = mu;
      out$phi[j-burn, ]   = phi;
      out$W[j-burn, ]     = W;
    }

    if (j %% verbose == 0) { cat("Dyn Logit PG: Iteration", j, "\n"); }
  }

  end.time = proc.time();
  diff.time = end.time - start.time;
  out$diff.time = diff.time
  
  ## Add new data to output.
  out$"y" = y;
  out$"X" = X;
  out$"n" = n;
  
  out
} ## logit.gibbs.R

## SYNTHETIC DATA ##
if (FALSE) {

  T = 100;
  P = 1;

  beta = array(0, dim=c(P, T+1));
  X = matrix(1, T, P);
  tpred = 1:T;

  N = nrow(X);

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

  psi = apply(t(X) * beta[,tpred+1], 2, sum); # log-odds
  p   = exp(psi) / (1 + exp(psi));
  y = rbinom(N, 1, p);
  
  ## Simulate
  out = dyn.logit.PG(y, X, tpred, samp=600, burn=0, verbose=100,
                     m.0=b.m0, C.0=b.C0,
                     mu.m0=NULL, mu.V0=NULL,
                     phi.m0=NULL, phi.V0=NULL,
                     W.a0=W.a0, W.b0=W.b0,
                     beta.true=NULL,
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
  
  points(1:T, y*(ymax-ymin) + ymin, cex=0.1);
  
  lines(0:T, out$beta[300,1,], col=3);
  
}

################################################################################
                                  ## TEST 2 ##
################################################################################

## SYNTHETIC DATA ##
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
  p   = exp(psi) / (1 + exp(psi));
  y = rbinom(N, 1, p);

  ## Simulate
  out = dyn.logit.PG(y, X, tpred, samp=600, burn=0, verbose=100,
                     m.0=b.m0, C.0=b.C0,
                     mu.m0=NULL, mu.V0=NULL,
                     phi.m0=NULL, phi.V0=NULL,
                     W.a0=W.a0, W.b0=W.b0,
                     beta.true=NULL,
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
  out = dyn.logit.PG(y, X, tpred, samp=500, burn=100, verbose=100,
                     m.0=b.m0, C.0=b.C0,
                     mu.m0=NULL, mu.V0=NULL,
                     phi.m0=NULL, phi.V0=NULL,
                     W.a0=W.a0, W.b0=W.b0,
                     beta.true=NULL,
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

  lambda = exp(out$beta[,1,]);
  prob = lambda / (1 + lambda);
  
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
  out = dyn.logit.PG(y, X, tpred, samp=100, burn=0, verbose=100,
                     m.0=b.m0, C.0=b.C0,
                     mu.m0=NULL, mu.V0=NULL,
                     phi.m0=NULL, phi.V0=NULL,
                     W.a0=W.a0, W.b0=W.b0,
                     beta.true=NULL,
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

if (FALSE) {

  psi = apply(out$beta[,,tpred+1], 1, function(z){apply(t(X) * z, 2, sum)})
  psi = t(psi)
  prob = exp(psi) / (1 + exp(psi));
  
  prob.mean = apply(prob, 2, mean);
  prob.95   = apply(prob, 2, function(x){quantile(x, 0.95)});
  prob.05   = apply(prob, 2, function(x){quantile(x, 0.05)});
  
  ymin = min(prob.mean, prob)
  ymax = max(prob.mean, prob)

  par(mfrow=c(1,1));
  
  plot(0:T, prob.mean, col=2, type="l", ylim=c(ymin,ymax));
  lines(0:T, prob.95, col="pink")
  lines(0:T, prob.05, col="pink")
  abline(h=mean(prob.mean), col=2, lty=c(2,2));

  points(jitter(tpred), y*(ymax-ymin) + ymin, cex=0.1);
  
}
