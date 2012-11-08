
## Independent componenent AR1.
source("Stationary.R");

## Define normal mixture.  FS&F p. 3511.
normal.mixture = list(
  w = c(0.00397, 0.0396, 0.168, 0.147, 0.125, 0.101, 0.104, 0.116, 0.107, 0.088),
  m = c(5.09, 3.29, 1.82, 1.24, 0.764, 0.391, 0.0431, -0.306, -0.673, -1.06),
  v = c(4.50, 2.02, 1.10, 0.422, 0.198, 0.107, 0.0778, 0.0766, 0.0947, 0.146)
  )
  normal.mixture$s = sqrt(normal.mixture$v)
  normal.mixture$N = length(normal.mixture$w)

ev.samp = -log(rexp(10000));
normal.mixture$mar.mean = mean(ev.samp);
normal.mixture$mar.var  = var(ev.samp);

NM = normal.mixture;

FFBS <- function(ell, phi, X, r, W, m0, C0)
{
  ## ell = - log(tau) : list
  ## phi : vector
  ## X : design matrix
  ## r : mixture indicators : list
  ## W : covariance matrix of innovations of beta.
  ## m0 : prior mean on (alpha, beta)
  ## C0 : prior var on (alpha, beta)

  T = length(ell);
  N = ncol(X);
  N.b  = length(phi);
  N.a = N - N.b;
  a.idc = 1:N.a;
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
    n.i = length(ell[[i.l]]);
    one = rep(1, n.i);

    a[,i]  = d * m[,i-1];
    R[,,i] = D %*% C[,,i-1] %*% D + big.W;

    tF.i = one %*% t(X[i.l,]);
    f.i  = tF.i %*% m[,i-1] + NM$m[r[[i.l]]];

    Q    = tF.i %*% R[,,i] %*% t(tF.i) + diag(NM$v[r[[i.l]]], n.i);
    QI = solve(Q);

    ## nm.p = 1 / NM$v[r[[i.l]]];
    Rtx = R[,,i] %*% X[i.l,];
    ## xRtx = X[i.l,] %*% Rtx;
    ## QI  = diag(nm.p, n.i) - (nm.p / (1/xRtx + sum(nm.p))) %*% t(nm.p);

    e.i = ell[[i.l]] - f.i;

    QI1 = QI %*% one;
    ## QI1 = nm.p - (nm.p / (1/xRtx + sum(nm.p))) * (t(nm.p) %*% one);

    ## A.i = R[,,i] %*% t(tF.i) %*% QI;
    A.i = Rtx %*% t(QI1);

    ## We could simplify further.
    m[,i] = a[,i] + A.i %*% e.i;
    ## C[,,i] = R[,,i] - A.i %*% Q %*% t(A.i);
    C[,,i] = R[,,i] - Rtx %*% t(Rtx) * (t(one) %*% QI1)[1];
    
  }

  ## Backward Sample
  ## L = t( chol(C[,,T+1]) );
  evd = eigen(C[,,T+1]);
  Rt = evd$vectors %*% diag(sqrt(evd$values)) %*% t(evd$vectors);
  theta = m[,T+1] + Rt %*% rnorm(N);
  alpha = theta[a.idc];
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
}

draw.tau <- function(y, lambda)
{
  T = length(y);
  
  tau = list();
  
  for (i in 1:T) {
    if (y[i] > 0) {
      U = runif(y[i]);
      R = c(0, sort(U));
      tau[[i]] = c(diff(R), 1 - R[y[i]+1] + rexp(1, lambda[i]))
    }
    else {
      tau[[i]] = 1 + rexp(1, lambda[i]);
    }
  }

  ## I am encountering a problem whereby tau==0 causes problems.  So I am doing
  ## the following.
  zero.idx = tau[[i]] <= 1e-20;
  num.zero = sum(zero.idx);
  if (num.zero > 0) {
    cat("Found", num.zero, "numbers <= 1e-20.\n");
    tau[[i]][zero.idx] = 1e-20;
  }

  tau
}

draw.r <- function(ell, lambda)
{
  T = length(ell);
  
  r = list();

  log.wds = log(NM$w) - log(NM$s);

  for (i in 1:T) {
    ## Check for ell==inf
    inf.idx = ell[[i]]==Inf;
    if (sum(inf.idx) > 0) {
      cat(i, "ell = inf\n");
    }
    
    ## Safer to work on log scale.  Columns correspond to outcome index i!
    res = ell[[i]] - log(lambda[i]);    
    log.post = -0.5 * outer(-1*NM$m, res, "+")^2 / NM$v + log.wds;    
    unnrm.post = exp(log.post);

    r[[i]] = apply(unnrm.post, 2, function(x){sample.int(NM$N, 1, prob=x)})

    ## r[[i]] = rep(0, length(ell[[i]]));
    ## for (j in 1:length(ell[[i]])) {
    ##   cat(i, j, res[j], "w:", unnrm.post[,j], "\n\n");
    ##   r[[i]][j] = sample.int(NM$N, 1, prob=unnrm.post[,j]);
    ## }
  }
  
  r
}

dyn.poisson <- function(y, Xa, Xb, samp=1000, burn=100, verbose=10000,
                        a.m0, a.C0, b.m0, b.C0,
                        phi.m0, phi.V0, W.a0, W.b0,
                        beta.known=NULL, alpha.known=NULL, phi.known=NULL, W.known=NULL,
                        r.known=NULL, tau.known=NULL)
{
  Xa = as.matrix(Xa);
  Xb = as.matrix(Xb);
  X = cbind(Xa, Xb);
  
  T = length(y)
  N = ncol(X)
  M = samp

  N.a = ncol(Xa)
  N.b = ncol(Xb)
  a.idc = 1:N.a;
  b.idc = 1:N.b+N.a;

  m0 = c(a.m0, b.m0);
  C0 = matrix(0, N, N);
  C0[a.idc, a.idc] = a.C0;
  C0[b.idc, b.idc] = b.C0;

  mu = rep(0, N);

  out = list(
    beta = array(0, dim=c(N.b, T+1, M)),
    alpha = array(0, dim=c(N.a, M)),
    phi  = array(0, dim=c(N.b, M)),
    W    = array(0, dim=c(N.b, M))
    ## r    = list(),
    ## tau  = list()
    )

  lambda = rep(mean(y), T);
  tau = draw.tau(y, lambda);
  ell = lapply(tau, function(x){min(-1 * log(x), 40)});
  r   = draw.r(ell, lambda)
  phi = rep(phi.m0, N.b);
  W   = W.b0 / W.a0;

  ## Set known values.  Comment out corresponding draws below.
  if (!is.null(beta.known))  beta = cbind(beta.known);
  if (!is.null(alpha.known)) alpha = alpha.known;
  if (!is.null(phi.known)) phi = phi.known;
  if (!is.null(W.known))   W   = W.known;
  if (!is.null(r.known))   r = r.known;
  if (!is.null(tau.known)) tau = tau.known;
  
  ## MCMC
  for (i in 1:(samp+burn)) {

    ## FFBS
    ab = FFBS(ell, phi, X, r, diag(W, N.b), m0, C0);
    alpha = ab$alpha;
    beta  = ab$beta;

    tbeta = t(beta);
    lambda = exp( Xa %*% alpha + apply(Xb * tbeta[-1,], 1, sum) );
    
    ## Param - Must change to make W not diagonal.
    ## phi = draw.phi.R(beta, mu, W, phi.m0, phi.V0, phi)
    W   = draw.W.R  (beta, mu, phi, W.a0, W.b0)
    
    ## Joint tau then indicators.

    tau = draw.tau(y, lambda);
    ell = lapply(tau, function(x){min(-1 * log(x), 40)});
    print(sum(ell>=40));
    
    r = draw.r(ell, lambda)
      
    if (i > burn) {
      out$beta[,,i-burn] = beta;
      out$alpha[,i-burn] = alpha;
      out$phi[,i-burn]   = phi;
      out$W[,i-burn]     = W;
      ## out$tau[[i-burn]]  = tau;
      ## out$r[[i-burn]]    = r;
    }

    if (i %% verbose == 0) cat("Iteration", i, "\n");
    
  }

  out
}

################################################################################
                                ## SYNTHETIC ##
################################################################################

generate.tau.r.y.t <- function(lambda)
{
  r   = c();
  tau = c();
  ep  = c();
  total = 0;

  j = 0;
  while (total < 1) {
    j = j + 1;
    r.j = sample.int(NM$N, 1, prob=NM$w);
    ep.j = rnorm(1, NM$m[r.j], NM$s[r.j]);
    ell.j = log(lambda) + ep.j;
    tau.j = exp(-1.0 * ell.j);
    r   = c(r, r.j);
    tau = c(tau, tau.j);
    ep  = c(ep.j);
    total = total + tau.j;
  }

  list("tau"=tau, "r"=r, "y"=j-1, "ep"=ep);
}

################################################################################
                                   ## MAIN ##
################################################################################

if (FALSE) {

  T = 1000;
  N.a = 1;
  N.b = 1;

  beta = array(0, dim=c(N.b, T+1));
  Xa = matrix(1, T, N.a);
  Xb = matrix(1, T, N.b);

  ## Parameters
  alpha = 1;
  W = 0.01;
  phi = 0.95;

  ## Prior
  a.m0 = alpha;
  a.C0 = 4;
  b.m0 = 0.0;
  b.C0 = 4;
  phi.m0 = 0.9
  phi.V0 = 0.1;
  W.a0   = 10;
  W.b0   = 0.1;

  ## Synthetic
  beta[,1] = 0.0;
  for (i in 2:(T+1)) {
    beta[,i] = phi* beta[,i-1] + sqrt(W) * rnorm(1);
  }

  lambda = exp( Xa %*% alpha + apply(Xb * t(beta)[-1,], 1, sum) );

  y = rep(0, T);
  tau = list();
  r = list();
  ep = list();

  for (i in 1:T) {
    x = generate.tau.r.y.t(lambda[i]);
    y[i]     = x$y;
    tau[[i]] = x$tau;
    r[[i]]   = x$r;
    ep[[i]] = x$ep;
  }

  ## Simulate
  out = dyn.poisson(y, Xa, Xb, samp=500, burn=0, verbose=1,
                    a.m0=a.m0, a.C0=a.C0, b.m0=b.m0, b.C0=b.C0,
                    phi.m0=phi.m0, phi.V0=phi.V0, W.a0=W.a0, W.b0=W.b0,
                    beta.known=NULL, alpha.known=NULL, phi.known=NULL, W.known=NULL,
                    r.known=NULL, tau.known=NULL)
  
}

if (FALSE) {

  tau.post = array(0, dim=c(samp, y[1]+1));
  for (i in 1:samp) {
    tau.post[i,] = out$tau[[i]][[1]];
  }
  
}

if (FALSE) {

  beta.mean = apply(out$beta, c(1,2), mean);
  
  adj.ell.1 = rep(0, T+1);
  ell.1 = rep(0, T+1);
  ep.1  = rep(0, T+1);
  for (i in 1:T) {
    ell.1[i+1] = -1 * log(tau[[1]][1]);
    ep.1[i+1]  = ep[[i]][1];
  }
  adj.ell.1 = ell.1 - ep.1 - alpha;
  
  ymin = min(beta.mean, beta, adj.ell.1);
  ymax = max(beta.mean, beta, adj.ell.1);
  
  plot(beta[1,], type="l", ylim=c(ymin,ymax));
  points(beta[1,])
  
  lines(adj.ell.1, col="gray")
  points(adj.ell.1, col="gray")
  
  abline(h=mean(adj.ell.1), col="gray");
  lines(beta.mean[1,], col=2);

  abline(h=mean(beta.mean[1,]), col=2, lty=c(2,2));
  lines(out$beta[1,,100], col=3);
  
}

if (FALSE) {

  flu = read.csv("DataSets/google-flu.txt");
  na.cols = is.na(flu[1,]);
  flu.clean = flu[,!na.cols]
  tx = flu.clean$Texas
  T = length(tx);
  
  Xa = rep(1.0, T);
  Xb = rep(1.0, T);

  Xa = as.matrix(Xa);
  Xb = as.matrix(Xb);
  
  ## Prior
  a.m0 = log(mean(tx)) - NM$mar.mean;
  a.C0 = 1e-4;
  b.m0 = 0.0;
  b.C0 = 1.0;
  phi.m0 = 0.95
  phi.V0 = 0.10;
  W.a0   = 10;
  W.b0   = 0.1;  
  
  ## Simulate
  out = dyn.poisson(tx, Xa, Xb, samp=100, burn=0, verbose=1,
                    a.m0=a.m0, a.C0=a.C0, b.m0=b.m0, b.C0=b.C0,
                    phi.m0=phi.m0, phi.V0=phi.V0, W.a0=W.a0, W.b0=W.b0,
                    beta.known=NULL, alpha.known=NULL, phi.known=1.0, W.known=NULL,
                    r.known=NULL, tau.known=NULL)

  beta.mean  = apply(out$beta, c(1,2), mean);
  alpha.mean = apply(out$alpha, 1, mean);
  tbeta = t(beta.mean);
  log.lambda = Xa %*% alpha.mean + apply(Xb * tbeta[-1,], 1, sum);

  plot(log(tx) - NM$mar.mean, ylim=c(0, 20));
  lines(log.lambda, col=2);
    
}

################################################################################
                                 ## APPENDIX ##
################################################################################

## When count data is large, you DO NOT want to save tau and r!

################################################################################
################################################################################
################################################################################

## FFBS <- function(ell, phi, X, r, W, m0, C0)
## {
##   ## ell = - log(tau) : list
##   ## phi : vector
##   ## X : design matrix
##   ## r : mixture indicators : list
##   ## W : covariance matrix of innovations of beta.
##   ## m0 : prior mean on (alpha, beta)
##   ## C0 : prior var on (alpha, beta)

##   T = length(ell);
##   N = ncol(X);
##   N.b  = length(phi);
##   N.a = N - N.b;
##   a.idc = 1:N.a;
##   b.idc = 1:N.b+N.a;
  
##   m = array(m0, dim=c(N, T+1));
##   C = array(C0, dim=c(N, N, T+1));
##   R = array(0., dim=c(N, N, T+1));
##   a = array(0., dim=c(N, T+1));

##   beta = array(0, dim=c(N.b, T+1));
  
##   d = c( rep(1, N.a), phi );
##   big.W = matrix(0, N, N); big.W[b.idc, b.idc] = W;

##   ## Feed Forward
##   for (i in 2:(T+1)) {
##     i.l = i-1;
##     n.i = length(ell[[i.l]]);
##     one = rep(1, n.i);

##     a[,i]  = d * m[,i-1];
##     R[,,i] = (C[,,i-1] * d) %*% diag(d, N) + big.W;

##     cat("R:\n");
##     print(R[,,i]);
    
##     f.i = X[i.l,] %*% m[,i-1] + NM$m[r[[i.l]]];
    
##     p.m = as.vector(1.0 / NM$v[r[[i.l]]]);

##     Rx  = R[,,i] %*% X[i.l,];
##     xRx = as.numeric(X[i.l,] %*% Rx);
##     temp = (p.m / (1/xRx - sum(p.m))) %*% t(p.m);
    
##     Q = (one * xRx) %*% t(one) + diag(NM$v[r[[i.l]]], n.i)

##     cat("Q:\n");
##     print(Q)
    
##     QI = solve(Q);
##     ## QI  = diag(p.m, n.i) - (p.m / (1/xRx + sum(p.m))) %*% t(p.m);

##     cat("QI:\n");
##     print(QI);
    
##     e.i = ell[[i.l]] - f.i;
##     QI1 = QI %*% one;

##     A.i = Rx %*% QI;

##     ## We could simplify further.
##     m[,i] = a[,i] + A.i %*% e.i;
##     C[,,i] = R[,,i] - (Rx / (one %*% QI1)[1]) %*% t(Rx);
    
##   }

##   print(R[,,T+1]);
##   print(C[,,T+1]);

##   ## Backward Sample
##   ## L = t( chol(C[,,T+1]) );
##   evd = eigen(C[,,T+1]);
##   Rt = evd$vectors %*% diag(sqrt(evd$values)) %*% t(evd$vectors);
##   theta = m[,T+1] + Rt %*% rnorm(N);
##   alpha = theta[a.idc];
##   beta[,T+1] = theta[b.idc];
  
##   for (i in (T+1):2) {

##     print(R[,,i]);
    
##     B = C[,,i-1] %*% (solve(R[,,i]) * d);

##     print(B);
    
##     theta.V = C[,,i-1] - B %*% R[,,i] %*% t(B);

##     print(theta.V);
    
##     L = t( chol(theta.V[b.idc, b.idc]) );
    
##     e = beta[,i] - a[b.idc,i];
##     beta.m = m[b.idc,i-1] + B[b.idc, b.idc] %*% e;

##     beta[,i-1] = beta.m + L %*% rnorm(N.b);
##   }

##   list("alpha"=alpha, "beta"=beta);
## }
