library("rootSolve")

## CUBS: conjugate updating backward sampling
## Use a MH step to to correct.

log.odds.mv <- function(rs)
{
  F1 = digamma(rs[1]) - digamma(rs[2])
  F2 = trigamma(rs[1]) + trigamma(rs[2])
  out = c("F1"=F1, "F2"=F2)
  out
}

rs.solve <- function(rs, fq)
{
  F1 = digamma(rs[1]) - digamma(rs[2]) - fq[1]
  F2 = trigamma(rs[1]) + trigamma(rs[2]) - fq[2]
  out = c("F1"=F1, "F2"=F2)
  out
}

cubs.draw <- function(y, X, n, mu, phi, W, m0, C0)
{
  ## When tracking (alpha, beta_t).  It may be the case that there is no alpha.
  ## z_t = x_t (alpha_t, beta_t) + ep_t, ep_t \sim N(0, V_t).
  ## beta_t = mu + phi * (beta_t - mu) + omega_t, omega_t \sim N(0,W).
  ## alpha_t = alpha_{t-1}
  
  ## y : vector of observations (T)
  ## X : design matrix (T x N)
  ## mu : mu (K)
  ## phi : vector (K)
  ## W : covariance MATRIX of innovations of beta (K x K)
  ## m0 : prior mean on (beta_0, alpha_0) (N)
  ## C0 : prior var on (beta_0, alpha_0) (N x N).

  T = length(y);
  N.b = length(mu);
  N   = ncol(X);
  N.a = N - N.b;
  a.idc = 1:N.a
  b.idc = (N.a+1):N;
  with.alpha = N.a > 0;
  
  m = array(m0, dim=c(N, T+1));
  C = array(C0, dim=c(N, N, T+1));
  R = array(0., dim=c(N, N, T+1));
  a = array(0., dim=c(N, T+1));

  beta = array(0, dim=c(N.b, T+1));
  
  d = c( rep(1, N.a), phi );
  D = diag(d, N);
  mu = c(rep(0, N.a), mu);
  big.W = matrix(0, N, N); big.W[b.idc, b.idc] = W;
  
  ## Filter Forward
  for (i in 2:(T+1)) {
    i.l = i-1;

    a[,i]  = d * m[,i-1] + (1-d) * mu;
    R[,,i] = D %*% C[,,i-1] %*% D  + big.W;

    x.i = t(X[i.l,])
    f.i = x.i %*% a[,i];
    q.i = ( x.i %*% R[,,i] %*% X[i.l,] )[1];
    
    rho.i = R[,,i] %*% X[i.l,];
    A.i = rho.i / q.i;

    ## Conjugate update
    rs.i = multiroot(rs.solve, start=c(1,1), fq=c(f.i, q.i));
    rstar.i = rs.i$root[1] + y[i.l];
    sstar.i = n[i.l] - y[i.l] + rs.i$root[2];
    fqstar.i = log.odds.mv(c(rstar.i, sstar.i));

    ## cat("f,q:", f.i, q.i, "root:", rs.i$root, "f.root:", rs.i$f.root, "\n");
    
    m[,i]  = a[,i] + A.i * (fqstar.i[1] - f.i);
    C[,,i] = R[,,i] + rho.i %*% t(rho.i) * ( (fqstar.i[2] / q.i - 1) / q.i );

  }

  ## Backwards sample
  L = t( chol(C[,,T+1]) );
  ## evd = eigen(C[,,T+1]);
  ## Rt = evd$vectors %*% diag(sqrt(evd$values), N) %*% t(evd$vectors);
  theta = m[,T+1] + L %*% rnorm(N);
  alpha = ifelse(with.alpha, theta[a.idc], 0);
  beta[,T+1] = theta[b.idc];
  
  for (i in (T+1):2) {

    A.bs = C[b.idc, b.idc, i-1] %*% (solve(R[b.idc, b.idc, i]) * phi);
    V.bs = C[b.idc, b.idc, i-1] - A.bs %*% R[b.idc, b.idc, i] %*% t(A.bs);
    m.bs = m[b.idc, i-1] + A.bs %*% (beta[,i] - a[b.idc,i]);

    L = t(chol(V.bs));
    beta[,i-1] = m.bs + L %*% rnorm(N.b);
  }
  
  out = list("alpha"=alpha, "beta"=beta, "m"=m, "C"=C);
  out
}

cubs.dens <- function(y, X, n, beta, alpha=NULL, m0, C0, log=TRUE)
{

  ba = c(beta, alpha);
  psi = X * 
  
  
}

################################################################################
                                   ## TEST ##
################################################################################

if (FALSE) {

  T = 100;
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
  y = rbinom(T, 1, p);
  w = rpg.devroye(T, 1, psi);
  n = rep(1, T)

  m0 = mu
  C0 = diag(b.C0, P);
  M = 1000
  verbose = 100
  out = list("beta"=array(0, dim=c(M, P, T+1)))
  for(i in 1:M) {
    if (i %% verbose == 0) cat("CUBS: iteration", i, "\n");
    draw = cubs.draw(y, X, n, mu, phi, W, m0, C0);
    out$beta[i,,] = draw$beta
  }
  
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
