## FFBS for generalized dynamic linear models.

if(!is.loaded("FFBS.so")) dyn.load("FFBS.so");

FFBS.R <- function(z, X, mu, phi, W, V, m0, C0)
{
  ## When tracking (alpha, beta_t).  It may be the case that there is no alpha.
  ## z_t = x_t (alpha_t, beta_t) + ep_t, ep_t \sim N(0, V_t).
  ## beta_t = mu + phi * (beta_t - mu) + omega_t, omega_t \sim N(0,W).
  ## alpha_t = alpha_{t-1}
  
  ## z : vector of observations (T)
  ## X : design matrix (T x N)
  ## mu : mu (K)
  ## phi : vector (K)
  ## W : covariance MATRIX of innovations of beta (K x K)
  ## V: time varying variances (T)
  ## m0 : prior mean on (beta_0, alpha_0) (N)
  ## C0 : prior var on (beta_0, alpha_0) (N x N).

  T = length(z);
  N.b = length(mu);
  N   = ncol(X);
  N.a = N - N.b;
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

    tF.i = t(X[i.l,])
    f.i  = tF.i %*% a[,i];

    Q    = tF.i %*% R[,,i] %*% t(tF.i) + diag(V[i.l], 1);
    QI   = solve(Q);

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
  L = t( chol(C[,,T+1]) );
  ## evd = eigen(C[,,T+1]);
  ## Rt = evd$vectors %*% diag(sqrt(evd$values), N) %*% t(evd$vectors);
  theta = m[,T+1] + L %*% rnorm(N);
  alpha = ifelse(with.alpha, theta[1:N.a], 0);
  beta[,T+1] = theta[b.idc];
  
  for (i in (T+1):2) {

    ## B = C[,,i-1] %*% (solve(R[,,i]) * d);
    ## theta.V = C[,,i-1] - B %*% R[,,i] %*% t(B);
    ## L = t( chol(theta.V[b.idc, b.idc]) );
    ## e = beta[,i] - a[b.idc,i];
    ## beta.m = m[b.idc,i-1] + B[b.idc, b.idc] %*% e;
    ## beta[,i-1] = beta.m + L %*% rnorm(N.b);

    A.bs = C[b.idc, b.idc, i-1] %*% (solve(R[b.idc, b.idc, i]) * phi);
    V.bs = C[b.idc, b.idc, i-1] - A.bs %*% R[b.idc, b.idc, i] %*% t(A.bs);
    m.bs = m[b.idc, i-1] + A.bs %*% (beta[,i] - a[b.idc,i]);

    L = t(chol(V.bs));
    beta[,i-1] = m.bs + L %*% rnorm(N.b);
  }
  
  list("alpha"=alpha, "beta"=beta);
} ## FFBS

FFBS.C <- function(z, X, mu, phi, W, V, m0, C0)
{
  T = length(z)

  T = length(z);
  N.b = length(mu);
  N   = ncol(X);
  N.a = N - N.b;

  ## Check
  ok = rep(6);
  ok[1] = length(mu)==length(phi)
  ok[2] = length(mu)==nrow(W) 
  ok[3] = nrow(W)==ncol(W) 
  ok[4] = T==nrow(X) 
  ok[5] = N==length(m0) 
  ok[6] = N==nrow(C0) && N==ncol(C0);

  if (!prod(ok)) {
    printf("FFBS.C: problem", ok, "\n");
    return(0)
  }    

  alpha = rep(0, max(N.a, 1));
  beta  = array(0, dim=c(N.b, T+1));
  
  OUT <- .C("ffbs", alpha, beta,
            z, X,
            mu, phi, W,
            V, m0, C0,
            as.integer(N.b), as.integer(N), as.integer(T),
            PACKAGE="FFBS");

  ## void ffbs(double *alpha_, double *beta_,
  ##         double *z_, double *X_, double *mu_, double *phi_, double *W_, double *V_,
  ##         double *m0_, double *C0_, int N_b, int N, int T)
     
  out = list("alpha"=OUT[[1]], "beta"=OUT[[2]]);
}
