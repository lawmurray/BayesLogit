library("BayesLogit")
source("NB-Shape.R");
source("ComputeMixture.R")
source("Benchmark-Utilities.R")

## Note: working on the log mean scale.

################################################################################

draw.GP <- function(z, KXX, KXF, KFF, v, return.all=FALSE)
{
  ## z: data
  ## X: covariate points
  ## F: forecast points
  ## v: vector of scalar variances
  
  N = ncol(KXX)
  U = chol(KXX + diag(v, N));

  ## Mean
  m.f = backsolve(U, z, transpose=TRUE);
  m.f = backsolve(U, m.f)
  m.f = t(KXF) %*% m.f;

  ## Var
  S = backsolve(U, KXF, transpose=TRUE);
  V.f = KFF - t(S) %*% S;

  f.draw = m.f + t(chol(V.f)) %*% rnorm(N);

  if (return.all) return(out = list("draw"=f.draw, "mean"=m.f, "var"=V.f))
  
  f.draw
}

draw.GP.post.mean <- function(z, IKXX, v, return.all=FALSE)
{
  N = ncol(IKXX)

  ## U = chol(Posterior Precisions)
  U = chol(IKXX + diag(1/v, N));
  
  kappa = z / v;
  m = backsolve(U, kappa, transpose=TRUE);
  m = backsolve(U, m);

  draw = m + backsolve(U, rnorm(N));
  
  draw
}

################################################################################

NB.PG.GP.gibbs <- function(y, X, F=NULL, K=NULL,
                           samp=1000, burn=500, verbose=500,
                           ups.true = NULL, w.true = NULL, d.true=NULL)
{
  ## y: response
  ## X: training predictors 
  ## F: testing predictors

  ignore.F = FALSE
  if (is.null(F)) { ignore.F = TRUE; F = matrix(0, nrow=1, ncol=ncol(X)); }
  
  X = as.matrix(X);
  F = as.matrix(F);
  y = as.matrix(y);
  
  P   = ncol(X)
  N.X = nrow(X)
  N.F = nrow(F)

  if (is.null(K)) K <- function(x,y) { ep = x - y; exp( -0.5 * (t(ep) %*% ep) ) }

  KXX = matrix(0, N.X, N.X);
  KXF = matrix(0, N.X, N.F);
  KFF = matrix(0, N.F, N.F);
  
  for (i in 1:N.X) for (j in 1:N.X) KXX[i,j] = K(X[i,], X[j,]);
  for (i in 1:N.X) for (j in 1:N.F) KXF[i,j] = K(X[i,], F[j,]);
  for (i in 1:N.F) for (j in 1:N.F) KFF[i,j] = K(F[i,], F[j,]);

  IKXX = chol2inv(chol(KXX))
  
  ## Initialize.
  ups  = rep(0.0, N.X)
  d    = 1;
  psi  = ups - log(d)
  zero = rep(0.0, N.X)

  ## Set known.
  if (!is.null(ups.true)) ups = ups.true; 
  if (!is.null(w.true))   w   = w.true;  
  if (!is.null(d.true))   d   = d.true;
  if (ignore.F)           f   = NA;

  ## Preprocess
  ymax = max(y);
  CDF = cumsum(hist(y, breaks=0:(ymax+1)-0.5, plot=FALSE)$counts)
  G   = N - CDF;
  
  out <- list(w   = matrix(0, nrow=samp, ncol=N.X),
              ups = matrix(0, nrow=samp, ncol=N.X),
              f   = matrix(0, nrow=samp, ncol=N.F),
              d   = rep(0, samp)
              )

  start.time = proc.time()
  
  ## Sample
  for ( j in 1:(samp+burn) )
  {
    if (j==burn+1) start.ess = proc.time()

    ## WARNING: JOINT DRAW.
    ## draw (d, w | beta)
    ## draw (d | beta)
    mu = exp(ups)
    d  = draw.df(d, mu, G, ymax);
    ## draw (w | d, beta)
    psi = ups - log(d);
    w = rpg.devroye(N, y+d, psi);

    ## draw beta
    kappa = 0.5 * (y-d)
    z     = kappa / w + log(d)
    ## ups   = draw.GP(z, KXX, KXX, KXX, 1/w)
    ups   = draw.GP.post.mean(z, IKXX, 1/w); 

    if (!ignore.F) f = draw.GP(z, KXX, KXF, KFF, zero)
    
    ## Record if we are past burn-in.
    if (j>burn) {
        out$w[j-burn,]    = w
        out$ups[j-burn,]  = ups
        out$d[j-burn]     = d
        out$f[j-burn,]    = f
    }

    if (j %% verbose == 0) { print(paste("NBPG GP logmean: Iteration", j)); }
  }

  end.time = proc.time()
  out$total.time = end.time - start.time
  out$ess.time   = end.time - start.ess
  
  out
} ## NB.PG.gibbs

################################################################################

NB.FS.GP.gibbs <- function(y, X, F=NULL, K=NULL,
                           samp=1000, burn=500, verbose=500,
                           ups.true = NULL, d.true=NULL, lambda.true=NULL, r.true=NULL)
{
  ## X: n by p matrix
  ## y: n by 1 vector, counts.

  ignore.F = FALSE
  if (is.null(F)) { ignore.F = TRUE; F = matrix(0, nrow=1, ncol=ncol(X)); }
  
  X = as.matrix(X);
  F = as.matrix(F);
  y = as.matrix(y);
  
  P   = ncol(X)
  N.X = nrow(X)
  N.F = nrow(F)

  if (is.null(K)) K <- function(x,y) { ep = x - y; exp( -0.5 * (t(ep) %*% ep) ) }

  KXX = matrix(0, N.X, N.X);
  KXF = matrix(0, N.X, N.F);
  KFF = matrix(0, N.F, N.F);

  for (i in 1:N.X) for (j in 1:N.X) KXX[i,j] = K(X[i,], X[j,]);
  for (i in 1:N.X) for (j in 1:N.F) KXF[i,j] = K(X[i,], F[j,]);
  for (i in 1:N.F) for (j in 1:N.F) KFF[i,j] = K(F[i,], F[j,]);

  IKXX = chol2inv(chol(KXX))
  
  ## Initialize.
  ups  = rep(0.0, N.X)
  d    = 1;
  psi  = ups - log(d)
  zero = rep(0.0, N.X)
  
  ## Set known.
  if (!is.null(ups.true))    ups = beta.true;
  if (!is.null(d.true))      d   = d.true;
  if (!is.null(lambda.true)) lambda = lambda.true;
  if (!is.null(r.true))      r   = r.true;
  if (ignore.F)              f   = NA;

  ## Preprocess
  ymax = max(y);
  CDF  = cumsum(hist(y, breaks=0:(ymax+1)-0.5, plot=FALSE)$counts)
  G    = N - CDF;
  
  cat("Finished Preprocessing.\n");
  
  out <- list(ups = matrix(nrow=samp, ncol=N.X),
              d = rep(0, samp),
              lambda = matrix(nrow=samp, ncol=N.X),
              f = matrix(nrow=samp, ncol=N.F),
              r = matrix(nrow=samp, ncol=N)
              )

  start.time = proc.time()
  
  ## Sample
  for ( j in 1:(samp+burn) )
  {
    if(j==burn+1) start.ess = proc.time()

    ## WARNING: JOINT DRAW.
    ## draw (d, lambda, r | beta)
    ## draw (d | beta)
    mu  = exp(ups)
    d = draw.df(d, mu, G, ymax);
    ## draw (lambda | d, beta)
    psi = ups - log(d);
    p = 1 / (1 + exp(-psi))
    lambda = rgamma(N, y+d, scale=p)
    ## draw (r | d, lambda, beta)
    nmix = compute.mixture(d);
    res  = psi - log(lambda)
    r    = draw.indicators.C(res, nmix);

    ## draw beta
    z   = log(lambda) + log(d) + nmix$m[r];
    ## ups = draw.GP(z, KXX, KXX, KXX, nmix$v[r])
    ups   = draw.GP.post.mean(z, IKXX, nmix$v[r]); 

    if (!ignore.F) f = draw.GP(z, KXX, KXF, KFF, zero)
    
    # Record if we are past burn-in.
    if (j>burn) {
        out$r[j-burn,]      = r
        out$ups[j-burn,]    = ups
        out$d[j-burn]       = d
        out$lambda[j-burn,] = lambda
        out$lambda[j-burn,] = f
    }
    
    if (j %% verbose == 0) { print(paste("NBFS-logmean: Iteration", j)); }
  }

  end.time = proc.time()
  out$total.time = end.time - start.time
  out$ess.time   = end.time - start.ess

  out
}

################################################################################
                                   ## TEST ##
################################################################################

if (FALSE)
{

  ## N = 1000;
  P = 2;

  ## Random
  X = matrix(runif(N*P), nrow=N, ncol=P);
  
  ## On a grid.  1/19 with len.scale 0.1 yields an invertible KXX.
  X1 = as.matrix(seq(0, 1, 1/15))
  L  = length(X1)
  M = matrix(rep(X1, L), nrow=L, ncol=L);
  X = cbind(as.numeric(M), as.numeric(t(M)));
  N = nrow(X)
  
  len.scale = 0.1
  nug = 0.0001
  ## Squared error covariance function.
  K <- function(x,y) {
    ep = (x - y) / len.scale ;
    nug * all(x==y) + exp( -0.5 * (t(ep) %*% ep) )
  }

  KXX = matrix(0, N, N);
  for (i in 1:N) for (j in 1:N) KXX[i,j] = K(X[i,], X[j,]);

  evd = eigen(KXX); root = evd$vec %*% diag(sqrt(abs(evd$values))) %*% t(evd$vec);
  ups = root %*% rnorm(N)
  ups = ups + 2.7
  z   = matrix(ups, L, L)

  mu = exp(ups)  
  d  = 4;
  nmix = compute.mixture(d)
  r = sample.int(nmix$nc, N, replace=TRUE, prob=nmix$p);
  lambda = exp(ups - log(d) - rnorm(N, nmix$m[r], sqrt(nmix$v[r])));
  ## lambda = (mu / d) * rgamma(N, d, 1);
  y = rpois(N, lambda);
  psi = ups - log(d);
  w = rpg.devroye(N, y+d, psi);

  ## library("rgl")
  ## points3d(X[,1], X[,2], y/8)

  wireframe(ups ~ X[,1] * X[,2])
  cloud(ups ~ X[,1] * X[,2])
}

if (FALSE) {
  
  ## source("NB-GaussianProcess.R")

  samp = 10000
  burn = 2000
  verbose = 100
  ntrials = 10
  bench = list();
  
  F = X

  sstat = list(ups=list());
  ess.time = rep(0, ntrials);

  ## USING PG
  for(i in 1:ntrials) {
  
    out.pg <- NB.PG.GP.gibbs(y, X, F=NULL, K=K,
                             samp=samp, burn=burn, verbose=verbose,
                             ups.true = NULL, w.true = NULL, d.true=NULL)

    ess.time[i] = out.pg$ess.time[3]
    sstat$ups[[i]] = sum.stat(out.pg$ups, out.pg$ess.time[3])

  }
  
  sstat$ups = simplify2array(sstat$ups)
  bench[["PG"]] = list("ess.time"=ess.time, "sstat"=sstat,"info"=list("ave.arate"=1))
  
  sstat = list(ups=list());
  ess.time = rep(0, ntrials);

  ## USING FS
  for(i in 1:ntrials) {
  
    out.fs <- NB.FS.GP.gibbs(y, X, F=NULL, K=K,
                             samp=samp, burn=burn, verbose=verbose,
                             ups.true = NULL, d.true=NULL, lambda.true=NULL, r.true=NULL)
    
    ess.time[i] = out.fs$ess.time[3]
    sstat$ups[[i]] = sum.stat(out.fs$ups, out.fs$ess.time[3])
    
  }
  
  sstat$ups = simplify2array(sstat$ups)
  bench[["FS"]] = list("ess.time"=ess.time, "sstat"=sstat, "info"=list("ave.arate"=1));

  setup.table(bench, "ups")
  
  plot3d(X[,1], X[,2], ups)
  plot3d(X[,1], X[,2], ups.pg.m, add=TRUE, col=2)
  plot3d(X[,1], X[,2], log(y+1e-2), add=TRUE, col=3)
  
  plot3d(X[,1], X[,2], ups)
  plot3d(X[,1], X[,2], ups.fs.m, add=TRUE, col=2)
  plot3d(X[,1], X[,2], log(y+1e-2), add=TRUE, col=3)

  cbind(ups, ups.pg.m, ups.fs.m, log(y))


  
}
