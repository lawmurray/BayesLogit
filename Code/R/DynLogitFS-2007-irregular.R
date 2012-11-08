
## Here we follow FS&F (2007)'s dynamic logit model.  There is one binary
## response at each time point.  This is in contrast to the possibility that you
## may get to see multiple responses at each time point.  Check out
## DynLogitMixture2.R for that.

## Independent AR(1)'s.  Maybe should change this.
source("Stationary.R");
source("LogitByMixture.R")

## Define normal mixture.  FS&F (2007) p. 3511.
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

FFBS <- function(y.u, phi, X, r, W, m0, C0)
{
  ## y.u : latent utility
  ## phi : vector
  ## X : design matrix
  ## r : mixture indicators
  ## W : covariance matrix of innovations of beta.
  ## m0 : prior mean on (alpha, beta)
  ## C0 : prior var on (alpha, beta)

  ## y.u = alpha X_1t + beta_t X_2t + Normal Mixture
  ## beta_t ~ AR(1).

  T = length(y.u);
  N = ncol(X);
  N.b  = length(phi);
  N.a = N - N.b;
  a.idc = 1:N.a;
  b.idc = 1:N.b+N.a;
  
  m = array(m0, dim=c(N, T+1));
  C = array(C0, dim=c(N, N, T+1));
  R = array(0., dim=c(N, N, T+1));
  a = array(0., dim=c(N, T+1));

  alpha = rep(0, N.a);
  beta = array(0, dim=c(N.b, T+1));
  
  d = c( rep(1, N.a), phi );
  D = diag(d, N);
  big.W = matrix(0, N, N); big.W[b.idc, b.idc] = W;

  ## Feed Forward
  for (i in 2:(T+1)) {
    i.y = i-1;
    n.i = length(y.u[i.y]);
    one = rep(1, n.i);

    a[,i]  = d * m[,i-1];
    R[,,i] = D %*% C[,,i-1] %*% D + big.W;

    tF.i = t(X[i.y,]);
    f.i  = tF.i %*% m[,i-1] + NM$m[r[i.y]];
    Q    = tF.i %*% R[,,i] %*% t(tF.i) + NM$v[r[i.y]];

    ## QI = solve(Q);
    ## QI  = diag(p.m, n.i) - (p.m / (1/xRx + sum(p.m))) %*% t(p.m);

    e.i = y.u[i.y] - f.i;
    ## QI1 = QI %*% one;

    A.i = R[,,i] %*% t(tF.i) / Q[1]

    ## We could simplify further.
    m[,i] = a[,i] + A.i %*% e.i;
    C[,,i] = R[,,i] - A.i %*% t(A.i) * Q[1];
    
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

dyn.logit.mix <- function(y, Xa, Xb, samp=1000, burn=100, verbose=10000,
                          a.m0, a.C0, b.m0, b.C0,
                          phi.m0, phi.V0, W.a0, W.b0,
                          alpha.known=NULL, beta.known=NULL, phi.known=NULL, W.known=NULL,
                          y.u.known=NULL, r.known=NULL)
{
  Xa = as.matrix(Xa);
  Xb = as.matrix(Xb);
  X  = cbind(Xa, Xb);
  
  T = nrow(X)
  P = ncol(X)
  M = samp

  N.a = ncol(Xa)
  N.b = ncol(Xb)
  a.idc = 1:N.a;
  b.idc = 1:N.b+N.a;

  ## Default prior parameters.
  m0 = c(a.m0, b.m0);
  C0 = matrix(0, N, N);
  C0[a.idc, a.idc] = a.C0;
  C0[b.idc, b.idc] = b.C0;

  mu = rep(0, N);
  
  out = list(
    alpha = array(0, dim=c(M, N.a)),
    beta  = array(0, dim=c(M, N.b, T+1)),
    phi  = array(0, dim=c(M, N.b)),
    W    = array(0, dim=c(M, N.b)),
    y.u   = array(0, dim=c(M, T)),
    r     = array(0, dim=c(M, T))
    )

  alpha = matrix(0.0, N.a);
  beta  = matrix(0.0, N.b, T+1);  # even odds.
  phi = rep(phi.m0, N.b);
  W   = W.b0 / W.a0;
  y.u   = matrix(0.0, N);
  r     = matrix(1.0, N);

  ## Always be careful with indicators.  You don't want to start with a very
  ## unlikely value (given the data).  It doesn't matter here because we can
  ## sample y.u without respect to r to start.  Nonetheless, seed randomly out
  ## of principle.
  r = sample.int(NM$N, N, prob=NM$w, replace=TRUE);

  ## In case we are doing testing.  You need to comment out specific sections below.
  if (!is.null(alpha.known)) alpha = alpha.known;
  if (!is.null(beta.known))  beta = beta.known;
  if (!is.null(phi.known)) phi = phi.known;
  if (!is.null(W.known))   W   = W.known;
  if (!is.null(y.u.known))  y.u  = y.u.known;
  if (!is.null(r.known))    r    = r.known;

  time.start = proc.time();
  
  for (i in 1:(samp+burn)) {

    tbeta = t(beta);
    lambda = exp( Xa %*% alpha + apply(Xb * tbeta[-1,], 1, sum) );

    ## (y.u, r | beta)  = (y.u | beta) (r | y.u, \beta)
    ## (y.u | beta, r) != (y.u | beta)
    y.u  = draw.utility(y, lambda)
    r    = draw.indicators(y.u, lambda)

    ## (alpha, beta | y, r)
    ab = FFBS(y.u, phi, X, r, diag(W, N.b), m0, C0);
    alpha = ab$alpha;
    beta  = ab$beta;

    ## AR(1) - phi, W assumed to be diagonal !!!
    ## phi = draw.phi.R(beta, mu, W, phi.m0, phi.V0, phi)
    W   = draw.W.R  (beta, mu, phi, W.a0, W.b0)

    ## cat("phi:", phi, "W:", W, "\n");
    
    if (i > burn) {
      out$beta[i-burn,,]  = beta;
      out$alpha[i-burn, ] = alpha;
      out$phi[i-burn, ]   = phi;
      out$W[i-burn, ]     = W;
      out$y.u[i-burn,]    = y.u;
      out$r[i-burn,]      = r;
    }

    if (i %% verbose == 0) cat("Iteration", i, "\n");
    
  }

  time.end = proc.time()
  out$time = time.end - time.start;

  out
}

################################################################################
                                   ## TEST ##
################################################################################

if (FALSE) {

  T = 1000;
  N.a = 1;
  N.b = 1;
  N = N.a + N.b;

  beta = array(0, dim=c(N.b, T+1));
  Xa = matrix(1, T, N.a);
  Xb = matrix(1, T, N.b);

  ## Parameters
  alpha = 0.0;
  W = 0.01;
  phi = 0.95;

  ## Prior
  a.m0 = alpha;
  a.C0 = 1;
  b.m0 = 0.0;
  b.C0 = 1;
  phi.m0 = 0.9
  phi.V0 = 0.1;
  W.a0   = 10;
  W.b0   = W.a0 * W;

  ## Synthetic
  beta[,1] = 0.0;
  for (i in 2:(T+1)) {
    beta[,i] = phi* beta[,i-1] + sqrt(W) * rnorm(1);
  }

  lambda = exp( Xa %*% alpha + apply(Xb * t(beta)[-1,], 1, sum) );

  r    = sample.int(NM$N, T, prob=NM$w, replace=TRUE);
  ep   = rnorm(T, NM$m[r], NM$s[r]);
  y.u  = log(lambda) + ep;
  y.u0 = -1 * log(rexp(T));
  y    = as.numeric(y.u > y.u0);

  ## Simulate
  out = dyn.logit.mix(y, Xa, Xb, samp=500, burn=100, verbose=100,
                      a.m0=alpha, a.C0=0.0001, b.m0=b.m0, b.C0=b.C0,
                      phi.m0=phi.m0, phi.V0=phi.V0, W.a0=W.a0, W.b0=W.b0,
                      beta.known=NULL, alpha.known=NULL, phi.known=1.0, W.known=NULL,
                      y.u.known=NULL, r.known=NULL)
  
}

if (FALSE) {

  beta.mean = apply(out$beta, c(2,3), mean);
  adj.y.u   = y.u - alpha - NM$m[r];

  y.u.upper = apply(out$y.u, 2, function(x){quantile(x, 0.95)});
  y.u.lower = apply(out$y.u, 2, function(x){quantile(x, 0.05)});
  
  ymin = min(beta.mean, beta, adj.y.u);
  ymax = max(beta.mean, beta, adj.y.u);
  
  plot(beta[1,], type="l", ylim=c(ymin,ymax));
  points(beta[1,])
  
  lines(y.u, col="gray")
  ## points(y.u, col="gray")

  lines(y.u.upper - alpha - NM$m[r], col="pink");
  lines(y.u.lower - alpha - NM$m[r], col="pink");
  
  abline(h=mean(adj.y.u), col="gray");
  lines(beta.mean[1,], col=2);
  abline(h=mean(beta.mean[1,]), col=2, lty=c(2,2));
  
  lines(out$beta[500,1,], col=3);
  
}

################################################################################
                                 ## APPENDIX ##
################################################################################

## Note: Big realization: there is a huge difference between phi = 0.95 and phi
## = 1.0 when tracking a noisy process.  Take the above example and GENERATE
## beta using phi = 0.95.  Then ESTIMATE (beta | y.u, phi) when phi = 0.95 and
## phi = 1.0.  You do a MUCH BETTER JOB tracking beta when you use phi = 1.0!

## What does this suggest?  You do a worse job estimating alpha is my guess.
## There is not stationary mean so it may not even make sense to do so.  In
## fact, I think things are misspecified in the case that Xa = 1.  In that case
## you can do a change of variables b = alpha + beta and end up with an
## identical model (when phi = 1).  So things aren't identified.  Nonetheless,
## the prior on alpha and y.u anchor things down.

## Okay.  Now try that same exercise to ESTIMATE (beta | phi) when phi = 0.95
## and phi = 1.0.  Things look crazy for phi = 1.0.  You know longer have y.u to
## anchor you down.

## But things are misspecified in that case because there should be no alpha if
## beta = 1, at least in the case that Xa = 1.  So try setting m.a0 = alpha and
## C.a0 = 1e-4, i.e. essentially fixing alpha.  In that case things look good.
## That is In the case that W = 1.0.  In the case that W = 0.01 then you have a
## very low signal-to-noise ratio.  In that case, things won't look so good.
## Basically, you are going to have a hard time "seeing" y.u (I am thinking) if
## W is small.  And really, when you thinking about it, all you get to see is
## y.u > y.u0, which is very little information.

## So I think the sampler is working properly, but you must think about the
## model you are using carefully.  Things might look horrible if you don't.
