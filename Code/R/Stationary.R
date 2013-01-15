
## library("tseries");

generate.data <- function(P, N, mu, phi, W)
{
  ## Set up parameters.
  mu  = array(mu, P);
  phi = array(phi, P);
  W   = array(W  , P);
  ## Create data.
  y = array(0, dim=c(P,N));
  y[,1] = rnorm(P, 0, sqrt(W / (1 - phi^2)));
  for (i in 2:N) {
    y[,i] = phi * y[,i-1] + rnorm(P, 0, sqrt(W));
  }
  y = apply(y, 2, function(x){x+mu});
  y
}

##------------------------------------------------------------------------------

a.phi.R <- function(ppsl, prev)
{
  alpha = (1-ppsl^2) / (1-prev^2);
  alpha = ifelse(ppsl >= 0 & ppsl < 1, alpha, 0);
  alpha
} # a.phi.R

draw.phi.R <- function(y, mu, W, m0, V0, phi.old, zero.mean=FALSE)
{
  ## Dim.
  P = dim(y)[1];
  N = dim(y)[2];

  ## To Return.
  phi = 0;

  ## Posterior Calculations.
  if (!zero.mean)  { y = matrix(apply(y, 2, function(x){x-mu}), P, N); }

  g.0 = rep(0, P);
  g.1 = y[,1] * y[,2];

  for (i in 2:(N-1)) {
    g.0 = g.0 + y[,i]^2;         # ~ n cov_0
    g.1 = g.1 + y[,i] * y[,i+1]; # ~ n cov_1
  }

  m1 = g.1 / g.0; # Likelihood mean
  V1 =   W / g.0; # Likelihood var

  pst.var  = V0 * V1 / (V0 + V1);
  pst.mean = (V1 * m0 + V0 * m1) / (V0 + V1);

  ## Propose and check if accept.
  phi = rnorm(P, pst.mean, sqrt(pst.var));
  phi = ifelse ( a.phi.R(phi, phi.old) > runif(P), phi, phi.old );

  phi
}

##------------------------------------------------------------------------------

draw.W.R <- function(y, mu, phi, a0, b0, zero.mean=FALSE)
{
  ## Assuming diagonal W.  phi is a vector.
  ## Dim.
  P = dim(y)[1];
  N = dim(y)[2];

  ## To Return.
  W = array(0, P);

  ## Posterior Calculations.
  if (!zero.mean) { y = matrix(apply(y, 2, function(x){x-mu}), nrow=P, ncol=N); }

  ## Sum of square residuals.
  res = matrix(y[,-1] - phi * y[,-N], P, N-1);

  S2 = (1-phi^2) * y[,1]^2;
  for(i in 1:(N-1)){
    S2 = S2 + res[,i]^2;
  }

  ## Posterior parameters.
  a1 = a0 + N + 1;
  b1 = b0 + S2;

  ## Draw -- Inverse Chi^2.
  W = 1 / rgamma(P, a1 * 0.5, rate = b1 * 0.5);

  W
}

##------------------------------------------------------------------------------

draw.mu.R <- function(y, phi, W, m0, V0)
{
  ## Dim.
  P = dim(y)[1];
  N = dim(y)[2];

  ## To Return.
  mu = array(0, P);

  ## Posterior Calculations.
  mu.hat = rowMeans(y);
  total  = (N * (1 - phi) + 1 + phi);
  wgt    = N * (1 - phi) / total;

  ## Likelihood
  m1 = wgt * mu.hat + (1 - wgt) * (y[,1] + y[,N]);
  V1 = W / ( (1 - phi) * total );

  ## Posterior
  pst.var  = V0 * V1 / (V0 + V1);
  pst.mean = (V1 * m0 + V0 * m1) / (V0 + V1);

  mu = rnorm(P, pst.mean, sqrt(pst.var));

  mu;
}

##------------------------------------------------------------------------------

ar1.gibbs <- function(y, samp=1000,
                      mu.m0=0, mu.V0=1e10, phi.m0=0.5, phi.V0=1e10, W.a0=0, W.b0=0,
                      burn=0, use.C=FALSE)
{
  if (!is.matrix(y)) { y = as.matrix(y, 1, length(y)); }

  P = dim(y)[1];
  N = dim(y)[2];

  ## To return.
  out = list(
    "mu"  = array(0, dim=c(P,samp)),
    "phi" = array(0, dim=c(P,samp)),
    "W"   = array(0, dim=c(P,samp))
    );

  ## Make sure priors are okay.
  mu.m0  = array(mu.m0 , P);
  mu.V0  = array(mu.V0 , P);
  phi.m0 = array(phi.m0, P);
  phi.V0 = array(phi.V0, P);
  W.a0   = array(W.a0  , P);
  W.b0   = array(W.b0  , P);

  ## Seed.
  mu  = array(0, P);
  phi = array(0, P);
  W   = array(0, P);

  for(i in 1:P){
    mu[i]  = mean(y[i,]);
    x      = y[i,] - mu[i];
    phi[i] = cor(x[-1], x[-N]);
    x      = x[-1] - phi[i] * x[-N];
    W[i]   = sd(x);
  }

  if (phi < 0 || phi >=1) { print("Warning: phi out of range."); }

  draw.phi = draw.phi.R;
  draw.W   = draw.W.R;
  draw.mu  = draw.mu.R;
  if (use.C) {
    draw.phi = draw.phi.C;
    draw.W   = draw.W.C;
    draw.mu  = draw.mu.C;
  }

  ## Gibbs sampling.
  for (i in 1:(samp+burn)){
    phi = draw.phi(y, mu, W, phi.m0, phi.V0, phi)
    W   = draw.W  (y, mu, phi, W.a0, W.b0)
    mu  = draw.mu (y, phi, W, mu.m0, mu.V0)
    if (i > burn) {
      out$phi[,i-burn] = phi;
      out$W  [,i-burn] = W;
      out$mu [,i-burn] = mu;
    }
  }
  out;
}

################################################################################

## "main" portion of script.

if (!exists("do.main"))  { do.main  = FALSE; }
if (!exists("gen.data")) { gen.data = FALSE; }

if (do.main) {
  ## Number of data points.
  P = 2;
  N = 1000;

  ## Generate Data.
  if (gen.data) { y  = generate.data(P, N, -2.0, 0.8, 0.2); }

  ## Number of samples
  samp = 1000;
  burn = 0;

  ## Prior
  mu.m0 = -2.0;
  mu.V0 = 10.0;
  phi.m0 = 0.8;
  phi.V0 = 2.0;
  W.a0 = 1;
  W.b0 = 1;

  ## Plot
  par(mfrow=c(4, 4));

  ## Estimate Posteriors
  gibbs = ar1.gibbs(y, samp = samp, burn = burn, W.a0 = 2, W.b0 = 1);

  for (i in 1:P) {
    ## Freq. point estimates.
    ar1   = arma(y[i,], c(1,0));
    phi.est = ar1$coef[1];
    mu.est  = ar1$coef[2] / (1 - phi.est);
    W.est   = mean(ar1$resid^2, na.rm=TRUE);

    plot(y[i,], type="l");
    plot(gibbs$mu[i,], pch=1, col="#00000030");
    abline(mean(gibbs$mu), 0, col=3);
    abline(mu.est, 0, col=4, lty=5);
    legend("topright", c("Gibbs", "Freq"), col=c(3,4), lty=c(1,1));

    plot(gibbs$phi[i,], pch=1, col="#00000030");
    abline(mean(gibbs$phi), 0, col=3);
    abline(phi.est, 0, col=4, lty=5);
    legend("topright", c("Gibbs", "Freq"), col=c(3,4), lty=c(1,1));

    plot(sqrt(gibbs$W[i,]), pch=1, col="#00000030");
    abline(sqrt(mean(gibbs$W)), 0, col=3);
    abline(sqrt(W.est), 0, col=4, lty=5);
    legend("topright", c("Gibbs", "Freq"), col=c(3,4), lty=c(1,1));
  }

  ## Estimate Posteriors
  gibbs = ar1.gibbs(y, samp = samp, burn = burn, W.a0 = 2, W.b0 = 1, use.C=TRUE);

  for (i in 1:P) {
    ## Freq. point estimates.
    ar1   = arma(y[i,], c(1,0));
    phi.est = ar1$coef[1];
    mu.est  = ar1$coef[2] / (1 - phi.est);
    W.est   = mean(ar1$resid^2, na.rm=TRUE);

    plot(y[i,], type="l");
    plot(gibbs$mu[i,], pch=1, col="#00000030");
    abline(mean(gibbs$mu[i,]), 0, col=3);
    abline(mu.est, 0, col=4, lty=5);
    legend("topright", c("Gibbs", "Freq"), col=c(3,4), lty=c(1,1));

    plot(gibbs$phi[i,], pch=1, col="#00000030");
    abline(mean(gibbs$phi[i,]), 0, col=3);
    abline(phi.est, 0, col=4, lty=5);
    legend("topright", c("Gibbs", "Freq"), col=c(3,4), lty=c(1,1));

    plot(sqrt(gibbs$W[i,]), pch=1, col="#00000030");
    abline(mean(sqrt(gibbs$W[i,])), 0, col=3);
    abline(sqrt(W.est), 0, col=4, lty=5);
    legend("topright", c("Gibbs", "Freq"), col=c(3,4), lty=c(1,1));
  }

}
