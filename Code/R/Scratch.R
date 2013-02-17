
digauss <- function(x, mu, lambda)
{
  out = abs(lambda / (2 * pi * x^3))^0.5 * exp(-lambda * (x-mu)^2 / (2 * mu^2 * x));
  out
}

dens.igamma <- function(x, a, b)
{
  out = a * log(b) - lgamma(a) - (a+1) * log(x) - b / x;
  out = exp(out)
  out
}

## Approximate density of J^*(h,z).
C.h <- function(X, z=0, h=1, N=100, m=1)
{
  ## X = xgrid
  ## N = number of terms to sum
  ## m = number of terms to return in proposal
  yy = X
  ppsl = X
  L = length(yy)
  for(i in 1:L){
    x = X[i]
    n = 0:N
    a = (-1)^n * gamma(n+h) / gamma(n+1) * (2*n+h) / (2 * pi * x^3)^0.5;
    e = exp( - (2*n+h)^2 / (2*x) )
    c = 2^h / gamma(h) * exp(-0.5 * x * z^2) / cosh(z/2)^h;
    s = c * sum(a*e);
    ppsl[i] = c * sum(a[1:m]*e[1:m]);
    yy[i] = s
  }

  out = list("y"=yy, "ppsl"=ppsl);

  out
}

## A guess at the upper bound of x for which a_n(x) > a_{n+1}(x) for all n.
ubound <- function(h, N=100)
{
  n = 0:N;
  
  top = 8 * n + 4 * h + 4;
  bot = lgamma(n+1+h) - lgamma(n+2) + log(2*n+h+2) - (lgamma(n+h) - lgamma(n+1) + log(2*n+h));

  out = 0.5 * top / bot;

  min(out[-1]);
}

################################################################################

## To test how things change when degrees of freedom increases.

if (FALSE) {

  ## The x-values (density scale).
  ## The h-values: parameter for J^*(h,z).
  x = seq(0.1, 6, 0.01);
  ## x = seq(0.1, 30, 1);
  h = seq(1,4,0.01);

  ## z in J^*(h,z).
  z = 0.0
  
  n = length(x);
  p = length(h);

  y = matrix(nrow=n, ncol=p)
  w = matrix(nrow=n, ncol=p)
  u = 1:p;
  w.2 = w;
  w.3 = w;
  gam = w
  ig  = w
  
  ## Calculate.
  for (i in 1:p) {
    ## N: number of terms to sum.
    ## m: number of terms in proposal.
    out = C.h(x, z=z, h=h[i], N=160, m=1);
    u[i] = ubound(h[i], N=100);
    y[,i] = out$y;
    w[,i] = out$ppsl;
    out = C.h(x, z=z, h=h[i], N=2, m=2);
    w.2[,i] = out$ppsl;
    out = C.h(x, z=z, h=h[i], N=2, m=3);
    #w.3[,i] = out$ppsl;
    ## IGamma near 0, Gamma tail?
    ## gam[,i] = dgamma(x, h[i], pi^2 / 8 + z^2 / 2)
    ## gam[,i] = (pi/2) * x^(h[i] - 1) * exp(- pi^2 / 8 * x) / gamma(h[i]);
    ## GOOD FOR (1,2) and kind of (0,1)
    gam[,i] = (pi / 2)^{h[i]} * x^(h[i] - 1) * exp(- pi^2 / 8 * x) / gamma(h[i]);
    ig [,i] = dens.igamma(x, 0.5, h[i]^2 / 2);
  }

  for (i in seq(1,p,10)) {
    par(mfrow=c(1,2))
    ## plot(x, y[,i], type="l")
    all = c(w[,i], y[,i], w.2[,i])
    ymax = max(all);
    ymin = min(all);
    ylim = c(ymin, ymax);
    
    plot(x, y[,i], type="l", col=1, ylim=ylim, main=paste("h:", h[i], "z:", z))
    lines(x, w[,i], col=2)
    #lines(x, w.2[,i], col=3);
    #lines(x, w.3[,i], col=4);
    lines(x, gam[,i], col=7);
    #lines(x, ig [,i], col=8);
    abline(v=u[i], col=5)
    legend("topright", lty=c(1,1,1,1), col=c(1,2,3,5), legend=c("density", "1-term", "2-term", "upper bound"));

    ## Ratio
    plot(x, y[,i] / gam[,i], type="l", main=paste("y/gam, h:", h[i]), col=7, ylim=c(0,2))
    lines(x, y[,i] / w[,i], type="l", col=2);
    abline(h=1.0, col=1)
    readline("<ENTER>")
  }

  ratio.gam = y / gam
  ratio.S0  = y / w

  dist.gam = abs(ratio.gam-1);
  dist.S0  = abs(ratio.S0-1)
  v.min = pmin(dist.gam, dist.S0)
  gam.win = v.min==dist.gam

  max.idx = apply(gam.win, 2, which.max)
  x.switch = x[max.idx]
  gam.max = max.idx; for(i in 1:p) gam.max[i] = dist.gam[max.idx[i],i]
  S0.max  = max.idx; for(i in 1:p) S0.max[i] = dist.S0[max.idx[i],i]

  plot(h, S0.max)
  plot(h, x.switch)

  lm.1 = lm(x.switch[200:301] ~ h[200:301])
  
}
