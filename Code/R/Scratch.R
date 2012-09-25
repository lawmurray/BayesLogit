
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
  x = seq(0.1, 10, 0.05);
  h = seq(1,20,1);

  ## z in J^*(h,z).
  z = 0.5
  
  n = length(x);
  p = length(h);

  y = matrix(nrow=n, ncol=p)
  w = matrix(nrow=n, ncol=p)
  u = 1:p;
  w.2 = w;

  ## Calculate.
  for (i in 1:p) {
    ## N: number of terms to sum.
    ## m: number of terms in proposal.
    out = C.h(x, z=z, h=h[i], N=100, m=1);
    u[i] = ubound(h[i], N=100);
    y[,i] = out$y;
    w[,i] = out$ppsl;
    out = C.h(x, z=z, h=h[i], N=2, m=2);
    w.2[,i] = out$ppsl;
  }

  for (i in 1:p) {
    ## plot(x, y[,i], type="l")
    all = c(w[,i], y[,i], w.2[,i])
    ymax = max(all);
    ymin = min(all);
    ylim = c(ymin, ymax);
    
    plot(x, y[,i], type="l", col=1, ylim=ylim, main=paste("h:", h[i], "z:", z))
    lines(x, w[,i], col=2)
    lines(x, w.2[,i], col=3);
    abline(v=u[i], col=4)
    legend("topright", lty=c(1,1,1,1), col=c(1,2,3,4), legend=c("density", "1-term", "2-term", "upper bound"));
    readline("<ENTER>")
  }
  
}
