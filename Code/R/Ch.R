
h14 = read.table("h1to4.txt")
t14 = read.table("t1to4.txt")
d14 = read.table("d1to4.txt")
e14 = 1/d14

p.left <- function(t, h)
{
  2^h * pgamma(1/t, shape=0.5, rate=0.5*h^2, lower.tail=FALSE)
}

p.right <- function(t, h)
{
  (4/pi)^h * pgamma(t, shape=h, rate=pi^2 * 0.125, lower.tail=FALSE)
}

a.coef <- function(n, x, h)
{
  ## a_n(x,h).
  ## You could precompute lgamma(h), log(2).
  d.n = (2 * n + h)
  l.out = h * log(2) - lgamma(h) + lgamma(n+h) - lgamma(n+1) + log(d.n) -
    0.5 * log(2 * pi * x^3) - 0.5 * d.n^2 / x;
  out = exp(l.out)
  out
}

c.2.coef <- function(n, x)
{
  c.n = (n+1/2) * pi
  out = (1 - 1/(x*c.n^2)) * c.n^2 * x * exp(-c.n^2 * x / 2)
  if (n > 0) out = 2 * out
  out
}

##------------------------------------------------------------------------------

rch.1 <- function(h)
{
  rate = pi^2 / 8
  idx  = floor((h-1) * 100) + 1;
  t    = t14[idx];
  
  num.trials = 0;
  total.iter = 0;

  while (TRUE)
    {
      num.trials = num.trials + 1;

      p.l = p.left (t, h);
      p.r = p.right(t, h);
      p   = p.l / (p.l + p.r);
      
      if ( runif(1) < p ) {
        ## Truncated Gamma
        X = rgamma(1, shape=h, rate=rate)
      }
      else {
        ## Truncated Inverse Gamma
        X = rtigauss(Z)
      }

      ## C = cosh(Z) * exp( -0.5 * Z^2 * X )

      ## Don't need to multiply everything by C, since it cancels in inequality.
      S = a.coef(0,X)
      Y = runif(1)*S
      n = 0

      while (TRUE)
        {
          n = n + 1
          total.iter = total.iter + 1;
          if ( n %% 2 == 1 )
            {
              S = S - a.coef(n,X)
              if ( Y<=S ) break
            }
          else
            {
              S = S + a.coef(n,X)
              if ( Y>S ) break
            }
        }

      if ( Y<=S ) break
    }

  ## 0.25 * X
  list("x"=0.25 * X, "n"=num.trials, "total.iter"=total.iter)
}


################################################################################
## 
################################################################################

if (FALSE)
{

  ## a_n is absolutely summable.  When h = 2 sum |a_n| = 0.5 as x -> infty.
  
  ## source("Ch.R")
  dx    = 0.1
  xgrid = seq(dx, 10, dx)
  y1    = xgrid
  y2    = xgrid
  n     = length(xgrid)

  for (i in 1:n) {
    y1[i] = 0
    y2[i] = 0
    for (j in 0:200) {
      y1[i] = y1[i] + (-1)^j * a.coef(j,xgrid[i],2)
      y2[i] = y2[i] + c.2.coef(j,xgrid[i])
    }
    
  }

  plot(xgrid, y1)
  points(xgrid, y2, col=2)
  
}
