
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
  exp(l.out)
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
