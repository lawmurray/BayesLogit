t = 0.64
p = 0.57810262346829443
q = 0.422599094

a <- function(n,x)
{
  if ( x>t )
    pi * (n+0.5) * exp( -(n+0.5)^2*pi^2*x/2 )
  else
    (2/pi/x)^1.5 * pi * (n+0.5) * exp( -2*(n+0.5)^2/x )
}

rJstar <- function()
{

  while (TRUE)
  {
    if ( runif(1) < p/(p+q) )
      X = t + 8 * rexp(1) / pi^2
    else
      while (TRUE)
      {
        E = rexp(2)
        if ( E[1]^2 <= 2*E[2]/t )
        {
          X = t/(1+t*E[1])^2
          break
        }
      }

    S = a(0,X)
    Y = runif(1)*S
    n = 0

    while (TRUE)
    {
      n = n + 1
      if ( n %% 2 == 1 )
      {
        S = S - a(n,X)
        if ( Y<S ) break
      }
      else
      {
        S = S + a(n,X)
        if ( Y>S ) break
      }
    }

    if ( Y<S ) break
  }

  X
}

## x2 = c()
## for ( i in 1:10000 ) x2[i] = rJstar()

## x2 = x2/4.0
