
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

pg.a.coef <- function(n, x, h, z=0)
{
  cosh(z/2)^h * exp(-0.5 * z^2 * x) * 4 * a.coef(n, 4 * x, h)
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
                            ## PLOTTING DENSITIES ##
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

##------------------------------------------------------------------------------
                           ## PLOTTING PG DENSITY ##
##------------------------------------------------------------------------------

if (FALSE)
{
  ## source("Ch.R")
  dx    = 0.01
  xgrid = seq(dx, 3, dx)
  y1    = xgrid
  y2    = xgrid
  y3    = xgrid
  y4    = xgrid
  y5    = xgrid
  n     = length(xgrid)

  for (i in 1:n) {
    y1[i] = 0
    y2[i] = 0
    y3[i] = 0
    y4[i] = 0
    y5[i] = 0
    for (j in 0:200) {
      y1[i] = y1[i] + (-1)^j * pg.a.coef(j,xgrid[i],1)
      y2[i] = y2[i] + (-1)^j * pg.a.coef(j,xgrid[i],2)
      y3[i] = y3[i] + (-1)^j * pg.a.coef(j,xgrid[i],3)
      y4[i] = y4[i] + (-1)^j * pg.a.coef(j,xgrid[i],1, 2.0)
      y5[i] = y5[i] + (-1)^j * pg.a.coef(j,xgrid[i],1, 4.0)
    }
    
  }

  ## png(filename="pg-dens.png", width=800, height=400)
  
  par(mfrow=c(1,2))
  ymax = max(c(y1,y2,y3))
  plot(xgrid, y1, type="l", ylim=c(0,ymax), main="Density of PG(b,0)", xlab="x", ylab="f(x|b,0)")
  lines(xgrid, y2, type="l", col=2, lty=2)
  lines(xgrid, y3, type="l", col=4, lty=4)

  legend("topright", legend=c("b=1", "b=2", "b=3"), col=c(1,2,4), lty=c(1,2,4))

  ## hist(rpg.devroye(10000, 1, 0), add=TRUE, prob=TRUE, breaks=100,
  ##      col=rgb(0, 0, 0, 16, maxColorValue=255))
  ## hist(rpg.devroye(10000, 2, 0), add=TRUE, prob=TRUE, breaks=100,
  ##      col=rgb(255, 0, 0, 16, maxColorValue=255))
  ## hist(rpg.devroye(10000, 3, 0), add=TRUE, prob=TRUE, breaks=100,
  ##      col=rgb(0, 0, 255, 16, maxColorValue=255))

  ymax = max(c(y1,y4,y5))
  plot(xgrid, y1, type="l", ylim=c(0,ymax), xlim=c(0,1), main="Density of PG(1,c)", xlab="x", ylab="f(x|1,c)")
  lines(xgrid, y4, type="l", col=2, lty=2)
  lines(xgrid, y5, type="l", col=4, lty=4)

  legend("topright", legend=c("c=0", "c=2", "c=4"), col=c(1,2,4), lty=c(1,2,4))

  ## hist(rpg.devroye(10000, 1, 0), add=TRUE, prob=TRUE, breaks=100,
  ##      col=rgb(0, 0, 0, 16, maxColorValue=255))
  ## hist(rpg.devroye(10000, 1, 2), add=TRUE, prob=TRUE, breaks=100,
  ##      col=rgb(255, 0, 0, 16, maxColorValue=255))
  ## hist(rpg.devroye(10000, 1, 4), add=TRUE, prob=TRUE, breaks=100,
  ##      col=rgb(0, 0, 255, 16, maxColorValue=255))

  dev.off()
  
}
