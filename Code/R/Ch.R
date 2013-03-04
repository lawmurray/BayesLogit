
h14 = scan("h1to4.txt")
t14 = scan("t1to4.txt")
d14 = scan("d1to4.txt")
e14 = 1/d14

################################################################################
                                ##    J^*    ##
################################################################################

p.ch.left <- function(t, h)
{
  2^h * pgamma(1/t, shape=0.5, rate=0.5*h^2, lower.tail=FALSE)
}

p.ch.right <- function(t, h)
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
                          ## SAMPLE TRUNCATED GAMMA ##
##------------------------------------------------------------------------------

## A RELATIVELY EASY METHOD TO IMPLEMENT
rltgamma.dagpunar.1 <- function(shape=1, rate=1, trnc=1)
{
  ## y ~ Ga(shape, rate, trnc)
  ## x = y/t
  ## x ~ Ga(shape, rate t, 1)
  a = shape
  b = rate * trnc

  if (shape <  1) return(NA)
  if (shape == 1) return(rexp(1) / rate + trnc);
  
  d1 = b-a
  d3 = a-1
  c0 = 0.5 * (d1 + sqrt(d1^2 + 4 * b)) / b

  x = 0
  accept = FALSE
  
  while (!accept) {
    x  = b + rexp(1) / c0
    u  = runif(1)

    l.rho = d3 * log(x) - x * (1-c0);
    l.M   = d3 * log(d3 / (1-c0)) - d3

    accept = log(u) <= (l.rho - l.M)
  }

  x = x / b
  y = trnc * x
  y
}

rltgamma.dagpunar <- function(num=1, shape=1, rate=1, trnc=1)
{
  shape = array(shape, dim=num)
  rate  = array(rate , dim=num)
  trnc = array(trnc, dim=num)

  y = rep(0, num)
  for (i in 1:num) y[i] = rltgamma.dagpunar.1(shape[i], rate[i], trnc[i]);

  y
}

rrtinvchi2.1 <- function(h, trnc)
{
  R = trnc / h^2
  E = rexp(2)
  while ( (E[1]^2) > (2 * E[2] / R)) {
    ## cat("E", E[1], E[2], E[1]^2, 2*E[2] / R, "\n")
    E = rexp(2)
  }
  ## cat("E", E[1], E[2], "\n")
  ## W^2 = (1 + R*E[1])^2 / R is left truncated chi^2(1) I_{(R,\infty)}.
  ## So X is right truncated inverse chi^2(1).
  X = R / (1 + R*E[1])^2
  X = h^2 * X
  X
}

rrtinvchi2 <- function(num, h, trnc)
{
  out = rep(0, num)
  for (i in 1:num)
    out[i] = rrtinvchi2.1(h, trnc)
  out
}

##------------------------------------------------------------------------------

a.coef.1 <- function(n,x, trnc)
{
  if ( x>trnc )
    pi * (n+0.5) * exp( -(n+0.5)^2*pi^2*x/2 )
  else
    (2/pi/x)^1.5 * pi * (n+0.5) * exp( -2*(n+0.5)^2/x )
}

g.tilde <- function(x, h, trnc)
{
  if (x > trnc)
    (0.5 * pi)^h * x^(h-1) / gamma(h) * exp(-pi^2 / 8 * x)
  else
    2^h * h * (2 * pi * x^3)^(-0.5) * exp(-0.5 * h^2 / x)
}

rch.1 <- function(h)
{
  if (h < 1) return(NA);
  
  rate  = pi^2 / 8
  idx   = floor((h-1) * 100) + 1;
  trnc  = t14[idx];
  
  num.trials = 0;
  total.iter = 0;

  p.l = p.ch.left (trnc, h);
  p.r = p.ch.right(trnc, h);
  p   = p.r / (p.l + p.r);

  max.inner = 100
      
  while (TRUE)
    {
      num.trials = num.trials + 1;

      if ( runif(1) < p ) {
        ## Left truncated gamma
        X = rltgamma.dagpunar.1(shape=h, rate=rate, trnc=trnc)
      }
      else {
        ## Right truncated inverse Chi^2
        X = rrtinvchi2.1(h, trnc)
      }

      ## C = cosh(Z) * exp( -0.5 * Z^2 * X )
      ## Don't need to multiply everything by C, since it cancels in inequality.

      S = a.coef(0,X,h)
      ## B = a.coef.1(0, X, trnc)
      D = g.tilde(X, h, trnc)
      Y = runif(1) * D
      n = 0

      ## cat("B,C,left", B, C, X<trnc, "\n")

      a.n = S
      decreasing = FALSE

      while (n < max.inner)
        {
          n = n + 1
          total.iter = total.iter + 1;
          
          a.prev = a.n
          a.n = a.coef(n, X, h)
          ## b.n = a.coef.1(n, X, trnc)
          decreasing = a.n < a.prev
          ## if (!decreasing) cat("n:", n, "; ");
          
          if ( n %% 2 == 1 )
            {
              S = S - a.n
              ## B = B - b.n
              if ( Y<=S && decreasing) break
            }
          else
            {
              S = S + a.n
              ## B = B + b.n
              if ( Y>S && decreasing) break
            }
        }

      ## cat("S,B =", S, B, "\n")

      if ( Y<=S ) break     
    }

  out = list("x"=X, "n"=num.trials, "total.iter"=total.iter)
  out
}

rch <- function(num, h)
{
  h = array(h, num)
  x = rep(0, num)
  for (i in 1:num) {
    out = rch.1(h[i])
    x[i] = out$x
  }
  x
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


################################################################################
                                ## TILTED J^* ##
################################################################################

## pigauss - cumulative distribution function for Inv-Gauss(mu, lambda).
##------------------------------------------------------------------------------
pigauss <- function(x, Z=1, lambda=1)
{
  ## I believe this works when Z = 0
  ## Z = 1/mu
  b = sqrt(lambda / x) * (x * Z - 1);
  a = sqrt(lambda / x) * (x * Z + 1) * -1.0;
  y = exp(pnorm(b, log.p=TRUE)) + exp(2 * lambda * Z + pnorm(a, log.p=TRUE));
                                        # y2 = 2 * pnorm(-1.0 / sqrt(x));
  y
}

p.tilt.left <- function(trnc, h, z)
{
  out = 0
  if (z == 0) {
    out = p.ch.left(trnc, h)
  }
  else {
    out = (2^h * exp(-z*h)) * pigauss(trnc, Z=z/h, h^2)
  }
  out
}

p.tilt.right <- function(trcn, h, z)
{
  ## Note: this works when z=0
  lambda.z = pi^2/8 + z^2/2
  (pi/2/lambda.z)^h * pgamma(trcn, shape=h, rate=lambda.z, lower.tail=FALSE)
}

## rigauss - sample from Inv-Gauss(mu, lambda).
##------------------------------------------------------------------------------
rigauss <- function(mu, lambda)
{
  nu = rnorm(1);
  y  = nu^2;
  x  = mu + 0.5 * mu^2 * y / lambda -
    0.5 * mu / lambda * sqrt(4 * mu * lambda * y + (mu*y)^2);
  if (runif(1) > mu / (mu + x)) {
    x = mu^2 / x;
  }
  x
}

rrtigauss <- function(h, z, R)
{
  ## R is truncation point
  z = abs(z);
  mu = h/z;
  X = R + 1;
  if (mu > R) {
    alpha = 0.0;
    while (runif(1) > alpha) {
      X = rrtinvchi2.1(h, R)
      alpha = exp(-0.5 * z^2 * X);
    }
    ## cat("rtigauss, part i:", X, "\n");
  }
  else {
    while (X > R) {
      lambda = h^2;
      Y = rnorm(1)^2;
      X = mu + 0.5 * mu^2 / lambda * Y -
        0.5 * mu / lambda * sqrt(4 * mu * lambda * Y + (mu * Y)^2);
      if ( runif(1) > mu / (mu + X) ) {
        X = mu^2 / X;
      }
    }
    ## cat("rtiguass, part ii:", X, "\n");
  }
  X;
}

rpg.1 <- function(h, z)
{
  z = z/2
  if (h < 1 || h > 4) return(NA);
  
  rate  = pi^2 / 8 + z^2 / 2
  idx   = floor((h-1) * 100) + 1;
  trnc  = t14[idx];
  
  num.trials = 0;
  total.iter = 0;

  p.l = p.tilt.left (trnc, h, z);
  p.r = p.tilt.right(trnc, h, z);
  p   = p.r / (p.l + p.r);

  ## cat("prob.right:", p, "\n")
  
  max.inner = 100
      
  while (num.trials < 10)
    {
      num.trials = num.trials + 1;

      uu = runif(1)
      if ( uu < p ) {
        ## Left truncated gamma
        X = rltgamma.dagpunar.1(shape=h, rate=rate, trnc=trnc)
      }
      else {
        ## Right truncated inverse Chi^2
        ## Note: this sampler works when z=0.
        X = rrtigauss(h, z, trnc)
      }

      ## C = cosh(Z) * exp( -0.5 * Z^2 * X )
      ## Don't need to multiply everything by C, since it cancels in inequality.

      S = a.coef(0,X,h)
      ## B = a.coef.1(0, X, trnc)
      D = g.tilde(X, h, trnc)
      Y = runif(1) * D
      n = 0

      ## cat("B,C,left", B, C, X<trnc, "\n")
      ## cat("test gt:", g.tilde(trnc * 0.1, h, trnc), "\n");
      ## cat("X, Y, S, gt:", X, Y, S, D,"\n");
      
      a.n = S
      decreasing = FALSE

      while (n < max.inner)
        {
          n = n + 1
          total.iter = total.iter + 1;
          
          a.prev = a.n
          a.n = a.coef(n, X, h)
          ## b.n = a.coef.1(n, X, trnc)
          decreasing = a.n < a.prev
          ## if (!decreasing) cat("n:", n, "; ");
          
          if ( n %% 2 == 1 )
            {
              S = S - a.n
              ## B = B - b.n
              if ( Y<=S && decreasing) break
            }
          else
            {
              S = S + a.n
              ## B = B + b.n
              if ( Y>S && decreasing) break
            }
        }

      ## cat("S,B =", S, B, "\n")

      ## Need to check max.outer
      if ( Y<=S ) break     
    }

  X = 0.25 * X
  out = list("x"=X, "num.trials"=num.trials, "total.iter"=total.iter)
  out
}

rpg.4 <- function(num, h, z)
{
  if (h < 1) return(NA)

  z = array(z, num)
  h = array(h, num)

  out = list(
    draw = rep(0, num),
    num.trials = rep(0, num),
    total.iter = rep(0, num)
    )
  
  for (i in 1:num) {
    draw = rpg.1(h[i], z[i])
    out$draw[i] = draw$x
    out$num.trials[i] = draw$num.trials
    out$total.iter[i] = draw$total.iter
  }

  out
}

################################################################################
                                ## CHECK rpg ##
################################################################################

if (FALSE)
{

  ## source("Ch.R")
  h = 2.3
  z = 1.1

  num = 20000

  samp.1 = rpg.4(num, h, z)
  samp.2 = rpg.gamma(num, h, z)

  mean(samp.1$draw)
  mean(samp.2)

  mean(samp.1$total.iter)

  hist(samp.1$draw, prob=TRUE, breaks=100)
  hist(samp.2, prob=TRUE, add=TRUE, col="#22000022", breaks=100)
  
}

if (FALSE) {

  ## source("Ch.R")
  source("ManualLoad.R")

  nsamp = 10000
  n = 1
  z = 0

  seed = sample.int(10000, 1)
  ## seed = 8922

  set.seed(seed)
  samp.a = rpg.alt(nsamp, n, z)
  ## samp.d = rpg.devroye(nsamp, n, z)
  set.seed(seed)
  samp.4 = rpg.4(nsamp, n, z)
  
  mean(samp.a)
  ## mean(samp.d)
  mean(samp.4$draw)

  hist(samp.a, prob=TRUE, breaks=100)
  hist(samp.4$draw, prob=TRUE, breaks=100, col="#99000022", add=TRUE)
  
  h = 1.5
  z = 0

  set.seed(seed)
  samp.a = rpg.alt(nsamp, h, z)
  ## samp.g = rpg.gamma(nsamp, h, z)
  set.seed(seed)
  samp.4 = rpg.4(nsamp, h, z)

  mean(samp.a)
  ## mean(samp.g)
  mean(samp.4$draw)

  hist(samp.a, prob=TRUE, breaks=100)
  hist(samp.4$draw, prob=TRUE, breaks=100, col="#99000022", add=TRUE)
  
}

if (FALSE)
{

  ## source("Ch.R")
  source("ManualLoad.R")

  reps   = 100
  nsamp = 10000
  n = 4
  z = 0

  seed = sample.int(10000, 1)
  ## seed = 8922

  set.seed(seed)
  start.a = proc.time()
  for (i in 1:reps) {
    samp.a = rpg.alt(nsamp, n, z)
  }
  end.a = proc.time();
  diff.a = end.a - start.a

  set.seed(seed)
  start.d = proc.time()
  for (i in 1:reps) {
    samp.d = rpg.devroye(nsamp, n, z)
  }
  end.d = proc.time()
  diff.d = end.d - start.d

  mean(samp.a)
  mean(samp.d)

  ## hist(samp.a, prob=TRUE, breaks=100)
  ## hist(samp.4$draw, prob=TRUE, breaks=100, col="#99000022", add=TRUE)
  
  h = 3.8
  z = 0

 set.seed(seed)
  start.a = proc.time()
  for (i in 1:reps) {
    samp.a = rpg.alt(nsamp, h, z)
  }
  end.a = proc.time();
  diff.a = end.a - start.a

  set.seed(seed)
  start.g = proc.time()
  for (i in 1:reps) {
    samp.g = rpg.gamma(nsamp, h, z)
  }
  end.g = proc.time()
  diff.g = end.g - start.g

  mean(samp.a)
  mean(samp.g)

  ## hist(samp.a, prob=TRUE, breaks=100)
  ## hist(samp.4$draw, prob=TRUE, breaks=100, col="#99000022", add=TRUE)
  
}
