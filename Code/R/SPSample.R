

################################################################################

y.func <- function(v)
{
  ## tan(sqrt(v))/sqrt(v))
  out = v
  gt.idc = v >  1e-6
  lt.idc = v < -1e-6
  md.idc = v <= 1e-6 & v >= -1e-6

  r = sqrt(abs(v))
  gt.val = r[gt.idc]
  out[gt.idc] = tan(gt.val)/gt.val
  lt.val = r[lt.idc]
  out[lt.idc] = tanh(lt.val)/ lt.val
  md.val = r[md.idc]
  out[md.idc] = 1 - (1/3) * md.val^2 + (2/15) * md.val^4 - (17/315) * md.val^6
  out
}

v.approx <- function(y)
{
  ifelse(y >= 1, atan(y * pi / 2)^2, -1 / y^2)
}

v.func.1 <- function(y)
{
  if (y==1) return(list(root=0.0, f.root=1.0, iter=0, estim.prec=1e-16))
  f <- function(v) { y - y.func(v) }
  lowerb = 0; upperb = 0;
  if (y > 1) upperb = (pi/2)^2
  if (y < 1) lowerb = min(-5, v.approx(y/2))
  out = uniroot(f, lower=lowerb, upper=upperb, maxiter=10000, tol=1e-8)
  out
}

v.func.1.alt <- function(y, lowerb, upperb)
{
  if (y==1) return(list(root=0.0, f.root=1.0, iter=0, estim.prec=1e-16))
  f <- function(v) { y - y.func(v) }
  out = uniroot(f, lower=lowerb, upper=upperb, maxiter=10000, tol=1e-8)
  out
}

v.secant.1 <- function(y, vb, va, tol=1e-8, maxiter=10)
{
  ## Assumes increasing and convex function.
  yb = y.func(vb)
  ya = y.func(va)
  if (yb > ya) return(NA)
  iter = 0
  ydiff = tol + 1
  while(abs(ydiff) > tol && iter < maxiter) {
    iter = iter + 1
    m = (ya - yb) / (va - vb)
    ydiff = y - yb
    vstar = ydiff / m + vb;
    ystar = y.func(vstar)
    if (ystar < y) {
      vb = vstar
      yb = ystar
    }
    else {
      va = vstar
      yb = ystar
    }
  }
  out = list(iter=iter, ydiff=ydiff, y=ystar, v=vstar)
  out
}

v.func <- function(y, lowerb=0, upperb=0)
{
  ## y inverse.
  N = length(y)
  a = rep(0,N)
  out = list(root=a, f.root=a, iter=a, estim.prec=a)
  for (i in 1:N) {
    temp = v.func.1(y[i], lowerb, upperb)
    out$root[i] = temp$root
    out$f.root[i] = temp$f.root
    out$iter[i] = temp$iter
    out$estim.prec[i] = temp$estim.prec
  }
  ## out = simplify2array(out)
  out$root
}

v.iterative <- function(y, ygrid, vgrid, dy)
{
  ## Assume ygrid is equally spaced.
  ## Assume ygrid is increasing.
  N   = length(ygrid)
  idx = floor((y-ygrid[1]) / dy) + 1
  if (idx >= N || idx<1) return(NA)
  dv   = vgrid[idx+1] - vgrid[idx]
  ## dy   = ygrid[idx+1] - ygrid[idx]
  ## Just a single secant approximation
  ## vapprox = (dv / dy) * (y - ygrid[idx]) + vgrid[idx]

  vb = vgrid[idx]
  va = vgrid[idx+1]
  ## vout = v.secant.1(y, vb, va)
  ## vapprox = vout$v

  vout = v.func.1.alt(y, vb, va)
  vapprox = vout$root

  print(vout)
  vapprox
}

################################################################################

if (FALSE) {
  dyl   = 0.01
  yleft = seq(0.1, 1, dyl)
  vleft = v.func(yleft)

  dyr    = 0.01
  yright = seq(1, 8, dyr)
  vright = v.func(yright)
}

v.table <- function(y)
{
  out = 0
  if (y <= 1)
    out = v.iterative(y, yleft, vleft, dyl)
  else
    out = v.iterative(y, yright, vright, dyr)

  if (is.na(out))
    out = v.approx(y)

  out
}

################################################################################
                   ## CALCULATE SADDLE POINT APPROXIMATION ##
################################################################################

x.left <- function(s)
{
  ## tanh(sqrt(2s))/sqrt(2s)
  if (any(s<0)) {
    print("x.left: s must be >= 0.")
    return(NA)
  }

  s = sqrt(2*s)
  tanh(s) / s
}

x.right <- function(s)
{
  ## tan(sqrt(2s))/sqrt(2s)
  if (any(s<0)) {
    print("x.left: s must be >= 0.")
    return(NA)
  }

  s = sqrt(2*s)
  tan(s) / s
}

cgf <- function(s, z)
{
  v = 2*s + z^2;
  v = sqrt(abs(v));
  out = log(cosh(v))
  out[s>0] = log(cos(v[s>0]));
  ## out = ifelse(s >= 0, log(cosh(u)), log(cos(u)))
  out = log(cosh(z)) - out
  out
}

sp.approx.1 <- function(x, n=1, z=0)
{
  ## v = v.table(x)
  v = v.func(x)
  u = v / 2
  t = u + z^2/2
  m = y.func(-z^2)

  temp = sqrt(abs(v))
  phi      = -log(cosh(temp))
  phi[v>0] = -log(cos (temp[v>0]))
  phi      = phi + log(cosh(z)) - t*x;
  
  K2 = x^2 + (1-x)/(2*u)
  K2[u<1e-5 & u>-1e-5] = x^2 - 1/3 + 2/15 * (2*u)

  spa = (0.5*n/pi)^0.5 * K2^-0.5 * exp(n*phi)

  out = list("spa"=spa, "phi"=phi, "K2"=K2, "t"=t, "x"=x, "u"=u, "m"=m);
  out
}

sp.approx.df <- function(x, n=1, z=0)
{
  N = length(x)
  a = rep(0, N)
  df = data.frame("spa"=a, "phi"=a, "K2"=a, "t"=a, "x"=a, "u"=a, "m"=a);
  for (i in 1:N) {
    temp = sp.approx.1(x[i], n, z);
    df[i,] = as.numeric(temp)
  }
  df
}

sp.approx <- function(x, n=1, z=0)
{
  N = length(x)
  spa = rep(0, N)
  for (i in 1:N) {
    temp = sp.approx.1(x[i], n, z)
    spa[i] = temp$spa
  }
  spa
}

##------------------------------------------------------------------------------
                                ## UNIT TEST ##
##------------------------------------------------------------------------------

if (FALSE) {

  ## source("SaddlePointApprox.R"); source("SPSample.R")
  n = 10
  z = 0
  
  dx = 0.01
  xgrid = seq(dx, 4, dx)
  spa = sp.approx(xgrid, n, z)

  plot(xgrid, spa, type="l")

  y1    = xgrid
  N     = length(xgrid)

  for (i in 1:N) {
    y1[i] = 0
    y2[i] = 0
    for (j in 0:200) {
      y1[i] = y1[i] + (-1)^j * a.coef(j,xgrid[i] * n ,n, z) * n
    }
  }
  
  lines(xgrid, y1, col=2)

  (1:N)[y1>spa]
  
}

################################################################################
                          ## POINTS OF INTERSECTION ##
################################################################################

delta.func1 <- function(x, m=1)
{
  val = ifelse(x >= m, log(x) - log(m), 0.5 * (1-1/x) - 0.5 * (1-1/m))
  der = ifelse(x >= m, 1/x, 0.5 / x^2)
  out = data.frame(val=val, der=der)
  out
}

delta.func2 <- function(x, m=1)
{
  val = ifelse(x >= m, log(x) - log(m), 0)
  der = ifelse(x >= m, 1/x, 0)
  out = data.frame(val=val, der=der)
  out
}

phi.func <- function(x, z)
{
  v = v.func(x)
  u = v / 2
  t = u + z^2/2
  
  temp = sqrt(abs(v))
  phi      = -log(cosh(temp))
  phi[v>0] = -log(cos (temp[v>0]))
  phi      = phi + log(cosh(z)) - t*x;

  val = phi
  der = -t

  out = data.frame(val=val, der=der)
  out
}

phi.eta.delta <- function(x, z=0, m=1)
{
  ## versions of phi, eta, delta
  phi = phi.func(x, z)
  ## phi$val = phi$val - phi.func(m,z)$val
  
  delta   = delta.func1(x,m)

  eta = phi - delta

  out <- data.frame(phi=phi$val, eta=eta$val, delta=delta$val,
              phi.d=phi$der, eta.d=eta$der, delta.d=delta$der)
  ## out <- cbind(phi, eta, delta)
  out
}

tangent.lines.eta <- function(x, z=0, m=1)
{
  phed = phi.eta.delta(x, z, m)
  slope = phed$eta.d
  icept = - phed$eta.d * x + phed$eta
  out = data.frame(slope=slope, icept=icept)
  out
}

## Maybe remove
find.ev.mode.1 <- function(z)
{
  z = abs(z)
  m = y.func(-z^2)
  f <- function(v) {
    (v + z^2 - 1/m^2) * y.func(v) + 1
  }
  check = z^2 * tanh(abs(z))^4;
  upperb = 0
  lowerb = -4
  if (check > 1) {
    lowerb = 0
    upperb = check
  }
  out = uniroot(f, lower=-100, upper=0, maxiter=10000, tol=1e-8)
  x = y.func(out$root)
  x
}

find.ev.mode <- function(z)
{
  N = length(z)
  x = rep(0, N)
  for (i in 1:N) {
    x[i] = find.ev.mode.1(z[i])
  }
  x
}

ig.mode <- function(mu, lambda)
{
  mu * ( (1+9*mu^2 / (4 * lambda))^0.5 - 1.5 * mu / lambda)
}

##------------------------------------------------------------------------------
                                ## UNIT TEST ##
##------------------------------------------------------------------------------

if (FALSE)
{

  n = 10
  z = 10
  
  dx = 0.001
  xgrid = seq(dx*5, 2, dx)
  spa = sp.approx(xgrid, n, z)
  spa.m = sp.approx.1(xgrid[1], n, z)
  spa.m = sp.approx.1(spa.m$m, n, z)
  
  m = spa.m$m
  ## m = 1
  ## m = spa.m$m + 0.1
  kappa = 1.1
  ## xl = m / kappa
  ## xl = find.ev.mode(z)
  ## xl = spa.m$m / 2
  ## xl = spa.m$m
  ## xr = m * kappa

  xl = spa.m$m
  m  = spa.m$m * 1.1
  xr = spa.m$m * 1.2
  
  phi.m = phi.func(m,z)$val
  
  ## plot(xgrid, sp.approx, type="l")

  phed = phi.eta.delta(xgrid, z, m)
  
  par(mfrow=c(1,2))

  ## idxp = xgrid >= 0
  idxp = xgrid<=0.2 & xgrid >=0.001
  
  plot(xgrid[idxp], phed$phi[idxp], type="l")
  lines(xgrid[idxp], phed$eta[idxp], col=2)
  lines(xgrid[idxp], -phed$delta[idxp] + phi.m, col=3)

  tl = tangent.lines.eta(c(xl, xr, m), z, m)

  lines(xgrid[idxp], tl$slope[1]*xgrid[idxp] + tl$icept[1], col=4)
  lines(xgrid[idxp], tl$slope[2]*xgrid[idxp] + tl$icept[2], col=5)

  lines(xgrid[idxp], -0.5*z^2 * (xgrid[idxp]-m) + phi.m, col=6, lty=4)

  plot(xgrid[idxp], phed$phi[idxp], type="l")

  phi.ev1 = tl$slope[1]*xgrid + tl$icept[1] + phed$delta
  phi.ev2 = tl$slope[2]*xgrid + tl$icept[2] + phed$delta
  phi.ev3 = rep(phi.m, length(xgrid))
  
  lines(xgrid[idxp], phi.ev1[idxp], col=4)
  lines(xgrid[idxp], phi.ev2[idxp], col=5)

  ##----------------------------------------------------------------------------

  par(mfrow=c(1,2))
  spa.df = sp.approx.df(xgrid, n, z)
  adj.phi = spa.df$phi - 0.5 * log(spa.df$K2) / n
  
  ra = c(spa.df$phi[idxp], adj.phi[idxp]); ymax = max(ra); ymin = min(ra);
  plot(xgrid[idxp], spa.df$phi[idxp], type="l", ylim=c(ymin, ymax))
  lines(xgrid[idxp], adj.phi[idxp], col=2)
  lines(xgrid[idxp], phi.ev1[idxp], col=4)
  lines(xgrid[idxp], phi.ev2[idxp], col=5)
  lines(xgrid[idxp], phi.ev3[idxp], col=6)
  ## lines(xgrid, phed$phi, col=7)

  plot(xgrid[idxp], spa.df$spa[idxp], type="l")
  lines(xgrid[idxp], spa.df$spa[idxp] * spa.df$K2[idxp]^0.5, col=2)

  ## ----------------------------------------
  
  plot(xgrid[idxp], exp(spa.df$phi[idxp]), type="l")
  lines(xgrid[idxp], exp(adj.phi[idxp]), col=2)
  lines(xgrid[idxp], exp(phi.ev1[idxp]), col=4, lty=2)
  lines(xgrid[idxp], exp(phi.ev2[idxp]), col=5, lty=2)
  lines(xgrid[idxp], exp(phi.ev3[idxp]), col=6)
  ## lines(xgrid[idxp], spa.df$spa[idxp] * (0.5*n/pi)^-0.5 * spa.df$K2[idxp]^0.5, col=1, lty=2)

  idr = xgrid >= m
  ## idx.alt = which.max(idr)
  idl = xgrid <= m
  ar  = max(xgrid[idr]^2 / spa.df$K2[idr])
  al  = max(xgrid[idl]^3 / spa.df$K2[idl])

  ev3 = al^0.5 * xgrid^(-1.5) * (0.5*n/pi)^0.5 
  ## ev3[xgrid >= m] = ar^0.5 * xgrid[xgrid>=m]^(-1) * (0.5*n/pi)^0.5
  ev3 = ev3 * exp(n * phi.m)

  ev1 = (al)^0.5 * exp(n*phi.ev1) * xgrid^-1.5 * (0.5*n/pi)^0.5
  ev2 = (ar)^0.5 * exp(n*phi.ev2) * xgrid^-1 * (0.5*n/pi)^0.5
  ev4 = exp(n*phi.ev1) * spa.df$K2^-0.5 * (0.5*n/pi)^0.5
  
  ra = c(spa.df$spa[idxp], ev1[xgrid<m]); ymax = max(ra); ymin = min(ra);
  plot(xgrid[idxp], spa.df$spa[idxp], type="l", ylim=c(ymin, ymax))
  lines(xgrid[idxp], (0.5*n/pi)^0.5 * exp(n*adj.phi[idxp]), col=3, lty=2)
  lines(xgrid[idxp], ev4[idxp], col=8, lty=2)
  lines(xgrid[idxp], ev1[idxp], col=4, lty=2)
  lines(xgrid[idxp], ev2[idxp], col=5, lty=2)
  lines(xgrid[idxp], ev3[idxp], col=6, lty=2)

  ## ----------------------------------------

  plot(xgrid[idxp], spa.df$phi[idxp], type="l")
  lines(xgrid[idxp], spa.df$phi[idxp] - 1.5/n*log(xgrid[idxp]), col=2, lty=2)

  plot(xgrid[idxp], spa.df$K2[idxp]^-0.5)  

}

################################################################################
                                ## rrtigauss ##
################################################################################

rrtinvch2.1 <- function(scale, trnc=1)
{
  R = trnc / scale
  ## E = rexp(2)
  ## while ( (E[1]^2) > (2 * E[2] / R)) {
  ##   ## cat("E", E[1], E[2], E[1]^2, 2*E[2] / R, "\n")
  ##   E = rexp(2)
  ## }  
  ## ## cat("E", E[1], E[2], "\n")
  ## ## W^2 = (1 + R*E[1])^2 / R is left truncated chi^2(1) I_{(R,\infty)}.
  ## ## So X is right truncated inverse chi^2(1).
  ## X = R / (1 + R*E[1])^2
  ## X = scale * X
  E = rtnorm(1, 0, 1, left=1/sqrt(R), right=Inf)
  X = scale / E^2
  X
}

rrtinvch2 <- function(num, scale, trnc)
{
  out = rep(0, num)
  for (i in 1:num)
    out[i] = rrtinvch2.1(scale, trnc)
  out
}

rrtigauss.ch2.1 <- function(mu, lambda, trnc)
{
  iter = 0
  alpha = 0.0;
  accept = FALSE
  while (!accept) {
    iter = iter + 1
    X = rrtinvch2.1(lambda, trnc)
    alpha = exp(- 0.5*lambda/mu^2 * X)
    ## l.alpha = 0.5*lambda/mu^2 * (X-trnc)
    ## cat(X, alpha, "\n")
    accept = runif(1) < alpha
  }
  out = list(x=X, iter=iter)
  out
}

rrtigauss.ch2 <- function(num, mu, lambda, trnc)
{
  x = rep(0, num)
  df = data.frame(x=x, iter=x)
  for (i in 1:num) {
    temp = rrtigauss.ch2.1(mu, lambda, trnc)
    df[i,] = as.numeric(temp)
  }
  df
}

rrtigauss.1 <- function(mu, lambda, trnc=1)
{
  ## trnc is truncation point
  accept = FALSE
  X = trnc + 1;
  if (trnc < mu) {
    alpha = 0.0;
    while (!accept) {
      X = rrtinvch2.1(lambda, trnc)
      ## alpha = exp(0.5*lambda/mu^2 * (X - trnc))
      l.alpha = - 0.5*lambda/mu^2 * X
      ## cat(X, alpha, "\n")
      accept = log(runif(1)) < l.alpha
    }
    ## cat("rtigauss.ch, part i:", X, "\n");
  }
  else { ## trnc >= mu
    while (X > trnc) {
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

rrtigauss <- function(num, mu, lambda, trnc=1)
{
  x = rep(0, num)
  for (i in 1:num)
    x[i] = rrtigauss.1(mu, lambda, trnc)
  x
}

rigauss.1 <- function(mu, lambda)
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

rigauss <- function(num, mu, lambda)
{
  x = rep(0, num)
  for (i in 1:num)
    x[i] = rigauss.1(mu, lambda)
  x
}

##------------------------------------------------------------------------------
                                ## UNIT TEST ##
##------------------------------------------------------------------------------

if (FALSE) {

  ## source("SPSample.R")
  lambda = 5
  mu = 2
  
  draw1 = rigauss(10000, mu, lambda)
  draw2 = rrtigauss(5000, mu, lambda, mu/2)
  draw3 = rrtigauss(5000, mu, lambda, mu*2)
  draw8 = rrtigauss.ch2(5000, mu, lambda, mu/2)

  par(mfrow=c(1,2))
  
  hist(draw1[draw1<=mu/2], prob=TRUE)
  hist(draw2, prob=TRUE, add=TRUE, col="#88000088")

  summary(draw1[draw1<=mu/2])
  summary(draw2)
  summary(draw8$x)
  
  hist(draw1[draw1<=mu*2], prob=TRUE)
  hist(draw3, prob=TRUE, add=TRUE, col="#88000088")

  summary(draw1[draw1<=mu*2])
  summary(draw3)

  trnc = 2
  lambda = 100
  draw4 = lambda / rchisq(20000, 1)
  draw5 = rrtinvch2(10000, lambda, trnc)

  hist(draw4[draw4<=trnc], prob=TRUE)
  hist(draw5, prob=TRUE, add=TRUE, col="#88000088")

  summary(draw4[draw4<=trnc])
  summary(draw5)

  draw6 = rrtigauss(5000, 1/sqrt(2*rl), n, 1)
  draw7 = rigauss(100000, 1/sqrt(2*rl), n)
  hist(draw6, prob=TRUE)
  hist(draw7[draw7<=1], prob=TRUE, add=TRUE, col="#88000088")
  
}

################################################################################
                                ## SP SAMPLER ##
################################################################################

if (FALSE)
{

  m = 1
  kappa = 1.1
  xl = m/kappa
  xr = m*kappa

  ul = 0.5 * v.func(xl)
  ur = 0.5 * v.func(xr)
  
  delta = delta.func1(c(xl,xr))  
  tl = tangent.lines.eta(c(xl,xr), 0)

  rlb = -tl$slope[1]
  rrb = -tl$slope[2]
  ilb = cgf(ul, 0) - delta$val[1] + delta$der[1] * xl
  irb = cgf(ur, 0) - delta$val[2] + delta$der[2] * xr
  
  ## alphal = 3/2
  ## alphar = 3/2
  alpha = 3/2
  inflate = alpha^0.5
  
}

sp.sampler.1 <- function(n=1, z=0, maxiter=100)
{
  if (n < 1) print("sp.sampler.1: n must be >= 1.")
  ## rl = 0.5 * z^2 - tl$slope[1]
  ## rr = 0.5 * z^2 - tl$slope[2]

  xl  = y.func(-z^2)
  mid = xl * 1.1
  xr  = xl * 1.2
  ## xl = mid / kappa
  ## xr = mid * kappa

  v.mid  = v.func(mid)
  K2.mid = mid^2 + (1-mid) / v.mid
  al     = mid^3 / K2.mid
  ar     = mid^2 / K2.mid

  tl = tangent.lines.eta(c(xl,xr), z, mid)
  rl = -tl$slope[1]
  rr = -tl$slope[2]
  il = tl$icept[1]
  ir = tl$icept[2]

  ## rl = rlb + 0.5 * z^2
  ## rr = rrb + 0.5 * z^2
  ## il = log(cosh(z)) + ilb
  ## ir = log(cosh(z)) + irb
  
  wl = al^0.5 * exp(-n*sqrt(2*rl) + n * il + 0.5*n - 0.5*n*(1-1/mid)) *
    pigauss(mid, Z=sqrt(2*rl), lambda=n)

  wr = ar^0.5 * (0.5*n/pi)^0.5 * exp(-n * log(n * rr) + n * ir - n * log(mid)) *
    gamma(n) * pgamma(mid, shape=n, rate=n * rr, lower=FALSE)
  
  w  = wl + wr

  go   = TRUE
  iter = 0
  X    = 2
  FX   = 0
  
  while (go && iter<maxiter) {
    iter = iter + 1
    
    if (w*runif(1) < wl) {
      ## sample left
      X  = rrtigauss.1(mu=1/sqrt(2*rl), lambda=n, trnc=1)
      ## while (X > 1) X = rigauss.1(1/sqrt(2*rl), n)
      phi.ev = n * (-rl * X + il) + 0.5 * n * ( (1-1/X) - (1-1/mid))
      FX = al^0.5 * (0.5 * n / pi)^0.5 * X^(-1.5) * exp(phi.ev);

    }
    else {
      ## sample right
      X  = rltgamma.dagpunar.1(shape=n, rate=n*rr, trnc=1)
      phi.ev = n * (-rr * X + ir) + n * (log(X) - log(mid))
      FX = ar^0.5 * (0.5 * n / pi)^0.5 * exp(phi.ev) / X;
    }

    spa = sp.approx.1(X, n, z)

    ## cat("FX, phi.ev, spa, phi", FX, phi.ev, spa$spa, spa$phi,"\n")

    if (FX*runif(1) < spa$spa) go = FALSE
    
  }
  
  out = list(x=X, iter=iter, wl=wl, wr=wr)
  out
}

sp.sampler <- function(num, n=1, z=0, return.df=FALSE)
{
  x  = rep(0, num)
  df = data.frame(x=x, iter=x)
  for (i in 1:num) {
    temp = sp.sampler.1(n,z)
    df$x[i] = temp$x
    df$iter[i] = temp$iter
  }
  if (!return.df) df = df$x
  df
}

##------------------------------------------------------------------------------
                                ## UNIT TEST ##
##------------------------------------------------------------------------------

if (FALSE)
{

  ## source("SPSample.R")
  n = 10
  z = 100

  dx = 0.001
  xgrid = seq(dx*5, 0.2, dx)
  spa = sp.approx(xgrid, n, z)

  spa.m = sp.approx.1(xgrid[1], n, z)
  spa.m = sp.approx.1(spa.m$m, n, z)
  xl  = spa.m$m
  mid = spa.m$m * 1.1
  xr  = spa.m$m * 1.2

  wla = sum(spa[xgrid<mid]) * dx
  wra = sum(spa[xgrid>=mid]) * dx

  d1 = sp.sampler.1(n,z)

  d1$wl / (d1$wl + d1$wr)
  wla   / (wla + wra)

  df = sp.sampler(5000, n, z, TRUE)
  
  tl = tl = tangent.lines.eta(c(xl,xr), z)
  rl = -tl$slope[1]
  rr = -tl$slope[2]
  draw.ig = rrtigauss(2000, mu=1/sqrt(2*rl), lambda=n, trnc=1)
  draw.ga = rltgamma.dagpunar(2000, shape=n, rate=n*rr, trnc=1)
  draw.jj = 4 * rpg(2000, n, 2*z) / n
  
  par(mfrow=c(1,3))
  
  plot(xgrid, spa)
  hist(df$x, prob=TRUE, add=TRUE, breaks=20)
  hist(draw.jj, prob=TRUE, add=TRUE, col="#22000022")

  plot(xgrid[xgrid<mid], spa[xgrid<mid] /  (wla / (wla+wra)))
  hist(df$x[df$x<mid], prob=TRUE, add=TRUE)
  hist(draw.ig, prob=TRUE, col="#00004444", add=TRUE)
  
  plot(xgrid[xgrid>=mid], spa[xgrid>=mid] / (wra / (wla+wra)))
  hist(df$x[df$x>=mid], prob=TRUE, add=TRUE)
  hist(draw.ga, prob=TRUE, col="#22000022", add=TRUE)

  c(mean(draw.jj), var(draw.jj))
  c(mean(df$x), var(df$x))
  
}
