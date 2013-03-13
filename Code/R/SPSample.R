
## if (!exists(ugrid) && !exists(xgrid)) {
  
## }

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

v.func <- function(y)
{
  ## y inverse.
  N = length(y)
  a = rep(0,N)
  out = list(root=a, f.root=a, iter=a, estim.prec=a)
  for (i in 1:N) {
    temp = v.func.1(y[i])
    out$root[i] = temp$root
    out$f.root[i] = temp$f.root
    out$iter[i] = temp$iter
    out$estim.prec[i] = temp$estim.prec
  }
  ## out = simplify2array(out)
  out$root
}

v.secant <- function(y, ygrid, vgrid, dy)
{
  ## Assume ygrid is equally spaced.
  N   = length(ygrid)
  idx = floor((y-ygrid[1]) / dy) + 1
  if (idx >= N || idx<1) return(NA)
  dv   = vgrid[idx+1] - vgrid[idx]
  ## dy   = ygrid[idx+1] - ygrid[idx]
  vapprox = (dv / dy) * (y - ygrid[idx]) + vgrid[idx]
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
    out = v.secant(y, yleft, vleft, dyl)
  else
    out = v.secant(y, yright, vright, dyr)

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
  out[s>0] = log(cos(v[s<0]));
  ## out = ifelse(s >= 0, log(cosh(u)), log(cos(u)))
  out = log(cosh(z)) - out
  out
}

sp.approx.1 <- function(x, z=0, n=1)
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

sp.approx <- function(x, z=0, n=1)
{
  N = length(x)
  spa = rep(0, N)
  for (i in 1:N) {
    temp = sp.approx.1(x[i], z, n)
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
  sp.approx = sp.dens(xgrid, z, n)

  plot(xgrid, sp.approx, type="l")

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

  (1:N)[y1>sp.approx]
  
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
}

phi.eta.delta <- function(x, z=0, m=1)
{
  ## normalized versions of phi, eta, delta
  phi = phi.func(x, z)
  phi$val = phi$val - phi.func(m,z)$val
  
  delta   = delta.func1(x,m)

  eta = phi - delta

  out <- data.frame(phi=phi$val, eta=eta$val, delta=delta$val,
              phi.d=phi$der, eta.d=eta$der, delta.d=delta$der)
  ## out <- cbind(phi, eta, delta)
  out
}

tangent.lines <- function(x, z=0, m=1)
{
  phed = phi.eta.delta(x, z, m)
  slope = phed$eta.d
  icept = - phed$eta.d * x + phed$eta
  out = data.frame(slope=slope, icept=icept)
  out
}

##------------------------------------------------------------------------------
                                ## UNIT TEST ##
##------------------------------------------------------------------------------

if (FALSE)
{

  n = 100
  z = 2
  
  dx = 0.01
  xgrid = seq(dx*5, 4, dx)
  spa = sp.approx(xgrid, z, n)
  spa.m = sp.approx.1(xgrid[1], z, n)
  m = spa.m$m
  phi.m = phi.func(m,z)$val
  
  ## plot(xgrid, sp.approx, type="l")

  phed = phi.eta.delta(xgrid, z, m)

  par(mfrow=c(1,2))
  
  plot(xgrid, phed$phi, type="l")
  lines(xgrid, phed$eta, col=2)
  lines(xgrid, -phed$delta, col=3)

  xl = m / 2
  xr = m * 2

  tl = tangent.lines(c(xl, xr), z, m)

  lines(xgrid, tl$slope[1]*xgrid + tl$icept[1], col=4)
  lines(xgrid, tl$slope[2]*xgrid + tl$icept[2], col=5)

  lines(xgrid, -0.5*z^2 * (xgrid-m) + phi.m, col=6, lty=4)

  plot(xgrid, phed$phi, type="l")

  lines(xgrid, tl$slope[1]*xgrid + tl$icept[1] + phed$delta, col=4)
  lines(xgrid, tl$slope[2]*xgrid + tl$icept[2] + phed$delta, col=5)
  
}
