
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

left.J <- function(xgrid, N, h=1)
{
  a = outer(0:N, xgrid, function(n,x){
    d.n = 2*n+h;
    (-1)^n * (gamma(n+h)/gamma(n+1)) * (d.n/sqrt(2*pi*x^3)) *
      exp(-d.n^2 * 0.5 / x)
  })
  a = (2^h/gamma(h)) * a
  s = apply(a, 2, cumsum)
  out = list("a"=a, "s"=s)
  out 
}

right.J1 <- function(xgrid, N)
{
  a = outer(0:N, xgrid, function(n,x){
    d.n = pi * (n + 0.5)
    (-1)^n * d.n * exp(-d.n^2 * 0.5 * x)
  })
  s = apply(a, 2, cumsum)
  out = list("a"=a, "s"=s)
  out
}

right.Jn <- function(x, h=1)
{
  c.0 = pi^2 / 8
  (0.5 * pi)^(h) * x^(h-1) * exp(-c.0*x) / (gamma(h))
}

##------------------------------------------------------------------------------

polar.cos.1 <- function(r, theta)
{
  rl =  cos(r * cos(theta)) * cosh(r * sin(theta))
  im =  -1 * sin(r * cos(theta)) * sinh(r * sin(theta))
  cbind(rl, im)
}

to.polar <- function(x,y) {
  r = sqrt(x^2+y^2)
  theta = atan(y/x)
  theta = ifelse(y==0 & x==0, 0, theta)
  theta = ifelse(x < 0, theta - pi, theta)
  cbind(r, theta)
}

complex.cos.rt.1 <- function(x, y)
{
  rt = to.polar(x,y)
  polar.cos.1(sqrt(rt[,1]), rt[,2]/2)
}

cos.rt.rl <- function(x, y)
{
  rl = outer(x,y, function(u,v){
    ## cos(u * cos(v)) * cosh(u * sin(v))
    complex.cos.rt.1(u,v)[,1]
  })
  rl
}

cos.rt.im <- function(x,y)
{
  im = outer(x, y, function(u,v){
    complex.cos.rt.1(u,v)[,2]
  })
  im
}

cos.rt <- function(x,y)
{
  rl = cos.rt.rl(x,y)
  im = cos.rt.im(x,y)
  list("rl"=rl, "im"=im)
}

log.cos.rt <- function(t.xgrid, t.ygrid)
{
  lx = length(t.xgrid)
  ly = length(t.ygrid)
  M  = cos.rt(t.xgrid, t.ygrid)
  rt = to.polar(as.numeric(M$rl), as.numeric(M$im))
  alpha = log(rt[,1])
  theta = rt[,2]
  alpha = array(alpha, dim=c(lx,ly)) ## 
  theta = array(theta, dim=c(lx,ly)) ##
  out = list("rl"=alpha, "im"=theta)
  out
}

sp.exponent <- function(txg, tyg, xval)
{
  lx = length(txg)
  ly = length(tyg)
  lcrt = log.cos.rt(txg, tyg)
  lcrt$rl = -1 * lcrt$rl - xval * outer(txg, rep(1, ly), "*")
  lcrt$im = -1 * lcrt$im -  xval * outer(rep(1, lx), tyg, "*")
  ## lcrt$rl = -1 * lcrt$rl
  ## lcrt$im = -1 * lcrt$im
  lcrt
}

polar.cos <- function(r, theta)
{
  rl = outer(r, theta, function(x,y){
    cos(x * cos(y)) * cosh(x * cos(y))
  })
  im = outer(r, theta, function(x,y){
    sin(x * sin(y)) * sinh(x * sin(y))
  })
  list("rl"=rl, "im"=im)
}

##------------------------------------------------------------------------------

tanu <- function(sgrid)
{
  u = sqrt(abs(sgrid))
  out = ifelse(sgrid < 0, tanh(u) / u, tan(u) / u)
  0.5 * out
}

################################################################################

## To test how things change when degrees of freedom increases.

if (FALSE) {

  ## The x-values (density scale).
  ## The h-values: parameter for J^*(h,z).
  xg = seq(0.1, 6, 0.01);
  ## x = seq(0.1, 30, 1);
  h = seq(1,4,0.01);

  ## z in J^*(h,z).
  z = 0.0
  
  n = length(xg);
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
    out = C.h(xg, z=z, h=h[i], N=160, m=1);
    u[i] = ubound(h[i], N=100);
    y[,i] = out$y;
    w[,i] = out$ppsl;
    out = C.h(xg, z=z, h=h[i], N=2, m=2);
    w.2[,i] = out$ppsl;
    out = C.h(xg, z=z, h=h[i], N=2, m=3);
    #w.3[,i] = out$ppsl;
    ## IGamma near 0, Gamma tail?
    ## gam[,i] = dgamma(x, h[i], pi^2 / 8 + z^2 / 2)
    ## gam[,i] = (pi/2) * x^(h[i] - 1) * exp(- pi^2 / 8 * x) / gamma(h[i]);
    ## GOOD FOR (1,2) and kind of (0,1)
    gam[,i] = (pi / 2)^{h[i]} * xg^(h[i] - 1) * exp(- pi^2 / 8 * xg) / gamma(h[i]);
    ig [,i] = dens.igamma(xg, 0.5, h[i]^2 / 2);
  }

  ## For pg corrections
  x = xg / 4

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

  ##----------------------------------------------------------------------------

  idc = seq(1,p,10)
  N   = length(idc)
  
  for (j in 1:N) {
    i = idc[j]
    ## Ratio
    col.l = "#FF000022"
    col.r = "#0000FF22"
    
    if (j==1) {
      col.l = "#FF0000FF"
      col.r = "#0000FFFF"
      plot(x, y[,i] / gam[,i], type="l", main=paste("f(x|h)/proposal; h=1..4"),
           col=col.l, ylim=c(0,1.2), ylab="ratio", xlab="x")
      legend("topleft", legend=c("l(x|h)", "r(x|h)"), col=c(col.l, col.r),
             lty=c(1,1))
    }
    else {
      lines(x, y[,i] / gam[,i], col=col.l)
    }
    lines(x, y[,i] / w[,i], type="l", col=col.r);
    abline(h=1.0, col=1)
    ## readline("<ENTER>")
  }
  
}

################################################################################
                                 ## PLOT J1 ##
################################################################################

if (FALSE) {

  x = seq(0.01, 4, 0.01);
  N = 3
  Nmax = 20
  
  out.r = right.J1(x, Nmax)
  out.l = left.J  (x, Nmax)

  ## png("left-right-partial-sums.png", width=600, height=300)
  w = x / 4
  par(mfrow=c(1,2))
  plot(w, 4*out.l$s[Nmax,], type="l", ylim=c(0,5), ylab="f(x)", xlab="x",
       main="(Left) Partial Sums", col=1)
  legend("topright", legend=paste("S", 1:N-1, sep=""), col=1+1:N, lty=1:N)
  for (i in 1:N) {
    lines(w, 4*out.l$s[i,], col=i+1, lty=i, lwd=2)
  }
  abline(v=0.64/4, lty=3)
  
  plot(w, 4*out.r$s[Nmax,], type="l", ylim=c(0,5), ylab="f(x)", xlab="x",
       main="(Right) Partial Sums", col=1)
  legend("topright", legend=paste("S", 1:N-1, sep=""), col=1+1:N, lty=1:N)
  for (i in 1:N) {
    lines(w, 4*out.r$s[i,], col=i+1, lty=i, lwd=2)
  }
  abline(v=0.64/4, lty=3)
  
}

################################################################################
                              ## RIGHT PROPOSAL ##
################################################################################

if (FALSE)
{

  par(mfrow=c(1,2))
  ## source("Scratch.R")
  ## png("pg-alternate-envelope.png", width=600, height=300);  par(mfrow=c(1,2))
  ## png("pg-left.png", width=600, height=300);  par(mfrow=c(1,2))
  xgrid = seq(0.01, 6, 0.01)
  N = 3
  Nmax = 20
  h = 3
  
  out.l = left.J  (xgrid, Nmax, h)
  out.r = right.Jn(xgrid, h)

  idint = which.min(abs(out.l$s[1,-1]-out.r[-1]))
  
  w = xgrid / 4
  plot(w, 4*out.l$s[1,], type="l", ylim=c(0,3.2), ylab="f(x)", xlab="x",
       main=paste("PG(", h, ")", sep=""), col=2)

  ## lines(w, 4*out.r, col=2, lty=2, lwd=2)
  ## legend("topright", legend=c(paste("Left S", c(1:N)-1, sep=""), "Right"), col=c(1+1:N,2), lty=c(1:N,2))
  legend("topright", legend=c(paste("Left S", c(1:N)-1, sep="")), col=c(1+1:N), lty=c(1:N))
  
  for (i in 2:N) {
    lines(w, 4*out.l$s[i,], col=i+1, lty=i, lwd=2)
  }
  lines(w, 4*out.l$s[Nmax,], col=1)
  ## abline(v=w[idint], lty=3)
  
}

################################################################################
                            ## Alternate Proposal ##
################################################################################

if (FALSE) {


  
}

################################################################################
                               ## SADDLE POINT ##
################################################################################

if (FALSE) {

  ## source("Scratch.R")

  xval = 1/2
  
  xgrid = seq(-4, (pi/2)^2+1, length.out=100)
  ygrid = seq(-3, 3, length.out=100)
  spe = sp.exponent(xgrid, ygrid, xval)

  ## png("sp-exponent-contours.png", width=600, height=300)
  par(mfrow=c(1,2))
  contour(2*xgrid, 2*ygrid, spe$rl, main="Real Contours K(s)-sx", ylab="Im", xlab="Re")
  contour(2*xgrid, 2*ygrid, spe$im)

  image(2*xgrid, 2*ygrid, spe$rl, main="Real Contours K(s)-sx", ylab="Im", xlab="Re")
  image(2*xgrid, 2*ygrid, spe$im)

  ## 1/2 <-> 0

  tval = 1.5
  xval = tanu(tval)
  
  spe.ra = sp.exponent(xgrid, 0, xval)$rl[,1]
  spe.ia = sp.exponent(tval, ygrid, xval)$rl[1,]
  spe.0  = sp.exponent(xgrid, 0, 0)$rl[,1]
  tu     = tanu(xgrid)
  
  ## png("saddlepoint-location.png", width=400, height=200)
  par(mfrow=c(1,3))
  idx = which(xgrid<(pi/2)^2)
  plot(2*xgrid[idx], tu[idx], type="l",
       main="K'(s) along Re Axis", xlab="Re", ylab="K'")
  abline(h=xval, col=2, lty=2, lwd=2)
  abline(v=2*tval, col=3, lty=3, lwd=2)
  plot(2*xgrid[idx], spe.ra[idx], type="l",
       main="K(s) - sx along Re Axis", ylab="", xlab="Re")
  abline(v=2*tval, col=3, lty=3, lwd=2)
  plot(2*ygrid, spe.ia, type="l",
       main="Re K(s) - sx along Im Axis", ylab="", xlab="Im")
}

################################################################################
################################################################################

if (FALSE)
{

  source("ManualLoad.R")

  num = 10000;
  n   = 2
  z   = 0
  
  samp.d = rpg.devroye(num, n, z)
  samp.s = rpg.sp(num, n, z)

  
}

################################################################################

################################################################################

if (FALSE) {

  N = 10
  
  A = diag(rnorm(N), N)
  S = t(A) %*% A
  h = rep(1/N, N)

  loss <- function(w, lambda, p) {
    ## assume sum(w) = 1
    d = w - h
    q = t(d) %*% S %*% d
    c = q[1] + lambda * sum(abs(w)^p)
  }

  lddir <- function(w, alpha) {
    d1 = sum( (alpha - 1) * log(w) )
    d2 = sum(lgamma(alpha)) - lgamma(sum(alpha))
    if (any(alpha==0)) print(alpha)
    d1 - d2
  }

  mdir1 <- function(w1, alpha1, w2, alpha2) {
    d1 = sum( (alpha1 - 1) * log(w1) )
    d2 = sum( (alpha2 - 1) * log(w2) )
    d3 = sum( lgamma(alpha1) ) - sum( lgamma(alpha2) )
    d1 - d2 - d3
  }

  lam = 10
  p   = 0.5
  c   = 40
  
  max.iter = 2000
  wold  = rep(1/N, N)
  ## wold = c(0.99, 0.01)
  ## wold = c(0.01, 0.99)
  loold = -1 * loss(wold, lam, p)
  wbes  = wold
  lobes = loold
  temper = 1
  thresh = 100

  M = max.iter
  naccept = 0
  out <- list("w"=matrix(0, nrow=M, ncol=N),
              "le"=rep(0, M))
  
  for (i in 1:max.iter) {
    anew  = c * wold
    anew  = ifelse(anew > thresh, thresh, anew)
    wnew  = as.numeric(rdirichlet(1, anew))
    lonew = -1 * loss(wnew, lam, p)
    aold  = c * wnew
    aold  = ifelse(aold > thresh, thresh, aold)
    plls  = lonew - loold
    plpr  = log(ddirichlet(wold, aold)) - log(ddirichlet(wnew, anew))
    if (any(wnew==0)) {
      ## print(wnew)
      plpr = -Inf
    }
    ## plpr  = lddir(wold, aold) - lddir(wnew, anew)
    if (log(runif(1)) < (temper * plls + plpr)) {
      wold  = wnew
      loold = lonew
      naccept = naccept + 1
    }
    if (lobes < loold) {
      wbes  = wold
      lobes = loold
    }

    out$w[i,]=wold
    out$le[i]=loold
    out$arate = naccept / max.iter
    
  }
  
}
