## In NB we need to estimate d from NB(d, p).  These routines sample d via
## random walk MH.

df.llh <- function(y, d, mu, G, ymax)
{
  p =  1 / (1 + d / mu)
  ## sum( log(d+0:(ymax-1)) * G[1:ymax] ) + d * sum(log(1-p)) + sum(y * log(p)) ;
  llh1 = sum( log(d+0:(ymax-1)) * G[1:ymax] )
  llh2 = d * sum(log(d) - log(mu+d))
  llh3 = sum(y * (log(mu) - log(mu+d)))
  ## cat(llh1, llh2, llh3, "\n")
  llh1 + llh2 + llh3
}

draw.df <- function(y, d.prev, mu, G, ymax)
{
  ## optim.out <- optim(d.prev, fn=df.llh, gr = NULL,
  ##                    mu, G, ymax,                              ## for minimization
  ##                    method="L-BFGS-B", lower=1, hessian=TRUE, control=list(fnscale=-1));
  
  ## mle = optim.out$par
  ## fim = -1 / optim.out$hessian

  ## cat("d.prev:", d.prev, "mle:", mle, "fim:", fim, "\n");

  d.new = d.prev
  
  ## Mixture of MH kernels.
  w = c(0.5, 0.5);
  ## k = sample.int(2, 1, prob=w);
  k = 1
  
  if (k==1) {
    ## Kernel 1: RW MH.
    rw.lower = max(d.prev - 1, 1);
    rw.upper = d.prev + 1;
    rw.grid  = rw.lower:rw.upper;
    rw.n     = length(rw.grid)
    rw.p     = rep(1/rw.n, rw.n);
    
    d.prop = sample(rw.grid, 1, prob=rw.p);

    ltarget = df.llh(y, d.prop, mu, G, ymax) - df.llh(y, d.prev, mu, G, ymax)
    lppsl = log(ifelse(d.prop==1, 1/2, 1/3)) - log(ifelse(d.prev==1, 1/2, 1/3));
    lratio = ltarget + lppsl
    
    if (runif(1) < exp(lratio)) d.new = d.prop
  }

  if (k==2) {
    ## Kernel 2: Ind MH.
    d.m = max(1, mle);
    d.prop = rpois(1, d.m);

    if (d.prop==0) return(d.prev);
    
    p.prop = 1 / (1 + d.prop / mu)
    p.prev = 1 / (1 + d.prev / mu)
    ltarget = df.llh(d.prop, mu, G, ymax) - df.llh(d.prev, mu, G, ymax)
    lppsl   = dpois(d.prev, d.m, log=TRUE) - dpois(d.prop, d.m, log=TRUE)
    lratio  = ltarget + lppsl

    if (runif(1) < exp(lratio)) d.new = d.prop
  }
    
  d.new
}
