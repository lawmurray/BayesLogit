library("coda")

## Functions needed for benchmarking.

################################################################################
                                 ## SUM STAT ##
################################################################################

## Generates summary statistics for static coefficient.
sum.stat <- function(beta, runtime, thin=1)
{
  beta = as.matrix(beta)
  n = nrow(beta);
  p = ncol(beta);
  beta = beta[seq(1, n, thin),];
  if (p==1) {beta = as.matrix(beta);}
  sstat = matrix(nrow=p, ncol=7);

  sstat[,1] = apply(beta, 2, mean);
  sstat[,2] = apply(beta, 2, sd);
  sstat[,3] = apply(beta, 2, effectiveSize);
  sstat[,4] = sstat[,3] / runtime;
  ## sstat[,5] = apply(beta, 2, ESS);
  ## sstat[,6] = sstat[,5] / runtime;
  sstat[,5] = sstat[,1] / sstat[,2];
  sstat[,6] = apply(beta, 2, function(x){quantile(x, 0.1)});
  sstat[,7] = apply(beta, 2, function(x){quantile(x, 0.9)});

  ## colnames(sstat) = c("mean", "sd", "ESS", "ESS.sec", "myESS", "myESS.sec", "t", "Q.1", "Q.9");
  colnames(sstat) = c("mean", "sd", "ESS", "ESR", "t", "Q.1", "Q.9");

  sstat
}

## Generates summary statistics for dynamic coefficient.
sum.stat.dyn <- function(beta, runtime, thin=1)
{
  beta.dim = dim(beta);
  n = nrow(beta);
  p = ncol(beta);
  if (length(beta.dim)==2) beta = as.array(beta, dim=c(n, 1, p));
                             
  n = nrow(beta);
  p = ncol(beta);
  T = dim(beta)[3];
  beta = beta[seq(1, n, thin),,];
  if (p==1) {beta = array(beta, dim=c(nrow(beta), 1, T));}
  sstat = array(0, dim=c(p, T, 7));

  sstat[,,1] = apply(beta, c(2,3), mean);
  sstat[,,2] = apply(beta, c(2,3), sd);
  sstat[,,3] = apply(beta, c(2,3), effectiveSize);
  sstat[,,4] = sstat[,,3] / runtime;
  ## sstat[,5] = apply(beta, 2, ESS);
  ## sstat[,6] = sstat[,5] / gbs$runtime;
  sstat[,,5] = sstat[,,1] / sstat[,,2];
  sstat[,,6] = apply(beta, c(2,3), function(x){quantile(x, 0.1)});
  sstat[,,7] = apply(beta, c(2,3), function(x){quantile(x, 0.9)});

  ## colnames(sstat) = c("mean", "sd", "ESS", "ESS.sec", "myESS", "myESS.sec", "t", "Q.1", "Q.9");
  dimnames(sstat)[[3]] = c("mean", "sd", "ESS", "ESR", "t", "Q.1", "Q.9");

  sstat
}

################################################################################
                                   ## PLOT ##
################################################################################

plot.bench <- function(bmark1, bmark2=NULL)
{
  P.b = dim(bmark1$gb$beta)[2];
  P.a = dim(bmark1$gb$iota)[2];
  T   = dim(bmark1$gb$beta)[3] - 1;

  num.plots = 1 + !is.null(bmark2);
  par(mfrow=c(1,num.plots));
  
  ## Static Coefficients.
  for(i in 1:P.a) {
    for (n in 1:num.plots) {
      if (n==2) bmark=bmark2
      else bmark=bmark1
      hist(bmark$gb$iota[,i]);
    }
    readline("Press <ENTER> for next...");
  }
  
  ## Dynamic Coefficients.
  for(i in 1:P.b) {
    for (n in 1:num.plots) {
      if (n==2) bmark=bmark2
      else bmark=bmark1
      
      beta.mean = bmark$sstat.beta[i,,"mean",1]
      beta.9    = bmark$sstat.beta[i,,"Q.9" ,1]
      beta.1    = bmark$sstat.beta[i,,"Q.1" ,1]
      
      ymin = min(beta.mean, bmark$gb$beta[,i,]);
      ymax = max(beta.mean, bmark$gb$beta[,i,]);
      
      plot(0:T, beta.mean, col=2, type="l", ylim=c(ymin,ymax));
      lines(0:T, beta.9, col="pink")
      lines(0:T, beta.1, col="pink")
      abline(h=mean(beta.mean), col=2, lty=c(2,2));

    }
    readline("Press <ENTER> for next...");
  }

} ## plot.bench

##------------------------------------------------------------------------------
plot.check.logit <- function(y, X.dyn, n, X.stc=NULL, bmark1, bmark2=NULL)
{
  X   = cbind(X.stc, X.dyn);
  P   = ncol(X);
  P.b = ncol(X.dyn);
  P.a = P - P.b;
  T   = nrow(X);
  
  num.plots = 1 + !is.null(bmark2);
  par(mfrow=c(1,num.plots));

  for (j in 1:num.plots) {
    if (j==2) bmark=bmark2
    else bmark=bmark1

    beta.mean = matrix(bmark$sstat.beta[,,"mean",1], nrow=P.b);
    
    lodds = apply(X.dyn * t(beta.mean)[-1,], 1, sum);
    if (P.a > 0) lodds = lodds + X.stc %*% iota.mean;
    
    ymin = min(lodds);
    ymax = max(lodds);
    
    plot(1:T, lodds, ylim=c(ymin, ymax), col=2, type="l", main="log odds");
    points(1:T, (y/n)*(ymax-ymin) + ymin, cex=0.1);
  }
  
} ## plot.check.logit

##------------------------------------------------------------------------------
plot.check.NB <- function(y, X.dyn, X.stc=NULL, bmark1, bmark2=NULL)
{
  X   = cbind(X.stc, X.dyn);
  P   = ncol(X);
  P.b = ncol(X.dyn);
  P.a = P - P.b;
  T   = nrow(X);
  
  num.plots = 1 + !is.null(bmark2);
  par(mfrow=c(1,num.plots));

  for (n in 1:num.plots) {
    if (n==2) bmark=bmark2
    else bmark=bmark1

    beta.mean = matrix(bmark$sstat.beta[,,"mean",1], nrow=P.b);
    
    lmean = apply(X.dyn * t(beta.mean)[-1,], 1, sum);
    if (P.a > 0) lmean = lmean + X.stc %*% iota.mean;
    
    ymin = min(lmean);
    ymax = max(lmean);
    
    plot(1:T, lmean, ylim=c(ymin, ymax), col=2, type="l", main="log mean and log(y)");
    points(1:T, log(y), cex=0.1);
  }
  
} ## plot.check.NB

##------------------------------------------------------------------------------

setup.table <- function(b.out, var.name="beta") {

  ave.sstat = lapply(b.out, function(x) apply(x$sstat[[var.name]], c(1,2), mean));
  ave.sstat = simplify2array(ave.sstat);

  ave.time  = lapply(b.out, function(x) mean(x$ess.time));
  ave.time  = simplify2array(ave.time)

  ## min, median, max
  mmm.ess = apply(ave.sstat[,3,], 2, function(x) quantile(x, c(0.0, 0.5, 1.0)) )
  mmm.esr = apply(ave.sstat[,4,], 2, function(x) quantile(x, c(0.0, 0.5, 1.0)) )

  ave.arate = simplify2array(lapply(b.out, function(x) x$info$ave.arate));
  
  the.table = rbind(ave.time, ave.arate, mmm.ess, mmm.esr);
  rownames(the.table) = c("time", "ARate", "ESS.min", "ESS.med", "ESS.max", "ESR.min", "ESR.med", "ESR.max");

  ## ave.alpha.indmh = mean(b.out[["IndMH"]]$sstat$alpha[1,1,])
  ## ave.alpha.ram = mean(b.out[["RAM"]]$sstat$alpha[1,1,])

  ## I don't always have IndMH defined.
  ## eval = eigen(-1*b.out[["IndMH"]]$gb$hess)$values;
  ## cn.P = eval[1] / eval[length(eval)];

  out <- list("ave.sstat"=ave.sstat,
              ## "eval.V"=eval, "cn.P"=cn.P, "cn.X"=b.out[["IndMH"]]$info$cn.X,
              "table"=t(the.table));

  out
}

##------------------------------------------------------------------------------

setup.table.dyn <- function(b.out, var.name="beta") {

  ave.sstat = lapply(b.out, function(x) apply(x$sstat[[var.name]], c(1,2,3), mean));
  ave.sstat = simplify2array(ave.sstat);

  ave.time  = lapply(b.out, function(x) mean(x$ess.time));
  ave.time  = simplify2array(ave.time)

  ## min, median, max
  mmm.ess = drop(apply(ave.sstat[,,3,,drop=FALSE], c(3,4), function(x) quantile(x, c(0.0, 0.5, 1.0)) ))
  mmm.esr = drop(apply(ave.sstat[,,4,,drop=FALSE], c(3,4), function(x) quantile(x, c(0.0, 0.5, 1.0)) ))

  ## ave.arate = simplify2array(lapply(b.out, function(x) x$info$ave.arate));
  ave.arate = 1
  
  the.table = rbind(ave.time, ave.arate, mmm.ess, mmm.esr);
  rownames(the.table) = c("time", "ARate", "ESS.min", "ESS.med", "ESS.max", "ESR.min", "ESR.med", "ESR.max");

  out <- list("ave.sstat"=ave.sstat,
              "table"=t(the.table));

  out
}
