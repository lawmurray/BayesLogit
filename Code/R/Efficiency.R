
update.distance <- function(beta)
{
  ## Assumes trials are by row.
  ## Folling Holmes and Held 2006.
  beta = as.matrix(beta);
  N = nrow(beta);

  d.beta = diff(beta);
  norm.d.beta = apply(d.beta, 1, function(x){ sqrt(x %*% x) });
  ave.dist = mean(norm.d.beta)

  ave.dist
}

update.distance.1d <- function(samp, relative=FALSE)
{
  ## samp is one dimensioanl.
  ## Folling Holmes and Held 2006.
  N = length(samp);

  d.samp = diff(samp);
  ave.dist = mean(abs(d.samp))

  if (relative) ave.dist = ave.dist / abs(mean(samp));
  
  ave.dist
}

## Just making sure that this is the same as acf in R.
naive.acf <- function(x, lag.max=30, type=c("correlation", "covariance"))
{
  m.x = mean(x)
  n = length(x);
  L = min(n-1, lag.max);
  ac = rep(0, L+1);
  
  for (i in 0:L) {
    idx = i + 1;
    ac[idx] = sum( (x[1:(n-i)] - m.x) * (x[(1+i):n] - m.x) ) / n;
  }

  if (type[1]=="correlation") ac = ac / ac[1];
  
  ac
}

ESS <- function(x, method=c("monotone", "positive"))
{
  ## x is one dimensional.
  ## Following Geyer 1992 and Holmes and Held 2006.
  
  acf.x = acf(x, lag.max=100, type="correlation", plot=FALSE)$acf;
  if (is.na(acf.x[1])) return(NA); # If x is constant return NA.
  
  n = length(acf.x);
  if (n %% 2 == 1) n = n - 1;
  m = n / 2;
  
  ## find last positive.
  Gamma.hat = acf.x[seq(1,n,2)] + acf.x[seq(2,n,2)];
  stay.pos = Gamma.hat > 0;
  m.pos = which.min(stay.pos) - 1;
  if (all(stay.pos)) m.pos = m;
  
  ## find last monotone
  ## d.Gamma_i = Gamma_{i+1} - Gamma_{i}, i=1..N
  d.Gamma = diff(Gamma.hat);
  stay.mon = d.Gamma < 0 ;  ## Check if decreasing.
  m.mon = which.min(stay.mon)
  if (all(stay.mon)) m.mon = m;

  ## monotone sequence estimator cutoff
  m.star = min(m.pos, m.mon);

  if (method[1]=="positive") m.star = last.pos;

  ## ## Bartlett cut-off for significant autocorrelations.
  ## bartlett = acf.x >= sqrt(2/length(x));
  ## stay.pos = bartlett > 0;
  ## n.pos = which.min(stay.pos) - 1;
  ## if (all(stay.pos)) n.pos = n;
  ## n.star = min(2*m.star, n.pos);
  
  ## ESS
  denom = 1 + 2 * sum(acf.x[2:(2*m.star)]);
  ## denom = -acf.x[1] + 2 * sum(Gamma.hat[1:m.star]);
  ## denom = 1 + 2 * sum(acf.x[2:n.star]);
  ESS = length(x) / denom;
  
  ## ## To check
  ## plot(Gamma)
  ## points(Gamma[1:m.star], col=2);
  ## plot(acf.x);
  ## points(acf.x[1:n.star], col=2);
  ## readline("ENTER")

  ESS
}

ESS.alt <- function(x)
{
  require("coda");
  effectiveSize(x)
}

convergence.heuristics <- function(samp, burn=1000)
{
  ## Assumes samples are by row.
  ## Coerce to matrix if 1-D.
  samp = as.matrix(samp);
  P = ncol(samp);
  N = nrow(samp) - burn;
  samp = samp[1:N+burn,];

  ess.nm = paste("ESS", 1:P, sep=".");
  dist.nm = paste("Dist.1d", 1:P, sep=".");
  rel.nm  = paste("Rel.Dist.1d", 1:P, sep=".");

  ess = apply(samp, 2, ESS);
  dst = update.distance(samp);
  dst.1d = apply(samp, 2, update.distance.1d);
  rel.dst.1d = apply(samp, 2, function(x){update.distance.1d(x,relative=TRUE)});

  ave.ess = mean(ess);
  ave.dst = mean(dst);
  ave.dst.1d = mean(dst.1d);
  ave.rel.dst.1d = mean(rel.dst.1d);

  ## sd.ess = sd(ess);
  ## sd.dst.1d = sd(dst.1d);
  ## sd.rel.dst.1d = sd(rel.dst.1d);

  everything = c(ave.ess, ave.dst, ave.dst.1d, ave.rel.dst.1d, ess, dst.1d, rel.dst.1d);
  e.name = c("Ave.ESS", "Ave.Ave.Dist", "Ave.Ave.Dist.1d", "Ave.Ave.Rel.Dist.1d", ess.nm, dist.nm, rel.nm);

  e.mat = matrix(everything, 1, 4 + 3*P);
  colnames(e.mat) = e.name;
  e.df = as.data.frame(e.mat);

  e.df
}

sstat <- function(x)
{
  the.names = c("Min", "Q1", "Med", "Q3", "Max", "Mean", "SD");
  out = c(min(x), quantile(x, c(0.25, 0.5, 0.75)), max(x), mean(x), sd(x));
  names(out) = the.names;
  out
}

convergence.summary <- function(lmcmc, ltime, burn=1000)
{
  ## lmcmc: list of samples.
  ## ltime: list of the times.
  trials = length(lmcmc);
  nsamp  = nrow(lmcmc[[1]]);
  nvar   = ncol(lmcmc[[1]]);
  ess    = matrix(nrow=trials, ncol=nvar);
  rtime  = rep(0, trials);
  samp.m = matrix(nrow=trials, ncol=nvar);
  samp.s = matrix(nrow=trials, ncol=nvar);
  
  ## Calculate ESS
  for (i in 1:trials) {
    the.samp = lmcmc[[i]]
    ess[i,] = apply(the.samp, 2, ESS);
    rtime[i] = ltime[[i]][1]
    samp.m[i,] = apply(the.samp, 2, mean);
    samp.s[i,] = apply(the.samp, 2, sd);
  }

  ## Calculate Summary Statistics.
  sstat.ess = t(apply(ess, 1, sstat));

  ess.sec = ess / rtime;
  sstat.sec = t(apply(ess.sec, 1, sstat));

  colnames(sstat.ess) = paste("ESS", colnames(sstat.ess), sep=".");
  colnames(sstat.sec) = paste("ESR", colnames(sstat.sec), sep=".");

  out <- list("ess"=ess, "ess.sec"=ess.sec, "sstat.ess"=sstat.ess, "sstat.sec"=sstat.sec,
              "rtime"=rtime, "nsamp"=nsamp, "nvar"=nvar, "samp.mean"=samp.m);

  out
}

convergence.table <- function(cs)
{
  sstat.ess.mean = apply(cs$sstat.ess, 2, mean);
  sstat.ess.sd   = apply(cs$sstat.ess, 2, sd);

  sstat.sec.mean = apply(cs$sstat.sec, 2, mean);
  sstat.sec.sd   = apply(cs$sstat.sec, 2, sd);

  rtime.mean = mean(cs$rtime);
  rtime.sd   = sd(cs$rtime);

  the.mean = c("rtime"=rtime.mean, sstat.ess.mean, sstat.sec.mean);
  the.sd   = c("rtime"=rtime.sd, sstat.ess.sd, sstat.sec.sd);

  out = rbind(the.mean, the.sd);
  rownames(out) = c("Mean", "SD");

  out
}

################################################################################

if (FALSE) {

  ## Setup
  
  prefix = c("BL", "HH", "FSF");
  ## prefix = c("BL", "Pro");
  data.name = c("diabetes", "heart", "australia", "germany");

  ## Storage
  ss = list(
    diabetes = list(beta.ave.ave=list(), beta.ave.sd=list()),
    heart = list(beta.ave.ave=list(), beta.ave.sd=list()),
    aus = list(beta.ave.ave=list(), beta.ave.sd=list()),
    ger = list(beta.ave.ave=list(), beta.ave.sd=list())
    );

  master.table = matrix(0, 4, 8);
  rownames(master.table) = data.name;
  colnames(master.table) = c("PG.ESS", "PG.Dist",
                             "HH.ESS", "HH.Dist",
                             "FSF.ESS", "FSF.Dist",
                             "BL.to.HH", "BL.to.FSF");

  time.table = matrix(0, 4, 8);
  rownames(time.table) = data.name;
  colnames(time.table) = c("PG.Time", "HH.Time", "FSF.Time",
                           "PG.ESS.per.Time", "HH.ESS.per.Time", "FSF.ESS.per.TIME",
                           "PG.HH.ratio", "PG.FSF.ratio");
  
}

## CALCULATE MASTER TABLE FOR ESS OF BETA using convergence.heuristics ##

if (FALSE) {

  ## Now calculate everything!

  for (j in c(1,2)) { ## method
    for (k in 1:4) { ## data

      cat("Working on", prefix[j], data.name[k], "\n");
      
      file.name = paste(prefix[j], data.name[k], "RData", sep=".");
      load(paste("BetaCompareData/Test02/", file.name, sep=""));
      ## load(paste("BetaCompareData/Test03/", file.name, sep=""));
      
      ch.list = lapply(the.beta, convergence.heuristics)
      beta.ave.list = lapply(the.beta, colMeans);
      
      ch.df = ch.list[[1]];
      beta.ave = beta.ave.list[[1]];
      time.df = the.time[[1]];
      
      N = length(ch.list);
      for (i in 2:N) {
        beta.ave = rbind(beta.ave, beta.ave.list[[i]]);
        ch.df = rbind(ch.df, ch.list[[i]]);
        time.df = rbind(time.df, the.time[[i]]);
      }

      cm.df = colMeans(ch.df);
      beta.ave.ave = apply(beta.ave, 2, mean);
      beta.ave.sd  = apply(beta.ave, 2, sd);

      master.table[k, 2*(j-1)+1:2] = cm.df[1:2];
      ss[[k]]$beta.ave.ave[[j]] = beta.ave.ave
      ss[[k]]$beta.ave.sd[[j]]  = beta.ave.sd

      ave.time = apply(time.df, 2, mean);
      time.table[k,j] = ave.time[1];
      time.table[k,j+3] = cm.df[1] / ave.time[1];

    }
  }

  master.table[,7] = master.table[,1] / master.table[,3];
  master.table[,8] = master.table[,1] / master.table[,5];

  time.table[,7] = time.table[,4] / time.table[,5];
  time.table[,8] = time.table[,4] / time.table[,6];

  ## ## Check coefficients
  ## for (k in 1:4) {
  ##   plot(ss[[k]]$beta.ave.ave[[1]], main=data.name[k], xlab="Index", ylab="Beta");
  ##   for(j in 2:3) points(jitter(ss[[k]]$beta.ave.ave[[j]], amount=0.01), col=j);
  ##   readline("Press <ENTER>...");
  ## }
      
}

################################################################################

## ESS FOR PROB ##

if (FALSE) {

  data.name = c("diabetes", "heart", "australia", "germany");
  
  load("~/RPackage/BayesLogit/Code/R/DataSets/australia.RData")
  load("~/RPackage/BayesLogit/Code/R/DataSets/germany.RData")
  load("~/RPackage/BayesLogit/Code/R/DataSets/heart.RData")
  load("~/RPackage/BayesLogit/Code/R/DataSets/diabetes.RData")
  
  all.data = list(
    diabetes = list("y"=y.diabetes, "X"=X.diabetes),
    heart = list("y"=y.heart, "X"=X.heart),
    aus = list("y"=y.aus, "X"=X.aus),
    ger = list("y"=y.ger, "X"=X.ger)
    )
  
  ## prefix = c("BL", "HH", "FSF");
  prefix = c("BL", "Pro");
  
  ## Now calculate everything!

  if (!exists("ess")) ess = matrix(nrow=4, ncol=2);
  colnames(ess) = prefix;
  rownames(ess) = data.name;
  
  for (j in c(1,2)) { ## method
    for (k in 1:4) {
      
      cat("Working on", prefix[j], data.name[k], "\n");
      
      file.name = paste(prefix[j], data.name[k], "RData", sep=".");
      ## load(paste("BetaCompareData/Test02/", file.name, sep=""));
      load(paste("BetaCompareData/Test03/", file.name, sep=""));
      
      ave.ess = rep(0, 10);

      ## If using log odds.
      for (i in 1:10) {
        X = all.data[[k]]$X;
        psi = the.beta[[i]][1001:10000,] %*% t(X);
        if (j==1) prob = exp(psi) / (1 + exp(psi))
        else prob = pnorm(psi, 0.0, 1.0);
        ess.vec = apply(prob, 2, ESS)
        num.na = sum(is.na(ess.vec))
        if (num.na > 0) cat("Warning: ", num.na, "NA\n");
        ave.ess[i] = mean(ess.vec, na.rm=TRUE);
      }
      
      ess[k,j] = mean(ave.ess);
    }
  }

}

################################################################################

if (FALSE) {

  ## Now calculate everything!
  prefix = c("BL", "FSF", "HH");
  data.name = c("diabetes", "heart", "australia", "germany");

  the.tables = list();
  means.betas = list();
  sd.means.betas = list();
  rtime.list = list();
  
  for (k in 1:4) { ## data

    temp.table = c();
    beta.means = c();
    rtime.all  = c();
    
    for (j in c(1,2,3)) { ## method

      cat("Working on", prefix[j], data.name[k], "\n");
      
      file.name = paste(prefix[j], data.name[k], "RData", sep=".");
      load(paste("BetaCompareData/Test04/", file.name, sep=""));

      cs = convergence.summary(the.beta, the.time);
      ct = convergence.table(cs);

      rownames(ct) = c(paste("Mean", prefix[j], sep="."), paste("SD", prefix[j], sep="."));

      temp.table = rbind(temp.table, ct);
      beta.means = rbind(beta.means, apply(cs$samp.mean, 2, mean));
      rtime.all  = rbind(rtime.all, cs$rtime)
      
    }

    rownames(beta.means) = prefix;
    rownames(rtime.all)  = prefix;

    the.tables[[data.name[k]]] = temp.table
    means.betas[[data.name[k]]] = beta.means
    sd.means.betas[[data.name[k]]] = apply(beta.means, 2, sd);
    rtime.list[[data.name[k]]] = rtime.all;
    
  }
      
}

################################################################################
                          ## USING UPDATED VERSION ##
################################################################################

if (FALSE) {

  meth = c("PG", "FS", "HH", "DFAE");
  dset  = c("australia", "germany", "diabetes", "heart");

  n.meth = length(meth);
  n.dset  = length(dset);

  mat.summary = c();
  cs.list = list();

  i = 4

  df = NULL;
  for (j in 1:2) {
    filename = paste("BetaCompareData/Test05/logit", meth[j], dset[i], "RData", sep=".");
    load(filename)
    ## cs = cs.list[[j]];
    cs = convergence.summary(the.beta, ess.time);
    cs.list[[j]] = cs;
    temp = c("method"=j, "nvar"=cs$nvar, "samp"=cs$nsamp, colMeans(cs$sstat.ess), colMeans(cs$sstat.sec))
    if (is.null(df)) {
      df = as.data.frame(t(temp))
    } else {
      df = rbind(df, t(temp))
    }
  }
  
}
