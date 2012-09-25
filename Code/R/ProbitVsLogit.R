
library("BayesLogit", lib.loc="~/RPackage/BayesLogit/Code/BLPackage/Test/");
source("~/RPackage/BayesLogit/Code/R/logit-MCMC.R")
source("~/RPackage/BayesLogit/Code/R/Probit.R")
source("~/RPackage/BayesLogit/Code/R/Efficiency.R")

y = y.heart
X = as.matrix(X.heart)

start.time = proc.time();
out.log = logit.gibbs.np.R(y, X, samp=10000, burn=0, verbose=1000);
log.time = proc.time() - start.time;

start.time = proc.time();
out.pro = probit(y, X, samp=10000, burn=0, verbose=1000);
pro.time = proc.time() - start.time;

ess.log = apply(out.log$beta, 2, ESS);
ess.pro = apply(out.pro$beta, 2, ESS);

cat("Logit: Time", log.time[1], "AVE ESS", mean(ess.log), "\n");
cat("Probit: Time", pro.time[1], "AVE ESS", mean(ess.pro), "\n");

ratio = (mean(ess.log) / log.time[1]) / (mean(ess.pro) / pro.time[1]);
cat("Log ESS/sec to Pro ESS/sec:", ratio, "\n");

par(mfrow=c(2,4));

for (i in 1:4) {
  hist(out.log$beta[,i], breaks=40, xlab=paste("beta", i), main=paste("Logit: pred", i));
}

for (i in 1:4) {
  hist(out.pro$beta[,i], breaks=40, xlab=paste("beta", i), main=paste("Probit: pred", i));
}

## Standardize?
x = out.pro$beta[,1];
x = (x - mean(x)) / sd(x);


