
if (is.loaded("BayesLogit.so"))  { dyn.unload("BayesLogit.so"); dyn.load("BayesLogit.so"); }
if (!is.loaded("BayesLogit.so")) { dyn.load("../C/BayesLogit.so"); }

source("LogitWrapper.R")
source("FFBS.R")
source("CUBS.R")
source("AR1.R")
source("Stationary.R")
source("Indicators.R");
source("Logit-Indicators.R")
