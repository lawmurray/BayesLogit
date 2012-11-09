
if (is.loaded("BayesLogit.so"))  { dyn.unload("BayesLogit.so"); dyn.load("BayesLogit.so"); }
if (!is.loaded("BayesLogit.so")) { dyn.load("../C/BayesLogit.so"); }

source("LogitWrapper.R")
source("FFBS.R")
source("NB-Indicators.R");
source("Logit-Indicators.R")
