## Copyright 2013 Nick Polson, James Scott, and Jesse Windle.

## This file is part of BayesLogit, distributed under the GNU General Public
## License version 3 or later and without ANY warranty, implied or otherwise.


if (is.loaded("BayesLogit.so"))  { dyn.unload("BayesLogit.so"); dyn.load("BayesLogit.so"); }
if (!is.loaded("BayesLogit.so")) { dyn.load("../C/BayesLogit.so"); }

source("LogitWrapper.R")
source("FFBS.R")
source("CUBS.R")
source("AR1.R")
source("Stationary.R")
source("Indicators.R");
source("Logit-Indicators.R")
source("DynExpFamMHWrapper.R")
