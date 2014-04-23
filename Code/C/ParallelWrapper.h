////////////////////////////////////////////////////////////////////////////////

// Copyright 2014 Nick Polson, James Scott, and Jesse Windle.

// This file is part of BayesLogit.

// BayesLogit is free software: you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or any later version.

// BayesLogit is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
// A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along with
// BayesLogit.  If not, see <http://www.gnu.org/licenses/>.

////////////////////////////////////////////////////////////////////////////////

#include "RNG.hpp"
#include "PolyaGamma.h"
#include "PolyaGammaAlt.h"
#include "PolyaGammaSP.h"
#include <vector>

using std::vector;

void rpg_hybrid_par(double *x, 
		    double *h, 
		    double *z, 
		    int* num, 
		    int* nthreads, 
		    vector<RNG> *r, 
		    vector<PolyaGamma> *dv, 
		    vector<PolyaGammaAlt> *alt, 
		    vector<PolyaGammaSP> *sp);
