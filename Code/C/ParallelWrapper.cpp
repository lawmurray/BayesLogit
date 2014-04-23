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
#include <exception>
#include <stdio.h>
#include <omp.h>

void rpg_hybrid_par(double *x, 
		    double *h, 
		    double *z, 
		    int* num, 
		    int* nthreads, 
		    vector<RNG> *r, 
		    vector<PolyaGamma> *dv, 
		    vector<PolyaGammaAlt> *alt, 
		    vector<PolyaGammaSP> *sp)
{
#ifdef USE_R
    printf("Currently, you must use GSL for parallel draw.\n");
#endif

    int i, tid;

#pragma omp parallel shared(num, h, z, r, dv, alt, sp) private(i, tid) num_threads(*nthreads)
    {
        // // Get thread number and write out.
        tid = omp_get_thread_num();
        // printf("Hello from thread %i.\n", tid);
        // // Check if master thread.
        // if (tid==0) {
        //   nthreads = omp_get_num_threads();
        //   printf("There are %i threads.\n", nthreads);
        // }

        // Now let's do this in parallel.

        #pragma omp for schedule(dynamic) nowait
        for(int i=0; i < *num; ++i){
            double b = h[i];
            if (b > 170) {
                double m = (*dv)[tid].pg_m1(b,z[i]);
                double v = (*dv)[tid].pg_m2(b,z[i]) - m*m;
                x[i] = (*r)[tid].norm(m, sqrt(v));
            }
            else if (b > 13) {
                (*sp)[tid].draw(x[i], b, z[i], (*r)[tid]);
            }
            else if (b==1 || b==2) {
                x[i] = (*dv)[tid].draw((int)b, z[i], (*r)[tid]);
            }
            else if (b > 1) {
                x[i] = (*alt)[tid].draw(b, z[i], (*r)[tid]);
            }
            else if (b > 0) {
                x[i] = (*dv)[tid].draw_sum_of_gammas(b, z[i], (*r)[tid]);
            }
            else {
                x[i] = 0.0;
            }
        }

    }

}
