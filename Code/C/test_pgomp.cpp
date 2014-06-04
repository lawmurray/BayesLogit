// -*- mode: c++; c-basic-offset: 4; -*-

#include "PolyaGammaOMP.h"
#include "RNG.h"
#include <algorithm>
#include <getpot>
#include <vector>
#include <sys/time.h>
#include <numeric>
#include <stdio.h>

using std::vector;
using std::fill;
using std::accumulate;

double calculateSeconds(const timeval &time1, const timeval &time2) {
    return time2.tv_sec - time1.tv_sec + (double)(time2.tv_usec - time1.tv_usec) / 1000000.0;
}

double addsq(double x, double y) {return x + y*y;}

int main(int argc, char** argv)
{
    GetPot args(argc, argv);

    int    nthread  = args.follow(1  , "--nthread");
    int    samp     = args.follow(100, "--samp");
    int    reps     = args.follow(1  , "--reps");

    double shape    = args.follow(1.0, "--shape");
    double tilt     = args.follow(0.0, "--tilt");

    vector<RNG> rngs(nthread);
    PolyaGammaOMP<double> pgomp(nthread);
    
    vector<double> x(samp);
    vector<double> b(samp);
    vector<double> z(samp);

    fill(b.begin(), b.end(), shape);
    fill(z.begin(), z.end(), tilt);

    struct timeval start, stop;
    gettimeofday(&start, NULL);

    for (int i = 0; i < reps; i++) {
	// fprintf(stderr, "Iteration %i\n", i);
	pgomp.draw_hybrid(&x[0], &b[0], &z[0], samp, rngs);
    }

    gettimeofday(&stop, NULL);
    double diff = calculateSeconds(start, stop);
    fprintf(stderr, "Time: %f sec. for %i samp (%i times) using %i threads.\n", diff, samp, reps, nthread);

    // Check output
    double m1_hat = accumulate(x.begin(), x.end(), 0.0) / samp;
    double m2_hat = accumulate(x.begin(), x.end(), 0.0, addsq) / samp;

    PolyaGamma dv;

    fprintf(stderr, "Sample moments: m1: %g; m2: %g\n", m1_hat, m2_hat);
    fprintf(stderr, "Actual moments: m1: %g, m2: %g\n", dv.pg_m1(shape, tilt), dv.pg_m2(shape, tilt));

    return 0;
}
