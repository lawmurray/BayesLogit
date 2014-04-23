#include <getpot>
#include <vector>
#include <sys/time.h>
#include "ParallelWrapper.h"
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

    int    nthreads = args.follow(1  , "--threads");
    int    samp     = args.follow(100, "--samp");
    int    reps     = args.follow(1  , "--reps");

    double shape    = args.follow(1.0, "--shape");
    double tilt     = args.follow(0.0, "--tilt");
    
    fprintf(stderr, "Will draw %i samples %i times using %i threads.\n", samp, reps, nthreads);
    fprintf(stderr, "Shape: %g; Tilt: %g\n", shape, tilt);

    vector<double> x(samp);
    vector<double> h(samp);
    vector<double> z(samp);

    fill(h.begin(), h.end(), shape);
    fill(z.begin(), z.end(), tilt);

    vector<RNG> r(nthreads);
    vector<PolyaGamma> dv(nthreads);
    vector<PolyaGammaAlt> alt(nthreads);
    vector<PolyaGammaSP> sp(nthreads);

    struct timeval start, stop;
    gettimeofday(&start, NULL);

    for (int i = 0; i < reps; i++)
	rpg_hybrid_par(&x[0], &h[0], &z[0], &samp, &nthreads, &r, &dv, &alt, &sp);

    gettimeofday(&stop, NULL);
    double diff = calculateSeconds(start, stop);
    fprintf(stderr, "Time: %f sec. for %i samp (%i times)\n", diff, samp, reps);

    // Check output
    double m1_hat = accumulate(x.begin(), x.end(), 0.0) / samp;
    double m2_hat = accumulate(x.begin(), x.end(), 0.0, addsq) / samp;

    fprintf(stderr, "Sample moments: m1: %g; m2: %g\n", m1_hat, m2_hat);
    fprintf(stderr, "Actual moments: m1: %g, m2: %g\n", dv[0].pg_m1(shape, tilt), dv[0].pg_m2(shape, tilt));

}
