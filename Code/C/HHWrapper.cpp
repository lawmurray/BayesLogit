
#include "HHWrapper.hpp"
#include "RNG.hpp"
// #include "/Users/jwindle/RV-Project/Code/C_Examples/MyLib/RNG/RNG.hpp"
#include <cmath>
#include <iostream>

#ifdef USE_R
#include "R.h"
#include "Rmath.h"
#include <R_ext/Utils.h>
#endif

const double __PI = 3.141592653589793238462643383279502884197;
const double MAX_ITER = 65536;

bool right(double U, double lambda, double r)
{
  double Z = 1.0;
  double X = exp(-0.5 * lambda);
  int j = 0;

  int iter = 0;
  while (iter < MAX_ITER) {
    j = j + 1;
    Z = Z - pow((j+1), 2) * pow(X, (pow(j+1, 2) -1));
    // You want to exit the loop if Z==U.
    if (Z >= U) return true;
    j = j + 1;
    Z = Z + pow((j+1), 2) * pow(X, (pow(j+1, 2) -1));
    if (Z < U) return false;

    // Break infinite loop.
    #ifdef USE_R
    if (++iter%10==0) R_CheckUserInterrupt();
    #endif

  }

  std::cout << "Right: failed to accept/reject after " << iter << " iterations.\n";
  std::cout << "U: " << U 
	    << " lambda: " << lambda
	    << " Z: " << Z 
	    << " r: " << r << "\n";

  return false;
}

bool left(double U, double lambda, double r)
{
  double H = 0.5 * log(2.0) + 2.5 * log(__PI) 
    - 2.5 * log(lambda) - pow(__PI, 2) / 2 / lambda 
    + 0.5 * lambda;
  
  // To guard against lambda = 0, which above produces H = -nan.
  if (lambda<=1e-100) {
    H = -4.93e+100;
    std::cout << "Warning: lambda = " << lambda << " <= 1e-100.\n";
  }

  double log_U = log(U);
  double Z = 1;
  double X = exp(-pow(__PI, 2) / (2 * lambda));
  double K = lambda / pow(__PI, 2);
  int j = 0;

  int iter = 0;
  while (iter < MAX_ITER) {
    j = j + 1;
    Z = Z - K * pow(X, (pow(j,2)-1));
    // You want to exit the loop if Z==U.
    if (H + log(Z) >= log_U) return true;
    j = j + 1;
    Z = Z + pow((j+1), 2) * pow(X, (pow(j+1,2) - 1 ));
    if (H + log(Z) < log_U) return false;

    // Break infinite loop.
    #ifdef USE_R
    if (++iter%10==0) R_CheckUserInterrupt();
    #endif
  }

  // There appears to be a problem when r = z - X beta = 0.
  // In that case, we draw extremely small lambda.

  std::cout << "Left: failed to accept/reject after " << iter << " iterations.\n";
  std::cout << "U: " << U 
	    << " lambda: " << lambda
	    << " Z: " << Z 
	    << " H: " << H
	    << " r: " << r << "\n";

  return false;
}

void hh_lambda(double *lambdap, double *rp, int *iterp)
{
  double lambda = 0.0;
  double r = fabs(*rp);
  int iter   = 0;

  bool ok = false;

  RNG rng;

  double Y, U, Z;

  #ifdef USE_R
  GetRNGstate();
  #endif

  while (!ok) {
    
    // Propose lambda ~ GIG(0.5, 1, r^2)
    Y = pow(rng.norm(0.0, 1.0), 2);
    // Unstable calculation -> You can produce negative numbers.
    // Y = 1 + (Y - sqrt(Y * (4*r + Y))) / (2 * r);
    Z = 1 + Y / (2 * r);
    Y = Z + sqrt(Z*Z-1);
    U = rng.unif();
    lambda = r * Y;
    if (U > 1/(1+Y)) lambda = r/Y;

    U = rng.unif();
    if (lambda > 4/3) {
      ok = right(U, lambda, r);
    }
    else {
      ok = left(U, lambda, r);
    }

    iter = iter + 1;

    #ifdef USE_R
    if (iter%10000==0) R_CheckUserInterrupt();
    #endif

  }

  #ifdef USE_R
  PutRNGstate();
  #endif

  *lambdap = lambda;
  *iterp   = iter;
}


double hh_lambda(double r, int *iterp)
{
  double lambda;
  hh_lambda(&lambda, &r, iterp);
  return lambda;
}

void hh_lambda_vec(double *lambda, double *r, int *n)
{
  int iter;

  for(int i = 0; i < *n; i++) {
    #ifdef USE_R
    if (iter%10==0) R_CheckUserInterrupt();
    #endif
    hh_lambda(lambda+i, r+i, &iter);
    // lambda[i] = hh_lambda(r[i], &iter);
  }

}

////////////////////////////////////////////////////////////////////////////////
				 // APPENDIX //
////////////////////////////////////////////////////////////////////////////////

// I am not sure if this is the absolute fastest way to do things.  I mean the
// following.  Currently, I am compiling using g++.  But I think that it may be
// better to use R CMD COMPILE and R CMD SHLIB.  I'm not sure if R does some
// sort of internal checking when you use g++.  Otherwise you could presumably
// do something nasty.  I did try compiling the above code using R CMD SHLIB
// (you need to define USE_R and then uncomment the long version of RNG.hpp).
// But that didn't seem to help that much.
