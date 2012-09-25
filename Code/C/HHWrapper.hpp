
#ifndef __HHWRAPPER__
#define __HHWRAPPER__

extern "C" {

  void hh_lambda(double *lambda, double *r, int *iter);

  void hh_lambda_vec(double *lambda, double *r, int *n);

}

double hh_lambda(double r, int *iter);

bool right(double U, double lambda, double r);

bool left(double U, double lambda, double r);

#endif
