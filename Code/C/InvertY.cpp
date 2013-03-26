#include "InvertY.hpp"

//------------------------------------------------------------------------------

double y_eval(double v)
{
  double y   = 0.0;
  double r   = sqrt(fabs(v));
  if (v > tol)
    y = tan(r) / r;
  else if (v < -1*tol)
    y = tanh(r) / r;
  else
    y = 1 + (1/3) * v + (2/15) * v * v + (17/315) * v * v * v;
  return y;
}

void ydy_eval(double v, double* yp, double* dyp)
{
  // double r   = sqrt(fabs(v));

  double y = y_eval(v);
  *yp = y;

  if (fabs(v) >= tol) 
    *dyp = 0.5 * (y*y + (1-y) / v);
  else
    *dyp = 0.5 * (y*y - 1/3 - (2/15) * v);

}

double f_eval(double v, void * params)
{
  double y = *((double*) params);
  return y_eval(v) - y;
}

void fdf_eval(double v, void* params, double* fp, double* dfp)
{
  double y = *((double*)params);
  ydy_eval(v, fp, dfp);
  *fp  -= y;
}

double df_eval(double v, void * params)
{
  double f, df;
  ydy_eval(v, &f, &df);
  return df;
}

//------------------------------------------------------------------------------

YV::YV() : T(gsl_root_fdfsolver_newton), s(NULL)
{
  s = gsl_root_fdfsolver_alloc (T);

  FDF.f   = &f_eval;
  FDF.df  = &df_eval;
  FDF.fdf = &fdf_eval;
  FDF.params = 0;
  
  ylower = y_eval(vlower);
  yupper = y_eval(vupper);
}

YV::~YV() { 
  if (s != NULL) { 
    gsl_root_fdfsolver_free(s); 
    s = NULL; 
  } 
}

double YV::y_func(double v) {
  return y_eval(v);
}

double YV::v_func(double y, int maxiter) 
{
  double v = 1.0;

  double ycopy = y;
  FDF.params = &ycopy;

  if (y < ylower) {
    return -1. / (y*y);
  } else if (y > yupper) {
    v = atan(0.5 * y * IYPI);
    return v*v;
  }
    
  double id = (log(y) / log(2) + 4.0) / 0.1;
  // printf("y, id, y[id], v[id]: %g, %g, %g, %g\n", y, id, ygrid[(int)id], vgrid[(int)id]);
  
  gsl_root_fdfsolver_set(s, &FDF, vgrid[(int)id]);

  int iter = 0;
  int status = 0;
  double vp = 0.0;
  // double fval, dfval;

  do {
    iter++;
    status = gsl_root_fdfsolver_iterate (s);
    vp = v;
    v  = gsl_root_fdfsolver_root(s);
    status = gsl_root_test_delta(v, vp, 0, 1e-8);

    // ydy_eval(v, &fval, &dfval);
    // printf("yval, dyval, v: %g, %g, %g\n", fval, dfval, v);

    // fdf_eval(v, FDF.params, &fval, &dfval);
    // printf("fval, dfval: %g, %g\n", fval, dfval);

  } while (status == GSL_CONTINUE && iter < maxiter);

  if (iter >= maxiter) fprintf(stderr, "YV: v reached maxiter.\n");

  return v;
}

double YV::upperIncompleteGamma(double x, double shape, double rate)
{
  double t = rate * x;
  return gsl_sf_gamma_inc(shape, t);
}
