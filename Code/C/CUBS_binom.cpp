#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_sf.h>
#include <stdio.h>
#include <cstring>

void binom_transform (const double* rs, const double* fq, double* out)
{
  double r=rs[0]; double s=rs[1];
  double f=fq[0]; double q=fq[1];
  double E = gsl_sf_psi(r) - gsl_sf_psi(s);
  double V = gsl_sf_psi_1(r) + gsl_sf_psi_1(s);
  out[0] = E - f;
  out[1] = V - q;
  // printf("r=%g, s=%g, f=%g, q=%g, out=%g, %g\n", r, s, f, q, out[0], out[1]);
}

void utest_binom_transform(double r, double s)
{
  double rs[2]; rs[0] = r; rs[1] = s;
  double fq[2] = {0, 0};
  double out[2];

  binom_transform(rs, fq, out);
  printf("r=%g, s=%g, f=%g, q=%g\n", rs[0], rs[1], out[0], out[1]);
}

int binom_transform_gsl (const gsl_vector* x, void* p, gsl_vector* f) {
  double rs[2];
  double out[2];
  double* fq = (double *)p;

  rs[0] = gsl_vector_get(x, 0);
  rs[1] = gsl_vector_get(x, 1);

  binom_transform(rs, fq, out);

  gsl_vector_set (f, 0, out[0]);
  gsl_vector_set (f, 1, out[1]);

  return GSL_SUCCESS;
}

void utest_binom_transform_gsl(double r, double s)
{
  gsl_vector* x = gsl_vector_alloc(2);
  x->data[0] = r;
  x->data[1] = s;

  double fq[2] = {0, 0};
  void* fqv = (void *)fq;

  gsl_vector* out = gsl_vector_alloc(2);
  binom_transform_gsl(x, fqv, out);

  printf("r=%g, s=%g, f=%g, q=%g\n", r, s, out->data[0], out->data[1]);
}

void binom_solver(const double* fq, double* rs, double epsabs, double epsrel, int max_iter)
{
  double params[2]; memmove(params, fq, 2 * sizeof(double));
  // fq[0] = prior[0]; fq[1] = prior[1];

  const gsl_multiroot_fsolver_type * T = gsl_multiroot_fsolver_hybrid;
  gsl_multiroot_fsolver            * s = gsl_multiroot_fsolver_alloc(T, 2);

  gsl_multiroot_function F;

  // Set up F.
  F.f = &binom_transform_gsl;
  F.n = 2;
  F.params = (void *)params;

  // Set up initial vector.
  gsl_vector* x = gsl_vector_alloc(2);
  gsl_vector_set_all(x, 1.0);

  gsl_multiroot_fsolver_set(s, &F, x);
  // printf("x: %g, %g \t f: %g, %g\n", s->x->data[0], s->x->data[1], s->f->data[0], s->f->data[0]);

  int i = 0;
  int msg = GSL_CONTINUE;
  for(i = 0; i < max_iter && msg != GSL_SUCCESS; i++) {
    gsl_multiroot_fsolver_iterate(s);
    // printf("x: %g, %g \t f: %g, %g\n", s->x->data[0], s->x->data[1], s->f->data[0], s->f->data[0]);
    // check |dx| < epsabs + epsrel * |x|
    msg = gsl_multiroot_test_delta(s->dx, s->x, epsabs, epsrel);
  }

  memmove(rs, s->x->data, 2 * sizeof(double));
}

void binom_post(const double* prior, double* post, double y, double n, double epsrel, int max_iter)
{
  double rs[2];
  double fq[2] = {0, 0};
  binom_solver(prior, rs, 0.0, epsrel, max_iter);
  rs[0] = rs[0] + y;
  rs[1] = rs[1] + n - y;
  binom_transform(rs, fq, post);
}

int main(int argc, char** argv)
{
  utest_binom_transform(1, 1);
  utest_binom_transform_gsl(1, 1);

  double prior[2] = {0.0, 3.3};
  double post[2];
  double rs[2];

  binom_solver(prior, rs, 0.0, 1e-8, 100);
  printf("rs = %g, %g\n", rs[0], rs[1]);

  binom_post(prior, post, 1, 1, 1e-8, 100);
  printf("post = %g, %g\n", post[0], post[1]);

  return 0;
}
