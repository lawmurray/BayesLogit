// -*- mode: c++; -*-

#include "Matrix/Matrix.h"
#include "RNG/RNG.hpp"
#include "PolyaGamma.hpp"
#include "Logit.hpp"
#include "MultLogit.hpp"

int test_mult()
{
  int P = 3;
  int J = 3;
  int N = 1000;

  RNG r;

  Matrix w, beta;

  Matrix y(N, J-1);
  Matrix X(N,   P);
  Matrix n(N);
  y.read("IO/y.dat", false);
  X.read("IO/X.dat", false);
  n.read("IO/n.dat", false);

  Matrix ty(J-1, N);
  Matrix tX(  P, N);
  ty.copy_transpose(y);
  tX.copy_transpose(X);

  Matrix m0(P, 1, J-1);
  Matrix P0(P, P, J-1);

  MultLogit multlogit(ty, tX, n);
  multlogit.gibbs(w, beta, m0, P0, 10000, 100, r); 

  return 0;
}

int main(int argc, char **argv)
{
  uint P = 1;
  uint N = 20;

  RNG r;

  Matrix tX(P, N);
  Matrix beta_known(P); beta_known.fill(0.5);
  Matrix y(N);
  Matrix n(N); n.fill(1.0);

  r.norm(tX, 0, 1.0);
  Matrix psi_known(tX, beta_known, 'T', 'N');

  Matrix p(N);
  for(uint i=0; i < N; i++){
    p(i) = 1 / ( 1 + exp(psi_known(i)) );
    double u = r.unif();
    if (u < p(i)) y(i) = 1.0;
  }

  Logit logit(y, tX, n, 0.5, Matrix(0.0), 1.0);

  Matrix w, beta;
  logit.gibbs(w, beta, 10000, 100, r);

  // test_mult();

  printf("End.\n");

  // // If I want to do clever Gibbs sampling.
  // int curr = 0;
  // int prev = 0;
  // int samp = 10;
  // int burn = 10;
  // for(int m = 0; m < (burn+samp); ++m){
  //   prev = m < burn+1 ? 0 : curr++;
  //   cout << prev << " " << curr << "\n";
  // }
  // cout << "\n";

  // PoylaGamma pg;
  // cout << pg.draw((int)1, 2.0, r) << "\n";

  return 0;
}
