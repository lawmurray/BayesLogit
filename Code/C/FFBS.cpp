
#include <Eigen/Dense>
#include <stdexcept>
#include "RNG.hpp"
#include "FFBS.h"

using Eigen::Map;

////////////////////////////////////////////////////////////////////////////////
				  // EXTERN //
////////////////////////////////////////////////////////////////////////////////

void ffbs(double *alpha_, double *beta_,
	  double *z_, double *X_,  double *V_,
	  double *mu_, double *phi_, double *W_,
	  double *m0_, double *C0_, int *N_b_, int *N_, int *T_)

{
  RNG r;
  
  int T = *T_;
  int N = *N_;
  int N_b = *N_b_;
  int N_a = N - N_b;

  if (T > 10000 || N > 1000) {
    fprintf(stderr, "ffbs: T=%i or N=%i is very large.  Aborting.\n", T, N);
    // throw std::runtime_error("T or N is very large.  Aborting.\n");
    return;
  }

  Map<VectorXd> alpha(alpha_, N_a > 1 ? N_a : 1);
  Map<MatrixXd> beta (beta_ , N_b, T+1);
  Map<VectorXd> z    (z_    , T       );
  Map<MatrixXd> X    (X_    , T  , N  );
  Map<VectorXd> mu   (mu_   , N_b     );
  Map<VectorXd> phi  (phi_  , N_b     );
  Map<MatrixXd> W    (W_    , N_b, N_b);
  Map<VectorXd> V    (V_    , T       );
  Map<VectorXd> m0   (m0_   , N       );
  Map<MatrixXd> C0   (C0_   , N  , N  );

  ffbs(alpha, beta, z, X, V, mu, phi, W, m0, C0, r);
}
