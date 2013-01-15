// -*- mode: c++; -*-

#ifndef __FFBS__
#define __FFBS__

#include <vector>
#include <Eigen/Dense>
#include <stdio.h>
#include "RNG.hpp"

#ifdef USE_R
#include <R_ext/Utils.h>
#endif

using std::vector;
using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::MatrixBase;
using Eigen::Matrix;
using Eigen::Dynamic;

////////////////////////////////////////////////////////////////////////////////

extern "C" {

  void ffbs(double *alpha_, double *beta_,
	    double *z_, double *X_, double *V_,
	    double *mu_, double *phi_, double *W_, 
	    double *m0_, double *C0_, int *N_b_, int *N_, int *T_);
  
}

////////////////////////////////////////////////////////////////////////////////
			   // TEMPLATED FUNCTIONS //
////////////////////////////////////////////////////////////////////////////////

template <typename dM, typename dV>
void ffbs(MatrixBase<dV> &alpha, MatrixBase<dM> &beta,
	  MatrixBase<dV> &z, MatrixBase<dM> &X, MatrixBase<dV> &V,
	  MatrixBase<dV> &mu, MatrixBase<dV> &phi, MatrixBase<dM> &W, 
	  MatrixBase<dV> &m0, MatrixBase<dM> &C0, RNG& r)
{
  // When tracking (alpha, beta_t).  It may be the case that there is no alpha.
  // z_t = x_t (alpha_t, beta_t) + ep_t, ep_t \sim N(0, V_t).
  // beta_t = mu + phi * (beta_t - mu) + omega_t, omega_t \sim N(0,W).
  // alpha_t = alpha_{t-1}
  
  // alpha: vector of static coefficients (N_a)
  // beta: matrix of dynamic coefficients (N_b x T+1)

  // z : vector of observations (T)
  // X : design matrix (T x N)
  // mu : vector (N_b)
  // phi : vector (N_b)
  // W : covariance MATRIX of innovations of beta (N_b x N_b)
  // V: vector, time varying variances (T)
  // m0 : prior mean on (beta_0, alpha_0) (N)
  // C0 : prior var on (beta_0, alpha_0) (N x N).

  int T   = z.size();
  int N   = X.cols();
  int N_b = mu.size();
  int N_a = N - N_b;

  // printf("T=%i, N=%i, N_b=%i, N_a=%i\n", T, N, N_b, N_a);

  bool with_alpha = N_a > 0;

  // Setup objects to track.
  vector<VectorXd> m(T+1, VectorXd(N)   ); m[0] = m0;
  vector<MatrixXd> C(T+1, MatrixXd(N, N)); C[0] = C0;
  vector<VectorXd> a(T+1, VectorXd(N)   );
  vector<MatrixXd> R(T+1, MatrixXd(N, N));

  // beta.resize(N_b, T+1);

  // Setup "big" evolution coefficients.
  VectorXd big_phi(N);
  VectorXd big_mu = VectorXd::Zero(N);
  MatrixXd big_W  = MatrixXd::Zero(N, N);

  if (with_alpha) big_phi.segment(0, N_a).setOnes();

  big_phi.segment(N_a, N_b) = phi;
  big_mu.segment(N_a, N_b)  = mu;
  big_W.block(N_a, N_a, N_b, N_b) = W;

  VectorXd _1m_big_phi = VectorXd::Ones(N) - big_phi; // 1 - Phi
  MatrixXd Phi = big_phi.asDiagonal();

  // Filter Forward
  for (int i=1; i<(T+1); i++) {
    int i_l = i-1;

    a[i] = big_phi.array() * m[i-1].array() + (_1m_big_phi).array() * big_mu.array();
    R[i] = Phi * C[i-1] * Phi + big_W;

    Matrix<double, 1, Dynamic> tF  = X.row(i_l);
    double                     f   = tF * a[i];

    Matrix<double, 1, 1> V_t;  V_t(0) = V(i_l);
    Matrix<double, 1, 1>       Q      = tF * R[i] * tF.transpose() + V_t;
    Matrix<double, 1, 1> QI;   QI(0)  = 1/Q(0);

    Matrix<double, Dynamic, 1> RF = R[i] * tF.transpose();
    double                     e  = z(i_l) - f;
    Matrix<double, Dynamic, 1> A  = RF * QI;

    m[i] = a[i] + A * e;
    C[i] = R[i] - RF * QI * RF.transpose();

    #ifdef USE_R
    R_CheckUserInterrupt();
    #endif
  }

  // Backwards sample.
  VectorXd draw(N); r.norm(draw, 1.0);

  VectorXd theta = m[T] + C[T].llt().matrixL() * draw;

  if (with_alpha) alpha = theta.segment(0, N_a);
  beta.col(T) = theta.segment(N_a, N_b);

  draw.resize(N_b);

  for (int i=T; i>0; i--) {
    MatrixXd Rsub = R[ i ].block(N_a, N_a, N_b, N_b); // Could use map.
    MatrixXd Csub = C[i-1].block(N_a, N_a, N_b, N_b); // Could use map.
    MatrixXd tA   = Rsub.llt().solve(Csub * phi.asDiagonal());

    VectorXd e    = beta.col(i) - a[i].segment(N_a, N_b);
    VectorXd m_bs = m[i-1].segment(N_a, N_b) + tA.transpose() * e;
    MatrixXd V_bs = Csub - tA.transpose() * Rsub * tA;

    r.norm(draw, 1.0);
    beta.col(i-1) = m_bs + V_bs.llt().matrixL() * draw;

    #ifdef USE_R
    R_CheckUserInterrupt();
    #endif
  }

} // ffbs

#endif

////////////////////////////////////////////////////////////////////////////////
				 // APPENDIX //
////////////////////////////////////////////////////////////////////////////////

// MatrixBase is base class for all dense matrices, vectors, and expressions.
// See the following about writing functions in Eigen.
// <http://eigen.tuxfamily.org/dox/classEigen_1_1MatrixBase.html#details>.
