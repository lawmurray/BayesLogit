#include "DynExpFamMH.hpp"
#include <iostream>
#include "RNG.hpp"

using std::cout;
using std::cerr;

void test_make_A(){

  VectorXd phi(2);
  phi << 0.8, 0.9;

  MatrixXd A;

  make_A(A, phi, 5);

  cout << "A:\n" << A << "\n";

}

void test_make_L()
{

  VectorXd phi(2);
  phi << 0.8, 0.9;

  MatrixXd L;

  make_L(L, phi, 5);

  cout << "L:\n" << L << "\n";

}

void test_A_and_L()
{

  int num = 5;
  int bsize = 2;
  int N   = num * bsize;

  VectorXd phi(bsize);
  phi << 0.8, 0.9;

  MatrixXd A(N, N); A.fill(0.0);
  MatrixXd L(N, N); L.fill(0.0);

  make_A(A, phi, 5);

  cout << "A:\n" << A << "\n";

  make_L(L, phi, 5);

  cout << "L:\n" << L << "\n";

  cout << "L * A\n" << L * A << "\n";

  cout << "A * L\n" << A * L << "\n";

}

void test_make_XiL()
{
  int num = 5;
  int bsize = 2;
  int nb   = num * bsize;

  VectorXd phi(bsize);
  phi << 0.8, 0.9;

  MatrixXd L(nb, nb); L.fill(0.0);

  make_L(L, phi, 5);

  MatrixXd tX(bsize, num); tX.fill(1.0);

  MatrixXd Xi(num, nb);
  
  make_Xi(Xi, tX, L);

  cout << "Xi:\n" << Xi << "\n";

}

void test_laplace_omega()
{
  int num = 5;
  int bsize = 2;
  int nb   = num * bsize;

  VectorXd phi(bsize);
  phi << 0.8, 0.9;

  MatrixXd L(nb, nb); L.fill(0.0);
  make_L(L, phi, 5);

  MatrixXd tX(bsize, num); tX.fill(1.0);
  VectorXd y(num); y.fill(1.0);
  
  MatrixXd Xi(num, nb);
  make_Xi(Xi, tX, L);

  int block_start = 1;
  int num_blocks  = 2;

  llh_struct llh(num);
  Gaussian nbp(bsize*num_blocks);

  MatrixXd prior_prec = VectorXd::Ones(nb).asDiagonal();
  VectorXd omega(bsize*num); omega.fill(1.0);

  llh.l0.fill(1.0);
  llh.l1.fill(1.0);
  llh.l2.fill(-1.0);
  llh.pl2.fill(1.0);

  cout << "Call laplace omega.\n";

  VectorXd omega_vec = omega.array();

  laplace_omega(nbp, omega_vec, llh,
				y, Xi, prior_prec,
				block_start, num_blocks);

  RNG r;

  VectorXd draw(bsize*num_blocks);
  rNorm(draw, nbp, r);

  cout << "m:\n" << nbp.m << "\n";
  cout << "d:\n" << draw << "\n";

}

void test_draw_omega_block()
{
  int num = 5;
  int bsize = 2;
  int nb   = num * bsize;

  int block_start = 0;
  int num_blocks  = num;

  VectorXd ntrials(num); ntrials.fill(1);
  // int ntrials = 1;

  VectorXd phi(bsize);
  phi << 0.8, 0.9;

  MatrixXd Phi = phi.asDiagonal();

  MatrixXd L(nb, nb); L.fill(0.0);
  make_L(L, phi, 5);

  MatrixXd tX(bsize, num); tX.fill(1.0);
  VectorXd y(num); y.fill(1.0);
  
  MatrixXd Xi(num, nb);
  make_Xi(Xi, tX, L);

  llh_struct llh(num);
  Gaussian nbp(bsize*num_blocks);

  MatrixXd prior_prec = VectorXd::Ones(nb).asDiagonal();

  MatrixXd omega(bsize, num); omega.fill(0.0);
  MatrixXd beta(bsize, num);

  omega_to_dyn_beta(beta, Phi, omega, 0);

  dyn_beta_to_psi(llh.psi_dyn, tX, beta, 0);
  llh.psi_stc.fill(0.0);
  log_logit_likelihood(&y(0), &ntrials(0), llh, 0);

  cout << "Call draw omega block.\n";

  RNG r;

  VectorXd offset(nb); offset.fill(0.0);
  
  draw_omega_block(omega, beta, llh,
				   y, tX, ntrials, offset,
				   Xi, L,
				   prior_prec, Phi,
				   block_start, num_blocks, 
				   r, &log_logit_likelihood);

  cout << "omega:\n" << omega << "\n";

}

void test_draw_omega(int reps=1)
{
  int num = 5;
  int bsize = 2;
  int nb   = num * bsize;

  VectorXd ntrials(num); ntrials.fill(1);

  VectorXd phi(bsize);
  phi << 0.8, 0.9;

  MatrixXd Phi = phi.asDiagonal();

  MatrixXd L(nb, nb); L.fill(0.0);
  make_L(L, phi, 5);

  MatrixXd tX(bsize, num); tX.fill(1.0);
  VectorXd y(num); y.fill(1.0);
  
  MatrixXd Xi(num, nb);
  make_Xi(Xi, tX, L);

  llh_struct llh(num);
  Gaussian nbp(bsize*num);

  MatrixXd prior_prec = VectorXd::Ones(nb).asDiagonal();

  MatrixXd omega(bsize, num); omega.fill(0.0);
  MatrixXd beta(bsize, num);

  omega_to_dyn_beta(beta, Phi, omega, 0);

  dyn_beta_to_psi(llh.psi_dyn, tX, beta, 0);
  llh.psi_stc.fill(0.0);
  log_logit_likelihood(&y(0), &ntrials(0), llh, 0);

  cout << "Call draw omega.\n";

  // cout << "Starting psi:\n" << llh.psi << "\n";

  RNG r;

  VectorXd offset(nb); offset.fill(0.0);

  MatrixXi starts(2,1);
  starts << 0, 3;

  for (int i=0; i<reps; i++) {
    draw_omega(omega, beta, llh,
			   y, tX, ntrials, offset,
			   Xi, L,
			   prior_prec, Phi,
			   starts, 
			   r, &log_logit_likelihood, true);
    cout << "omega:\n" << omega << "\n";
  }

  cout << "omega:\n" << omega << "\n";

  cout << "beta:\n" << beta << "\n";

  cout << "psi:\n" << llh.psi.transpose() << "\n";

}

void test_draw_omega_wrapper(int reps=1)
{
  int num = 500;
  int bsize = 2;
  int nb   = num * bsize;

  VectorXd ntrials(num); ntrials.fill(1.0);
  // int ntrials = 1;

  VectorXd phi(bsize);
  phi << 0.8, 0.9;

  MatrixXd Phi = phi.asDiagonal();

  MatrixXd L(nb, nb); L.fill(0.0);
  make_L(L, phi, num);

  MatrixXd tX(bsize, num); tX.fill(1.0);
  VectorXd y(num); y.fill(1.0);
  
  llh_struct llh(num);

  MatrixXd prior_prec = VectorXd::Ones(nb).asDiagonal();

  MatrixXd omega(bsize, num); omega.fill(0.0);
  MatrixXd beta(bsize, num);

  omega_to_dyn_beta(beta, Phi, omega, 0);

  dyn_beta_to_psi(llh.psi_dyn, tX, beta, 0);
  llh.psi_stc.fill(0.0);
  log_logit_likelihood(&y(0), &ntrials(0), llh, 0);

  cout << "Call draw omega wrapper.\n";

  // cout << "Starting psi:\n" << llh.psi << "\n";

  RNG r;

  VectorXd offset(nb); offset.fill(0.0);

  int nstarts = 5;
  MatrixXi starts(nstarts,1);
  starts << 0, 99, 199, 299, 399;

  int  naccept  = 0;
  bool just_max = true;
  int type      = 0;

  for (int i=0; i<reps; i++) {
    draw_omega(&omega(0), &beta(0), 
	       &llh.psi_dyn(0), &llh.psi_stc(0),
	       &llh.psi(0), &llh.l0(0), &llh.l1(0), &llh.l2(0), &llh.pl2(0),
	       &y(0), &tX(0), &ntrials(0), &offset(0),
	       &prior_prec(0), &phi(0),
	       &starts(0), &num, &bsize, &nstarts, &naccept, &just_max, &type);
    // cout << "omega:\n" << omega << "\n";
  }

  cout << "omega:\n" << omega << "\n";

  cout << "beta:\n" << beta << "\n";

  cout << "psi:\n" << llh.psi.transpose() << "\n";

}

void test_draw_omega_real_data(int reps=1)
{
  int num = 10;
  int bsize = 1;
  int nb   = num * bsize;

  VectorXd ntrials(num); ntrials.fill(20.0);
  // int ntrials = 20;

  VectorXd phi(bsize);
  phi << 0.9;

  MatrixXd Phi = phi.asDiagonal();

  MatrixXd L(nb, nb); L.fill(0.0);
  make_L(L, phi, num);

  MatrixXd tX(bsize, num); tX.fill(1.0);
  VectorXd y(num); 
  y << 4, 13, 14, 15, 13, 15, 19, 17, 14, 16;
  
  MatrixXd Xi(num, nb);
  make_Xi(Xi, tX, L);

  llh_struct llh(num);
  Gaussian nbp(bsize*num);

  MatrixXd prior_prec = VectorXd::Ones(nb).asDiagonal();
  prior_prec(0,0) = 0.5524862;

  MatrixXd omega(bsize, num); omega.fill(0.0);
  MatrixXd beta(bsize, num);

  omega_to_dyn_beta(beta, Phi, omega, 0);

  dyn_beta_to_psi(llh.psi_dyn, tX, beta, 0);
  llh.psi_stc.fill(0.0);
  log_logit_likelihood(&y(0), &ntrials(0), llh, 0);

  cout << "Call draw omega with real data.\n";

  // cout << "Starting psi:\n" << llh.psi << "\n";

  RNG r;
  VectorXd offset(nb); offset.fill(0.0);

  int nstarts = 2;
  MatrixXi starts(nstarts,1);
  starts << 0, 3;

  for (int i=0; i<reps; i++) {
    cout << "omega: " << omega << "\n";
    // cout << "beta: " << beta << "\n";
    draw_omega(omega, beta, llh,
			   y, tX, ntrials, offset,
			   Xi, L,
			   prior_prec, Phi,
			   starts, 
			   r, &log_logit_likelihood, true);
  }

  cout << "omega: " << omega << "\n";
  // cout << "beta: " << beta << "\n";
  // cout << "psi:\n" << llh.psi.transpose() << "\n";

}

void test_draw_omega_wrapper_real_data(int reps=1)
{
  int num = 10;
  int bsize = 1;
  int nb   = num * bsize;

  VectorXd ntrials(num); ntrials.fill(20.0);
  // int ntrials = 20;

  VectorXd phi(bsize);
  phi << 0.9;

  MatrixXd Phi = phi.asDiagonal();

  MatrixXd L(nb, nb); L.fill(0.0);
  make_L(L, phi, num);

  MatrixXd tX(bsize, num); tX.fill(1.0);
  VectorXd y(num); 
  y << 4, 13, 14, 15, 13, 15, 19, 17, 14, 16;
  
  MatrixXd Xi(num, nb);
  make_Xi(Xi, tX, L);

  llh_struct llh(num);
  Gaussian nbp(bsize*num);

  MatrixXd prior_prec = VectorXd::Ones(nb).asDiagonal();
  prior_prec(0,0) = 0.5524862;

  MatrixXd omega(bsize, num); omega.fill(0.0);
  MatrixXd beta(bsize, num);

  omega_to_dyn_beta(beta, Phi, omega, 0);

  dyn_beta_to_psi(llh.psi_dyn, tX, beta, 0);
  llh.psi_stc.fill(0.0);
  log_logit_likelihood(&y(0), &ntrials(0), llh, 0);

  cout << "Call draw omega wrapper with real data.\n";

  // cout << "Starting psi:\n" << llh.psi << "\n";

  RNG r;
  VectorXd offset(nb); offset.fill(0.0);

  int nstarts = 2;
  MatrixXi starts(nstarts,1);
  starts << 0, 5;

  int  naccept  = 0;
  bool just_max = true;

  for (int i=0; i<reps; i++) {
    cout << "omega:\n" << omega << "\n";
    draw_omega(&omega(0), &beta(0), 
	       &llh.psi_dyn(0), &llh.psi_stc(0),
	       &llh.psi(0), &llh.l0(0), &llh.l1(0), &llh.l2(0), &llh.pl2(0),
	       &y(0), &tX(0), &ntrials(0), &offset(0),
	       &prior_prec(0), &phi(0),
	       &starts(0), &num, &bsize, &nstarts, &naccept, &just_max, 0);
  }

  cout << "omega:\n" << omega << "\n";
  // cout << "beta:\n" << beta << "\n";
  // cout << "psi:\n" << llh.psi.transpose() << "\n";

}

int main(int argc, char** argv)
{
  // test_make_A();
  // test_make_L();
  // test_A_and_L();
  // test_make_XiL();
  // test_laplace_omega();
  // test_draw_omega_block();
  // test_draw_omega(5);
  test_draw_omega_wrapper(100);
  // test_draw_omega_real_data(2);
  // test_draw_omega_wrapper_real_data(2);
}
