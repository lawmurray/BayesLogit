
#include <Eigen/dense>
#include "RNG.hpp"

using std::cout;
using std::endl;

using namespace Eigen;

template<typename Mat>
void fun1(Mat &a){
  a.array() += 1;
}

template<typename M1>
void fun2(MatrixBase<M1> &a){
  a.array() += 1;
}

template<typename V1>
void fun3(MatrixBase<V1> &a){
  a.array() += 1;
}

int main(int argc, char** argv)
{
  RNG r;

  // MatrixXd A(4,1);
  // Matrix<double, 1, 1> mu;  mu(0)  = 0.0;
  // Matrix<double, 1, 1> sig; sig(0) = 1.0;

  MatrixXd A(4,1);
  VectorXd v(4);

  r.norm(A, 1.0);

  Map<MatrixXd> M(&A(0), 4, 1);
  fun2(M);

  fun3(v);

  cout << "A:\n" << A << endl;

  return 0;
}
