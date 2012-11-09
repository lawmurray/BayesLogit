// -*- mode: c++; fill-column: 70; -*-

//////////////////////////////////////////////////////////////////////

// Copyright 2012 Jesse Windle - jwindle@ices.utexas.edu

// This program is free software: you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License
// as published by the Free Software Foundation, either version 3 of
// the License, or (at your option) any later version.

// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this program.  If not, see
// <http://www.gnu.org/licenses/>.

//////////////////////////////////////////////////////////////////////

// This class "implements" the MatrixFrame class.  You can think the
// Matrix class as a data container and the MatrixFrame class as a
// view of that container.  See MatrixFrame.h for documentation of
// MatrixFrame class.

// Everything here is COLUMN MAJOR for compatibility with Fortran.

// When compiling include -lblas -llapack.

//////////////////////////////////////////////////////////////////////

#ifndef __MATRIX__
#define __MATRIX__

#include <vector>
#include <iostream>
#include <string>
#include <sstream>
#include <cmath>
#include <algorithm> // For min,max.
#include <stdio.h>

#ifndef DISABLE_FIO
#include <fstream>
using std::ofstream;
using std::ifstream;
#endif

#include "MatrixFrame.h"

using std::vector;
using std::string;
using std::ostream;
using std::istream;
using std::stringstream;
using std::min;
using std::max;

#ifndef vdouble
typedef vector<double> vdouble;
#endif

//////////////////////////////////////////////////////////////////////
			     // Matrix //
//////////////////////////////////////////////////////////////////////

class Matrix : public MatrixFrame
{
 protected:
  vdouble v;   // vector to store data.

 public:

  // Constructors.
  Matrix() : MatrixFrame(), v(1)
    { p = &v[0]; nr = 1; nc = 1; nm = 1; }
  Matrix(uint r, uint c=1, uint n=1) : MatrixFrame(), v(r*c*n)
    { p = &v[0]; nr = r; nc = c; nm = n; }
  Matrix( int r,  int c=1,  int n=1) : MatrixFrame(), v(r*c*n)
    { p = &v[0]; nr = (uint)r; nc = (uint)c; nm = (int)n; }
  Matrix(uint r, uint c, uint n, double f): MatrixFrame(), v(r*c*n, f)
    { p = &v[0]; nr = r; nc = c; nm = n; }
  Matrix(double d, const Matrix& M) : MatrixFrame(), v(M.nr*M.nc*M.nm, d)
    { p = &v[0]; nr = M.nr; nc = M.nc; nm = M.nm; }
  Matrix(const Matrix& M) : MatrixFrame(), v(M.v)
    { p = &v[0]; nr = M.nr; nc = M.nc; nm = M.nm; }
  Matrix(const MatrixFrame& M, uint n=1) : MatrixFrame(), v(M.vol()*n)
    { p = &v[0]; nr = M.rows(); nc = M.cols(); nm = M.mats()*n;
      for(uint i = 0; i < nr*nc*nm; i++) v[i] = M(i % (nr*nc*M.mats()) ); }
  Matrix(double d) : MatrixFrame(), v(1, d)
  { p = &v[0]; nr = 1; nc = 1; nm = 1; }

  Matrix(const double *ptr, uint r, uint c, uint n=1) : MatrixFrame(), v(r*c*n)
  { p = &v[0]; nr = r; nc = c; nm = n;
    for(uint i = 0; i < nr*nc*nm; i++) v[i] = ptr[i]; }
  Matrix(const double *ptr,  int r,  int c,  int n=1) : MatrixFrame(), v(r*c*n)
  { p = &v[0]; nr = (uint)r; nc = (uint)c; nm = (uint)n;
    for(uint i = 0; i < nr*nc*nm; i++) v[i] = ptr[i]; }

  // For predefined types of matrices.
  Matrix(const string& s, uint r, uint n=1);

  // To instatiate from matrix multiplication.
  Matrix(const MatrixFrame& a, const MatrixFrame& b, char ta='N', char tb='N', double alpha = 1.0);

  // ~Matrix() {};

  // Test equality and assign equality.
  Matrix& operator= (const Matrix &M);      // Makes copy.
  Matrix& operator= (const MatrixFrame &M); // Makes copy.
  // bool    operator==(const Matrix &M) const;

  // Iterators...
  vector<double>::iterator begin()
  { return v.begin(); }
  vector<double>::iterator end()
  { return v.end(); }
  vector<double>::const_iterator begin() const
  { return v.begin(); }
  vector<double>::const_iterator end() const
  { return v.end(); }

  // Utility functions.
  void resize(uint r, uint c=1, uint n=1)
  { v.resize(r*c*n); p = &v[0]; nr = r; nc = c; nm = n; }
  void clone(const MatrixFrame& M);
  void clone(const MatrixFrame& M, const MatrixFrame& rs, const MatrixFrame& cs);
  void clone(const MatrixFrame& M, const MatrixFrame& rs, uint c);
  void clone(const MatrixFrame& M, uint r, const MatrixFrame& cs);
  // void copy(const Matrix& M);
  //void cbind(const MatrixFrame& M);
  //void rbind(const MatrixFrame& M);

  // Read //
  uint read(      istream&  is, bool header=0, bool binary=0);
  uint readstring(const string& s, bool header=0);
  #ifndef DISABLE_FIO
  uint read(const string& file, bool header=0, bool binary=0);
  #endif

  // Writing is taken care of in MatrixFrame.h.
  // Matrix operations are taken care of in MatrixFrame.h

  
}; // Matrix

//////////////////////////////////////////////////////////////////////
			  // Constructors //
//////////////////////////////////////////////////////////////////////

Matrix::Matrix(const string& s, uint r, uint n) : MatrixFrame(), v(1)
{
  switch (s[0])
    {
    case 'I': // The identity matrix.
      resize(r, r, n);
      for(uint k = 0; k < nm; k++)
	for(uint i = 0; i < nr; i++)
	  operator()(i,i,k) = 1;
      break;
    case '1': // The "unity" column vectors.
      resize(r, 1, n);
      for(uint k = 0; k < nm; k++)
	for(uint i = 0; i < nr; i++)
	  operator()(i,0,k) = 1;
      break;
    case 'N': // The Natural numbers.
      resize(r, 1, 1);
      for(uint i = 0; i < nr; i++)
	operator()(i,0,0) = (double)(i+1);
      break;
    case 'W': // The Whole numbers.
      resize(r, 1, 1);
      for(uint i = 0; i < nr; i++)
	operator()(i,0,0) = (double)i;
      break;
    // case 'Z': // A sequence of integers.
    //   int diff = (int)n - (int)r;
    //   int dist = abs(diff);
    //   resize(dist+1, 1, 1);
    //   int sgn = diff / dist;
    //   int idx = 0;
    //   while(idx <= dist){
    // 	operator()(idx++) = r;
    // 	r = r + sgn;
    //   }
    //   break;
    default: // Set to scalar zero.
      resize(1,1,1);
    }
}

Matrix::Matrix(const MatrixFrame& a, const MatrixFrame& b, char ta, char tb, double alpha) : MatrixFrame(), v(1)
{
  uint opa_rows = ta=='T' ? a.cols() : a.rows();
  uint opb_cols = tb=='T' ? b.rows() : b.cols();
  resize(opa_rows, opb_cols, 1);

  // We trick g++ here.
  MatrixFrame c(&v[0], opa_rows, opb_cols, (uint)1);
  gemm(c, a, b, ta, tb, alpha);
}

//////////////////////////////////////////////////////////////////////
		       // Utility Functions //
//////////////////////////////////////////////////////////////////////

void Matrix::clone(const MatrixFrame& M)
{
  resize(M.rows(), M.cols(), M.mats());
  copy(M);
} // copy

void Matrix::clone(const MatrixFrame& M, const MatrixFrame& rs, const MatrixFrame& cs)
{
  resize(rs.area(), cs.area(), 1);
  copy(M, rs, cs);
}

void Matrix::clone(const MatrixFrame& M, const MatrixFrame& rs, uint c)
{
  resize(rs.area(), 1);
  copy(M, rs, c);
}

void Matrix::clone(const MatrixFrame& M, uint r, const MatrixFrame& cs)
{
  resize(1, cs.area());
  copy(M, r, cs);
}

// void Matrix::cbind(const MatrixFrame& M)
// {
//   sizecheck(mats()==M.mats() && rows()==M.rows());
//   Matrix temp(*this);
//   resize(rows(), cols() + M.cols(), mats());
//   for(uint m = 0; m < mats(); ++m){
//     copy(temp[m], 0, 0);
//     col(nc, M.cols()).copy(M[m], 0, 0);
//   }
// }

// Is it a bad idea to overload a function found in MatrixFrame?
// According to Effective C++ it is, but this makes things mroe
// intuitive.

// void Matrix::copy(const Matrix& M)
// {
//   resize(M.rows(), M.cols(), M.mats());
//   for(uint i = 0; i < vol(); i++) v[i] = M.vec(i);
// } // copy

//////////////////////////////////////////////////////////////////////
		  // Assgiment and Test Equality //
//////////////////////////////////////////////////////////////////////

Matrix& Matrix::operator= (const Matrix &M)
{
  clone(M);
  return *this;
} // operator=

// May not need both of these operators.

Matrix& Matrix::operator= (const MatrixFrame &M)
{
  clone(M);
  return *this;
} // operator=

// bool Matrix::operator==(const Matrix &M) const
// {
//   if(p==&M(0) && nr==M.rows() && nc==M.cols() && nm==M.mats()) return true;
//   if(vol() != M.vol()) return false;
//   for(uint i = 0; i < vol(); i++) if(vec(i)!=M.vec(i)) return false;
//   return true;
// }

//////////////////////////////////////////////////////////////////////
			      // READ //
//////////////////////////////////////////////////////////////////////

// Read in a matrix from a stream.  The header contains the
// dimensionality information, i.e. rows, cols, mats.  If header is
// set to true these values are read from the stream and used to set
// the dimensions of this matrix.

uint Matrix::read( std::istream& is, bool header, bool binary)
{
  // Tell us if something is wrong.
  if (!is || is.eof())  return 0;

  uint n = 0; // The number of items read.

  // Read binary.
  if(binary){
    if(header){
      uint r,c,m;
      is.read((char*) &r, sizeof(nr));
      is.read((char*) &c, sizeof(nc));
      is.read((char*) &m, sizeof(nm));
      resize(r, c, m);
    }
    n = scan(is, false, true);
  }
  // Write human.
  if(!binary){
    if(header){
      uint r,c,m;
      is >> r;
      is >> c;
      is >> m;
      resize(r, c, m);
    }
    n = scan(is, false, false);
  }
  return n;
} // read

#ifndef DISABLE_FIO
uint Matrix::read(const string& file, bool header, bool binary)
{
  std::ifstream ifs(file.c_str());
  if(!ifs){
    fprintf(stderr, "Cannot read file %s.\n", file.c_str());
    return 0;
  }
  return read(ifs, header, binary);
} // read
#endif

uint Matrix::readstring(const string& s, bool header)
{
  stringstream ss(s);
  return read(ss, header, false);
} // readstring

//////////////////////////////////////////////////////////////////////
		      // END OF CLASS METHODS //
//////////////////////////////////////////////////////////////////////

#ifndef Mat
typedef Matrix Mat;
#endif

#ifndef Matrices
typedef Matrix Matrices;
#endif

//////////////////////////////////////////////////////////////////////
			    // CASTING //
//////////////////////////////////////////////////////////////////////

// This seems like it might be unnecessary.
Matrix cast(double d)
{
  return Matrix(1, 1, 1, d);
}

//////////////////////////////////////////////////////////////////////
	  // Hadamard operations by OVERLOADED OPERATORS //
//////////////////////////////////////////////////////////////////////

// These form of multiplication, addition, division, and subtraction
// will not be as fast as hprodeq, hsumeq, etc since it involves
// returning a copy of a Matrix.  However, it will make it easier to
// read code.  There are matrix packages out there, e.g. eigen, that
// will cleverly reduce the amount of overhead for concatenated
// operations.  But I found it to be lacking in documentation given
// its complexity.

#define MHOP(NAME, OP, OPEQ, TYPE)				\
  Matrix operator OP(const TYPE& a, const TYPE& b)		\
  {								\
    MatrixFrame small = a.area() < b.area() ? a : b;		\
    MatrixFrame big   = a.area() < b.area() ? b : a;		\
    sizecheck(hconform(big,small));				\
    Matrix c(big);						\
    NAME(c, small, 1.0);					\
    return c;							\
  }								\
  Matrix operator OP(const TYPE& a, double b)			\
  {								\
    Matrix c(a);						\
    NAME(c, b);							\
    return c;							\
  }								\
  Matrix operator OP(double b, const TYPE& a)			\
  {								\
    Matrix c(a);						\
    NAME(c, b);							\
    return c;							\
  }								\

//MHOP(hprodeq, *, *=, Matrix) MHOP(hsumeq, +, +=, Matrix)
//MHOP(hdiveq,  /, /=, Matrix) MHOP(hsubeq, -, -=, Matrix)

MHOP(hprodeq, *, *=, MF) MHOP(hsumeq, +, +=, MF)
MHOP(hdiveq,  /, /=, MF) MHOP(hsubeq, -, -=, MF)

#undef MHOP

//////////////////////////////////////////////////////////////////////
			// Basic Functions //
//////////////////////////////////////////////////////////////////////

#ifndef sq
double sq(double x){return x * x;}
#endif

#define UNARY(FUNC, FUNCEQ)					\
  MatrixFrame FUNC(MF a, MF b)					\
  {								\
    sizecheck(a.vol()==b.vol());					\
    for(uint l = 0; l < a.vol(); l++) a(l) = FUNC(b(l));	\
    return a;							\
  }								\
  Matrix FUNC(MF a)						\
  {								\
    Matrix c(a);						\
    FUNC(c, a);							\
    return c;							\
  }								\
  MatrixFrame FUNCEQ(MF a)					\
  {								\
    for(uint l = 0; l < a.vol(); l++) a(l) = FUNC(a(l));	\
    return a;							\
  }								\

UNARY(log, logq)   UNARY(exp, expq)   UNARY(sqrt, sqrtq)
UNARY(sin, sinq)   UNARY(cos, cosq)   UNARY(tan , tanq)
UNARY(asin, asinq) UNARY(acos, acosq) UNARY(atan, atanq)
UNARY(sinh, sinhq) UNARY(cosh, coshq) UNARY(tanh, tanhq)
UNARY(fabs, fabsq) UNARY(ceil, ceilq) UNARY(floor, floorq)
UNARY(sq, sqq)     UNARY(log10, log10q)

#undef UNARY

// NEED TO FIX THIS.

#define BINARY(FUNC, NAME)					\
  MatrixFrame NAME(MF c, MF a, MF b)				\
  {								\
    sizecheck(c.vol() == a.vol() && a.vol()==b.vol());		\
    for(uint l = 0; l < a.vol(); l++) c(l) = FUNC(a(l), b(l));	\
    return c;							\
  }								\
  MatrixFrame NAME(MF c, double a, MF b)			\
  {								\
    sizecheck(c.vol()==b.vol());					\
    for(uint l = 0; l < b.vol(); l++) c(l) = FUNC(a, b(l));	\
    return c;							\
  }								\
  MatrixFrame NAME(MF c, MF a, double b)			\
  {								\
    sizecheck(c.vol()==a.vol());					\
    for(uint l = 0; l < a.vol(); l++) c(l) = FUNC(a(l), b);	\
    return c;							\
  }								\
  Matrix NAME(MF a, MF b)					\
  {								\
    Matrix c(a);						\
    NAME(c, a, b);						\
    return c;							\
  }								\
  Matrix NAME(double a, MF b)					\
  {								\
    Matrix c(b);						\
    NAME(c, a, b);						\
    return c;							\
  }								\
  Matrix NAME(MF a, double b)					\
  {								\
    Matrix c(a);						\
    NAME(c, a, b);						\
    return c;							\
  }								\

BINARY(min, hmin) BINARY(max, hmax) BINARY(pow, pow)

#undef BINARY

//////////////////////////////////////////////////////////////////////
			       // R //
//////////////////////////////////////////////////////////////////////

Matrix seq(double start, double end, double delta=1.0)
{
  double sign = end - start < 0 ? -1.0 : 1.0;
  delta = sign * fabs(delta);
  double N_double = (end-start) / delta;
  if (N_double < 0) fprintf(stderr, "Problem in seq: N_double < 0");
  int N = floor(N_double) + 1;

  Matrix a(N);
  a(0) = start;
  for (int i=1; i<N; i++)
    a(i) = a(i-1) + delta;

  return a;
}

Matrix rowSums(MF M)
{
  uint nc = M.cols();
  uint nr = M.rows();
  Matrix a(nr);

  for (uint j=0; j < nc; j++)
    for (uint i=0; i < nr; i++)
      a(i) += M(i,j);

  return a;
}

void rowSums(Matrix& a, MF M)
{
  uint nc = M.cols();
  uint nr = M.rows();
  a = M.col(0);

  for (uint j=1; j < nc; j++)
    for (uint i=0; i < nr; i++)
      a(i) += M(i,j);
}

Matrix colSums(MF M)
{
  uint nc = M.cols();
  uint nr = M.rows();
  Matrix a(nc);

  for (uint j=0; j < nc; j++)
    for (uint i=0; i < nr; i++)
      a(j) += M(i,j);

  return a;
}

double maxAll(MF M)
{
  uint N = M.size();
  double mx = M(0);
  for(uint i=1; i<N; i++) {
    mx = M(i) > mx ? M(i) : mx;
  }
  return mx;
}

double minAll(MF M)
{
  uint N = M.size();
  double mx = M(0);
  for(uint i=1; i<N; i++) {
    mx = M(i) < mx ? M(i) : mx;
  }
  return mx;
}

//////////////////////////////////////////////////////////////////////
			   // transpose //
//////////////////////////////////////////////////////////////////////

// A = t(B)
void trans(Matrix& A, Matrix& B)
{
  uint T = B.mats();
  uint R = B.rows();
  uint C = B.cols();
  if(A.rows() != C || A.cols() != R || A.mats() != T) A.resize(C, R, T);
  for(uint t=0; t < T; ++t)
    for(uint i=0; i < R; ++i)
      for(uint j=0; j < C; ++j)
	A(j,i,t) = B(i,j,t);
}

//////////////////////////////////////////////////////////////////////
		   // SIMPLE CONJUGATE GRADIENTS //
//////////////////////////////////////////////////////////////////////

// http://www.cs.cmu.edu/~quake-papers/painless-conjugate-gradient.pdf
// by Jonathan Richard Shewchuk

// Assumes A is positive definite.

int cg(MF& x, const MF& A, const MF& b, double tol, int max_iter)
{
  uint P = b.rows();

  int iter = 0;

  Matrix r(b);
  gemm(r, A, x, 'N', 'N', -1.0, 1.0);

  Matrix d(r);

  double delta_new = dot(r, r);
  // See note below...
  // double delta_0   = delta_new;
  double delta_old = 0;

  Matrix q(P);
  double alpha;
  double beta;

  //while(iter < max_iter && delta_new > tol*delta_0){
  //The above condition is suggested by Shewchuk, but we want an absolute tolerance.
  while(iter < max_iter && delta_new > tol){
    gemm(q, A, d);
    alpha = delta_new / dot(d, q);
    axpy(alpha, d, x);

    // You might want to discard this step.
    if (iter % 50 == 0) {
      r.clone(b);
      gemm(r, A, x, 'N', 'N', -1.0, 1.0);
    }
    else{
      axpy(-1.0 * alpha, q, r);
    }

    delta_old = delta_new;
    delta_new = dot(r, r);

    beta = delta_new / delta_old;
    hsum(d, d, r, beta, 0.0);
    //for(uint j = 0; j < P; j++)
    //  d(j) = beta * d(j) + r(j);

    iter++;
  }

  return iter;
}

//////////////////////////////////////////////////////////////////////
		     // BLAS / LAPACK WRAPPERS //
//////////////////////////////////////////////////////////////////////

// Solve a symmetric, positive definite system of equations.
// ax = b -> b := x;
int symsolve(MF a, Matrix& b, char uplo='L')
{
  // b.fill(0.0);
  // for(uint i = 0; i < b.cols(); ++i)
  //   b(i,i) = 1.0;

  Matrix temp(a);
  int info = posv(temp, b, uplo);

  if (info) {
    fprintf(stderr, "Problem with symsolve: ");
    if (info < 0)
      fprintf(stderr, "%i th argument had illegal value.\n", info);
    if (info > 0)
      fprintf(stderr, "leading minor order %i is not pos. def.\n", info);

    throw std::runtime_error("potrf failed\n");
  }

  return info;
}

int syminv(MF a, Matrix& ainv, char uplo='L')
{
  ainv.resize(a.rows(), a.cols());
  ainv.fill(0.0);
  for(uint i = 0; i < ainv.cols(); ++i)
     ainv(i,i) = 1.0;

  return symsolve(a, ainv, uplo);
}

// Get the Cholesky decomposition of a matrix a.
int chol(Matrix& c, MF a, char uplo='L')
{
  c.clone(a);
  int info = chol(c, uplo);

  // FIX FIX Set other entries to zero.
  // Do I want to have an option for letting a person not do this?

  if (uplo=='L') {
    for (uint j = 1; j < c.cols(); ++j)
      for (uint i = 0; i < j; ++i)
	c(i,j) = 0.0;
  }
  else {
    for (uint j = 0; j < c.cols(); ++j)
      for (uint i = j+1; i < c.cols(); ++i)
	c(i,j) = 0.0;
  }

  if (info) {
    fprintf(stderr, "Problem with chol: ");
    if (info < 0)
      fprintf(stderr, "%i th argument had illegal value.\n", info);
    if (info > 0)
      fprintf(stderr, "leading minor order %i is not pos. def.\n", info);

    throw std::runtime_error("potrf failed\n");
  }

  return info;
}

extern "C" void dsyevd_(char* JOBZ, char* UPLO, int* N, double* A, int* LDA, double* W, double* WORK, int* LWORK, int* IWORK, int* LIWORK, int* INFO);

void dsyevd(char jobz, char uplo, int n, double* a, int lda, double* w, double* work, int lwork, int* iwork, int liwork, int* info)
{ return dsyevd_(&jobz, &uplo, &n, a, &lda, w, work, &lwork, iwork, &liwork, info); }

int symeigen(Matrix& evec, Matrix& eval, MF symmat)
{
  sizecheck(symmat.rows()==symmat.cols());
  int info = 0;
  int N = symmat.rows();
  
  evec.clone(symmat);
  eval.resize(N);

  // int lwork  = 1 + 6 * N + 2*N*N;
  // int liwork = 3 + 5*N;

  int lwork = -1;
  int liwork = -1;

  std::vector<double> work(1);
  std::vector<int>    iwork(1);

  dsyevd('V', 'U', N, &evec(0), N, &eval(0), &work[0], lwork, &iwork[0], liwork, &info);

  lwork = (int) work[0];
  liwork = (int) iwork[0];

  work.resize(lwork);
  iwork.resize(liwork);

  dsyevd('V', 'U', N, &evec(0), N, &eval(0), &work[0], lwork, &iwork[0], liwork, &info);
  
  if (info != 0) {
    fprintf(stderr, "problem in symeigen; info=%i.\n", info);
    throw std::runtime_error("symeigen failed\n");
  }

  return info;
}

int symsqrt(Matrix& rt, MF symmat)
{
  int N = symmat.rows();

  Matrix evec;
  Matrix eval;
  int info = symeigen(evec, eval, symmat);

  Matrix d(N);
  for(int i=0; i<N; i++) d(i) = sqrt(sqrt(eval(i)));

  rt.resize(N, N);
  prodonrow(evec, d);
  gemm(rt, evec, evec, 'N', 'T');

  return info;
}

int syminvsqrt(Matrix& rt, MF symmat)
{
  int N = symmat.rows();

  Matrix evec;
  Matrix eval;
  int info = symeigen(evec, eval, symmat);

  Matrix d(N);
  for(int i=0; i<N; i++) d(i) = sqrt(sqrt(1.0 / eval(i)));

  rt.resize(N, N);
  prodonrow(evec, d);
  gemm(rt, evec, evec, 'N', 'T');

  return info;
}

//--------------------------------------------------------------------

extern "C" void dgesvd_(char* JOBU, char* JOBVT, int* M, int* N, double* A, int* LDA, double* S, double* U, int* LDU, double* VT, int* LDVT, double* WORK, int* LWORK, int* INFO);

void dgesvd(char jobu, char jobvt, int m, int n, double* a, int lda, double* s, double* u, int ldu, double* vt, int ldvt, double*work, int lwork, int* info)
{ return dgesvd_(&jobu, &jobvt, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work, &lwork, info); }

int svd(Matrix& U, Matrix& S, Matrix& tV, Matrix& X)
{
  // Always calculate a thinned U and a non-thinned V.  Pad out the
  // singular values by 0.0 when needed.
  Matrix A(X);

  int m = A.rows();
  int n = A.cols();

  tV.resize(n,n);
  S.resize(n); S.fill(0.0); // The singular values.
  // Should be of size min(m,n), but I want to pad things out.
  
  if (m >= n) U.resize(m,n);
  else U.resize(m,m);
  int ldu = U.rows();

  char jobu  = 'S';
  char jobvt = 'A';

  int maxmn = m > n ? m : n;
  int minmn = m < n ? m : n;

  vector<double> work(1);
  int lwork = -1;

  int info;

  // Workspace query.
  dgesvd(jobu, jobvt, m, n, &A(0), m, &S(0), &U(0), ldu, &tV(0), n, &work[0], lwork, &info);
  printf("lwork: %g\n", work[0]);
  
  // SVD.
  lwork = (int)work[0];
  lwork = lwork > 1 ? lwork : 5 * minmn + maxmn;

  work.resize(lwork);
  dgesvd(jobu, jobvt, m, n, &A(0), m, &S(0), &U(0), ldu, &tV(0), n, &work[0], lwork, &info);

  if (info != 0) {
    fprintf(stderr, "problem in svd; info=%i.\n", info);
    throw std::runtime_error("svd failed\n");
  }

  return info;
}

//////////////////////////////////////////////////////////////////////
			  // END OF CLASS //
//////////////////////////////////////////////////////////////////////

#endif // MATRIX

/*
  class block_iterator
  {

  private:

    int row;
    int col;
    int stride;

    double *begin;
    double *end;
    int     offset;

  public:

    block_iterator(double *p, int r, int c)
      : row(r)
      , col(c)
      , stride(0)
      , begin(p)
      , end(p + row*col)
      , offset(0) {};

    virtual inline double& operator()(int i, int j){
      offset = (row + stride) * j + i;
      return *(begin + offset);
      // ptr = begin + (row + stride) * j + i;
      // return *ptr;
    }

    virtual inline double& operator[](int i){
      offset = (row + stride) * (i / col) + i % col;
      return *(begin + offset);
      // ptr = begin + (row + stride) * (i / col) + i % col;
      // return *ptr;
    }

    virtual inline double& operator++(){
      offset += offset % stride != 0 ? 1 : stride;
      return *(begin + offset);
      // ptr += (ptr - begin) % stride != 0 ? 1 : stride;
      // return *ptr;
    }

  };
 */
