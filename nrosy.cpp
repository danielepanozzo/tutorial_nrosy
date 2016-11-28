#include "nrosy.h"
#include <Eigen/Geometry>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <Eigen/Eigenvalues>
#include <iostream>

using namespace std;
using namespace Eigen;

MatrixXd nrosy
        (
        const MatrixXd& V,          // Vertices of the mesh
        const MatrixXi& F,          // Faces
        const MatrixXi& TT,         // Adjacency triangle-triangle
        const VectorXi& soft_id,    // Soft constraints face ids
        const MatrixXd& soft_value, // Soft constraints 3d vectors
        const int n                 // Degree of the n-rosy field
        )
{
  assert(soft_id.size() > 0); // One constraint is necessary to make the solution unique

  Matrix<double,Dynamic,3> T1(F.rows(),3), T2(F.rows(),3);

  // Compute the local reference systems for each face
  for (unsigned i=0;i<F.rows();++i)
  {
    Vector3d e1 =  V.row(F(i, 1)) - V.row(F(i, 0));
    Vector3d e2 =  V.row(F(i, 2)) - V.row(F(i, 0));
    T1.row(i) = e1.normalized();
    T2.row(i) = T1.row(i).cross(T1.row(i).cross(e2)).normalized();
  }

  // Build the sparse matrix, with an energy term for each edge
  std::vector< Triplet<std::complex<double> > > t;
  std::vector< Triplet<std::complex<double> > > tb;

  unsigned count = 0;
  for (unsigned fi=0;fi<F.rows();++fi)
  {
    for (unsigned ei=0;ei<F.cols();++ei)
    {
      // Look up the opposite face
      int fj = TT(fi,ei);

      // If it is a boundary edge, it does not contribute to the energy
      if (fj == -1) continue;

      // Avoid to count every edge twice
      if (fi > fj) continue;

      // Compute the complex representation of the common edge
      Vector3d v  = (V.row(F(fi,(ei+1)%3)) - V.row(F(fi,ei)));
      Vector2d vi = Vector2d(v.dot(T1.row(fi)),v.dot(T2.row(fi))).normalized();
      std::complex<double> ci(vi(0),vi(1));
      Vector2d vj = Vector2d(v.dot(T1.row(fj)),v.dot(T2.row(fj))).normalized();
      std::complex<double> cj(vj(0),vj(1));

      // Add the term conj(fi)^n*xi - conj(fj)^n*xj to the energy matrix
      t.push_back(Triplet<std::complex<double> >(count,fi,    std::pow(std::conj(ci),n)));
      t.push_back(Triplet<std::complex<double> >(count,fj,-1.*std::pow(std::conj(cj),n)));

      ++count;
    }
  }

  // Convert the constraints into the complex polynomial coefficients and add them as soft constraints
  for (unsigned r=0; r<soft_id.size(); ++r)
  {
    int f = soft_id(r);
    Vector3d v = soft_value.row(r);
    std::complex<double> c(v.dot(T1.row(f)),v.dot(T2.row(f)));
    t.push_back(Triplet<std::complex<double> >(count,f, 1000));
    tb.push_back(Triplet<std::complex<double> >(count,0, std::pow(c,n) * std::complex<double>(1000,0)));
    ++count;
  }

  // Solve the linear system
  typedef SparseMatrix<std::complex<double>> SparseMatrixXcd;
  SparseMatrixXcd A(count,F.rows());
  A.setFromTriplets(t.begin(), t.end());
  SparseMatrixXcd b(count,1);
  b.setFromTriplets(tb.begin(), tb.end());
  SparseLU< SparseMatrixXcd > solver;
  solver.compute(A.adjoint()*A);
  assert(solver.info()==Success);
  MatrixXcd x = solver.solve(A.adjoint()*MatrixXcd(b));
  assert(solver.info()==Success);

  // Convert the interpolated polyvector into Euclidean vectors
  MatrixXd R(F.rows(),3);
  for (int f=0; f<F.rows(); ++f)
  {
    // Find the roots of p(t) = (t - c0) using
    // https://en.wikipedia.org/wiki/Companion_matrix
    Eigen::MatrixXcd M = Eigen::MatrixXcd::Zero(n,n);
    for (int i=1;i<n;++i)
      M(i,i-1) = std::complex<double>(1,0);
    M(0,n-1) = x(f);
    std::complex<double> root = M.eigenvalues()(0);
    R.row(f) = T1.row(f) * root.real() + T2.row(f) * root.imag();
  }

  return R;
}
